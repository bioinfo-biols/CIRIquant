#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import re
import argparse


def main():
    from version import __version__
    import circ
    import pipeline
    from logger import get_logger
    from utils import check_file, check_dir, check_config, get_thread_num
    from utils import CIRCparser, TOOLS

    # Init argparser
    parser = argparse.ArgumentParser(prog='CIRIquant')

    # required arguments
    parser.add_argument('--config', dest='config_file', metavar='FILE',
                        help='Config file in YAML format', )
    parser.add_argument('-1', '--read1', dest='mate1', metavar='MATE1',
                        help='Input mate1 reads (for paired-end data)', )
    parser.add_argument('-2', '--read2', dest='mate2', metavar='MATE2',
                        help='Input mate2 reads (for paired-end data)', )

    # optional arguments
    parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                        help='Output directory, default: ./', )
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default=None,
                        help='Output sample prefix, default: input sample name', )
    parser.add_argument('-t', '--threads', dest='cpu_threads', default=4, metavar='INT',
                        help='Number of CPU threads, default: 4', )
    parser.add_argument('-a', '--anchor', dest='anchor', default=5, metavar='INT',
                        help='Minimum anchor length for junction alignment, default: 5', )
    parser.add_argument('-l', '--libary-type', dest='library_type', metavar='INT', default=0,
                        help='Library type, 0: unstranded, 1: read1 match the sense strand,'
                             '2: read1 match the antisense strand, default: 0', )

    parser.add_argument('-v', '--verbose', dest='verbosity', default=False, action='store_true',
                        help='Run in debugging mode', )
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('-e', '--log', dest='log_file', default=None, metavar='LOG',
                        help='Log file, default: out_dir/prefix.log', )

    # provide pre-defined list of circRNAs
    parser.add_argument('--bed', dest='bed', metavar='FILE', default=None,
                        help='bed file for putative circRNAs (optional)', )
    parser.add_argument('--circ', dest='circ', metavar='FILE', default=None,
                        help='circRNA prediction results from other softwares', )
    parser.add_argument('--tool', dest='tool', metavar='TOOL', default=None,
                        help='circRNA prediction tool, required if --circ is provided', )

    # when provide RNase R result, do RNase R correction
    parser.add_argument('--RNaseR', dest='rnaser', metavar='FILE', default=None,
                        help='CIRIquant result of RNase R sample', )

    # skip hisat2 alignment for RNA-seq data
    parser.add_argument('--bam', dest='bam', metavar='BAM', default=None,
                        help='hisat2 alignment to reference genome', )

    # skip stringtie prediction
    parser.add_argument('--no-gene', dest='gene_exp', default=False, action='store_true',
                        help='Skip stringtie estimation for gene abundance', )

    args = parser.parse_args()

    """Check required parameters"""
    # check input reads
    if args.mate1 and args.mate2:
        reads = [check_file(args.mate1), check_file(args.mate2)]
    else:
        sys.exit('No input files specified, please see manual for detailed information')

    try:
        lib_type = int(args.library_type)
    except ValueError:
        sys.exit('Wrong library type, please check your command.\nSupported types:\n0 - unstranded;\n'
                 '1 - read1 match the sense strand;\n2 - read1 match the antisense strand;')

    if lib_type not in [0, 1, 2]:
        sys.exit('Wrong library type, please check your command.\nSupported types:\n0 - unstranded;\n'
                 '1 - read1 match the sense strand;\n2 - read1 match the antisense strand;')

    # check configuration
    if args.config_file:
        config = check_config(check_file(args.config_file))
    else:
        sys.exit('A config file is needed, please see manual for detailed information.')

    """Check optional parameters"""
    # use circRNA bed file if provided
    bed_file = check_file(args.bed) if args.bed else None
    circ_file = check_file(args.circ) if args.circ else None
    circ_tool = args.tool

    # user provided RNase R CIRIquant results
    rnaser_file = check_file(args.rnaser) if args.rnaser else None

    # pre aligned hisat2 bam
    hisat_bam = check_file(args.bam) if args.bam else None

    # Output prefix
    if args.prefix is None:
        try:
            prefix = re.search(r'(\S+)[_/-][12]', os.path.basename(reads[0])).group(1)
        except AttributeError:
            sys.exit('Ambiguous sample name, please manually select output prefix')
    else:
        prefix = args.prefix

    # check output dir
    outdir = './' + prefix if args.output is None else args.output
    outdir = check_dir(outdir)

    # Parse arguments
    log_file = os.path.abspath(args.log_file) if args.log_file else '{}/{}.log'.format(outdir, prefix)
    verbosity = args.verbosity
    logger = get_logger('CIRIquant', log_file, verbosity)

    # Add lib to PATH
    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    os.chmod(lib_path + '/CIRI2.pl', 0o755)

    """Start Running"""
    os.chdir(outdir)
    logger.info('Input reads: ' + ','.join([os.path.basename(args.mate1), os.path.basename(args.mate2)]))

    if lib_type == 0:
        lib_name = 'unstranded'
    elif lib_type == 1:
        lib_name = 'ScriptSeq'
    elif lib_type == 2:
        lib_name = 'TAKARA SMARTer'
    else:
        sys.exit('Unsupported library type, please check the manual for instructions.')

    logger.info('Library type: {}'.format(lib_name))
    logger.info('Output directory: {}, Output prefix: {}'.format(outdir, prefix))
    logger.info('Config: {} Loaded'.format(config))

    thread = get_thread_num(int(args.cpu_threads))
    anchor = int(args.anchor)

    # Step1: Data Preparation
    # Step1.1: HISAT2 mapping
    if hisat_bam is None:
        logger.info('Align RNA-seq reads to reference genome ..')
        hisat_bam = pipeline.align_genome(log_file, thread, reads, outdir, prefix)
    else:
        logger.info('HISAT2 alignment bam provided, skipping alignment step ..')
    logger.debug('HISAT2 bam: {}'.format(os.path.basename(hisat_bam)))

    # Step1.2: Estimate Gene Abundance
    if args.gene_exp:
        logger.info('Skipping gene abundance estimation')
    else:
        pipeline.gene_abundance(log_file, thread, outdir, prefix, hisat_bam)

    # Step3: run CIRI2
    if bed_file:
        logger.info('Using user-provided circRNA bed file: {}'.format(bed_file))
    else:
        if circ_file or circ_tool:
            if circ_file and circ_tool:
                logger.info('Using predicted circRNA results from {}: {}'.format(circ_tool, circ_file))
                circ_parser = CIRCparser(circ_file, circ_tool)
            else:
                sys.exit('--circ and --tool must be provided in the same time!')
        else:
            logger.info('No circRNA information provided, run CIRI2 for junction site prediction ..')
            bwa_sam = pipeline.run_bwa(log_file, thread, reads, outdir, prefix)
            ciri_file = pipeline.run_ciri(log_file, thread, bwa_sam, outdir, prefix)
            circ_parser = CIRCparser(ciri_file, 'CIRI2')

        bed_file = '{}/{}.bed'.format(outdir, prefix)
        circ_parser.convert(bed_file)

    # Step4: estimate circRNA expression level
    out_file = circ.proc(log_file, thread, bed_file, hisat_bam, rnaser_file, reads, outdir, prefix, anchor, lib_type)

    # Remove temporary files
    pipeline.clean_tmp(outdir, prefix)

    logger.info('circRNA Expression profile: {}'.format(os.path.basename(out_file)))

    logger.info('Finished!')


if __name__ == '__main__':
    main()
