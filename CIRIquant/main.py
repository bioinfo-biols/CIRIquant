#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import re
import argparse


def main():
    import circ
    import pipeline
    from logger import get_logger
    from utils import check_file, check_dir, check_config, get_thread_num

    # Init argparser
    parser = argparse.ArgumentParser()

    parser.add_argument('--config', dest='config_file', metavar='CONFIG',
                        help='Config file', )
    parser.add_argument('-1', '--read1', dest='mate1', metavar='MATE1',
                        help='Input mate1 reads (for paired-end data)', )
    parser.add_argument('-2', '--read2', dest='mate2', metavar='MATE2',
                        help='Input mate2 reads (for paired-end data)', )

    parser.add_argument('-v', '--verbose', dest='verbosity', default=False, action='store_true',
                        help='Run in debugging mode', )
    parser.add_argument('--no-gene', dest='gene_exp', default=False, action='store_true',
                        help='Skip stringtie estimation for gene abundance', )
    parser.add_argument('-e', '--log', dest='log_file', default=None, metavar='LOG',
                        help='Log file', )
    parser.add_argument('--circ', dest='circ', metavar='FILE', default=None,
                        help='bed file for putative circRNAs (optional)', )
    parser.add_argument('--bam', dest='bam', metavar='BAM', default=None,
                        help='hisat2 alignment to reference genome', )
    parser.add_argument('-o', '--out', dest='output', metavar='DIR', default=None,
                        help='Output directory (default: current directory)', )
    parser.add_argument('-p', '--prefix', dest='prefix', metavar='PREFIX', default=None,
                        help='Output sample prefix (default: input sample name)', )
    parser.add_argument('-t', '--threads', dest='cpu_threads', default=4, metavar='INT',
                        help='Number of CPU threads', )
    parser.add_argument('-a', '--anchor', dest='anchor', default=5, metavar='INT',
                        help='Minimum anchor length for junction alignment', )
    parser.add_argument('--RNaseR', dest='rnaser', metavar='FILE', default=None,
                        help='CIRIquant result of RNase R sample', )

    args = parser.parse_args()

    # Add lib to PATH
    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    os.chmod(lib_path + '/CIRI2.pl', 0o755)

    # check input reads
    # reads containing input mate pair
    if args.mate1 and args.mate2 and check_file(args.mate1) and check_file(args.mate2):
        reads = [os.path.abspath(args.mate1), os.path.abspath(args.mate2)]
    else:
        sys.exit('No input files specified')

    # optional input: user provided circRNA junction sites
    if args.circ and check_file(args.circ):
        circ_file = os.path.abspath(args.circ)
    else:
        circ_file = None

    # optional input: user provided RNase R CIRIquant results
    if args.rnaser and check_file(args.rnaser):
        rnaser_file = os.path.abspath(args.rnaser)
    else:
        rnaser_file = None

    # hisat2 alignment results
    if args.bam and check_file(args.bam):
        hisat_bam = os.path.abspath(args.bam)
    else:
        hisat_bam = None

    # Output prefix
    if args.prefix is None:
        try:
            prefix = re.search(r'(\S+)[_/-][12]', os.path.basename(reads[0])).group(1)
        except AttributeError:
            sys.exit('Ambiguous sample name, please manually select output prefix')
    else:
        prefix = args.prefix

    # check output dir
    if args.output is None:
        outdir = os.path.abspath('./' + prefix)
    else:
        outdir = os.path.abspath(args.output)
    check_dir(outdir)

    # Parse arguments
    log_file = os.path.abspath(args.log_file) if args.log_file else '{}/{}.log'.format(outdir, prefix)
    verbosity = args.verbosity
    logger = get_logger('CIRIquant', log_file, verbosity)

    if args.config_file and check_file(args.config_file):
        config_file = os.path.abspath(args.config_file)
    else:
        sys.exit('A config file is needed, please see manual for detailed information.')
    config = check_config(config_file)

    # Start Running
    os.chdir(outdir)
    logger.info('Input reads: ' + ','.join([os.path.basename(args.mate1), os.path.basename(args.mate2)]))
    logger.info('Output directory: {}, Output prefix: {}'.format(outdir, prefix))
    logger.info('Config: {} Loaded'.format(config))

    thread = get_thread_num(int(args.cpu_threads))
    anchor = args.anchor

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
    if circ_file is None:
        logger.info('No circRNA information provided, run CIRI2 for junction site prediction ..')
        bwa_sam = pipeline.run_bwa(log_file, thread, reads, outdir, prefix)
        ciri_file = pipeline.run_ciri(log_file, thread, bwa_sam, outdir, prefix)
        circ_file = pipeline.convert_bed(ciri_file)
    else:
        logger.info('Using putative circRNA bed file: {}'.format(os.path.basename(circ_file)))

    # Step4: estimate circRNA expression level
    out_file = circ.proc(log_file, thread, circ_file, hisat_bam, rnaser_file, reads, outdir, prefix, anchor)

    logger.info('circRNA Expression profile: {}'.format(os.path.basename(out_file)))
    logger.info('Finished!')


if __name__ == '__main__':
    main()
