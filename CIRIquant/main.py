#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import re
import argparse
import simplejson as json


def main():
    import align
    import pipeline
    import circ
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

    args = parser.parse_args()

    # Add lib to PATH
    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    os.chmod(lib_path + '/CIRI2.pl', 0o755)

    # check input reads
    if args.mate1 and args.mate2 and check_file(args.mate1) and check_file(args.mate2):
        reads = [os.path.abspath(args.mate1), os.path.abspath(args.mate2)]
    else:
        sys.exit('No input files specified')

    # optional input
    if args.circ and check_file(args.circ):
        circ_file = os.path.abspath(args.circ)
    else:
        circ_file = None

    if args.bam and check_file(args.bam):
        hisat_bam = os.path.abspath(args.circ)
    else:
        hisat_bam = None

    # check output dir
    if args.prefix is None:
        try:
            prefix = re.search(r'(\S+)[_/-][12]', os.path.basename(reads[0])).group(1)
        except AttributeError:
            sys.exit('Ambiguous input file name, please manually select output prefix')
    else:
        prefix = args.prefix

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

    # Start Running
    os.chdir(outdir)
    logger.info('Input reads: ' + ','.join([os.path.basename(args.mate1), os.path.basename(args.mate2)]))
    logger.info('Output directory: ' + outdir)
    logger.info('Sample prefix: ' + prefix)

    config = check_config(config_file)
    thread = get_thread_num(int(args.cpu_threads))

    stat = {}
    # Step1: Data Preparation
    # Step1.1: HISAT2 mapping
    if hisat_bam is None:
        logger.info('Align RNA-seq reads to reference genome ..')
        hisat_bam = pipeline.align_genome(log_file, thread, reads, outdir, prefix, config)
    else:
        logger.info('HISAT2 alignment bam provided, skipping alignment step ..')
    logger.debug('HISAT2 bam: {}'.format(os.path.basename(hisat_bam)))

    # Step1.2: Estimate Gene Abundance
    if args.gene_exp:
        logger.info('Skipping gene abundance estimation')
    else:
        pipeline.gene_abundance(log_file, thread, outdir, prefix, hisat_bam, config)

    # Step2: candidate reads
    if args.bam:
        logger.info('Candidate reads provided, skip scanning step ..')
        cand_reads = reads
    else:
        logger.info('Extracting circular candidate reads ..')
        cand_reads, tmp = align.unmapped_reads(log_file, thread, outdir, prefix, hisat_bam, config)
        stat.update(tmp)
    logger.debug('circRNA reads: ' + ','.join(cand_reads))

    # Step3: run CIRI2
    if circ_file is None:
        logger.info('No circRNA information provided, run prediction step ..')
        bwa_sam = pipeline.run_bwa(log_file, thread, cand_reads, outdir, prefix, config)
        ciri_file = pipeline.run_ciri(log_file, thread, bwa_sam, outdir, prefix, config)
        circ_file = pipeline.convert_bed(ciri_file)
    else:
        logger.info('Using putative circRNA bed file: {}'.format(os.path.basename(circ_file)))

    # Step4: estimate circRNA expression level
    output_file, tmp = circ.proc(log_file, thread, circ_file, hisat_bam, cand_reads, outdir, prefix, config)
    stat.update(tmp)

    stat_file = '{}/{}.stat'.format(outdir, prefix)
    with open(stat_file, 'w') as out:
        json.dump(stat, out)

    logger.info('Finished, see {} for circRNA expression profile'.format(output_file))


if __name__ == '__main__':
    main()
