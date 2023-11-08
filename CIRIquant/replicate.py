#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import argparse
import logging
import subprocess
from version import __version__

LOGGER = logging.getLogger('CIRI_DE')


def main():
    global LOGGER
    from logger import get_logger
    from utils import check_file

    # Init argparser
    parser = argparse.ArgumentParser(prog="CIRIquant_DE_replicate")

    parser.add_argument('--lib', dest='lib', metavar='FILE', required=True,
                        help='library information')
    parser.add_argument('--bsj', dest='bsj', metavar='FILE', required=True,
                        help='circRNA expression matrix', )
    parser.add_argument('--gene', dest='gene', metavar='FILE', required=True,
                        help='gene expression matrix', )
    parser.add_argument('--out', dest='out', metavar='FILE', required=True,
                        help='output result of circRNA differential expression analysis', )
    parser.add_argument('--out2', dest='out2', metavar='FILE', required=True,
                        help='output result of gene differential expression analysis', )
    # Optional parameters
    args = parser.parse_args()
    LOGGER = get_logger('CIRI_DE', None, False)

    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']

    # Load GTF
    lib_file = check_file(args.lib)
    bsj_file = check_file(args.bsj)
    gene_file = check_file(args.gene)
    out_file = os.path.abspath(args.out)
    out_file2 = os.path.abspath(args.out2)

    LOGGER.info('Library information: {}'.format(lib_file))
    LOGGER.info('circRNA expression matrix: {}'.format(bsj_file))
    LOGGER.info('gene expression matrix: {}'.format(gene_file))
    LOGGER.info('Output circ DE results: {}'.format(out_file))
    LOGGER.info('Output gene DE results: {}'.format(out_file2))

    de_cmd = 'Rscript {}/CIRI_DE.R --lib={} --bsj={} --gene={} --out={} --out2={}'.format(lib_path, lib_file, bsj_file, gene_file, out_file, out_file2)
    subprocess.call(de_cmd, shell=True)
    LOGGER.info('Finished!')


if __name__ == '__main__':
    main()
