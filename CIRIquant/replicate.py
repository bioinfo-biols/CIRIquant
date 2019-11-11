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
    parser = argparse.ArgumentParser(prog="prep_CIRIquant")

    parser.add_argument('--lib', dest='lib', metavar='FILE', required=True,
                        help='library information')
    parser.add_argument('--bsj', dest='bsj', metavar='FILE', required=True,
                        help='circRNA expression matrix', )
    parser.add_argument('--gene', dest='gene', metavar='FILE', required=True,
                        help='gene expression matrix', )
    parser.add_argument('--out', dest='out', metavar='FILE', required=True,
                        help='output result of differential expression analysis', )

    # Optional parameters
    args = parser.parse_args()
    LOGGER = get_logger('CIRI_DE', None, False)

    lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    os.chmod(lib_path + '/CIRI_DE.R', 0o755)

    # Load GTF
    lib_file = check_file(args.lib)
    bsj_file = check_file(args.bsj)
    gene_file = check_file(args.gene)
    out_file = os.path.abspath(args.out)

    LOGGER.info('Library information: {}'.format(lib_file))
    LOGGER.info('circRNA expression matrix: {}'.format(bsj_file))
    LOGGER.info('gene expression matrix: {}'.format(gene_file))
    LOGGER.info('Output DE results: {}'.format(out_file))

    de_cmd = 'CIRI_DE.R --lib={} --bsj={} --gene={} --out={}'.format(lib_file, bsj_file, gene_file, out_file)
    subprocess.call(de_cmd, shell=True)
    LOGGER.info('Finished!')


if __name__ == '__main__':
    main()
