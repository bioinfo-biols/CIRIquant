#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import argparse
import logging
from version import __version__

LOGGER = logging.getLogger('prep_CIRIquant')


def main():
    global LOGGER
    from logger import get_logger
    from de import load_gtf
    from utils import check_file

    # Init argparser
    parser = argparse.ArgumentParser(prog="CIRIquant")

    parser.add_argument('-i', dest='input', metavar='FILE', required=True,
                        help='input sample list', )
    parser.add_argument('--bsj', dest='bsj', metavar='FILE', required=True,
                        help='output file of circRNA bsj reads number', )
    parser.add_argument('--ratio', dest='ratio', metavar='FILE', required=True,
                        help='output file of circRNA junction ratio matrix', )

    # Optional parameters
    args = parser.parse_args()
    LOGGER = get_logger('prep_CIRIquant', None, False)

    # Load GTF
    sample_lst = check_file(args.input)
    bsj_file, ratio_file = os.path.abspath(args.bsj), os.path.abspath(args.ratio)

    LOGGER.info('Input file name: {}'.format(os.path.basename(sample_lst)))
    LOGGER.info('Output BSJ matrix: {}'.format(bsj_file))
    LOGGER.info('Output Ratio matrix: {}'.format(ratio_file))

    all_circ = {}
    all_sample = []
    all_data = {}

    with open(sample_lst, 'r') as f:
        for line in f:
            sample, sample_file = line.rstrip().split()
            all_sample.append(sample)
            sample_header, sample_data = load_gtf(sample_file)
            all_circ.update({i: 1 for i in sample_data})
            all_data[sample] = sample_data

    with open(bsj_file, 'w') as bsj_out, open(ratio_file, 'w') as ratio_out:
        tmp_header = ["", ] + all_sample
        bsj_out.write(','.join(tmp_header) + '\n')
        ratio_out.write(','.join(tmp_header) + '\n')

        for circ_id in all_circ:
            tmp_bsj, tmp_ratio = [circ_id, ], [circ_id, ]
            for sample in all_sample:
                if circ_id in all_data[sample]:
                    tmp_bsj.append(all_data[sample][circ_id].bsj)
                    tmp_ratio.append(all_data[sample][circ_id].ratio)
                else:
                    tmp_bsj.append(0)
                    tmp_ratio.append(0)
            bsj_out.write(','.join([str(x) for x in tmp_bsj]) + '\n')
            ratio_out.write(','.join([str(x) for x in tmp_ratio]) + '\n')

    LOGGER.info('Finished!')


if __name__ == '__main__':
    main()
