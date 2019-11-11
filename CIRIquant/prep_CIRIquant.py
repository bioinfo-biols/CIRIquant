#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import argparse
import logging
from collections import namedtuple
from version import __version__

LOGGER = logging.getLogger('prep_CIRIquant')
CIRC = namedtuple('CIRC', 'bsj fsj ratio rnaser_bsj rnaser_fsj')
INFO = namedtuple('INFO', 'circ_type gene_id gene_name gene_type')


def load_gtf(in_file):
    from circ import GTFParser
    LOGGER.info('Loading CIRIquant result: {}'.format(in_file))

    circ_data = {}
    circ_info = {}
    with open(in_file, 'r') as f:
        header = {}
        for line in f:
            if line.startswith('##'):
                key, value = line.rstrip().strip('#').split(':')
                header.update({key.strip(): value.strip()})
                continue

            content = line.rstrip().split('\t')
            tmp_parser = GTFParser(content)
            circ_data[tmp_parser.attr['circ_id']] = CIRC(
                float(tmp_parser.attr['bsj']),
                float(tmp_parser.attr['fsj']),
                float(tmp_parser.attr['junc_ratio']),
                float(tmp_parser.attr['rnaser_bsj']) if 'rnaser_bsj' in tmp_parser.attr else None,
                float(tmp_parser.attr['rnaser_fsj']) if 'rnaser_fsj' in tmp_parser.attr else None,
            )
            circ_info[tmp_parser.attr['circ_id']] = INFO(
                tmp_parser.attr['circ_type'],
                tmp_parser.attr['gene_id'] if 'gene_id' in tmp_parser.attr else 'NA',
                tmp_parser.attr['gene_name'] if 'gene_name' in tmp_parser.attr else 'NA',
                tmp_parser.attr['gene_type'] if 'gene_type' in tmp_parser.attr else 'NA',
            )
    return header, circ_data, circ_info


def main():
    global LOGGER
    from logger import get_logger
    from utils import check_file

    # Init argparser
    parser = argparse.ArgumentParser(prog="prep_CIRIquant")

    parser.add_argument('-i', dest='input', metavar='FILE', required=True,
                        help='input sample list', )
    parser.add_argument('--lib', dest='lib', metavar='FILE', required=True,
                        help='output file of library information')
    parser.add_argument('--circ', dest='circ', metavar='FILE', required=True,
                        help='output file of circRNA annotation information')
    parser.add_argument('--bsj', dest='bsj', metavar='FILE', required=True,
                        help='output file of circRNA bsj reads number', )
    parser.add_argument('--ratio', dest='ratio', metavar='FILE', required=True,
                        help='output file of circRNA junction ratio matrix', )

    # Optional parameters
    args = parser.parse_args()
    LOGGER = get_logger('prep_CIRIquant', None, False)

    # Load GTF
    sample_lst = check_file(args.input)
    lib_file, info_file = os.path.abspath(args.lib), os.path.abspath(args.circ)
    bsj_file, ratio_file = os.path.abspath(args.bsj), os.path.abspath(args.ratio)

    LOGGER.info('Input file name: {}'.format(os.path.basename(sample_lst)))
    LOGGER.info('Output library information: {}'.format(lib_file))
    LOGGER.info('Output circRNA annotation: {}'.format(info_file))
    LOGGER.info('Output BSJ matrix: {}'.format(bsj_file))
    LOGGER.info('Output Ratio matrix: {}'.format(ratio_file))

    all_circ = {}
    all_sample = {}
    all_data = {}

    is_paired = 0
    with open(sample_lst, 'r') as f:
        for line in f:
            content = line.rstrip().split()
            sample, sample_file, group = content[0], content[1], content[2]
            sample_header, sample_data, sample_info = load_gtf(sample_file)
            all_sample[sample] = sample_header
            all_sample[sample]['Group'] = group
            if len(content) > 3:
                all_sample[sample]['Subject'] = content[3]
                is_paired = 1
            all_circ.update(sample_info)
            all_data[sample] = sample_data

    sample_names = sorted(all_sample.keys())
    circ_ids = sorted(all_circ.keys())

    # library information and circRNA annotation
    with open(lib_file, 'w') as lib_out, open(info_file, 'w') as info_out:
        info_out.write('circ_id,circ_type,gene_id,gene_name,gene_type\n')
        for circ_id in circ_ids:
            tmp_circ = all_circ[circ_id]
            tmp_line = [circ_id, tmp_circ.circ_type, tmp_circ.gene_id, tmp_circ.gene_name, tmp_circ.gene_type]
            info_out.write(','.join(['"{}"'.format(x) for x in tmp_line]) + '\n')

        if is_paired == 0:
            lib_out.write('Sample,Total,Mapped,Circular,Group\n')
        else:
            lib_out.write('Sample,Total,Mapped,Circular,Group,Subject\n')

        for sample in sample_names:
            tmp_sample = all_sample[sample]
            tmp_line = [sample, tmp_sample['Total_Reads'], tmp_sample['Mapped_Reads'],
                        tmp_sample['Circular_Reads'], tmp_sample['Group'], ]
            if is_paired != 0:
                if 'Subject' in tmp_sample:
                    tmp_line.append(tmp_sample['Subject'])
                else:
                    LOGGER.error('No subject ID found for {}, please check your input'.format(sample))
                    sys.exit()
            lib_out.write(','.join(tmp_line) + '\n')

    # BSJ and junction ratio
    with open(bsj_file, 'w') as bsj_out, open(ratio_file, 'w') as ratio_out:
        tmp_header = ["circ_id", ] + sample_names
        bsj_out.write(','.join(tmp_header) + '\n')
        ratio_out.write(','.join(tmp_header) + '\n')

        for circ_id in circ_ids:
            tmp_bsj, tmp_ratio = [circ_id, ], [circ_id, ]
            for sample in sample_names:
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
