#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import re
import argparse
from collections import defaultdict


def main():
    from logger import get_logger

    # Init argparser
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', dest='sample_list', metavar='SAMPLE', default=None,
                        help='Sample list', )
    parser.add_argument('-n', dest='control', metavar='CONTROL', default=None,
                        help='control gtf file', )
    parser.add_argument('-c', dest='case', metavar='CASE', default=None,
                        help='case gtf file', )
    parser.add_argument('-l', dest='read_len', metavar='LEN', default=None,
                        help='read length')
    parser.add_argument('-o', '--out-file', dest='output', metavar='OUT', default='CIRI_DE.list',
                        help='Output file (default: CIRI_DE.list)', )

    args = parser.parse_args()

    read_len = int(args.read_len)
    samples = []
    if args.sample_list:
        sample_list = os.path.abspath(args.sample_list)
        with_replicate(sample_list, read_len)
    elif args.control and args.case:
        pass
    else:
        sys.exit('Please provide input files')


def with_replicate(sample_lst, read_len):
    from math import ceil
    # lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
    # os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
    # os.chmod(lib_path + '/CIRI_DE.Rscript', 0o755)
    samples = []
    with open(sample_lst, 'r') as f:
        for line in f:
            samples.append(tuple(line.rstrip().split()))

    # ID,Type,Mapped_reads
    design_matrix = []

    g_matrix = defaultdict(dict)
    t_matrix = defaultdict(dict)
    gene_IDs = {}

    c_matrix = []
    j_matrix = []
    for sample, sample_type, prefix in samples:
        gene_f = '{}/gene/{}_out.gtf'.format(prefix, sample)
        circ_f = '{}/{}.gtf'.format(prefix, sample)

        with open(gene_f, 'r') as f:
            t_len = 0
            for line in f:
                if line.starswith('#'):
                    continue
                content = line.rstrip().split('\t')
                if content[2] == 'transcript':
                    if t_len > 0:
                        t_matrix[t_id].setdefault(sample, int(ceil(cov * t_len / read_len)))
                    t_id = RE_TRANSCRIPT_ID.search(content[-1]).group(1)
                    g_id = get_gene_ID(content[-1], chrom, t_id)
                    gene_IDs[t_id] = g_id
                    cov = get_cov(content[-1])
                    t_len = 0
                if content[2] == 'exon':
                    t_len += int(content[4]) - int(content[3]) + 1
            t_matrix[t_id].setdefault(sample, int(ceil(cov * t_len / read_len)))

        with open(circ_f, 'r') as f:
            header = {}
            for line in f:
                if line.startswith('##'):
                    key, value = line.strip('#').split(':')
                    header.update({key: value})
                    continue
                content = line.rstrip().split('\t')
                if content[2] == 'circRNA':
                    continue


    for t_id, t_exp in t_matrix.iteritems():
        for sample in t_exp:
            g_matrix[gene_IDs[t_id]][sample] = g_matrix[gene_IDs[t_id]].setdefault(sample, 0) + t_exp[sample]


RE_GENE_ID=re.compile('gene_id "([^"]+)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
RE_STRING=re.compile(re.escape(opts.string))
# t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)

def get_gene_ID(attr, chrom, t_id):
    m = RE_GENE_ID.search(attr)
    if m:
        return m.group(1)
    m = RE_GENE_NAME.search(attr)
    if m:
        return chrom + '|' + m.group(1)
    return t_id


def get_cov(attr):
    m = RE_COVERAGE.search(s)
    cov = max(float(m.group(1)), 0.0) if m else 0.0
    return cov
