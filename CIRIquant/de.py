#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import argparse
from multiprocessing import Pool
from collections import namedtuple

import numpy as np
np.random.seed(5)
import numexpr as ne
ne.set_num_threads(4)

from logger import ProgressBar

LOGGER = None
CIRC = namedtuple('CIRC', 'bsj fsj ratio rnaser_bsj rnaser_fsj')


def main():
    global LOGGER
    from circ import grouper
    from logger import get_logger
    from utils import check_file, check_dir

    # Init argparser
    parser = argparse.ArgumentParser()

    # Sample with replicates
    parser.add_argument('-i', dest='sample_list', metavar='SAMPLE', default=None,
                        help='Sample list', )

    # Sample without replicate
    parser.add_argument('-n', dest='control', metavar='CONTROL', default=None,
                        help='control gtf file', )
    parser.add_argument('-c', dest='case', metavar='CASE', default=None,
                        help='case gtf file', )

    # Optional parameters
    parser.add_argument('-l', dest='read_len', metavar='INT', default=None,
                        help='read length')
    args = parser.parse_args()

    thread = 4

    out_dir = check_dir('./CIRI_DE')
    de_file = out_dir + '/differential_expression.csv'
    ds_file = out_dir + '/differential_splicing.csv'
    de_mtx = out_dir + '/bsj_reads.mtx'
    ds_mtx = out_dir + '/junction_ratio.mtx'

    LOGGER = get_logger('CIRIquant', './CIRI_DE/run.log', True)

    read_len = int(args.read_len)
    size = 100000

    # Load GTF
    case_gtf, ctrl_gtf = check_file(args.case), check_file(args.control)

    case_header, case_data = load_gtf(case_gtf)
    ctrl_header, ctrl_data = load_gtf(ctrl_gtf)

    factor = float(case_header['Mapped_Reads']) / float(ctrl_header['Mapped_Reads'])

    # Multiprocessing
    all_circ = set(case_data.keys() + ctrl_data.keys())

    pool = Pool(thread, score_initializer, (case_data, ctrl_data, case_header, ctrl_header, size))
    jobs = []
    chunk_size = max(1000, len(all_circ) / thread + 1)
    for circ_chunk in grouper(all_circ, chunk_size):
        jobs.append(pool.apply_async(score_worker, (circ_chunk, factor, )))
    pool.close()
    pool.join()

    circ_scores = dict()
    for job in jobs:
        tmp_score = job.get()
        circ_scores.update(tmp_score)

    with open('./CIRI_DE/CIRI_DE.results.csv', 'w') as out:
        out.write('circRNA_ID\tDE_score\tDS_score\n')
        for circ_id in all_circ:
            case_bsj = int(case_data[circ_id].bsj) if circ_id in case_data else 0
            ctrl_bsj = int(ctrl_data[circ_id].bsj) if circ_id in ctrl_data else 0

            case_fsj = int(case_data[circ_id].fsj) if circ_id in case_data else 0
            ctrl_fsj = int(ctrl_data[circ_id].fsj) if circ_id in ctrl_data else 0

            tmp_de, tmp_ds = circ_scores[circ_id]

            out.write('{}\t{}\t{}\n'.format(circ_id, tmp_de, tmp_ds))

    LOGGER.info('Finished!')


CASE = None
CTRL = None
CASE_PRIOR = []
CTRL_PRIOR = []
SIZE = 10000


def score_initializer(case_data, case_header, ctrl_data, ctrl_header, size):
    global CASE, CTRL, CASE_PRIOR, CTRL_PRIOR
    global SIZE
    CASE, CTRL = case_data, ctrl_data
    CASE_PRIOR = gmm_sampling(case_header)
    CTRL_PRIOR = gmm_sampling(ctrl_header)
    SIZE = size


def score_worker(circ_ids, factor):
    score_data = {}
    for circ_id in circ_ids:
        if circ_id in CASE:
            if CASE[circ_id].rnaser_fsj and CASE[circ_id].fsj != 0:
                case_dis = prior_sampling(CASE_PRIOR, CASE[circ_id].rnaser_bsj, CASE[circ_id].rnaser_fsj, CASE[circ_id].fsj)
            else:
                case_dis = np.random.gamma(shape=int(CASE[circ_id].bsj) + 1 + 1, size=SIZE)
        else:
            case_dis = np.random.gamma(shape=1, size=SIZE)

        if circ_id in CTRL:
            if CTRL[circ_id].rnaser_fsj and CTRL[circ_id].fsj != 0:
                ctrl_dis = prior_sampling(CTRL_PRIOR, CTRL[circ_id].rnaser_bsj, CTRL[circ_id].rnaser_fsj, CTRL[circ_id].fsj)
            else:
                ctrl_dis = np.random.gamma(shape=int(CTRL[circ_id].bsj) + 1 + 1, size=SIZE)
        else:
            ctrl_dis = np.random.gamma(shape=1, size=SIZE)

        case_bsj = int(CASE[circ_id].bsj) if circ_id in CASE else 0
        ctrl_bsj = int(CTRL[circ_id].bsj) if circ_id in CTRL else 0



        tmp_de = de_score(max(case_bsj, 1), max(ctrl_bsj, 1), 1.0 / factor)

        case_fsj = int(CASE[circ_id].fsj) if circ_id in CASE else 0
        ctrl_fsj = int(CTRL[circ_id].fsj) if circ_id in CTRL else 0
        tmp_ds = ds_score(max(case_bsj, 1), max(case_fsj, 1), max(ctrl_bsj, 1), max(ctrl_fsj, 1))

        score_data[circ_id] = (tmp_de, tmp_ds)
    return score_data


def load_gtf(in_file):
    from circ import GTFParser

    LOGGER.info('Loading GTF {}'.format(in_file))

    circ_data = {}
    with open(in_file, 'r') as f:
        header = {}
        for line in f:
            if line.startswith('##'):
                key, value = line.strip('#').split(':')
                header.update({key: value})
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
    return header, circ_data


def de_score(a, b, factor=1, size=SIZE, pval=0.05):
    x1 = np.random.gamma(shape=a + 1, size=size)
    x2 = np.random.gamma(shape=b + 1, size=size)

    z = ne.evaluate("x1 / x2")
    z.sort()
    if np.mean(z) * factor >= 1:
        score = max(np.log2(z[int(round(size * pval))] * factor), 0)
    else:
        score = min(np.log2(z[int(round(size * (1 - pval)))] * factor), 0)
    return score


def de_score_corrected(dis_a, dis_b, factor=1, size=SIZE, pvalue=0.05):
    # from itertools import izip

    fc = ne.evaluate("dis_a / dis_b")
    # fc = sorted([np.log(j / i) for i, j in izip(dis_a, dis_b)])
    miu = np.mean(fc) + factor
    assert len(dis_a) > 0 and len(dis_b) > 0
    if miu >= 0:
        score = max((fc[int(size * pvalue)] + fc[int(size * pvalue) + 1]) / 2 + factor, 0)
    else:
        score = min((fc[-int(size * pvalue)] + fc[-int(size * pvalue) + 1]) / 2 + factor, 0)
    return score


def ds_score(a1, b1, a2, b2, size=SIZE, pval=0.05):
    x1 = np.random.beta(a1, b1, size=size)
    x2 = np.random.beta(a2, b2, size=size)

    z = ne.evaluate("x1 / x2")
    z.sort()
    if np.mean(z) >= 1:
        score = max(np.log2(z[int(round(size * pval))]), 0)
    else:
        score = min(np.log2(z[int(round(size * (1 - pval)))]), 0)
    return score


def gmm_sampling(header, size=SIZE):
    n = int(header['N'])
    w = header['W'].split(',')
    m = header['M'].split(',')
    sd = header['SD'].split(',')

    sample = np.array([])
    for i in range(n):
        sample = np.append(sample, np.random.normal(float(m[i]), float(sd[i]), int(size * float(w[i]))))
    return sample


def prior_sampling(sample, bsj, fsj, total_fsj):
    factor = depth_factor(bsj / (bsj + 0.5 * fsj))
    corrected_read = sample * 0.5 * factor * total_fsj
    return [np.random.gamma(shape=x + 1) for x in corrected_read if x > 0]


def depth_factor(r):
    return r / (1.0 - r)


# def with_replicate(sample_lst, read_len):
#     from math import ceil
#     # lib_path = os.path.dirname(os.path.split(os.path.realpath(__file__))[0]) + '/libs'
#     # os.environ['PATH'] = lib_path + ':' + os.environ['PATH']
#     # os.chmod(lib_path + '/CIRI_DE.Rscript', 0o755)
#     samples = []
#     with open(sample_lst, 'r') as f:
#         for line in f:
#             samples.append(tuple(line.rstrip().split()))
#
#     # ID,Type,Mapped_reads
#     design_matrix = []
#
#     g_matrix = defaultdict(dict)
#     t_matrix = defaultdict(dict)
#     gene_IDs = {}
#
#     c_matrix = []
#     j_matrix = []
#     for sample, sample_type, prefix in samples:
#         gene_f = '{}/gene/{}_out.gtf'.format(prefix, sample)
#         circ_f = '{}/{}.gtf'.format(prefix, sample)
#
#         with open(gene_f, 'r') as f:
#             t_len = 0
#             for line in f:
#                 if line.starswith('#'):
#                     continue
#                 content = line.rstrip().split('\t')
#                 if content[2] == 'transcript':
#                     if t_len > 0:
#                         t_matrix[t_id].setdefault(sample, int(ceil(cov * t_len / read_len)))
#                     t_id = RE_TRANSCRIPT_ID.search(content[-1]).group(1)
#                     g_id = get_gene_ID(content[-1], chrom, t_id)
#                     gene_IDs[t_id] = g_id
#                     cov = get_cov(content[-1])
#                     t_len = 0
#                 if content[2] == 'exon':
#                     t_len += int(content[4]) - int(content[3]) + 1
#             t_matrix[t_id].setdefault(sample, int(ceil(cov * t_len / read_len)))
#
#         with open(circ_f, 'r') as f:
#             header = {}
#             for line in f:
#                 if line.startswith('##'):
#                     key, value = line.strip('#').split(':')
#                     header.update({key: value})
#                     continue
#                 content = line.rstrip().split('\t')
#                 if content[2] == 'circRNA':
#                     continue
#
#
#     for t_id, t_exp in t_matrix.iteritems():
#         for sample in t_exp:
#             g_matrix[gene_IDs[t_id]][sample] = g_matrix[gene_IDs[t_id]].setdefault(sample, 0) + t_exp[sample]
#
#
# RE_GENE_ID=re.compile('gene_id "([^"]+)"')
# RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
# RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
# RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
# RE_STRING=re.compile(re.escape(opts.string))
# # t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
#
# def get_gene_ID(attr, chrom, t_id):
#     m = RE_GENE_ID.search(attr)
#     if m:
#         return m.group(1)
#     m = RE_GENE_NAME.search(attr)
#     if m:
#         return chrom + '|' + m.group(1)
#     return t_id
#
#
# def get_cov(attr):
#     m = RE_COVERAGE.search(s)
#     cov = max(float(m.group(1)), 0.0) if m else 0.0
#     return cov


if __name__ == '__main__':
    main()
