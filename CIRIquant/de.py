#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import argparse
import logging
from multiprocessing import Pool
from collections import namedtuple

import numpy as np
np.random.seed(5)
import numexpr as ne
ne.set_num_threads(4)

from version import __version__

LOGGER = logging.getLogger('CIRI_DE')
CIRC = namedtuple('CIRC', 'bsj fsj ratio rnaser_bsj rnaser_fsj')


def main():
    global LOGGER
    from circ import grouper
    from logger import get_logger
    from utils import check_file, get_thread_num

    # Init argparser
    parser = argparse.ArgumentParser(prog="CIRIquant")

    parser.add_argument('-n', dest='control', metavar='CONTROL', default=None,
                        help='control gtf file', )
    parser.add_argument('-c', dest='case', metavar='CASE', default=None,
                        help='case gtf file', )
    parser.add_argument('-o', dest='out', metavar='FILE', default='CIRI_DE.csv',
                        help='output file name', )
    parser.add_argument('-p', dest='pval', metavar='FLOAT', default=0.05,
                        help='P value threshold for DE and DS score calculation', )
    parser.add_argument('-t', '--threads', dest='cpu_threads', default=4, metavar='INT',
                        help='Number of CPU threads, default: 4', )

    # Optional parameters
    args = parser.parse_args()

    thread = get_thread_num(int(args.cpu_threads))
    LOGGER = get_logger('CIRI_DE', None, False)

    size = 100000

    pval = float(args.pval)
    assert 0 < pval < 0.5

    # Load GTF
    case_gtf, ctrl_gtf = check_file(args.case), check_file(args.control)
    out_csv = os.path.abspath(args.out)

    LOGGER.info('Experiment: {}'.format(os.path.basename(case_gtf)))
    LOGGER.info('Control: {}'.format(os.path.basename(ctrl_gtf)))
    LOGGER.info('Threads: {}'.format(thread))
    LOGGER.info('Threshold: {}'.format(args.pval))

    case_header, case_data = load_gtf(case_gtf)
    ctrl_header, ctrl_data = load_gtf(ctrl_gtf)

    factor = float(case_header['Mapped_Reads']) / float(ctrl_header['Mapped_Reads'])

    # Multiprocessing
    all_circ = set(case_data.keys() + ctrl_data.keys())

    jobs = []
    chunk_size = max(1000, len(all_circ) / thread + 1)

    if 'N' in case_header and 'N' in ctrl_header:
        pool = Pool(thread, correction_initializer, (case_data, case_header, ctrl_data, ctrl_header, size, pval))
        for circ_chunk in grouper(all_circ, chunk_size):
            jobs.append(pool.apply_async(correction_worker, (circ_chunk, factor, )))
    else:
        pool = Pool(thread, score_initializer, (case_data, ctrl_data, size, pval))
        for circ_chunk in grouper(all_circ, chunk_size):
            jobs.append(pool.apply_async(score_worker, (circ_chunk, factor, )))

    pool.close()
    pool.join()

    circ_scores = dict()
    for job in jobs:
        tmp_score = job.get()
        circ_scores.update(tmp_score)

    LOGGER.info('Output csv: {}'.format(out_csv))
    with open(out_csv, 'w') as out:
        out.write('circRNA_ID\tCase_BSJ\tCase_FSJ\tCase_Ratio\tCtrl_BSJ\tCtrl_FSJ\tCtrl_Ratio\tDE_score\tDS_score\n')
        for circ_id in all_circ:
            case_bsj = int(case_data[circ_id].bsj) if circ_id in case_data else 0
            case_fsj = int(case_data[circ_id].fsj) if circ_id in case_data else 0
            case_ratio = float(case_data[circ_id].ratio) if circ_id in case_data else 0

            ctrl_bsj = int(ctrl_data[circ_id].bsj) if circ_id in ctrl_data else 0
            ctrl_fsj = int(ctrl_data[circ_id].fsj) if circ_id in ctrl_data else 0
            ctrl_ratio = float(ctrl_data[circ_id].ratio) if circ_id in ctrl_data else 0

            tmp_de, tmp_ds = circ_scores[circ_id]

            tmp_line = [circ_id, case_bsj, case_fsj, case_ratio, ctrl_bsj, ctrl_fsj, ctrl_ratio, tmp_de, tmp_ds]

            out.write('\t'.join([str(x) for x in tmp_line]) + '\n')

    LOGGER.info('Finished!')


CASE = None
CTRL = None
CASE_PRIOR = []
CTRL_PRIOR = []
SIZE = 100000
PVAL = 0.05


def score_initializer(case_data, ctrl_data, size, pval):
    global CASE, CTRL
    global SIZE, PVAL
    CASE, CTRL = case_data, ctrl_data
    SIZE = size
    PVAL = pval


def correction_initializer(case_data, case_header, ctrl_data, ctrl_header, size, pval):
    global CASE, CTRL, CASE_PRIOR, CTRL_PRIOR
    global SIZE, PVAL
    CASE, CTRL = case_data, ctrl_data
    CASE_PRIOR = gmm_sampling(case_header)
    CTRL_PRIOR = gmm_sampling(ctrl_header)
    SIZE = size
    PVAL = pval


def score_worker(circ_ids, factor):
    score_data = {}
    for circ_id in circ_ids:
        case_bsj = int(CASE[circ_id].bsj) if circ_id in CASE else 0
        ctrl_bsj = int(CTRL[circ_id].bsj) if circ_id in CTRL else 0

        tmp_de = de_score(max(case_bsj, 1), max(ctrl_bsj, 1), 1.0 / factor, size=SIZE, pval=PVAL)

        case_fsj = int(CASE[circ_id].fsj) if circ_id in CASE else 0
        ctrl_fsj = int(CTRL[circ_id].fsj) if circ_id in CTRL else 0
        tmp_ds = ds_score(max(case_bsj, 1), max(case_fsj, 1),
                          max(ctrl_bsj, 1), max(ctrl_fsj, 1), size=SIZE, pval=PVAL)

        score_data[circ_id] = (tmp_de, tmp_ds)
    return score_data


def correction_worker(circ_ids, factor):
    score_data = {}
    for circ_id in circ_ids:
        if circ_id in CASE:
            if CASE[circ_id].rnaser_fsj and CASE[circ_id].fsj != 0:
                case_exp = prior_exp_sampling(CASE_PRIOR, CASE[circ_id].rnaser_bsj, CASE[circ_id].rnaser_fsj, CASE[circ_id].fsj)
            else:
                case_exp = np.random.gamma(shape=int(CASE[circ_id].bsj) + 1 + 1, size=SIZE)
        else:
            case_exp = np.random.gamma(shape=1, size=SIZE)

        if circ_id in CTRL:
            if CTRL[circ_id].rnaser_fsj and CTRL[circ_id].fsj != 0:
                ctrl_exp = prior_exp_sampling(CTRL_PRIOR, CTRL[circ_id].rnaser_bsj, CTRL[circ_id].rnaser_fsj, CTRL[circ_id].fsj)
            else:
                ctrl_exp = np.random.gamma(shape=int(CTRL[circ_id].bsj) + 1 + 1, size=SIZE)
        else:
            ctrl_exp = np.random.gamma(shape=1, size=SIZE)

        tmp_de = corrected_score(case_exp, ctrl_exp, 1.0 / factor, size=SIZE, pvalue=PVAL)
        score_data[circ_id] = (tmp_de, None)
    return score_data


def load_gtf(in_file):
    from circ import GTFParser

    LOGGER.info('Loading CIRIquant result: {}'.format(in_file))

    circ_data = {}
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


def corrected_score(dis_a, dis_b, factor=1, size=SIZE, pvalue=0.05):
    from itertools import izip

    # fc = ne.evaluate("dis_a / dis_b")
    fc = sorted([np.log(j / i) for i, j in izip(dis_a, dis_b)])
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


def prior_exp_sampling(sample, bsj, fsj, total_fsj):
    factor = depth_factor(bsj / (bsj + 0.5 * fsj))
    corrected_read = ne.evaluate("sample * 0.5 * factor * total_fsj")
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
