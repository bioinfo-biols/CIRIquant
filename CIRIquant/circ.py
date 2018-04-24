#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
import re
import logging
import subprocess
from collections import defaultdict

from align import SamParser
logger = logging.getLogger('CIRIquant')


class DenovoParser(SamParser):
    @property
    def is_linear(self):
        if self.cigar == '*':
            return 0
        elif self.__is_linear(self.cigar_stats[0]) and self.__is_linear(self.cigar_stats[-1]):
            return 1
        else:
            return 0

    @staticmethod
    def __is_linear(segment):
        if segment[0] == 'M' or segment[1] <= 10:
            return 1
        else:
            return 0


class BedParser(object):
    def __init__(self, content):
        self.chr = content[0]
        self.start = int(content[1])
        self.end = int(content[2])
        self.id = '{}:{}|{}'.format(self.chr, self.start, self.end)


def load_bed(infile):
    # load circRNA position
    circ_info = defaultdict(dict)
    with open(infile, 'r') as f:
        for line in f:
            content = line.rstrip().split('\t')
            parser = BedParser(content)
            circ_id = '{}:{}|{}'.format(parser.chr, parser.start, parser.end)
            circ_info[parser.chr][circ_id] = parser
    return circ_info


def load_fai(infile):
    # load fasta index
    faidx = {}
    with open(infile, 'r') as f:
        for line in f:
            content = line.rstrip().split('\t')
            chrom, length, start = content[:3]
            faidx[chrom] = [int(start), int(length)]
    return faidx


def extract_seq(fasta, start, length):
    with open(fasta, 'r') as f:
        f.seek(start, 1)
        seq = f.read(length)
        seq = re.sub('\n', '', seq)
    return seq


def generate_index(log_file, circ_info, config, circ_fasta):
    # generate circular reference
    import subprocess
    from logger import ProgressBar

    fai = config['genome'] + '.fai'
    if not os.path.exists(fai):
        logger.debug('Indexing FASTA')
        index_cmd = '{} faidx {}'.format(config['samtools'], config['genome'])
        with open(log_file, 'a') as log:
            subprocess.call(index_cmd, shell=True, stderr=log, stdout=log)
    fasta_index = load_fai(fai)

    logger.info('Extract circular sequence')
    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    with open(circ_fasta, 'w') as out:
        for chrom in sorted(circ_info.keys()):
            prog.update(100 * cnt / len(circ_info))
            cnt += 1
            if chrom not in fasta_index:
                sys.exit('Unconsistent chromosome id: {}'.format(chrom))
            # logger.debug('Processing {}'.format(chrom))
            chrom_start, chrom_length = fasta_index[chrom]
            chrom_seq = extract_seq(config['genome'], chrom_start, chrom_length)

            chrom_circ = circ_info[chrom]
            for circ_id in chrom_circ:
                parser = chrom_circ[circ_id]
                circ_seq = chrom_seq[parser.start:parser.end + 1] * 2
                if circ_seq.count('N') > len(circ_seq) * 0.5 or len(circ_seq) == 0:
                    continue
                out.write('>{}'.format(parser.id) + '\n')
                out.write(circ_seq + '\n')
    prog.update(100)

    return fasta_index


def build_index(log_file, thread, pseudo_fasta, outdir, prefix, config):
    logger.info('Building circular index ..')
    denovo_index = '{}/circ/{}_index'.format(outdir, prefix)
    with open(log_file, 'a') as log:
        build_cmd = '{}-build -p {} -f {} {}'.format(
            config['hisat2'],
            thread,
            pseudo_fasta,
            denovo_index
        )
        subprocess.call(build_cmd, shell=True, stderr=log, stdout=log)
    return denovo_index


def denovo_alignment(log_file, thread, reads, outdir, prefix, config):
    logger.info('De novo alignment for circular RNAs ..')
    denovo_bam = '{}/circ/{}_denovo.bam'.format(outdir, prefix)
    with open(log_file, 'a') as log:
        align_cmd = '{} -p {} -q -x {}/circ/{}_index -1 {} -2 {} | {} view -bS > {}'.format(
            config['hisat2'],
            thread,
            outdir,
            prefix,
            reads[0],
            reads[1],
            config['samtools'],
            denovo_bam,
        )
        subprocess.call(align_cmd, shell=True, stderr=log, stdout=log)
    return denovo_bam


def mapped_bsj_reads(log_file, denovo_bam, config, threshold=5):
    circ_data = defaultdict(dict)
    stat = {'circular_reads': 0}
    logger.info('Parsing de-novo mapped BSJ reads')

    log = open(log_file, 'a')
    p = subprocess.Popen('{} view {}'.format(config['samtools'], denovo_bam),
                         shell=True, stderr=log, stdout=subprocess.PIPE)
    for line in p.stdout:
        content = line.rstrip().split('\t')
        parser = DenovoParser(content)
        if parser.cigar == '*':
            continue
        circ_chrom, circ_pos = parser.rname.split(':')
        circ_start, circ_end = int(circ_pos.split('|')[0]), int(circ_pos.split('|')[1])
        is_bsj = 0

        for i in parser.alignment:
            if i[0] + threshold <= circ_end - circ_start <= i[1] - threshold:
                is_bsj = 1
        if is_bsj and parser.is_linear:
            stat['circular_reads'] += 1
            circ_data[parser.prefix][parser.pair] = parser.rname
    p.communicate()
    log.close()

    return circ_data, stat


def mapped_fsj_reads(log_file, circ_info, bsj_reads, hisat_bam, config):
    logger.info('Generate circRNA index')
    circ_index = defaultdict(dict)
    for chrom in circ_info:
        for circ_id in circ_info[chrom]:
            parser = circ_info[chrom][circ_id]
            start_div, end_div = parser.start / 500, parser.end / 500
            for i in [start_div, end_div]:
                circ_index[chrom].setdefault(i, {}).update({circ_id: 1})

    logger.info('Parsing Forward-Spliced Junction reads')
    log = open(log_file, 'a')
    p = subprocess.Popen('{} view {}'.format(config['samtools'], hisat_bam),
                         shell=True, stderr=log, stdout=subprocess.PIPE)

    start_info = defaultdict(dict)
    end_info = defaultdict(dict)
    for line in p.stdout:
        content = line.rstrip().split('\t')
        parser = DenovoParser(content)
        if parser.prefix in bsj_reads or parser.cigar == '*':
            continue
        if parser.rname not in circ_index:
            continue
        if not parser.is_linear:
            continue

        for r_start, r_end, q_start, q_end in parser.alignment:
            if r_end - r_start < 5:
                continue
            start_div = r_start / 500
            end_div = r_end / 500
            for i in range(start_div, end_div + 1):
                if i not in circ_index[parser.rname]:
                    continue

                for circ_id in circ_index[parser.rname][i]:
                    if r_start <= circ_info[parser.rname][circ_id].start <= r_end:
                        start_info[circ_id].setdefault(parser.prefix, {}).update({parser.pair: 1})
                    if r_start <= circ_info[parser.rname][circ_id].end <= r_end:
                        end_info[circ_id].setdefault(parser.prefix, {}).update({parser.pair: 1})
    p.communicate()
    log.close()
    return start_info, end_info


def proc(log_file, thread, circ_file, hisat_bam, reads, outdir, prefix, config):
    from utils import check_dir
    circ_dir = '{}/circ'.format(outdir)
    check_dir(circ_dir)

    circ_fasta = '{}/circ/{}_index.fa'.format(outdir, prefix)
    circ_info = load_bed(circ_file)

    # extract fasta file for reads alignment
    generate_index(log_file, circ_info, config, circ_fasta)

    # hisat2-build index
    denovo_index = build_index(log_file, thread, circ_fasta, outdir, prefix, config)
    logger.debug('De-novo index: {}'.format(denovo_index))

    # hisat2 de novo alignment for candidate reads
    denovo_bam = denovo_alignment(log_file, thread, reads, outdir, prefix, config)
    logger.debug('De-novo bam: {}'.format(denovo_bam))

    bsj_info = defaultdict(list)
    coverage = {}
    bsj_reads, stat = mapped_bsj_reads(log_file, denovo_bam, config)
    for read_id in bsj_reads:
        for pair, circ_id in bsj_reads[read_id].iteritems():
            circ_id = bsj_reads[read_id][pair]
            coverage[circ_id] = coverage.setdefault(circ_id, 0) + 1
        bsj_info[circ_id].append(read_id)

    start_info, end_info = mapped_fsj_reads(log_file, circ_info, bsj_reads, hisat_bam, config)

    fsj_info = defaultdict(dict)
    for circ_id in start_info:
        for read_id in start_info[circ_id]:
            fsj_info[circ_id]['start'] = fsj_info[circ_id].setdefault('start', 0) + len(start_info[circ_id][read_id])
    for circ_id in end_info:
        for read_id in end_info[circ_id]:
            fsj_info[circ_id]['end'] = fsj_info[circ_id].setdefault('end', 0) + len(end_info[circ_id][read_id])

    output_file = '{}/{}.CIRIquant'.format(outdir, prefix)
    header = ['chr', 'start', 'end', 'name', 'coverage', 'bsj', 'start_fsj' 'end_fsj', 'reads_id']
    with open(output_file, 'w') as out:
        out.write('\t'.join(header) + '\n')
        for chrom in sorted(circ_info.keys()):
            for circ_id in sorted(circ_info[chrom].keys(), key=lambda x: circ_info[chrom][x].start):
                parser = circ_info[chrom][circ_id]
                tmp_out = [parser.chr, parser.start, parser.end, circ_id]
                if circ_id in coverage:
                    tmp_out.append(coverage[circ_id])
                else:
                    continue

                if circ_id in bsj_info:
                    tmp_out.append(len(bsj_info[circ_id]))
                else:
                    continue

                if circ_id in fsj_info and 'start' in fsj_info[circ_id]:
                    tmp_out.append(fsj_info[circ_id]['start'])
                else:
                    tmp_out.append(0)

                if circ_id in fsj_info and 'end' in fsj_info[circ_id]:
                    tmp_out.append(fsj_info[circ_id]['end'])
                else:
                    tmp_out.append(0)

                if circ_id in bsj_info:
                    tmp_out.append(','.join(bsj_info[circ_id]))

                if tmp_out[4] > 0:
                    out.write('\t'.join([str(x) for x in tmp_out]) + '\n')

    return output_file, stat