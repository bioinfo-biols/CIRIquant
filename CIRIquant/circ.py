#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
import re
import logging
import subprocess
from collections import defaultdict

import pdb

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
                circ_seq = chrom_seq[parser.start:parser.end + 1]
                if circ_seq.count('N') > len(circ_seq) * 0.5:
                    continue
                out.write('>{}'.format(parser.id) + '\n')
                out.write(circ_seq + '\n')
    prog.update(100)

    return fasta_index


def build_index(log_file, thread, pseudo_fasta, outdir, prefix, config):
    logger.info('Building circular index ..')
    with open(log_file, 'a') as log:
        build_cmd = '{}-build -p {} -f {} {}/circ/{}_index'.format(
            config['hisat2'],
            thread,
            pseudo_fasta,
            outdir,
            prefix,
        )
        subprocess.call(build_cmd, shell=True, stderr=log, stdout=log)
    return 1


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
    return 1


def proc(log_file, thread, circ_file, reads, outdir, prefix, config):
    import subprocess
    circ_fasta = '{}/circ/{}_index.fa'.format(outdir, prefix)
    circ_info = load_bed(circ_file)

    # extract fasta file for reads alignment
    # generate_index(log_file, circ_info, config, circ_fasta)

    # hisat2-build index
    # build_index(log_file, thread, circ_fasta, outdir, prefix, config)

    # hisat2 de novo alignment for candidate reads
    denovo_alignment(log_file, thread, reads, outdir, prefix, config)

