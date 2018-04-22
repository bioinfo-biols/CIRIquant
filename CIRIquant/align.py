#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import re
import logging
import subprocess
from string import maketrans
logger = logging.getLogger('CIRIquant')

trantab = maketrans('ATCGN', 'TAGCN')
CIGAR = re.compile(r'\d+[MIDNSHPX]')
PREFIX = re.compile(r'(.+)[/_-][12]')


class SamParser(object):
    """ Parse sam format """
    def __init__(self, content):
        self.qname = content[0]
        self.flag = int(content[1])
        self.rname = content[2]
        self.pos = int(content[3])
        self.mapq = int(content[4])
        self.cigar = content[5]

        prefix_m = PREFIX.search(self.qname)
        self.prefix = prefix_m.group(1) if prefix_m else self.qname
        self.seq = content[9]
        self.qual = content[10]

    @property
    def flag_stats(self):
        tmp = self.flag
        info = []
        for i in range(12):
            if tmp is 0:
                info.append(0)
            else:
                tmp, rem = divmod(tmp, 2)
                info.append(rem)
        return info

    @property
    def cigar_stats(self):
        if self.cigar == '*':
            return None

        tmp = CIGAR.findall(self.cigar)
        align = []
        for i in tmp:
            frag = [i[-1], int(i[:-1])]
            align.append(frag)
        return align

    @property
    def pair(self):
        return 1 - self.flag_stats[6]

    @property
    def strand(self):
        return 1 - self.flag_stats[4]

    @property
    def alignment(self):
        if self.cigar_stats is None:
            return None
        cov = []

        rstart = self.pos
        qstart = 0
        for match_type, length in self.cigar_stats:
            if match_type == 'M':
                cov.append([rstart, rstart + length - 1,
                            qstart, qstart + length - 1])
                qstart += length
                rstart += length
            elif match_type == 'S' or match_type == 'H':
                qstart += length
            elif match_type == 'D' or match_type == 'N':
                rstart += length
        return cov

    @property
    def is_linear(self):
        if self.cigar == '*':
            return 0
        elif self.__is_linear(self.cigar_stats[0]) and self.__is_linear(self.cigar_stats[-1]) and self.mapq >= 10:
            return 1
        else:
            return 0

    @staticmethod
    def __is_linear(segment):
        if segment[0] == 'M' and segment[1] > 5:
            return 1
        else:
            return 0

    def output(self):
        if self.strand:
            return self.seq, self.qual
        else:
            return reverse_comp(self.seq), reverse(self.qual)


def reverse(seq):
    return seq[::-1]


def reverse_comp(seq):
    return ''.join(reverse(seq.translate(trantab)))


def write_fastq(rname, seq, qual, fh):
    fh.write('@{}\n'.format(rname))
    fh.write('{}\n'.format(seq))
    fh.write('+\n')
    fh.write('{}\n'.format(qual))


def unmapped_reads(log_file, thread, outdir, prefix, hisat_bam, config):
    from utils import check_dir, subprocess_setup
    logger.info('Extracting circular candidate reads ..')
    unmapped_dir = outdir + '/circ'
    check_dir(unmapped_dir)

    # open BAM file pipe
    log = open(log_file, 'a')
    fpipe = subprocess.Popen('{} view -h {}'.format(config['samtools'], hisat_bam), shell=True,
                             stdout=subprocess.PIPE, stderr=log,
                             preexec_fn=subprocess_setup)
    infile = fpipe.stdout

    cand_reads = ['{}/{}_unmapped_1.fq'.format(unmapped_dir, prefix),
                  '{}/{}_unmapped_2.fq'.format(unmapped_dir, prefix), ]
    with open(cand_reads[0], 'w') as r1, open(cand_reads[1], 'w') as r2:
        stat = proc_unmapped(infile, r1, r2)

    fpipe.communicate()
    log.close()
    logger.debug('circRNA reads: {0}_unmapped_1.fq, {0}_unmapped_2.fq'.format(prefix))
    return cand_reads, stat


def proc_unmapped(infile, r1, r2):
    read_pair = {}
    tmp_read = None
    stat = {'cleaned_reads': 0, 'mapped_reads': 0}

    for line in infile:
        if line.startswith('@'):
            continue
        content = line.rstrip().split('\t')
        parser = SamParser(content)
        if parser.flag_stats[8]:
            continue

        if tmp_read and parser.prefix != tmp_read:
            stat['cleaned_reads'] += 1
            if not read_pair[0].flag_stats[2] and not read_pair[1].flag_stats[2]:
                stat['mapped_reads'] += 1

            if read_pair[0].is_linear and read_pair[1].is_linear and read_pair[0].strand != read_pair[1].strand:
                pass
            else:
                read1, read2 = read_pair[0], read_pair[1]
                read1_seq, read1_qual = read1.output()
                write_fastq(read1.prefix, read1_seq, read1_qual, r1)
                read2_seq, read2_qual = read2.output()
                write_fastq(read2.prefix, read2_seq, read2_qual, r2)
            read_pair = {}
        tmp_read = parser.prefix
        read_pair[parser.pair] = parser

    # Recover last read pair
    stat['cleaned_reads'] += 1
    if not read_pair[0].flag_stats[2] and not read_pair[1].flag_stats[2]:
        stat['mapped_reads'] += 1
    if read_pair[0].is_linear and read_pair[1].is_linear and read_pair[0].strand != read_pair[1].strand:
        read1, read2 = read_pair[0], read_pair[1]
        read1_seq, read1_qual = read1.output()
        write_fastq(read1.prefix, read1_seq, read1_qual, r1)
        read2_seq, read2_qual = read2.output()
        write_fastq(read2.prefix, read2_seq, read2_qual, r2)
    return stat
