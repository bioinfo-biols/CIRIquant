#!/usr/bin/env python
import os
import sys
import re
import time
import pysam
import gzip
import subprocess
# import multiprocessing
from collections import defaultdict
from itertools import izip_longest
from string import maketrans

# LOCK = multiprocessing.Lock()
PREFIX=re.compile(r'(.+)[/_-][12]')
trantab = maketrans('ATCGN', 'TAGCN')


# def paral_write(info):
#     LOCK.acquire()
#     print info


def flag_stats(flag):
    info = []
    for i in range(12):
        if flag is 0:
            info.append(0)
        else:
            flag, rem = divmod(flag, 2)
            info.append(rem)
    return info


def linear_segment(segment):
    if segment[0] == 0 and segment[1] > 5:
        return 1
    else:
        return 0


def is_linear(query):
    if query.cigar and linear_segment(query.cigartuples[0]) and linear_segment(query.cigartuples[-1]):
        return 1
    else:
        return 0


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=None)


def sam_to_fastq(query):
    seq = query.query_sequence
    qual = ''.join([chr(i + 33) for i in query.query_qualities])
    if query.is_reverse:
        seq, qual = reverse_comp(seq), reverse(qual)
    return seq, qual


def reverse(seq):
    return seq[::-1]


def reverse_comp(seq):
    return ''.join(reverse(seq.translate(trantab)))


def unmapped_reads(bam, r1, r2):
    stat = defaultdict(int)

    args = ['samtools', 'view', '-u', '-F', '256', bam]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE)
    fd_child = proc.stdout.fileno()

    sam = pysam.AlignmentFile(fd_child, 'rb')
    for pair_num, (read1, read2) in enumerate(grouper(sam, 2)):
        read_m = PREFIX.search(read1.query_name)
        read_prefix = read_m.group(1) if read_m else read_m.query_name

        if read1.is_read2:
            read1, read2 = read2, read1

        # print read1, read1.is_reverse
        # print read2, read2.is_reverse
        # raw_input()

        stat['cleaned_reads'] += 1
        if not read1.is_unmapped and not read2.is_unmapped:
            stat['mapped_reads'] += 1

        if is_linear(read1) and is_linear(read2) and read1.is_reverse != read2.is_reverse:
            pass
        else:
            r1.write('@{}\n{}\n+\n{}\n'.format(read_prefix, *sam_to_fastq(read1)))
            r2.write('@{}\n{}\n+\n{}\n'.format(read_prefix, *sam_to_fastq(read2)))

    sam.close()
    proc.communicate()
    return stat


def main():
    start = time.time()
    chroms = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    cpu = 4

    bam_file = '/histor/zhao/zhangjy/dev/CIRIquant2/test_data/test/align/test.bam'
    # bam_file = '/histor/zhao/zhangjy/01.circExp/FDR_Compare/hs68_rep1/align/hs68_rep1.bam'
    read1 = '/histor/zhao/zhangjy/dev/CIRIquant2/test_data/test_pysam/read1.fq.gz'
    read2 = '/histor/zhao/zhangjy/dev/CIRIquant2/test_data/test_pysam/read2.fq.gz'

    with gzip.open(read1, 'wb') as r1, gzip.open(read2, 'wb') as r2:
        stat = unmapped_reads(bam_file, r1, r2)
    print stat
    # pool = multiprocessing.Pool(processes=cpu)
    # for x in chroms:
    #     pool.apply_async(unmapped_reads, (x, bam_file))
    # pool.close()
    # pool.join()

    end = time.time()
    print (end - start)


if __name__ == '__main__':
    main()
