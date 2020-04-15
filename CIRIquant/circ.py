#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
import re
import logging
import pysam
import time
import subprocess

from multiprocessing import Pool
from collections import defaultdict
from itertools import izip_longest

import utils

LOGGER = logging.getLogger('CIRIquant')
PREFIX = re.compile(r'(.+)[/_-][12]')


class BedParser(object):
    """
    Class for parsing circRNA information in bed file
    """

    def __init__(self, content):
        self.chrom = content[0]
        self.start = int(content[1])
        self.end = int(content[2])
        self.circ_id = content[3]
        self.strand = content[5]
        self.length = self.end - self.start + 1


def load_bed(fname):
    """
    Load Back-Spliced Junction Sites in bed file

    Parameters
    -----
    fname : str
        input file name

    Returns
    -----
    dict
        divide all circRNAs into different chromsomes

    """
    circ_info = defaultdict(dict)
    with open(fname, 'r') as f:
        for line in f:
            content = line.rstrip().split('\t')
            parser = BedParser(content)
            circ_info[parser.chrom][parser.circ_id] = parser
    return circ_info


def update_info(circ_info, rnaser_file):
    """
    Add information of RNase R circRNAs
    """
    circ_exp = {}
    LOGGER.info('RNase R file: {}'.format(rnaser_file))
    with open(rnaser_file, 'r') as f:
        header = {}
        for line in f:
            if line.startswith('##'):
                key, value = line.rstrip().strip('#').split(':')
                header.update({key: value})
                continue

            content = line.rstrip().split('\t')
            tmp_parser = GTFParser(content)
            circ_exp[tmp_parser.attr['circ_id']] = {
                'bsj': float(tmp_parser.attr['bsj']),
                'fsj': float(tmp_parser.attr['fsj']),
                'ratio': 2 * float(tmp_parser.attr['bsj']) / (2 * float(tmp_parser.attr['bsj']) + float(tmp_parser.attr['fsj']))
            }
            parser = BedParser([
                tmp_parser.chrom,
                tmp_parser.start,
                tmp_parser.end,
                tmp_parser.attr['circ_id'],
                '.',
                tmp_parser.strand,
            ])
            if parser.chrom in circ_info and parser.circ_id in circ_info[parser.chrom]:
                continue
            circ_info[parser.chrom].update({parser.circ_id: parser})

    return circ_exp, (int(header['Total_Reads']), int(header['Mapped_Reads']), int(header['Circular_Reads']))


def load_fai(fname):
    """
    Load fai index of fasta created by samtools

    Parameters
    -----
    fname : str
        input fai file name

    Returns
    -----
    dict
         chromsome and start / end position in file

    """
    faidx = {}
    with open(fname, 'r') as f:
        for line in f:
            content = line.rstrip().split('\t')
            chrom, length, start, eff_length, line_length = content
            shift_length = int(length) * int(line_length) / int(eff_length)
            faidx[chrom] = [int(start), shift_length]
    return faidx


def extract_seq(fasta, start, length):
    """
    Extract sequence from fasta according to given position

    Parameters
    -----
    fasta : str
        file name of fasta
    start : int
        offset of chrosome
    length : int
        length of chromsome sequence

    Returns
    -----
    str
        sequence from start to start + length

    """
    with open(fasta, 'r') as f:
        f.seek(start, 0)
        seq = f.read(length)
        seq = re.sub('\n', '', seq)
    return seq


def generate_index(log_file, circ_info, circ_fasta):
    """
    Generate pseudo circular index

    Parameters
    -----
    log_file : str
        file name of log file, used for subprocess.PIPE
    circ_info : dict
        back-spliced junction sites
    circ_fasta :
        output fasta file name

    Returns
    -----
    dict
        chromosomes used in BSJ sites

    """

    from logger import ProgressBar

    fai = utils.FASTA + '.fai'
    if not os.path.exists(fai):
        LOGGER.debug('Indexing FASTA')
        index_cmd = '{} faidx {}'.format(utils.SAMTOOLS, utils.FASTA)
        with open(log_file, 'a') as log:
            LOGGER.debug(index_cmd)
            subprocess.call(index_cmd, shell=True, stderr=log, stdout=log)
    fasta_index = load_fai(fai)

    LOGGER.info('Extract circular sequence')
    prog = ProgressBar()
    prog.update(0)
    cnt = 0
    with open(circ_fasta, 'w') as out:
        for chrom in sorted(circ_info.keys()):
            prog.update(100 * cnt / len(circ_info))
            cnt += 1
            if chrom not in fasta_index:
                sys.exit('Unconsistent chromosome id: {}'.format(chrom))

            chrom_start, chrom_length = fasta_index[chrom]
            chrom_seq = extract_seq(utils.FASTA, chrom_start, chrom_length)

            chrom_circ = circ_info[chrom]
            for circ_id in chrom_circ:
                parser = chrom_circ[circ_id]
                circ_seq = chrom_seq[parser.start - 1:parser.end] * 2
                if circ_seq.count('N') > len(circ_seq) * 0.5 or len(circ_seq) == 0:
                    continue
                out.write('>{}'.format(parser.circ_id) + '\n')
                out.write(circ_seq + '\n')
    prog.update(100)

    return fasta_index


def build_index(log_file, thread, pseudo_fasta, outdir, prefix):
    """
    Build hisat2 index for pseudo circular index

    Returns
    -----
    str
        index file name used in denovo mapping

    """
    LOGGER.info('Building circular index ..')
    denovo_index = '{}/circ/{}_index'.format(outdir, prefix)

    build_cmd = '{}-build -p {} -f {} {}'.format(
        utils.HISAT2,
        thread,
        pseudo_fasta,
        denovo_index
    )

    with open(log_file, 'a') as log:
        LOGGER.debug(build_cmd)
        subprocess.call(build_cmd, shell=True, stderr=log, stdout=log)

    return denovo_index


def denovo_alignment(log_file, thread, reads, outdir, prefix):
    """
    Call hisat2 to read re-alignment

    Returns
    -----
    str
        Output bam file

    """
    LOGGER.info('De novo alignment for circular RNAs ..')
    denovo_bam = '{}/circ/{}_denovo.bam'.format(outdir, prefix)
    sorted_bam = '{}/circ/{}_denovo.sorted.bam'.format(outdir, prefix)

    align_cmd = '{} -p {} --dta -q -x {}/circ/{}_index -1 {} -2 {} | {} view -bS > {}'.format(
        utils.HISAT2,
        thread,
        outdir,
        prefix,
        reads[0],
        reads[1],
        utils.SAMTOOLS,
        denovo_bam,
    )

    sort_cmd = '{} sort --threads {} -o {} {}'.format(
        utils.SAMTOOLS,
        thread,
        sorted_bam,
        denovo_bam,
    )

    index_cmd = '{} index -@ {} {}'.format(
        utils.SAMTOOLS,
        thread,
        sorted_bam,
    )

    with open(log_file, 'a') as log:
        LOGGER.debug(align_cmd)
        subprocess.call(align_cmd, shell=True, stderr=log, stdout=log)

        LOGGER.debug(sort_cmd)
        subprocess.call(sort_cmd, shell=True, stderr=log, stdout=log)

        LOGGER.debug(index_cmd)
        subprocess.call(index_cmd, shell=True, stderr=log, stdout=log)
    # wait for samtools index to correct time stamp
    time.sleep(5)
    return sorted_bam


def grouper(iterable, n, fillvalue=None):
    """
    Collect data info fixed-length chunks or blocks
    grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    """

    args = [iter(iterable)] * n
    return izip_longest(*args, fillvalue=None)


def proc_denovo_bam(bam_file, thread, circ_info, threshold, lib_type):
    """
    Extract BSJ reads in denovo mapped bam file

    Returns
    -----
    dict
        query_name -> mate_id -> pysam.AlignSegment

    """

    LOGGER.info('Detecting reads containing Back-splicing signals')

    sam = pysam.AlignmentFile(bam_file, 'rb')
    header = sam.header['SQ']
    sam.close()

    pool = Pool(thread, denovo_initializer, (bam_file, circ_info, threshold, ))
    jobs = []
    chunk_size = max(500, len(header) / thread + 1)
    for circ_chunk in grouper(header, chunk_size):
        jobs.append(pool.apply_async(denovo_worker, (circ_chunk, lib_type, )))
    pool.close()
    pool.join()

    cand_reads = defaultdict(dict)
    for job in jobs:
        tmp_cand = job.get()
        for read_id, mate_id, circ_id, blocks, cigartuples in tmp_cand:
            cand_reads[read_id][mate_id] = (circ_id, blocks, cigartuples)

    return cand_reads


BAM = None
THRESHOLD = None
CIRC = None


def denovo_initializer(infile, circ_info, threshold):
    """
    Initializer for passing bam file name
    """
    global BAM, CIRC, THRESHOLD
    BAM, CIRC, THRESHOLD = infile, circ_info, threshold


def denovo_worker(circ_chunk, lib_type):
    """
    Find candidate reads with junction signal

    Parameters
    -----
    circ_chunk : list
        list of Pysam header to process

    Returns
    -----
    list
        pysam.AlignedSegment, candidate reads with junction signal

    """
    #TODO: add support for stranded specific RNA_seq libraries
    sam = pysam.AlignmentFile(BAM, 'rb')
    cand_reads = []
    for d in circ_chunk:
        if d is None:
            continue
        circ_id, junc_site = d['SN'], int(d['LN']) / 2
        contig = circ_id.split(':')[0]
        parser = CIRC[contig][circ_id]

        tmp_cand = []
        circ_reads = defaultdict(dict)
        for read in sam.fetch(circ_id, multiple_iterators=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            read_strand = '-' if read.is_reverse else '+'
            if lib_type == 0:
                pass
            elif lib_type == 1:
                if read.is_read1 - read.is_read2 > 0:
                    if read_strand != parser.strand:
                        continue
                else:
                    if read_strand == parser.strand:
                        continue
            elif lib_type == 2:
                if read.is_read1 - read.is_read2 > 0:
                    if read_strand == parser.strand:
                        continue
                else:
                    if read_strand != parser.strand:
                        continue

            circ_reads[query_prefix(read.query_name)][read.is_read1 - read.is_read2] = 1
            if read.mapping_quality <= 10:
                continue
            if read.get_overlap(junc_site - THRESHOLD, junc_site + THRESHOLD) >= THRESHOLD * 2:
                tmp_cand.append((read.query_name, read.is_read1 - read.is_read2, circ_id, read.get_blocks(), read.cigartuples))

        for qname, mate_id, circ_id, blocks, cigartuples in tmp_cand:
            if -1 * mate_id in circ_reads[query_prefix(qname)]:
                cand_reads.append((qname, mate_id, circ_id, blocks, cigartuples))
    sam.close()

    return cand_reads


def proc_genome_bam(bam_file, thread, circ_info, cand_reads, threshold, tmp_dir):
    """
    Extract FSJ reads and check BSJ reads alignment information

    Returns
    -----
    dict
        bsj reads of circRNAs, pair_id -> mate_id -> circ_id
    dict
        fsj reads of circRNAs, pair_id -> mate_id -> circ_id

    """
    import cPickle
    LOGGER.info('Detecting FSJ reads from genome alignment file')

    sam = pysam.AlignmentFile(bam_file, 'rb')
    header = sam.header['SQ']
    sam.close()

    pool = Pool(thread, genome_initializer, (bam_file, circ_info, cand_reads, threshold))
    jobs = []
    for chrom_info in header:
        jobs.append(pool.apply_async(genome_worker, (chrom_info['SN'], tmp_dir, )))
    pool.close()
    pool.join()

    fp_bsj = defaultdict(dict)
    fsj_reads = defaultdict(dict)
    cand_to_genome = []

    for job in jobs:
        tmp = job.get()
        if tmp is None:
            continue
        res = tmp if isinstance(tmp, dict) else cPickle.load(open(tmp, 'rb'))
        chrom_fp_bsj, chrom_fsj, chrom_cand = res['fp_bsj'], res['fsj_reads'], res['cand_to_genome']
        for pair_id, mate_id in chrom_fp_bsj:
            fp_bsj[pair_id][mate_id] = 1
        for pair_id, mate_id, circ_id in chrom_fsj:
            fsj_reads[pair_id][mate_id] = circ_id
        cand_to_genome += chrom_cand

    circ_bsj = defaultdict(dict)
    circ_fsj = defaultdict(dict)
    for pair_id in cand_reads:
        for mate_id, (circ_id, blocks, cigartuples) in cand_reads[pair_id].iteritems():
            if pair_id in fp_bsj and mate_id in fp_bsj[pair_id]:
                continue
            circ_bsj[circ_id].update({query_prefix(pair_id): 1})

    for pair_id in fsj_reads:
        for mate_id, circ_id in fsj_reads[pair_id].iteritems():
            if pair_id in cand_reads and mate_id in cand_reads[pair_id] and not (pair_id in fp_bsj and mate_id in fp_bsj[pair_id]):
                continue
            circ_fsj[circ_id].update({query_prefix(pair_id): 1})

    # sam = pysam.AlignmentFile(bam_file, 'rb')
    # cand_sam = pysam.AlignmentFile(bam_file + '.fsj', 'w', template=sam)
    # cand_sam.close()
    # sam.close()
    #
    # with open(bam_file + '.fsj', 'w') as out:
    #     for read in cand_to_genome:
    #         out.write(read + '\n')
    return circ_bsj, circ_fsj


BSJ = None


def genome_initializer(bam_file, circ_info, cand_bsj, threshold):
    """
    Initializer for passing bam file name and circRNA_info
    """
    global BAM, CIRC, THRESHOLD, BSJ
    BAM, CIRC, BSJ, THRESHOLD = bam_file, circ_info, cand_bsj, threshold


def genome_worker(chrom, tmp_dir):
    """
    Find FSJ reads and re-check BSJ reads

    Parameters
    -----
    chrom : str
        chromosme or scaffold name for process

    Returns
    -----
    list
        false positive reads information,  (query_name, mate_id)
    list
        fsj_reads of circRNAs, (query_name, mate_id, circ_id)

    """
    import cPickle

    if chrom not in CIRC:
        return None

    sam = pysam.AlignmentFile(BAM, 'rb')
    cand_to_genome = []

    fp_bsj = []
    for read in sam.fetch(chrom, multiple_iterators=True):
        # If Reads is bsj candidate
        if read.is_unmapped:
            continue
        if read.query_name not in BSJ or read.is_read1 - read.is_read2 not in BSJ[read.query_name]:
            continue
        circ_id, blocks, cigartuples = BSJ[read.query_name][read.is_read1 - read.is_read2]
        cand_to_genome.append(read.to_string())

        # check alignment against reference genome
        qual_filter = 1 if mapping_quality(blocks) <= mapping_quality(read.get_blocks()) + 5 else 0
        linear_filter = 1 if is_linear(read.cigartuples[0]) and is_linear(read.cigartuples[-1]) else 0
        align_filter = 1 if read.mapping_quality <= 10 or read.is_secondary or read.is_supplementary else 0
        if qual_filter or linear_filter or align_filter:
            fp_bsj.append((read.query_name, read.is_read1 - read.is_read2))

    fsj_reads = []
    for circ_id, parser in CIRC[chrom].iteritems():
        # FSJ across start site
        for read in sam.fetch(region='{0}:{1}-{1}'.format(chrom, parser.start)):
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.mapping_quality <= 10:
                continue
            if not read.get_overlap(parser.start - 1, parser.start + THRESHOLD - 1) >= THRESHOLD:
                continue
            if is_mapped(read.cigartuples[0]) and is_mapped(read.cigartuples[-1]):
                fsj_reads.append((read.query_name, read.is_read1 - read.is_read2, circ_id))

        for read in sam.fetch(region='{0}:{1}-{1}'.format(chrom, parser.end)):
            if read.is_unmapped or read.is_supplementary:
                continue
            if read.mapping_quality <= 10:
                continue
            if not read.get_overlap(parser.end - THRESHOLD, parser.end) >= THRESHOLD:
                continue
            if is_mapped(read.cigartuples[0]) and is_mapped(read.cigartuples[-1]):
                fsj_reads.append((read.query_name, read.is_read1 - read.is_read2, circ_id))

    sam.close()

    res = {'fp_bsj': fp_bsj, 'fsj_reads': fsj_reads, 'cand_to_genome': cand_to_genome}

    res_to_string = cPickle.dumps(res, 0)
    if sys.getsizeof(res_to_string) > 1024 * 1024 * 1024:
        pkl_file = "{}/{}.pkl".format(tmp_dir, chrom)
        cPickle.dump(res, open(pkl_file, "wb"), -1)
        return pkl_file

    return res


def mapping_quality(blocks):
    return sum([j - i for i, j in blocks])


def is_mapped(cigar_tuple):
    """
    Whether end of alignment segment is a mapped end

    Parameters
    -----
    cigar_tuple : tuple of cigar

    Returns
    -----
    int
        1 for linear end, 0 for ambiguous end

    """
    if cigar_tuple[0] == 0 or cigar_tuple[1] <= 10:
        return 1
    else:
        return 0


def is_linear(cigar_tuple):
    """
    Whether end of alignment segment is a linear end

    Parameters
    -----
    cigar_tuple : tuple of cigar

    Returns
    -----
    int
        1 for linear end, 0 for ambiguous end

    """
    if cigar_tuple[0] == 0 and cigar_tuple[1] >= 5:
        return 1
    else:
        return 0


def query_prefix(query_name):
    """
    Get pair id without mate id marker

    Paramters
    -----
    read : pysam.AlignedSegment

    Returns
    -----
    str
        mate id of segment

    """
    prefix_m = PREFIX.search(query_name)
    prefix = prefix_m.group(1) if prefix_m else query_name
    return prefix


def proc(log_file, thread, circ_file, hisat_bam, rnaser_file, reads, outdir, prefix, anchor, lib_type):
    """
    Build pseudo circular reference index and perform reads re-alignment
    Extract BSJ and FSJ reads from alignment results

    Returns
    -----
    str
        output file name

    """
    from utils import check_dir
    circ_dir = '{}/circ'.format(outdir)
    check_dir(circ_dir)

    circ_fasta = '{}/circ/{}_index.fa'.format(outdir, prefix)
    circ_info = load_bed(circ_file)
    if rnaser_file:
        LOGGER.info('Loading RNase R results')
        rnaser_exp, rnaser_stat = update_info(circ_info, rnaser_file)

    # extract fasta file for reads alignment
    generate_index(log_file, circ_info, circ_fasta)

    # hisat2-build index
    denovo_index = build_index(log_file, thread, circ_fasta, outdir, prefix)
    LOGGER.debug('De-novo index: {}'.format(denovo_index))

    # hisat2 de novo alignment for candidate reads
    denovo_bam = denovo_alignment(log_file, thread, reads, outdir, prefix)
    LOGGER.debug('De-novo bam: {}'.format(denovo_bam))

    # Find BSJ and FSJ informations
    cand_bsj = proc_denovo_bam(denovo_bam, thread, circ_info, anchor, lib_type)
    bsj_reads, fsj_reads = proc_genome_bam(hisat_bam, thread, circ_info, cand_bsj, anchor, circ_dir)

    total_reads, mapped_reads = bam_stat(hisat_bam)
    circ_reads = sum([len(bsj_reads[i]) for i in bsj_reads]) * 2
    sample_stat = (total_reads, mapped_reads, circ_reads)

    sample_exp = expression_level(circ_info, bsj_reads, fsj_reads)

    # circRNA annotation
    header = [
        'Sample: {}'.format(prefix),
        'Total_Reads: {}'.format(total_reads),
        'Mapped_Reads: {}'.format(mapped_reads),
        'Circular_Reads: {}'.format(circ_reads),
    ]
    out_file = '{}/{}.gtf'.format(outdir, prefix)

    if rnaser_file:
        import coeff
        tmp_header, circ_exp = coeff.correction(sample_exp, sample_stat, rnaser_exp, rnaser_stat)
        header += tmp_header
    else:
        circ_exp = sample_exp

    from version import __version__
    header += ['version: {}'.format(__version__), ]
    gtf_info = index_annotation(utils.GTF)
    format_output(circ_info, circ_exp, sample_stat, header, gtf_info, out_file)

    return out_file


def expression_level(circ_info, bsj_reads, fsj_reads):
    LOGGER.info('Merge bsj and fsj results')
    circ_exp = defaultdict(dict)
    bsj_ids = {}
    for chrom in circ_info:
        for circ_id in circ_info[chrom]:
            bsj = len(bsj_reads[circ_id]) if circ_id in bsj_reads else 0
            fsj = len(fsj_reads[circ_id]) if circ_id in fsj_reads else 0
            if bsj == 0 and fsj == 0:
                continue
            bsj_ids[circ_id] = bsj_reads[circ_id].keys()
            junc = 2.0 * bsj / (2.0 * bsj + fsj)
            circ_exp[circ_id] = {'bsj': bsj, 'fsj': fsj, 'ratio': junc}

    return circ_exp


def bam_stat(bam_file):
    """
    Stat of bam file

    Returns
    -----
    int
        number of total reads
    int
        number of mapped reads

    """
    sam = pysam.AlignmentFile(bam_file, 'rb')
    total = sam.count(read_callback=total_callback, until_eof=True)
    sam.close()

    sam = pysam.AlignmentFile(bam_file, 'rb')
    unmapped = sam.count(read_callback=unmapped_callback, until_eof=True)
    sam.close()

    return total, total - unmapped


def unmapped_callback(read):
    """
    callback for counting unmapped reads
    """
    return read.is_unmapped or read.mate_is_unmapped and not read.is_supplementary and not read.is_secondary


def total_callback(read):
    """
    callback for counting total reads
    """
    return not read.is_supplementary and not read.is_secondary


def format_output(circ_info, circ_exp, sample_stat, header, gtf_index, outfile):
    """
    Output bsj information of circRNA expression levels
    """
    LOGGER.info('Output circRNA expression values')

    with open(outfile, 'w') as out:
        for h in header:
            out.write('##' + h + '\n')
        for chrom in sorted(circ_info.keys(), key=by_chrom):
            for circ_id in sorted(circ_info[chrom].keys(), cmp=by_circ, key=lambda x:circ_info[chrom][x]):
                if circ_id not in circ_exp or circ_exp[circ_id]['bsj'] == 0:
                    continue
                parser = circ_info[chrom][circ_id]
                strand = parser.strand
                tmp_line = [
                    chrom,
                    'CIRIquant',
                    'circRNA',
                    parser.start,
                    parser.end,
                    '{:.4f}'.format(2 * 1000.0 * 1000.0 * circ_exp[circ_id]['bsj'] / sample_stat[1]),
                    strand,
                    '.',
                ]

                field = circRNA_attr(gtf_index, parser)
                tmp_attr = 'circ_id "{}"; circ_type "{}"; bsj {:.3f}; fsj {:.3f}; junc_ratio {:.3f};'.format(
                    circ_id,
                    field['circ_type'] if field else "Unknown",
                    circ_exp[circ_id]['bsj'],
                    circ_exp[circ_id]['fsj'],
                    circ_exp[circ_id]['ratio'],
                )
                for key in 'rnaser_bsj', 'rnaser_fsj':
                    if key in circ_exp[circ_id]:
                        tmp_attr += ' {} {:.3f};'.format(key, circ_exp[circ_id][key])

                for key in 'gene_id', 'gene_name', 'gene_type':
                    if key in field:
                        tmp_attr += ' {} "{}";'.format(key, field[key])

                tmp_line.append(tmp_attr)

                out.write('\t'.join([str(x) for x in tmp_line]) + '\n')
    return 1


def by_chrom(x):
    """
    Sort by chromosomes
    """
    chrom = x
    if x.startswith('chr'):
        chrom = chrom.strip('chr')
    try:
        chrom = int(chrom)
    except ValueError as e:
        pass
    return chrom


def by_circ(x, y):
    """
    Sort circRNAs by the start and end position
    """
    return x.end - y.end if x.start == y.start else x.start - y.start


class GTFParser(object):
    """
    Class for parsing annotation gtf
    """

    def __init__(self, content):
        self.chrom = content[0]
        self.source = content[1]
        self.type = content[2]
        self.start, self.end = int(content[3]), int(content[4])
        self.strand = content[6]
        self.attr_string = content[8]


    @property
    def attr(self):
        """
        Parsing attribute column in gtf file
        """
        field = {}
        for attr_values in [re.split(r'[\s=]+', i.strip()) for i in self.attr_string.split(';')[:-1]]:
            key, value = attr_values[0], attr_values[1:]
            field[key] = ' '.join(value).strip('"')
        return field


def index_annotation(gtf):
    """
    Generate binned index for element in gtf
    """

    LOGGER.info('Loading annotation gtf ..')
    gtf_index = defaultdict(dict)
    with open(gtf, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            content = line.rstrip().split('\t')
            # only include gene and exon feature for now
            if content[2] not in ['gene', 'exon']:
                continue
            parser = GTFParser(content)
            # if 'gene_biotype' in parser.attr and parser.attr['gene_biotype'] in ['lincRNA', 'pseudogene']:
            #     continue
            # if 'gene_type' in parser.attr and parser.attr['gene_type'] in ['lincRNA', 'pseudogene']:
            #     continue

            start_div, end_div = parser.start / 500, parser.end / 500
            for i in xrange(start_div, end_div + 1):
                gtf_index[parser.chrom].setdefault(i, []).append(parser)
    return gtf_index


def circRNA_attr(gtf_index, circ):
    """
    annotate circRNA information
    """
    if circ.chrom not in gtf_index:
        LOGGER.warn('chrom of contig "{}" not in annotation gtf, please check'.format(circ.chrom))
        return {}
    start_div, end_div = circ.start / 500, circ.end / 500

    host_gene = {}
    start_element = defaultdict(list)
    end_element = defaultdict(list)

    for x in xrange(start_div, end_div + 1):
        if x not in gtf_index[circ.chrom]:
            continue
        for element in gtf_index[circ.chrom][x]:
            # start site
            if element.start <= circ.start <= element.end and element.strand == circ.strand:
                start_element[element.type].append(element)
            # end site
            if element.start <= circ.end <= element.end and element.strand == circ.strand:
                end_element[element.type].append(element)
            # if element.type != 'gene':
            #     continue
            if element.end < circ.start or circ.end < element.start:
                continue
            if element.attr['gene_id'] not in host_gene:
                host_gene[element.attr['gene_id']] = element

    circ_type = {}
    forward_host_gene = []
    antisense_host_gene = []

    if len(host_gene) > 0:
        for gene_id in host_gene:
            if host_gene[gene_id].strand == circ.strand:
                forward_host_gene.append(host_gene[gene_id])
                if 'exon' in start_element and 'exon' in end_element:
                    circ_type['exon'] = 1
                else:
                    circ_type['intron'] = 1
            else:
                antisense_host_gene.append(host_gene[gene_id])
                circ_type['antisense'] = 1
    else:
        circ_type['intergenic'] = 1

    if len(forward_host_gene) > 1:
        circ_type['gene_intergenic'] = 1

    field = {}
    if 'exon' in circ_type:
        field['circ_type'] = 'exon'
    elif 'intron' in circ_type:
        field['circ_type'] = 'intron'
    # elif 'gene_intergenic' in circ_type:
    #     field['circ_type'] = 'gene_intergenic'
    elif 'antisense' in circ_type:
        field['circ_type'] = 'antisense'
    else:
        field['circ_type'] = 'intergenic'

    if len(forward_host_gene) == 1:
        if 'gene_id' in forward_host_gene[0].attr:
            field['gene_id'] = forward_host_gene[0].attr['gene_id']
        if 'gene_name' in forward_host_gene[0].attr:
            field['gene_name'] = forward_host_gene[0].attr['gene_name']
        if 'gene_type' in forward_host_gene[0].attr:
            field['gene_type'] = forward_host_gene[0].attr['gene_type']
        elif 'gene_biotype' in forward_host_gene[0].attr:
            field['gene_type'] = forward_host_gene[0].attr['gene_biotype']
        else:
            pass
    elif len(forward_host_gene) > 1:
        tmp_gene_id = []
        tmp_gene_name = []
        tmp_gene_type = []
        for x in forward_host_gene:
            if 'gene_id' in x.attr:
                tmp_gene_id.append(x.attr['gene_id'])
            if 'gene_name' in x.attr:
                tmp_gene_name.append(x.attr['gene_name'])
            if 'gene_type' in x.attr:
                tmp_gene_type.append(x.attr['gene_type'])
            elif 'gene_biotype' in x.attr:
                tmp_gene_type.append(x.attr['gene_biotype'])
            else:
                pass
        if len(tmp_gene_id) > 0:
            field['gene_id'] = ','.join(tmp_gene_id)
        if len(tmp_gene_name) > 0:
            field['gene_name'] = ','.join(tmp_gene_name)
        if len(tmp_gene_type) > 0:
            field['gene_type'] = ','.join(tmp_gene_type)
    elif field['circ_type'] == 'antisense' and len(antisense_host_gene) > 0:
        tmp_gene_id = []
        tmp_gene_name = []
        tmp_gene_type = []
        for x in antisense_host_gene:
            if 'gene_id' in x.attr:
                tmp_gene_id.append(x.attr['gene_id'])
            if 'gene_name' in x.attr:
                tmp_gene_name.append(x.attr['gene_name'])
            if 'gene_type' in x.attr:
                tmp_gene_type.append(x.attr['gene_type'])
            elif 'gene_biotype' in x.attr:
                tmp_gene_type.append(x.attr['gene_biotype'])
            else:
                pass
        if len(tmp_gene_id) > 0:
            field['gene_id'] = ','.join(tmp_gene_id)
        if len(tmp_gene_name) > 0:
            field['gene_name'] = ','.join(tmp_gene_name)
        if len(tmp_gene_type) > 0:
            field['gene_type'] = ','.join(tmp_gene_type)
    else:
        pass

    return field
