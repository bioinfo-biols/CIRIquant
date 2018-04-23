#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import logging
import subprocess
from utils import check_dir
logger = logging.getLogger('CIRIquant')


def align_genome(log_file, thread, reads, outdir, prefix, config):
    # mapping to reference genome
    logger.info('Align reads to reference genome ..')
    align_dir = outdir + '/align'
    check_dir(align_dir)
    hisat_bam = '{}/{}.bam'.format(align_dir, prefix)
    hisat_cmd = '{} -p {} -q -x {} -1 {} -2 {} -t | {} view -bS > {}'.format(
        config['hisat2'],
        thread,
        config['hisat_index'],
        reads[0],
        reads[1],
        config['samtools'],
        hisat_bam
    )
    logger.debug(hisat_cmd)

    with open(log_file, 'a') as log:
        subprocess.call(hisat_cmd, shell=True, stderr=log)
    logger.debug('HISAT2 bam: ' + hisat_bam)
    return hisat_bam


def gene_abundance(log_file, thread, outdir, prefix, hisat_bam, config):
    align_dir = outdir + '/align'
    check_dir(align_dir)

    # sort hisat2 bam
    logger.info('Sorting hisat2 bam ..')
    sorted_bam = '{}/{}.sorted.bam'.format(align_dir, prefix)
    sort_cmd = '{} sort --threads {} -o {} {}'.format(
        config['samtools'],
        thread,
        sorted_bam,
        hisat_bam,
    )
    logger.debug(sort_cmd)
    with open(log_file, 'a') as log:
        subprocess.call(sort_cmd, shell=True, stderr=log)

    # estimate gene expression
    logger.info('Estimate gene abundance ..')
    gene_dir = '{}/gene'.format(outdir)
    check_dir(gene_dir)
    stringtie_cmd = '{0} {1} -e -G {2} -C {3}/{4}_cov.gtf -p {5} -o {3}/{4}_out.gtf -A {3}/{4}_genes.list'.format(
        config['stringtie'],
        sorted_bam,
        config['gtf'],
        gene_dir,
        prefix,
        thread,
    )
    logger.debug(stringtie_cmd)
    with open(log_file, 'a') as log:
        subprocess.call(stringtie_cmd, shell=True, stderr=log)
    logger.debug('Gene expression profile: {}/{}_genes.list'.format(gene_dir, prefix))
    return 1


def run_bwa(log_file, thread, cand_reads, outdir, prefix, config):
    logger.info('Running BWA-mem mapping candidate reads ..')
    bwa_sam = '{}/circ/{}_unmapped.sam'.format(outdir, prefix)
    logger.debug('BWA-mem sam: ' + bwa_sam)

    bwa_cmd = '{} mem -t {} -T 19 {} {} > {}'.format(
        config['bwa'],
        thread,
        config['bwa_index'],
        ' '.join(cand_reads),
        bwa_sam,
    )
    logger.debug(bwa_cmd)
    with open(log_file, 'a') as log:
        subprocess.call(bwa_cmd, shell=True, stderr=log, stdout=log)
    return bwa_sam


def run_ciri(log_file, thread, bwa_sam, outdir, prefix, config):
    logger.info('Running CIRI2 for circRNA detection ..')
    ciri_file = '{}/circ/{}.ciri'.format(outdir, prefix)
    ciri_cmd = 'CIRI2.pl -I {} -O {} -F {} -A {} -0 -T {} -G {}'.format(
        bwa_sam,
        ciri_file,
        config['genome'],
        config['gtf'],
        thread,
        log_file,
    )
    logger.debug(ciri_cmd)
    with open(log_file, 'a') as log:
        subprocess.call(ciri_cmd, shell=True, stderr=log, stdout=log)
    return ciri_file


def convert_bed(ciri_file):
    out_file = ciri_file + '.bed'
    with open(ciri_file, 'r') as f, open(out_file, 'w') as out:
        f.readline()
        for line in f:
            content = line.rstrip().split('\t')
            chrom, start, end = content[1:4]
            out.write('\t'.join([chrom, start, end]) + '\n')
    return out_file
