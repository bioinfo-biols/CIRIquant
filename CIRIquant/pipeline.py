#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import logging
import subprocess
import utils
LOGGER = logging.getLogger('CIRIquant')


def align_genome(log_file, thread, reads, outdir, prefix):
    # mapping to reference genome
    align_dir = outdir + '/align'
    utils.check_dir(align_dir)
    hisat_bam = '{}/{}.bam'.format(align_dir, prefix)
    hisat_cmd = '{} -p {} --dta -q -x {} -1 {} -2 {} -t | {} view -bS > {}'.format(
        utils.HISAT2,
        thread,
        utils.HISAT_INDEX,
        reads[0],
        reads[1],
        utils.SAMTOOLS,
        hisat_bam
    )

    # sort hisat2 bam
    sorted_bam = '{}/{}.sorted.bam'.format(align_dir, prefix)
    sort_cmd = '{} sort --threads {} -o {} {}'.format(
        utils.SAMTOOLS,
        thread,
        sorted_bam,
        hisat_bam,
    )

    index_cmd = '{} index -@ {} {}'.format(
        utils.SAMTOOLS,
        thread,
        sorted_bam,
    )

    with open(log_file, 'a') as log:
        LOGGER.debug(hisat_cmd)
        subprocess.call(hisat_cmd, shell=True, stderr=log)

        LOGGER.debug(sort_cmd)
        subprocess.call(sort_cmd, shell=True, stderr=log)

        LOGGER.debug(index_cmd)
        subprocess.call(index_cmd, shell=True, stderr=log)
    
    if os.path.getsize(sorted_bam + '.bai') <= 16:
        raise utils.PipelineError('Empty hisat2 bam generated, please re-run CIRIquant with -v and check the fastq and hisat2-index.')

    return sorted_bam


def gene_abundance(log_file, thread, outdir, prefix, hisat_bam):
    align_dir = outdir + '/align'
    utils.check_dir(align_dir)

    # estimate gene expression
    LOGGER.info('Estimate gene abundance ..')
    gene_dir = '{}/gene'.format(outdir)
    utils.check_dir(gene_dir)

    stringtie_cmd = '{0} {1} -e -G {2} -C {3}/{4}_cov.gtf -p {5} -o {3}/{4}_out.gtf -A {3}/{4}_genes.list'.format(
        utils.STRINGTIE,
        hisat_bam,
        utils.GTF,
        gene_dir,
        prefix,
        thread,
    )

    with open(log_file, 'a') as log:
        LOGGER.debug(stringtie_cmd)
        subprocess.call(stringtie_cmd, shell=True, stderr=log)

    LOGGER.debug('Gene expression profile: {}/{}_genes.list'.format(gene_dir, prefix))
    return 1


def run_bwa(log_file, thread, cand_reads, outdir, prefix):
    LOGGER.info('Running BWA-mem mapping candidate reads ..')
    utils.check_dir('{}/circ'.format(outdir))

    bwa_sam = '{}/circ/{}_unmapped.sam'.format(outdir, prefix)

    bwa_cmd = '{} mem -t {} -T 19 {} {} > {}'.format(
        utils.BWA,
        thread,
        utils.BWA_INDEX,
        ' '.join(cand_reads),
        bwa_sam,
    )
    with open(log_file, 'a') as log:
        LOGGER.debug(bwa_cmd)
        subprocess.call(bwa_cmd, shell=True, stderr=log, stdout=log)

    LOGGER.debug('BWA-mem sam: ' + bwa_sam)
    return bwa_sam


def run_ciri(log_file, thread, bwa_sam, outdir, prefix):
    LOGGER.info('Running CIRI2 for circRNA detection ..')
    ciri_file = '{}/circ/{}.ciri'.format(outdir, prefix)
    ciri_cmd = 'CIRI2.pl -I {} -O {} -F {} -A {} -0 -T {} -G {}'.format(
        bwa_sam,
        ciri_file,
        utils.FASTA,
        utils.GTF,
        thread,
        log_file,
    )

    with open(log_file, 'a') as log:
        LOGGER.debug(ciri_cmd)
        subprocess.call(ciri_cmd, shell=True, stderr=log, stdout=log)

    return ciri_file


def convert_bed(ciri_file):
    out_file = ciri_file + '.bed'
    with open(ciri_file, 'r') as f, open(out_file, 'w') as out:
        f.readline()
        for line in f:
            content = line.rstrip().split('\t')
            chrom, start, end = content[1:4]
            strand = content[10]
            out.write('\t'.join([chrom, start, end, '{}:{}|{}'.format(chrom, start, end), '.', strand]) + '\n')
    return out_file


def clean_tmp(outdir, prefix):
    tmp_files = [
        '{}/circ/{}_unmapped.sam'.format(outdir, prefix),
        '{}/circ/{}_denovo.bam'.format(outdir, prefix),
        '{}/align/{}.bam'.format(outdir, prefix),
    ]
    for f in tmp_files:
        if os.path.exists(f) and os.path.isfile(f):
            os.remove(f)
    return 0
