#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import logging
logger = logging.getLogger('CIRIquant')


def align_genome(log_file, thread, reads, outdir, prefix, config):
    logger.info('Align reads to reference genome ..')
