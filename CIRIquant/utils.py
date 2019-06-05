#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import logging
LOGGER = logging.getLogger('CIRIquant')

BWA = None
HISAT2 = None
STRINGTIE = None
SAMTOOLS = None

FASTA = None
GTF = None
BWA_INDEX = None
HISAT_INDEX = None


def subprocess_setup():
    """
    Python installs a SIGPIPE handler by default. This is usually not what
    non-Python subprocesses expect
    """
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def get_thread_num(thread):
    from multiprocessing import cpu_count
    cores = cpu_count()
    workers = min(cores, int(thread))
    LOGGER.info('{} CPU cores availble, using {}'.format(cores, workers))
    return workers


def check_file(file_name):
    if os.path.exists(file_name) and os.path.isfile(file_name):
        return 1
    else:
        sys.exit('File: {}, not found'.format(file_name))


def check_dir(dir_name):
    if os.path.exists(dir_name):
        if not os.path.isdir(dir_name):
            sys.exit('Directory: {}, clashed with existed files'.format(dir_name))
    else:
        os.mkdir(dir_name)
    return 1


def check_config(config_file):
    import yaml
    global BWA, HISAT2, STRINGTIE, SAMTOOLS
    global FASTA, GTF, BWA_INDEX, HISAT_INDEX

    # check config reliability
    LOGGER.info('Loading config file: {}'.format(os.path.basename(config_file)))
    with open(config_file, 'r') as infile:
        config = yaml.load(infile, Loader=yaml.BaseLoader)

    # check all tools
    if 'tools' not in config:
        sys.exit('Path of required software must be provided!')

    for i in 'bwa', 'hisat2', 'stringtie', 'samtools':
        if i not in config['tools']:
            sys.exit('Tool: {} need to be specificed'.format(i))
        globals()[i.upper()] = config['tools'][i]

    # check required software version
    check_samtools_version(config['tools']['samtools'])

    # check reference and index
    for i in 'fasta', 'gtf', 'bwa_index', 'hisat_index':
        if i not in config['reference']:
            sys.exit('Reference {} need to be specificed'.format(i))
        globals()[i.upper()] = config['reference'][i]

    return config['name']


def check_samtools_version(samtools):
    from commands import getoutput
    from distutils.version import LooseVersion
    version = getoutput('{} --version'.format(samtools).split('\n')[0].split(' ')[1])
    if version and cmp(LooseVersion(version), LooseVersion('1.3.1')) < 0:
        sys.exit('samtools version too low, 1.3.1 required')
    return 1
