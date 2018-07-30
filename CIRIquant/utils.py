#! /usr/bin/env python
# -*- encoding:utf-8 -*=
import os
import sys
import logging
LOGGER = logging.getLogger('CIRIquant')


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
    import simplejson as json
    # check config reliability
    LOGGER.info('Loading config file: {}'.format(os.path.basename(config_file)))
    config = {i: i for i in ['samtools', 'bwa', 'hisat2', 'stringtie']}
    with open(config_file, 'r') as infile:
        config.update(json.load(infile))

    # check required software version
    check_samtools_version(config['samtools'])
    if 'genome' not in config:
        sys.exit('Please provide reference genome')
    if 'hisat_index' not in config:
        sys.exit('Please provide HISAT2 index for reference genome')
    return config


def check_samtools_version(samtools):
    from commands import getoutput
    from distutils.version import LooseVersion
    version = getoutput('{} --version'.format(samtools).split('\n')[0].split(' ')[1])
    if version and cmp(LooseVersion(version), LooseVersion('1.3.1')) < 0:
        sys.exit('samtools version too low, 1.3.1 required')
    return 1
