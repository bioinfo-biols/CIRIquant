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
        return os.path.abspath(file_name)
    else:
        sys.exit('File: {}, not found'.format(file_name))


def check_dir(dir_name):
    if os.path.exists(dir_name):
        if os.path.isdir(dir_name):
            pass # Output directory already exists
        else:
            sys.exit('Directory: {}, clashed with existed files'.format(dir_name))
    else:
        os.mkdir(dir_name)
    return os.path.abspath(dir_name)


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
        globals()[i.upper()] = check_file(config['tools'][i])

    # check required software version
    check_samtools_version(config['tools']['samtools'])

    # check reference and index
    for i in 'fasta', 'gtf':
        if i not in config['reference']:
            sys.exit('Reference {} need to be specified'.format(i))
        globals()[i.upper()] = check_file(config['reference'][i])

    if 'bwa_index' not in config['reference']:
        sys.exit('BWA index not found')
    BWA_INDEX = os.path.splitext(check_file(config['reference']['bwa_index'] + '.bwt'))[0]

    if 'hisat_index' not in config['reference']:
        sys.exit('HISAT2 index not found')
    HISAT_INDEX = os.path.splitext(os.path.splitext(check_file(config['reference']['hisat_index'] + '.1.ht2'))[0])[0]
    
    return config['name']


def check_samtools_version(samtools):
    from commands import getoutput
    from distutils.version import LooseVersion
    version = getoutput('{} --version'.format(samtools).split('\n')[0].split(' ')[1])
    if version and cmp(LooseVersion(version), LooseVersion('1.9')) < 0:
        sys.exit('samtools version too low, 1.9 required')
    return 1


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


TOOLS = ['CIRI2']


class CIRCparser(object):
    """convert circRNA prediction results to bed file"""
    def __init__(self, circ, tool):
        self.circ = circ
        if tool not in TOOLS:
            sys.exit('Tool: {} not in supported list, please manually convert the circRNA coordinates to bed file.')
        self.tool = tool

    def _ciri2(self):
        circ_data = []
        with open(self.circ, 'r') as f:
            f.readline()
            for line in f:
                content = line.rstrip().split('\t')
                chrom, start, end, strand = content[1], int(content[2]), int(content[3]), content[10]
                circ_data.append([chrom, start, end, strand])
        return circ_data

    def convert(self, out_file):
        circ_data = getattr(self, '_' + self.tool.lower())()
        with open(out_file, 'w') as out:
            for chrom, start, end, strand in circ_data:
                tmp_line = [chrom, start, end, '{}:{}|{}'.format(chrom, start, end), '.', strand]
                out.write('\t'.join([str(x) for x in tmp_line]) + '\n')
        return out_file
