#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys
import time
import logging


def get_logger(logger_name='logger', fname=None, verbosity=False):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG)
    level = logging.DEBUG if verbosity else logging.INFO

    # LOG format
    fmt = "%(asctime)-15s [%(levelname)-5s] %(message)s"
    datefmt = "[%a %Y-%m-%d %H:%M:%S]"
    formatter = logging.Formatter(fmt, datefmt)

    # LOG file
    if fname is not None:
        fh = open(fname, 'w')
        fh.close()
        file_handler = logging.FileHandler(fname)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

    # LOG console
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def find_logger_basefilename(logger):
    log_file = None
    parent = logger.__dict__['parent']
    if parent.__class__.__name__ == 'RootLogger':
        for h in logger.__dict__['handlers']:
            if h.__class__.__name__ == 'FileHandler':
                log_file = h.baseFilename
    else:
        log_file = find_logger_basefilename(parent)
    return log_file


class ProgressBar(object):

    def __init__(self, width=50):
        self.last_x = -1
        self.width = width

    def update(self, x):
        assert 0 <= x <= 100
        if self.last_x == int(x):
            return
        self.last_x = int(x)
        p = int(self.width * (x / 100.0))
        time_stamp = time.strftime("[%a %Y-%m-%d %H:%M:%S]", time.localtime())
        sys.stderr.write('\r%s [%-5s] [%s]' % (time_stamp, str(int(x)) + '%', '#' * p + '.' * (self.width - p)))
        sys.stderr.flush()
        if x == 100:
            sys.stderr.write('\n')
