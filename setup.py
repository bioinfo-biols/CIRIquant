#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import codecs

from setuptools import setup, find_packages
from CIRIquant.version import __version__


def read(infile):
    return codecs.open(os.path.join(os.path.dirname(__file__), infile)).read()


setup(
    name='CIRI_Quant',
    version=__version__,
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'CIRIquant=CIRIquant.main:main',
        ]
    },
    url='https://github.com/Kevinzjy/CIRI_Quant2',
    license='MIT',
    author='Jinyang Zhang',
    author_email='zhangjinyang@biols.ac.cn',
    description='circular RNA quantification pipeline',
    long_description=read('README.md'),
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
    ],
    keywords='circRNA',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'argparse', 'simplejson',
    ],
)
