#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import codecs

from setuptools import setup, find_packages
from CIRIquant.version import __version__


def read(infile):
    return codecs.open(os.path.join(os.path.dirname(__file__), infile)).read()


setup(
    name='CIRIquant',
    version=__version__,
    author='Jinyang Zhang',
    author_email='zhangjinyang@biols.ac.cn',
    maintainer='Jinyang Zhang',
    maintainer_email='zhangjinyang@biols.ac.cn',
    url='https://github.com/Kevinzjy/CIRIquant',
    description='circular RNA quantification pipeline',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
    ],
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'CIRIquant=CIRIquant.main:main',
            'CIRI_DE=CIRIquant.de:main',
            'CIRI_DE_replicate=CIRIquant.replicate:main',
            'prep_CIRIquant=CIRIquant.prep_CIRIquant:main'
        ]
    },
    keywords='circRNA',
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'argparse==1.2.1', 'PyYAML==5.1.1', 'pysam==0.15.2', 'numpy==1.16.4',
        'scipy==1.2.2', 'scikit-learn==0.20.3', 'numexpr==2.6.9',
    ],
)
