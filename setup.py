#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import gbrs
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = []
on_rtd = os.environ.get('READTHEDOCS', None)
if not on_rtd:
    # requirements.append('numpy==1.8.2')
    requirements.append('g2gtools')
    requirements.append('emase')
    requirements.append('matplotlib')
    requirements.append('bx-python>=0.7.2')
    requirements.append('pysam>=0.8.1')
    requirements.append('biopython>=1.63')
    requirements.append('pysqlite>=2.6.3')
    requirements.append('tables>=3.1.0')

test_requirements = [
    'pytest'
]

setup(
    name='gbrs',
    version=gbrs.__version__,
    description='A suite of tools for Reconstructing Genomes and Quantifying Allele Specific Expression from RNA-Seq data',
    long_description=readme + '\n\n' + history,
    author='Kwangbom \"KB\" Choi, Ph.D.',
    author_email='kb.choi@jax.org',
    url='https://github.com/churchill-lab/gbrs',
    packages=[
        'gbrs',
    ],
    package_dir={'gbrs':
                 'gbrs'},
    include_package_data=True,
    install_requires=requirements,
    license="GPLv3",
    zip_safe=False,
    keywords='gbrs',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ],
    scripts=['scripts/gbrs',
             'scripts/convert-kallisto-result-for-gbrs',
             'scripts/convert-salmon-result-for-gbrs'],
    test_suite='tests',
    tests_require=test_requirements
)
