#!/usr/bin/env python
import os
import sys

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


def get_gbrs_version():
    sys.path.insert(0, "gbrs")
    import version
    return version.__version__

requirements = []
on_rtd = os.environ.get('READTHEDOCS', None)
if not on_rtd:
    with open("requirements.txt") as requirements_file:
        requirements_lines = requirements_file.readlines()
        for line in requirements_lines:
            requirements.append(line)

test_requirements = [
    'pytest'
]

setup(
    name='gbrs',
    version=get_gbrs_version(),
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
    license="MIT License",
    zip_safe=False,
    keywords='gbrs',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.10',
    ],
    scripts=['scripts/gbrs',
             'scripts/export-genoprob-file',
             'scripts/convert-kallisto-result-for-gbrs',
             'scripts/convert-salmon-result-for-gbrs',
             'scripts/run_gbrs_on_cluster.sh'],
    test_suite='tests',
    tests_require=test_requirements
)
