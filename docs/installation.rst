.. highlight:: shell

============
Installation
============

We recommend using conda virtual enviroment::

    $ git clone https://github.com/churchill-lab/gbrs.git
    $ cd gbrs
    $ conda create -n gbrs jupyter scipy=0.13.3 cython matplotlib biopython
    $ source activate gbrs
    (gbrs) $ conda install pytables=3.1.0
    (gbrs) $ conda install -c https://conda.binstar.org/bcbio pysam
    (gbrs) $ pip install pysqlite
    (gbrs) $ pip install bx-python
    (gbrs) $ pip install emase
    (gbrs) $ pip install g2gtools
    (gbrs) $ python setup.py install

Or if you have virtualenvwrapper installed::

    $ mkvirtualenv gbrs
    $ pip install gbrs

Or at the command line::

    $ easy_install gbrs

Then, make a folder to store GBRS specific data and set the following environment variable. You may want to add the second line (export) to your shell rc file (e.g., .bashrc or .bash_profile). For example,::

    $ mkdir /home/gbrs
    $ export GBRS_DATA=/home/gbrs

**(For DO, CC, or CCRIX)** Download data files to $GBRS_DATA folder::

    $ cd $GBRS_DATA
    $ wget ftp://churchill-lab.jax.org/pub/software/GBRS/R75-REL1410/\* .
    $ tar xzf gbrs.hybridized.targets.bowtie-index.tar.gz

