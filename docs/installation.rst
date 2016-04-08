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

**(For DO, CC, or CCRIX)** Make a GBRS_DATA folder, download data files to the folder, and make it visible from your shell. For example,::

    $ mkdir /home/myspace/gbrs
    $ cd /home/myspace/gbrs
    $ wget ftp://churchill-lab.jax.org/pub/software/GBRS/R75-REL1410/\* .
    $ tar xzf 8-way.bowtie1.index.tar.gz
    $ export GBRS_DATA=/home/myspace/gbrs
