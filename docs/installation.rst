.. highlight:: shell

============
Installation
============

We recommend using conda virtual enviroment::

    $ git clone https://github.com/churchill-lab/gbrs.git
    $ cd gbrs
    $ conda create -n gbrs scipy=0.13.3 cython matplotlib hdf5 ipython jupyter biopython
    $ source activate gbrs
    (gbrs) $ conda install pytables=3.1.0
    (gbrs) $ conda install -c https://conda.binstar.org/bcbio pysam
    (gbrs) $ pip install pysqlite
    (gbrs) $ pip install bx-python
    (gbrs) $ pip install emase
    (gbrs) $ pip install g2gtools
    (gbrs) $ python setup.py install

If you have virtualenvwrapper installed::

    $ mkvirtualenv gbrs
    $ pip install gbrs

Or, at the command line::

    $ easy_install gbrs

