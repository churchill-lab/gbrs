.. highlight:: shell

============
Installation
============

To avoid conflicts among dependencies, we highly recommend using conda virtual environment::

    $ conda create -n gbrs jupyter
    $ source activate gbrs
    (gbrs) $ conda install -c kbchoi gbrs

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

