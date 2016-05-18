.. highlight:: shell

============
Installation
============

Although GBRS is available at PyPI (https://pypi.python.org/pypi/gbrs/) for 'pip install' or 'easy_install', we highly recommend using Anaconda distribution of python (https://www.continuum.io/downloads) to install all the dependencies without issues. GBRS is also available on Anaconda Cloud, so add the following channels if you have not already::

    $ conda config --add channels r
    $ conda config --add channels bioconda

To avoid conflicts among dependencies, we highly recommend using conda virtual environment::

    $ conda create -n gbrs jupyter
    $ source activate gbrs

Once GBRS virtual environment is created and activated, your shell prompt will show '(gbrs)' at the beginning to specify what virtual environment you are currently in. Now please type the following and install GBRS::

    (gbrs) $ conda install -c kbchoi gbrs

Then, make a folder to store GBRS specific data and set the following environment variable. You may want to add the second line (export command) to your shell rc file (e.g., .bashrc or .bash_profile). For example, ::

    $ mkdir /home/<user>/gbrs
    $ export GBRS_DATA=/home/<user>/gbrs

**(For DO, CC, or CCRIX)** Download data files to $GBRS_DATA folder::

    $ cd $GBRS_DATA
    $ wget ftp://churchill-lab.jax.org/pub/software/GBRS/R84-REL1505/\* .
    $ tar xzf gbrs.hybridized.targets.bowtie-index.tar.gz

That's all! We note that you can go out from GBRS virtual environment anytime once you are done using GBRS::

    (gbrs) $ source deactivate
