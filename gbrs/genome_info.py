import os
import sys
import numpy as np

data_dir = os.getenv('GBRS_DATA')
faifile = os.path.join(data_dir, 'ref.fa.fai')

try:
    CHRS = np.loadtxt(faifile, usecols=(0,), dtype='string')
    clens = np.loadtxt(faifile, usecols=(1,), dtype=np.int)
    CHRLENS = dict(zip(CHRS, clens))
    NUM_CHRS = len(CHRS)

except (OSError, IOError) as e:
    print >> sys.stderr, "[%s] %s file not found." % (e, faifile)
