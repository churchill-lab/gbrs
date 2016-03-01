import os
import sys
import numpy as np
from collections import OrderedDict

data_dir = os.getenv('GBRS_DATA', '.')
faifile = os.path.join(data_dir, 'ref.fa.fai')
gposfile = os.path.join(data_dir, 'geneIDs.ordered.npz')  #

try:
    CHRLENS = OrderedDict(np.loadtxt(faifile, usecols=(0, 1), dtype='|S8,<i4'))
except:
    print >> sys.stderr, "Please make sure if $GBRS_DATA is set correctly: %s" % data_dir
    raise
else:
    CHRS = CHRLENS.keys()
    NUM_CHRS = len(CHRS)

try:
    GENEPOS = np.load(gposfile)
except:
    print >> sys.stderr, "Please make sure if $GBRS_DATA is set correctly: %s" % data_dir
    raise
else:
    pass
