import os
import numpy as np
import emase

DATA_DIR = os.getenv('GBRS_DATA', '.')

def hybridize(**kwargs):
    pass


def align(**kwargs):
    pass


def bam2emase(**kwargs):
    pass


def mask(**kwargs):
    alnfile = kwargs.get('alnfile')
    grpfile = kwargs.get('grpfile')
    gtypefile = kwargs.get('gtypefile')

    alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile)
    gtype_chr = np.load(gtypefile)
    gidfile = os.path.join(DATA_DIR, 'gene_ids.npz')
    gname_chr = np.load(gidfile)

    g2t = dict()
    with open(grpfile) as fh:
        for curline in fh:
            item = curline.rstrip().split()
            g2t[item[0]] = item[1:]

    gtmask = np.zeros((alnmat.num_haplotypes, alnmat.num_loci))

    hid = dict(zip(alnmat.hname, np.arange(alnmat.num_haplotypes)))
    tid = dict(zip(alnmat.lname, np.arange(alnmat.num_loci)))

    for chro in gtype_chr.keys():
        gtype = gtype_chr[chro]
        gname = gname_chr[chro]
        for gid, g in enumerate(gname):
            gt = gtype[gid]
            hid2set = np.array([hid[gt[0]], hid[gt[1]]])
            tid2set = np.array([tid[t] for t in g2t[g]])
            gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
    alnmat.multiply(gtmask, axis=2)

def quantify(**kwargs):
    alnfile = kwargs.get('alnfile')
    gtypefile = kwargs.get('gtypefile')
    if gtypefile is not None:
        alnmat = mask(alnmat, **kwargs)
    else:
        alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile)
    pass



