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
    """
    Applying genotype calls to multi-way alignment incidence matrix

    :param alnfile: alignment incidence file (h5),
    :param gtypefile: genotype calls by GBRS (npz),
    :param grpfile: gene ID to isoform ID mapping (tsv)
    :return: masked version of alignment incidence file (h5)
    """
    alnfile = kwargs.get('alnfile')
    alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile)

    gtypefile = kwargs.get('gtypefile')
    gtype_chr = np.load(gtypefile)

    gidfile = os.path.join(DATA_DIR, 'gene_ids.ordered.npz')
    gname_chr = np.load(gidfile)

    g2t = dict()
    grpfile = kwargs.get('grpfile')
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

    outfile = kwargs.get('outfile')
    alnmat.save(h5file=outfile)

def quantify(**kwargs):
    """
    Quantify expected read counts

    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :param:
    :return:
    """
    alnfile = kwargs.get('alnfile')
    alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile)



