import os
import numpy as np
from emase.EMfactory import EMfactory
from emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix as APM


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
    :param grpfile: gene ID to isoform ID mapping info (tsv)
    :return: masked version of alignment incidence file (h5)
    """
    gidfile = os.path.join(DATA_DIR, 'gene_ids.ordered.npz')
    gname_chr = np.load(gidfile)

    gtypefile = kwargs.get('gtypefile')
    gtype_chr = np.load(gtypefile)

    grpfile = kwargs.get('grpfile')
    if grpfile is not None:
        g2t = dict()
        with open(grpfile) as fh:
            for curline in fh:
                item = curline.rstrip().split()
                g2t[item[0]] = item[1:]

    alnfile = kwargs.get('alnfile')
    alnmat = APM(h5file=alnfile)

    gtmask = np.zeros((alnmat.num_haplotypes, alnmat.num_loci))
    hid = dict(zip(alnmat.hname, np.arange(alnmat.num_haplotypes)))
    lid = dict(zip(alnmat.lname, np.arange(alnmat.num_loci)))
    if grpfile is not None:
        for chro in gtype_chr.keys():
            gtype = gtype_chr[chro]
            gname = gname_chr[chro]
            for gid, g in enumerate(gname):
                gt = gtype[gid]
                hid2set = np.array([hid[gt[0]], hid[gt[1]]])
                tid2set = np.array([lid[t] for t in g2t[g]])
                gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
    else:
        for chro in gtype_chr.keys():
            gtype = gtype_chr[chro]
            gname = gname_chr[chro]
            for gid, g in enumerate(gname):
                gt = gtype[gid]
                hid2set = np.array([hid[gt[0]], hid[gt[1]]])
                gtmask[np.meshgrid(hid2set, lid[g])] = 1.0

    alnmat.multiply(gtmask, axis=2)
    outfile = kwargs.get('outfile')
    alnmat.save(h5file=outfile)


def quantify(**kwargs):
    """
    Quantify expected read counts

    :param alnfile: alignment incidence file (h5)
    :param grpfile: gene ID to isoform ID mapping info (tsv)
    :param lenfile: transcript lengths (tsv)
    :param multiread_model: emase model (default: 4)
    :param read_length: read length (default: 100)
    :param pseudocount: prior read count (default: 0.0)
    :param tolerance: tolerance for EM termination (default: 0.0001 in TPM)
    :param max_iters: maximum iterations for EM iteration
    :param report_alignment_counts: whether to report alignment counts (default: False)
    :param report_posterior:
    :return: Expected read counts (tsv)
    """
    alnfile = kwargs.get('alnfile')
    grpfile = kwargs.get('grpfile')
    outbase = kwargs.get('outbase')
    pseudocount = kwargs.get('pseudocount')
    lenfile = kwargs.get('lenfile')
    read_length = kwargs.get('read_length')
    multiread_model = kwargs.get('multiread_model')
    tolerance = kwargs.get('tolerance')
    max_iters = kwargs.get('max_iters')
    report_gene_counts = grpfile is not None
    report_alignment_counts = kwargs.get('report_alignment_counts')
    report_posterior = kwargs.get('report_posterior')

    em_factory = EMfactory(APM(h5file=alnfile, grpfile=grpfile))
    em_factory.prepare(pseudocount=pseudocount, lenfile=lenfile, read_length=read_length)
    em_factory.run(model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True)
    em_factory.report_depths(filename="%s.isoforms.tpm" % outbase, tpm=True)
    em_factory.report_effective_read_counts(filename="%s.isoforms.effective_read_counts" % outbase)
    if report_posterior:
        em_factory.export_posterior_probability(filename="%s.posterior.h5" % outbase)
    if report_gene_counts:
        em_factory.report_depths(filename="%s.genes.tpm" % outbase, tpm=True, grp_wise=True)
        em_factory.report_effective_read_counts(filename="%s.genes.effective_read_counts" % outbase, grp_wise=True)
    if report_alignment_counts:
        alnmat = APM(h5file=alnfile, grpfile=grpfile)
        alnmat.report_alignment_counts(filename="%s.isoforms.alignment_counts" % outbase)
        if report_gene_counts:
            alnmat._bundle_inline(reset=True)
            alnmat.report_alignment_counts(filename="%s.genes.alignment_counts" % outbase)
