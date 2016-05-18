import os
import sys
import numpy as np
import subprocess
from itertools import dropwhile
from collections import defaultdict
import emase


DATA_DIR = os.getenv('GBRS_DATA', '.')


def is_comment(s):
    return s.startswith('#')


def get_names(idfile):
    ids = dict()
    master_id = 0
    with open(idfile) as fh:
        for curline in fh:
            item = curline.rstrip().split("\t")
            g = item[0]
            if not ids.has_key(g):
                ids[g] = master_id
                master_id += 1
    num_ids = len(ids)
    names = {index:name for name, index in ids.iteritems()}
    return [names[k] for k in xrange(num_ids)]


def hybridize(**kwargs):
    outfile = kwargs.get('outfile')
    outdir = os.path.dirname(outfile)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Get pooled transcriptome
    outbase = os.path.splitext(outfile)[0]
    fastalist = kwargs.get('fastalist').split(',')
    haplist = kwargs.get('haplist').split(',')
    num_haps = len(fastalist)
    lenfile = outbase + '.info'
    seqout = open(outfile, 'w')
    lenout = open(lenfile, 'w')
    for hid in xrange(num_haps):
        fasta = fastalist[hid]
        hapname = haplist[hid]
        print >> sys.stderr, "[gbrs::hybridize] Adding suffix \'_%s\' to the sequence ID's of %s..." % (hapname, fasta)
        fh = open(fasta)
        curline = fh.next()  # The first fasta header
        curline = curline.rstrip().split()[0] + '_' + hapname
        seqout.write('%s\n' % curline)
        lenout.write('%s\t' % curline[1:])
        seqlen = 0
        for curline in fh:
            if curline[0] == '>':
                curline = curline.rstrip().split()[0] + '_' + hapname + "\n"
                lenout.write('%d\n%s\t' % (seqlen, curline[1:].rstrip()))
                seqlen = 0
            else:
                seqlen += len(curline.rstrip())
            seqout.write(curline)
        fh.close()
        lenout.write('%d\n' % seqlen)
    seqout.close()
    lenout.close()

    # Build bowtie index for the pooled transcriptome
    build_bowtie_index = kwargs.get('build_bowtie_index')
    if build_bowtie_index:
        out_index = outbase + '.bowtie1'
        print >> sys.stderr, "[gbrs::hybridize] Building bowtie1 index..."
        status = subprocess.call("bowtie-build %s %s" % (outfile, out_index), shell=True)


def align(**kwargs):
    raise NotImplementedError


def bam2emase(**kwargs):
    alnfile = kwargs.get('alnfile')
    lidfile = kwargs.get('lidfile')
    if lidfile is None:
        lidfile2chk = os.path.join(DATA_DIR, 'ref.transcripts.info')
        if os.path.exists(lidfile2chk):
            lidfile = lidfile2chk
        else:
            print >> sys.stderr, '[gbrs::bam2emase] Cannot find a locus id file.'
            return 2
    outfile = kwargs.get('outfile')
    if outfile is None:
        outfile = 'gbrs.bam2emased.' + os.path.splitext(os.path.basename(alnfile))[0] + '.h5'

    haplotypes = tuple(kwargs.get('haplogypes').split(','))
    index_dtype = kwargs.get('index_dtype')
    data_dtype = kwargs.get('data_dtype')

    loci = get_names(lidfile)
    alignmat_factory = emase.AlignmentMatrixFactory(alnfile)
    alignmat_factory.prepare(haplotypes, loci, outdir=os.path.dirname(outfile))
    alignmat_factory.produce(outfile, index_dtype=index_dtype, data_dtype=data_dtype)
    alignmat_factory.cleanup()


def intersect(**kwargs):
    flist = kwargs.get('emasefiles').split(',')
    outfile = kwargs.get('outfile')
    if outfile is None:
        outfile = 'gbrs.intersected.' + os.path.basename(flist[0])
    complib = kwargs.get('complib')
    alnmat = emase.AlignmentPropertyMatrix(h5file=flist[0])
    for f in flist[1:]:
        alnmat_next = emase.AlignmentPropertyMatrix(h5file=f)
        if np.all(alnmat.rname == alnmat_next.rname):
            alnmat = alnmat * alnmat_next
        else:
            print >> sys.stderr, "[gbrs::intersect] The read ID's are not compatible."
            return 2
    alnmat.save(h5file=outfile, complib=complib)


def compress(**kwargs):
    flist = kwargs.get('emasefiles').split(',')
    outfile = kwargs.get('outfile')
    complib = kwargs.get('complib')

    ec = defaultdict(int)
    for alnfile in flist:
        alnmat_rd = emase.AlignmentPropertyMatrix(h5file=alnfile)
        for h in xrange(alnmat_rd.num_haplotypes):
            alnmat_rd.data[h] = alnmat_rd.data[h].tocsr()
        if alnmat_rd.count is None:
            alnmat_rd.count = np.ones(alnmat_rd.num_reads)
        for curind in xrange(alnmat_rd.num_reads):
            ec_key = []
            for h in xrange(alnmat_rd.num_haplotypes):
                i0 = alnmat_rd.data[h].indptr[curind]
                i1 = alnmat_rd.data[h].indptr[curind+1]
                ec_key.append(','.join(map(str, sorted(alnmat_rd.data[h].indices[i0:i1]))))
            ec[':'.join(ec_key)] += alnmat_rd.count[curind]
    ec = dict(ec)
    num_ecs = len(ec)
    alnmat_ec = emase.AlignmentPropertyMatrix(shape=(alnmat_rd.num_loci, alnmat_rd.num_haplotypes, num_ecs))
    alnmat_ec.hname = alnmat_rd.hname
    alnmat_ec.lname = alnmat_rd.lname
    alnmat_ec.count = np.zeros(num_ecs)
    for row_id, ec_key in enumerate(ec):
        alnmat_ec.count[row_id] = ec[ec_key]
        nzlocs = ec_key.split(':')
        for h in xrange(alnmat_ec.num_haplotypes):
            nzlocs_h = nzlocs[h]
            if nzlocs_h != '':
                nzinds = np.array(map(np.int, nzlocs_h.split(',')))
                alnmat_ec.data[h][row_id, nzinds] = 1
    alnmat_ec.finalize()
    alnmat_ec.save(h5file=outfile, complib=complib)


def stencil(**kwargs):
    """
    Applying genotype calls to multi-way alignment incidence matrix

    :param alnfile: alignment incidence file (h5),
    :param gtypefile: genotype calls by GBRS (tsv),
    :param grpfile: gene ID to isoform ID mapping info (tsv)
    :return: genotyped version of alignment incidence file (h5)
    """
    alnfile = kwargs.get('alnfile')
    gtypefile = kwargs.get('gtypefile')
    grpfile = kwargs.get('grpfile')
    if grpfile is None:
        grpfile2chk = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if os.path.exists(grpfile2chk):
            grpfile = grpfile2chk
        else:
            print >> sys.stderr, '[gbrs::stencil] A group file is *not* given. Genotype will be stenciled as is.'

    # Load alignment incidence matrix ('alnfile' is assumed to be in multiway transcriptome)
    alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile, grpfile=grpfile)

    # Load genotype calls
    hid = dict(zip(alnmat.hname, np.arange(alnmat.num_haplotypes)))
    gid = dict(zip(alnmat.gname, np.arange(len(alnmat.gname))))
    gtmask = np.zeros((alnmat.num_haplotypes, alnmat.num_loci))
    gtcall_g = dict.fromkeys(alnmat.gname)
    with open(gtypefile) as fh:
        if grpfile is not None:
            gtcall_t = dict.fromkeys(alnmat.lname)
            for curline in dropwhile(is_comment, fh):
                item = curline.rstrip().split("\t")
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(alnmat.groups[gid[g]])
                gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
                for t in tid2set:
                    gtcall_t[alnmat.lname[t]] = gt
        else:
            for curline in dropwhile(is_comment, fh):
                item = curline.rstrip().split("\t")
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                gtmask[np.meshgrid(hid2set, gid[g])] = 1.0

    alnmat.multiply(gtmask, axis=2)
    for h in xrange(alnmat.num_haplotypes):
        alnmat.data[h].eliminate_zeros()
    outfile = kwargs.get('outfile')
    if outfile is None:
        outfile = 'gbrs.stenciled.' + os.path.basename(alnfile)
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
    if grpfile is None:
        grpfile2chk = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if os.path.exists(grpfile2chk):
            grpfile = grpfile2chk
        else:
            print >> sys.stderr, '[gbrs::quantify] A group file is not given. Group-level results will not be reported.'

    outbase = kwargs.get('outbase')
    gtypefile = kwargs.get('gtypefile')
    pseudocount = kwargs.get('pseudocount')
    lenfile = kwargs.get('lenfile')
    if lenfile is None:
        lenfile2chk = os.path.join(DATA_DIR, 'gbrs.hybridized.targets.info')
        if os.path.exists(lenfile2chk):
            lenfile = lenfile2chk
        else:
            print >> sys.stderr, '[gbrs::quantify] A length file is not given. Transcript length adjustment will *not* be performed.'

    read_length = kwargs.get('read_length')
    multiread_model = kwargs.get('multiread_model')
    tolerance = kwargs.get('tolerance')
    max_iters = kwargs.get('max_iters')
    report_group_counts = grpfile is not None  # If grpfile exist, always report groupwise results too
    report_alignment_counts = kwargs.get('report_alignment_counts')
    report_posterior = kwargs.get('report_posterior')

    # Load alignment incidence matrix ('alnfile' is assumed to be in multiway transcriptome)
    alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile, grpfile=grpfile)

    # Load genotype calls
    if gtypefile is not None:  # Genotype calls are at the gene level
        outbase = outbase + '.diploid'
        hid = dict(zip(alnmat.hname, np.arange(alnmat.num_haplotypes)))
        gid = dict(zip(alnmat.gname, np.arange(len(alnmat.gname))))
        gtmask = np.zeros((alnmat.num_haplotypes, alnmat.num_loci))
        gtcall_g = dict.fromkeys(alnmat.gname)
        gtcall_t = dict.fromkeys(alnmat.lname)
        with open(gtypefile) as fh:
            for curline in dropwhile(is_comment, fh):
                item = curline.rstrip().split("\t")
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(alnmat.groups[gid[g]])
                gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
                for t in tid2set:
                    gtcall_t[alnmat.lname[t]] = gt
        alnmat.multiply(gtmask, axis=2)
        for h in xrange(alnmat.num_haplotypes):
            alnmat.data[h].eliminate_zeros()
    else:
        outbase = outbase + ".multiway"
        gtcall_g = None
        gtcall_t = None

    # Run emase
    em_factory = emase.EMfactory(alnmat)
    em_factory.prepare(pseudocount=pseudocount, lenfile=lenfile, read_length=read_length)
    em_factory.run(model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True)
    em_factory.report_depths(filename="%s.isoforms.tpm" % outbase, tpm=True, notes=gtcall_t)
    em_factory.report_read_counts(filename="%s.isoforms.expected_read_counts" % outbase, notes=gtcall_t)
    if report_posterior:
        em_factory.export_posterior_probability(filename="%s.posterior.h5" % outbase)
    if report_group_counts:
        em_factory.report_depths(filename="%s.genes.tpm" % outbase, tpm=True, grp_wise=True, notes=gtcall_g)
        em_factory.report_read_counts(filename="%s.genes.expected_read_counts" % outbase, grp_wise=True, notes=gtcall_g)
    if report_alignment_counts:
        alnmat = emase.AlignmentPropertyMatrix(h5file=alnfile, grpfile=grpfile)
        alnmat.report_alignment_counts(filename="%s.isoforms.alignment_counts" % outbase)
        if report_group_counts:
            alnmat._bundle_inline(reset=True)
            alnmat.report_alignment_counts(filename="%s.genes.alignment_counts" % outbase)
