# standard library imports
from collections import defaultdict
from itertools import dropwhile
from pathlib import Path
import logging
import os
import subprocess
import sys

# 3rd party library imports
import emase
import numpy as np

# local library imports
#

DATA_DIR = os.getenv('GBRS_DATA', '.')
logger = logging.getLogger(__name__)


def is_comment(s):
    return s.startswith('#')


def get_names(idfile):
    ids = dict()
    master_id = 0
    with open(idfile) as fh:
        for curline in fh:
            item = curline.rstrip().split('\t')
            g = item[0]
            if g not in ids:
                ids[g] = master_id
                master_id += 1
    num_ids = len(ids)
    names = {index: name for name, index in ids.items()}
    return [names[k] for k in range(num_ids)]


def hybridize(
    fasta_list: list[Path | str],
    hap_list: str,
    out_file: Path | str = 'gbrs.hybridized.targets.fa',
    build_bowtie_index: bool = False
):
    out_dir = os.path.dirname(out_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for x in fasta_list:
        logger.info(f'Fasta File: {x}')
    logger.info(f'Haplotype List: {hap_list}')
    logger.info(f'Output File: {out_file}')

    # Get pooled transcriptome
    outbase = os.path.splitext(out_file)[0]
    hap_list = hap_list.split(',')
    num_haps = len(fasta_list)
    lenfile = f'{outbase}.info'
    seqout = open(out_file, 'w')
    lenout = open(lenfile, 'w')
    for hid in range(num_haps):
        fasta = fasta_list[hid]
        hapname = hap_list[hid]
        logger.warning(
            f"[gbrs::hybridize] Adding suffix '_{hapname}' to the sequence ID's of {fasta}"
        )
        fh = open(fasta)
        line = fh.readline()  # The first fasta header
        line = f'{line.rstrip().split()[0]}_{hapname}'
        seqout.write(f'{line}\n')
        lenout.write(f'{line[1:]}\t')
        seq_len = 0
        for line in fh:
            if line[0] == '>':
                line = f'{line.rstrip().split()[0]}_{hapname}\n'
                lenout.write(f'{seq_len}\n{line[1:].rstrip()}\t')
                seq_len = 0
            else:
                seq_len += len(line.rstrip())
            seqout.write(line)
        fh.close()
        lenout.write(f'{seq_len}\n')
    seqout.close()
    lenout.close()

    # Build bowtie index for the pooled transcriptome
    if build_bowtie_index:
        out_index = f'{outbase}.bowtie1'
        logger.warning('[gbrs::hybridize] Building bowtie1 index')
        status = subprocess.call(f'bowtie-build outfile out_index', shell=True)


def bam2emase(
    alignment_file: Path | str,
    haplotypes: str,
    locusid_file: Path | str = None,
    out_file: Path | str = None,
    delim: str = '_',
    index_dtype: str = 'uint32',
    data_dtype: str = 'uint8'
) -> None:

    if locusid_file is None:
        locusid_file = os.path.join(DATA_DIR, 'ref.transcripts.info')
        if not os.path.exists(locusid_file):
            logger.error('Cannot find a locus id file.')

    if out_file is None:
        out_file = (
            'gbrs.bam2emased.'
            f'{os.path.splitext(os.path.basename(alignment_file))[0]}.h5'
        )

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Locus ID File: {locusid_file}')
    logger.info(f'Output File: {out_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Delimiter: {delim}')
    logger.info(f'Index dtype: {index_dtype}')
    logger.info(f'Data dtype: {data_dtype}')

    haplotypes = tuple(haplotypes.split(','))
    loci = get_names(locusid_file)

    logger.debug('Constructing AMF')
    alignmat_factory = emase.AlignmentMatrixFactory(alignment_file)
    alignmat_factory.prepare(
        haplotypes, loci, delim=delim, outdir=os.path.dirname(out_file)
    )
    alignmat_factory.produce(
        out_file, index_dtype=index_dtype, data_dtype=data_dtype
    )
    alignmat_factory.cleanup()

    logger.debug('Done')


def intersect(
    emase_files: list[Path | str],
    out_file: Path | str = None,
    comp_lib: str = 'zlib'
) -> None:
    if out_file is None:
        out_file = 'gbrs.intersected.' + os.path.basename(emase_files[0])

    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {out_file}')
    logger.info(f'Compression Library: {comp_lib}')

    logger.info(f'Loading {emase_files[0]}')
    aln_mat = emase.AlignmentPropertyMatrix(h5file=emase_files[0])
    for f in emase_files[1:]:
        logger.info(f'Loading {f}')
        alnmat_next = emase.AlignmentPropertyMatrix(h5file=f)

        if np.all(aln_mat.rname == alnmat_next.rname):
            aln_mat = aln_mat * alnmat_next
        else:
            logger.error('The read ID\'s are not compatible.')
            raise ValueError('The read ID\'s are not compatible.')

    logger.info(f'Saving {out_file}')
    aln_mat.save(h5file=out_file, complib=comp_lib)
    logger.info('Done')


def compress(
    emase_files: list[Path | str],
    out_file: Path | str,
    comp_lib: str = 'zlib'
) -> None:
    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {out_file}')
    logger.info(f'Compression Library: {comp_lib}')

    num_loci = None
    num_haplotypes = None
    names_loci = None
    names_haplotypes = None

    ec = defaultdict(int)
    for aln_file in emase_files:
        logger.info(f'Loading {aln_file}')

        alnmat_rd = emase.AlignmentPropertyMatrix(h5file=aln_file)

        logger.debug(f'Number Loci: {alnmat_rd.num_loci}')
        logger.debug(f'Number Haplotypes: {alnmat_rd.num_haplotypes}')
        logger.debug(f'Number Reads: {alnmat_rd.num_reads}')

        # each file should be the same
        num_loci = alnmat_rd.num_loci
        num_haplotypes = alnmat_rd.num_haplotypes
        names_loci = alnmat_rd.lname
        names_haplotypes = alnmat_rd.hname

        for h in range(alnmat_rd.num_haplotypes):
            alnmat_rd.data[h] = alnmat_rd.data[h].tocsr()

        if alnmat_rd.count is None:
            alnmat_rd.count = np.ones(alnmat_rd.num_reads)

        for cur_ind in range(alnmat_rd.num_reads):
            ec_key = []
            for h in range(alnmat_rd.num_haplotypes):
                i0 = alnmat_rd.data[h].indptr[cur_ind]
                i1 = alnmat_rd.data[h].indptr[cur_ind + 1]
                ec_key.append(
                    ','.join(
                        map(str, sorted(alnmat_rd.data[h].indices[i0:i1]))
                    )
                )
            ec[':'.join(ec_key)] += alnmat_rd.count[cur_ind]

    ec = dict(ec)
    num_ecs = len(ec)

    logger.debug('Constructing APM')
    logger.debug(f'Number Loci: {num_loci}')
    logger.debug(f'Number Haplotypes: {num_haplotypes}')
    logger.debug(f'Number ECs: {num_ecs}')

    alnmat_ec = emase.AlignmentPropertyMatrix(
        shape=(num_loci, num_haplotypes, num_ecs)
    )
    alnmat_ec.hname = names_haplotypes
    alnmat_ec.lname = names_loci
    alnmat_ec.count = np.zeros(num_ecs)

    for row_id, ec_key in enumerate(ec):
        alnmat_ec.count[row_id] = ec[ec_key]
        nzlocs = ec_key.split(':')
        for h in range(alnmat_ec.num_haplotypes):
            nzlocs_h = nzlocs[h]
            if nzlocs_h != '':
                nzinds = np.array(list(map(int, nzlocs_h.split(','))))
                alnmat_ec.data[h][row_id, nzinds] = 1
    alnmat_ec.finalize()
    alnmat_ec.save(h5file=out_file, complib=comp_lib)

    logger.debug('Done')


def stencil(
    alignment_file: Path | str,
    genotype_file: Path | str,
    group_file: Path | str,
    out_file: Path | str
):
    """
    Applying genotype calls to multi-way alignment incidence matrix

    :param alnfile: alignment incidence file (h5),
    :param gtypefile: genotype calls by GBRS (tsv),
    :param grpfile: gene ID to isoform ID mapping info (tsv)
    :return: genotyped version of alignment incidence file (h5)
    """

    if group_file is None:
        group_file = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if not os.path.exists(group_file):
            logger.info('A group file is *not* given. Genotype will be stenciled as is.')

    if out_file is None:
        out_file = f'gbrs.stenciled.{os.path.basename(alignment_file)}'

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Genotype File: {genotype_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Output File: {out_file}')

    # Load alignment incidence matrix ('alnfile' is assumed to be in multiway transcriptome)
    logger.info(f'Loading {alignment_file}')
    aln_mat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)

    # Load genotype calls
    logger.info(f'Loading and processing genotype calls {genotype_file}')
    hid = dict(zip(aln_mat.hname, np.arange(aln_mat.num_haplotypes)))
    gid = dict(zip(aln_mat.gname, np.arange(len(aln_mat.gname))))
    gtmask = np.zeros((aln_mat.num_haplotypes, aln_mat.num_loci))
    gtcall_g = dict.fromkeys(aln_mat.gname)
    with open(genotype_file) as fh:
        if group_file is not None:
            gtcall_t = dict.fromkeys(aln_mat.lname)
            for line in dropwhile(is_comment, fh):
                item = line.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(aln_mat.groups[gid[g]])
                gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
                for t in tid2set:
                    gtcall_t[aln_mat.lname[t]] = gt
        else:
            for line in dropwhile(is_comment, fh):
                item = line.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                gtmask[np.meshgrid(hid2set, gid[g])] = 1.0

    aln_mat.multiply(gtmask, axis=2)
    for h in range(aln_mat.num_haplotypes):
        aln_mat.data[h].eliminate_zeros()

    logger.info(f'Saving {out_file}')
    aln_mat.save(h5file=out_file)
    logger.info('Done')


def quantify(
    alignment_file: Path | str,
    group_file: Path | str = None,
    length_file: Path | str = None,
    genotype_file: Path | str = None,
    outbase: str = 'gbrs.quantified',
    multiread_model: int = 4,
    pseudocount: float = 0.0,
    max_iters: int = 999,
    tolerance: float = 0.0001,
    report_alignment_counts: bool = False,
    report_posterior: bool = False
) -> None:
    """
    Quantify expected read counts

    :param alnfile: alignment incidence file (h5)
    :param grpfile: gene ID to isoform ID mapping info (tsv)
    :param lenfile: transcript lengths (tsv)
    :param multiread_model: emase model (default: 4)
    :param pseudocount: prior read count (default: 0.0)
    :param tolerance: tolerance for EM termination (default: 0.0001 in TPM)
    :param max_iters: maximum iterations for EM iteration
    :param report_alignment_counts: whether to report alignment counts (default: False)
    :param report_posterior:
    :return: Expected read counts (tsv)
    """
    if group_file is None:
        group_file = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if not os.path.exists(group_file):
            logger.warning('A group file is not given. Group-level results will not be reported.')

    if length_file is None:
        length_file = os.path.join(DATA_DIR, 'gbrs.hybridized.targets.info')
        if not os.path.exists(length_file):
            logger.warning('A length file is not given. Transcript length adjustment will *not* be performed.',)

    # If group_file exist, always report groupwise results too
    report_group_counts = (
        group_file is not None
    )

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Length File: {length_file}')
    logger.info(f'Genotype File: {genotype_file}')
    logger.info(f'Outbase: {outbase}')
    logger.info(f'Multiread Model: {multiread_model}')
    logger.info(f'Pseudocount: {pseudocount}')
    logger.info(f'Tolerance: {tolerance}')
    logger.info(f'Report Alignment Counts: {report_alignment_counts}')
    logger.info(f'Report Posterior: {report_posterior}')

    # load alignment incidence matrix ('alignment_file' is assumed to be in multiway transcriptome)
    logger.debug('Loading APM')
    aln_mat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)

    # Load genotype calls
    if genotype_file is not None:
        # Genotype calls are at the gene level
        logger.debug('Loading genotype calls')
        outbase = f'{outbase}.diploid'
        logger.debug(f'Outbase now: {outbase}')
        # haplotype as key and number as value
        hid = dict(zip(aln_mat.hname, np.arange(aln_mat.num_haplotypes)))
        # gene id as key and number as value
        gid = dict(zip(aln_mat.gname, np.arange(len(aln_mat.gname))))
        gtmask = np.zeros((aln_mat.num_haplotypes, aln_mat.num_loci))
        # genome genotype calls, gene_id as key
        gtcall_g = dict.fromkeys(aln_mat.gname)
        # transcript genotype calls, transcript_id as key
        gtcall_t = dict.fromkeys(aln_mat.lname)
        # count = 0
        # for k, v in gtcall_t.items():
        #    print(k, v)
        #    if count > 10:
        #        break
        #    count += 1
        with open(genotype_file) as fh:
            for curline in dropwhile(is_comment, fh):
                item = curline.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(aln_mat.groups[gid[g]])
                gtmask[tuple(np.meshgrid(hid2set, tid2set))] = 1.0
                for t in tid2set:
                    gtcall_t[aln_mat.lname[t]] = gt

        aln_mat.multiply(gtmask, axis=2)
        for h in range(aln_mat.num_haplotypes):
            aln_mat.data[h].eliminate_zeros()
    else:
        outbase = f'{outbase}.multiway'
        logger.debug(f'Outbase now: {outbase}')
        gtcall_g = None
        gtcall_t = None

    # Run emase
    logger.debug('Running emase')
    em_factory = emase.EMfactory(aln_mat)
    em_factory.prepare(pseudocount=pseudocount, lenfile=length_file)
    em_factory.run(
        model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True
    )
    em_factory.report_depths(
        filename=f'{outbase}.isoforms.tpm', tpm=True, notes=gtcall_t
    )
    em_factory.report_read_counts(
        filename=f'{outbase}.isoforms.expected_read_counts', notes=gtcall_t
    )
    if report_posterior:
        em_factory.export_posterior_probability(
            filename=f'{outbase}.posterior.h5'
        )
    if report_group_counts:
        em_factory.report_depths(
            filename=f'{outbase}.genes.tpm',
            tpm=True,
            grp_wise=True,
            notes=gtcall_g,
        )
        em_factory.report_read_counts(
            filename=f'{outbase}.genes.expected_read_counts',
            grp_wise=True,
            notes=gtcall_g,
        )
    if report_alignment_counts:
        alnmat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
        alnmat.report_alignment_counts(
            filename=f'{outbase}.isoforms.alignment_counts'
        )
        if report_group_counts:
            alnmat._bundle_inline(reset=True)
            alnmat.report_alignment_counts(
                filename=f'{outbase}.genes.alignment_counts'
            )
    logger.debug('Done')
