# standard library imports
from collections import defaultdict
from itertools import dropwhile
from pathlib import Path
import logging
import os
import subprocess

# 3rd party library imports
import emase
import numpy as np

# local library imports
#

DATA_DIR = os.getenv('GBRS_DATA', '.')
logger = logging.getLogger(__name__)


def is_comment(s: str) -> bool:
    return s.startswith('#')


def get_names(id_file: str) -> list[str]:
    ids = dict()
    master_id = 0
    with open(id_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            g = item[0]
            if g not in ids:
                ids[g] = master_id
                master_id += 1
    num_ids = len(ids)
    names = {index: name for name, index in ids.items()}
    return [names[k] for k in range(num_ids)]


def compress(
    emase_files: list[Path | str],
    out_file: Path | str,
    comp_lib: str = 'zlib'
) -> None:
    """
    Compress EMASE file to alignment incidence matrix

    Args:
        emase_files: list of EMASE files to compress
        out_file: name of the compressed EMASE file
        comp_lib: compression library to use
    """
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
        logger.info(f'Loading EMASE file: {aln_file}')
        aln_mat_rd = emase.AlignmentPropertyMatrix(h5file=aln_file)

        logger.debug(f'Number Loci: {aln_mat_rd.num_loci}')
        logger.debug(f'Number Haplotypes: {aln_mat_rd.num_haplotypes}')
        logger.debug(f'Number Reads: {aln_mat_rd.num_reads}')

        # each file should be the same
        num_loci = aln_mat_rd.num_loci
        num_haplotypes = aln_mat_rd.num_haplotypes
        names_loci = aln_mat_rd.lname
        names_haplotypes = aln_mat_rd.hname

        for h in range(aln_mat_rd.num_haplotypes):
            aln_mat_rd.data[h] = aln_mat_rd.data[h].tocsr()

        if aln_mat_rd.count is None:
            aln_mat_rd.count = np.ones(aln_mat_rd.num_reads)

        logger.debug('Creating unique ECs')
        for cur_ind in range(aln_mat_rd.num_reads):
            ec_key = []
            for h in range(aln_mat_rd.num_haplotypes):
                i0 = aln_mat_rd.data[h].indptr[cur_ind]
                i1 = aln_mat_rd.data[h].indptr[cur_ind + 1]
                ec_key.append(
                    ','.join(
                        map(str, sorted(aln_mat_rd.data[h].indices[i0:i1]))
                    )
                )
            ec[':'.join(ec_key)] += aln_mat_rd.count[cur_ind]

    ec = dict(ec)
    num_ecs = len(ec)

    logger.info('Constructing APM')
    logger.debug(f'Number Loci: {num_loci}')
    logger.debug(f'Number Haplotypes: {num_haplotypes}')
    logger.debug(f'Number ECs: {num_ecs}')

    aln_mat_ec = emase.AlignmentPropertyMatrix(
        shape=(num_loci, num_haplotypes, num_ecs)
    )
    aln_mat_ec.hname = names_haplotypes
    aln_mat_ec.lname = names_loci
    aln_mat_ec.count = np.zeros(num_ecs)

    logger.debug('Adding data to APM')
    for row_id, ec_key in enumerate(ec):
        aln_mat_ec.count[row_id] = ec[ec_key]
        nzlocs = ec_key.split(':')
        for h in range(aln_mat_ec.num_haplotypes):
            nzlocs_h = nzlocs[h]
            if nzlocs_h != '':
                nzinds = np.array(list(map(int, nzlocs_h.split(','))))
                aln_mat_ec.data[h][row_id, nzinds] = 1
    aln_mat_ec.finalize()

    logger.info(f'Saving EMASE Formatted File: {out_file}')
    aln_mat_ec.save(h5file=out_file, complib=comp_lib)
    logger.info('Done')


def stencil(
    alignment_file: Path | str,
    genotype_file: Path | str,
    group_file: Path | str = None,
    out_file: Path | str = None
) -> None:
    """
    Applying genotype calls to multi-way alignment incidence matrix.

    Args:
        alignment_file: alignment incidence file (h5)
        genotype_file: genotype calls by GBRS (tsv)
        group_file: gene ID to isoform ID mapping info (tsv)
        out_file: genotyped version of alignment incidence file (h5)
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

    # Load alignment incidence matrix ('alignment_file' is assumed to be in multiway transcriptome)
    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    # Load genotype calls
    logger.info(f'Loading and processing genotype calls from: {genotype_file}')
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

    logger.info(f'Saving EMASE Formatted File: {out_file}')
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
    Quantify expected read counts.

    Args:
        alignment_file: EMASE alignment incidence file (in hdf5 format)
        group_file: gene ID to isoform ID mapping info (tsv)
        length_file: transcript lengths (tsv)
        genotype_file: tab delimited file of locus and diplotype
        outbase: basename of all the generated output files
        multiread_model: EMASE model (default: 4)
        pseudocount: prior read count (default: 0.0)
        max_iters: maximum iterations for EM iteration
        tolerance: tolerance for EM termination (default: 0.0001 in TPM)
        report_alignment_counts: whether to report alignment counts
        report_posterior: whether to report posterior probabilities
    """
    if group_file is None:
        group_file = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if not os.path.exists(group_file):
            logger.warning('A group file is not given. Group-level results will not be reported.')

    if length_file is None:
        length_file = os.path.join(DATA_DIR, 'gbrs.hybridized.targets.info')
        if not os.path.exists(length_file):
            logger.warning('A length file is not given. Transcript length adjustment will *not* be performed.')

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
    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)

    # Load genotype calls
    if genotype_file is not None:
        # Genotype calls are at the gene level
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
        logger.info(f'Loading and processing genotype calls from: {genotype_file}')
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
    logger.info('Running EMASE')
    em_factory = emase.EMfactory(aln_mat)
    em_factory.prepare(pseudocount=pseudocount, lenfile=length_file)
    em_factory.run(
        model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True
    )

    logger.info(f'Generating isoform TPMs: {outbase}.isoforms.tpm')
    em_factory.report_depths(
        filename=f'{outbase}.isoforms.tpm', tpm=True, notes=gtcall_t
    )

    logger.info(f'Generating isoform Read Counts: {outbase}.isoforms.expected_read_counts')
    em_factory.report_read_counts(
        filename=f'{outbase}.isoforms.expected_read_counts', notes=gtcall_t
    )
    if report_posterior:
        logger.info(f'Generating Posterior Probabilities: {outbase}.posterior.h5')
        em_factory.export_posterior_probability(
            filename=f'{outbase}.posterior.h5'
        )
    if report_group_counts:
        logger.info(f'Generating gene TPMs: {outbase}.genes.tpm')
        em_factory.report_depths(
            filename=f'{outbase}.genes.tpm',
            tpm=True,
            grp_wise=True,
            notes=gtcall_g,
        )

        logger.info(f'Generating gene Read Counts: {outbase}.genes.expected_read_counts')
        em_factory.report_read_counts(
            filename=f'{outbase}.genes.expected_read_counts',
            grp_wise=True,
            notes=gtcall_g,
        )

    if report_alignment_counts:
        alnmat = emase.AlignmentPropertyMatrix(h5file=alignment_file, grpfile=group_file)

        logger.info(f'Generating isoform Alignment Counts: {outbase}.isoforms.alignment_counts')
        alnmat.report_alignment_counts(
            filename=f'{outbase}.isoforms.alignment_counts'
        )

        if report_group_counts:
            logger.info(f'Generating gene Alignment Counts: {outbase}.genes.alignment_counts')
            alnmat._bundle_inline(reset=True)
            alnmat.report_alignment_counts(
                filename=f'{outbase}.genes.alignment_counts'
            )
    logger.debug('Done')
