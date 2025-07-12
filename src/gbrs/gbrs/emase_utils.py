# standard library imports
import os
from collections import defaultdict
from itertools import dropwhile

# 3rd party library imports
import numpy as np

# local library imports
from gbrs import utils
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix
from gbrs.emase.EMfactory import EMfactory

DATA_DIR = os.getenv('GBRS_DATA', '.')
logger = utils.get_logger('gbrs')


def compress(
        emase_files: list[str],
        output_file: str,
        comp_lib: str = 'zlib'
) -> None:
    """
    Compress EMASE files by creating equivalence classes of identical alignment patterns.

    This function groups reads with identical alignment patterns across all haplotypes
    into equivalence classes (ECs), significantly reducing file size while preserving
    all alignment information. It is used for storage efficiency and downstream analysis.

    ALGORITHM:
    1. Load each EMASE file
    2. For each read in each file:
       - Extract alignment pattern across all haplotypes
       - Create a unique key representing the alignment pattern
       - Group reads with identical patterns into equivalence classes
    3. Build output matrix where each row represents one equivalence class
    4. Save compressed representation

    WHAT THIS DOES:
    - Takes one or more EMASE files (can be different samples or replicates)
    - Groups reads with IDENTICAL alignment patterns across all haplotypes
    - Creates equivalence classes where each EC represents a unique alignment pattern
    - Result: Compressed representation where identical patterns are merged
    - Output: EMASE file with fewer rows (ECs instead of individual reads)

    USE CASES:
    - Storage compression: Reduce file size by merging identical patterns
    - Downstream analysis: Work with equivalence classes instead of individual reads
    - Data merging: Combine multiple samples or replicates into single file
    - Memory efficiency: Reduce memory usage for large datasets

    DIFFERENCE FROM GET_COMMON_ALIGNMENTS:
    - This function: Groups identical alignment patterns (equivalence classes)
    - Get_common_alignments: Finds common alignments across files (intersection)
    - This is less restrictive - keeps all unique alignment patterns

    EQUIVALENCE CLASS DEFINITION:
    An equivalence class contains all reads that have the exact same alignment
    pattern across all haplotypes. For example, if reads A, B, and C all align
    to locus 1 in haplotype 0 and locus 5 in haplotype 1, they form one EC.

    Args:
        emase_files: List of EMASE files to compress. Files can have different reads.
        output_file: Name of the compressed EMASE file
        comp_lib: Compression library to use for output file

    Requirements:
        - All files must have same number of loci and haplotypes
        - Files can have different reads (unlike get_common_alignments)
        - Files can represent different samples or replicates

    Output Format:
        - Shape: (num_loci, num_haplotypes, num_equivalence_classes)
        - Each row represents one equivalence class
        - /count array contains number of reads in each EC
        - Sparse matrix format for memory efficiency
    """
    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    num_loci = None
    num_haplotypes = None
    names_loci = None
    names_haplotypes = None

    ec = defaultdict(int)
    for aln_file in emase_files:
        logger.info(f'Loading EMASE file: {aln_file}')
        aln_mat_rd = AlignmentPropertyMatrix(h5_file=aln_file)

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

        # Example dense matrix from file:
        #
        # Haplotype: A
        # 1 1 0 0 1
        # 0 0 1 0 0
        # 0 0 0 1 0
        # 0 0 0 1 0
        #
        # Haplotype: B
        # 0 0 0 0 0
        # 0 1 0 0 0
        # 0 0 0 1 0
        # 0 0 0 1 0
        #
        # Read 1, Haplotype A, ec_key='0,1,4'
        # Read 1, Haplotype B, ec_key=''
        # Read 1, ec_key = ['0,1,4','']
        # ec = {'0,1,4:': 1.0}
        #
        # Read 2, Haplotype A, ec_key='2'
        # Read 2, Haplotype B, ec_key='1'
        # Read 2, ec_key = ['2','1']
        # ec = {'0,1,4:': 1.0, '2:1': 1.0}
        #
        # Read 3, Haplotype A, ec_key='3'
        # Read 3, Haplotype B, ec_key='3'
        # Read 3, ec_key = ['3','3']
        # ec = {'0,1,4:': 1.0, '2:1': 1.0, '3:3': 1.0}
        #
        # Read 4, Haplotype A, ec_key='3'
        # Read 4, Haplotype B, ec_key='3'
        # Read 4, ec_key = ['3','3']
        # ec = {'0,1,4:': 1.0, '2:1': 1.0, '3:3': 2.0}
        #

        logger.debug('Creating unique ECs')
        for cur_ind in range(aln_mat_rd.num_reads):
            ec_key = []
            for h in range(aln_mat_rd.num_haplotypes):
                i0 = aln_mat_rd.data[h].indptr[cur_ind]
                i1 = aln_mat_rd.data[h].indptr[cur_ind + 1]
                # logger.debug(f'Read {cur_ind}, haplotype {h}: sparse indices {i0}-{i1}')

                # ec_key is a comma separated list of read indices, one entry per haplotype
                ec_key.append(
                    ','.join(
                        map(str, sorted(aln_mat_rd.data[h].indices[i0:i1]))
                    )
                )
                # logger.debug(f'{ec_key=}')
            # ec is a dictionary of equivalence classes with the value being
            # the number of occurrences
            ec[':'.join(ec_key)] += aln_mat_rd.count[cur_ind]
            # logger.debug(f'{ec}')

    ec = dict(ec)
    num_ecs = len(ec)

    logger.info('Constructing APM')
    logger.debug(f'Number Loci: {num_loci}')
    logger.debug(f'Number Haplotypes: {num_haplotypes}')
    logger.debug(f'Number ECs: {num_ecs}')

    aln_mat_ec = AlignmentPropertyMatrix(shape=(num_loci, num_haplotypes, num_ecs))
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

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    aln_mat_ec.save(h5_file=output_file, complib=comp_lib)
    logger.info('Done')


def compress_optimized(
        emase_files: list[str],
        output_file: str,
        comp_lib: str = 'zlib'
) -> None:
    """
    Compress EMASE files by creating equivalence classes of identical alignment patterns.

    This function is identical to the original compress function but uses tuple-based
    EC keys instead of string operations for better performance.

    Args:
        emase_files: List of EMASE files to compress. Files can have different reads.
        output_file: Name of the compressed EMASE file
        comp_lib: Compression library to use for output file
    """
    for x in emase_files:
        logger.info(f'EMASE file: {x}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    num_loci = None
    num_haplotypes = None
    names_loci = None
    names_haplotypes = None

    # Use tuple-based EC dictionary instead of string-based
    ec = defaultdict(int)
    for aln_file in emase_files:
        logger.info(f'Loading EMASE file: {aln_file}')
        aln_mat_rd = AlignmentPropertyMatrix(h5_file=aln_file)

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

        # Example dense matrix from file:
        #
        # Haplotype: A
        # 1 1 0 0 1
        # 0 0 1 0 0
        # 0 0 0 1 0
        # 0 0 0 1 0
        #
        # Haplotype: B
        # 0 0 0 0 0
        # 0 1 0 0 0
        # 0 0 0 1 0
        # 0 0 0 1 0
        #
        # Read 1, Haplotype A, ec_key=(0,1,4)
        # Read 1, Haplotype B, ec_key=()
        # Read 1, ec_key = [(0,1,4),()]
        # ec = {((0, 1, 4), ()): 1.0})
        #
        # Read 2, Haplotype A, ec_key=(2)
        # Read 2, Haplotype B, ec_key=(1)
        # Read 2, ec_key = [(2, ), (1, )]
        # ec = {((0, 1, 4), ()): 1.0, ((2, ), (1, )): 1.0})
        #
        # Read 3, Haplotype A, ec_key=(3)
        # Read 3, Haplotype B, ec_key=(3)
        # Read 3, ec_key = [(3, ), (3, )]
        # ec = {((0, 1, 4), ()): 1.0, ((2, ), (1, )): 1.0, ((3, ), (3, )): 1.0})
        #
        # Read 4, Haplotype A, ec_key=(3)
        # Read 4, Haplotype B, ec_key=(3)
        # Read 4, ec_key = [(3, ), (3, )]
        # ec = {((0, 1, 4), ()): 1.0, ((2, ), (1, )): 1.0, ((3, ), (3, )): 2.0})
        #

        logger.debug('Creating unique ECs with tuple keys')
        for cur_ind in range(aln_mat_rd.num_reads):
            ec_key_parts = []
            for h in range(aln_mat_rd.num_haplotypes):
                i0 = aln_mat_rd.data[h].indptr[cur_ind]
                i1 = aln_mat_rd.data[h].indptr[cur_ind + 1]
                logger.debug(f'Read {cur_ind}, haplotype {h}: sparse indices {i0}-{i1}')
                indices = aln_mat_rd.data[h].indices[i0:i1]

                # use tuple instead of string for better performance
                # ec_key_parts a list of tuples (one entry per haplotype)
                # each tuples values are read indices
                ec_key_parts.append(tuple(sorted(indices)))
                logger.debug(f'{ec_key_parts=}')

            # create tuple key instead of string
            ec_key = tuple(ec_key_parts)
            # ec is a dictionary of equivalence classes with the value being
            # the number of occurrences
            ec[ec_key] += aln_mat_rd.count[cur_ind]
            logger.debug(f'{ec=}')

    logger.debug('ec conversion')
    ec = dict(ec)
    logger.debug('ec conversion done')
    num_ecs = len(ec)

    logger.info('Constructing APM')
    logger.debug(f'Number Loci: {num_loci}')
    logger.debug(f'Number Haplotypes: {num_haplotypes}')
    logger.debug(f'Number ECs: {num_ecs}')

    aln_mat_ec = AlignmentPropertyMatrix(shape=(num_loci, num_haplotypes, num_ecs))
    aln_mat_ec.hname = names_haplotypes
    aln_mat_ec.lname = names_loci
    aln_mat_ec.count = np.zeros(num_ecs)

    logger.debug('Adding data to APM')
    for row_id, (ec_key, count) in enumerate(ec.items()):
        aln_mat_ec.count[row_id] = count

        # Process tuple-based EC key
        for h in range(aln_mat_ec.num_haplotypes):
            indices = ec_key[h]
            if len(indices) > 0:
                # Convert tuple back to array for indexing
                aln_mat_ec.data[h][row_id, list(indices)] = 1

    aln_mat_ec.finalize()

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    aln_mat_ec.save(h5_file=output_file, complib=comp_lib)
    logger.info('Done')


def stencil(
        alignment_file: str,
        genotype_file: str,
        group_file: str = None,
        output_file: str = None
) -> None:
    """
    This function transforms multi-way alignment data (containing all founder
    haplotypes) into a diploid representation based on individual genotype
    calls. It filters alignment data to retain only the haplotypes that are
    present in the individual's genotype.

    - Load EMASE alignment matrix and optional group information
    - Load genotype calls from GBRS analysis
    - Create genotype mask matrix based on called haplotypes
    - Apply mask to alignment matrix (element-wise multiplication)
    - Remove zero entries to maintain sparse matrix efficiency
    - Save filtered alignment matrix

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing
            multi-way alignment data.

        genotype_file: Path to the genotype calls file generated by GBRS
            analysis.
            Format: tab-separated, gene ID than called haplotypes

        group_file: Path to the group file containing transcript-to-gene
            mapping.
            Format: tab-separated transcript ID than gene ID.

        output_file: Path for the output stenciled EMASE file, which will
            contain diploid alignment data filtered by genotype calls.

    Raises:
        FileNotFoundError: If alignment_file, genotype_file, or group_file
            does not exist.
        ValueError: If genotype file format is invalid or incompatible with
            alignment data.
        RuntimeError: If the AlignmentPropertyMatrix cannot be loaded or
            processed.
    """
    if group_file is None:
        group_file = os.path.join(DATA_DIR, 'ref.gene2transcripts.tsv')
        if not os.path.exists(group_file):
            logger.info('A group file is *not* given. Genotype will be stenciled as is.')

    if output_file is None:
        output_file = f'gbrs.stenciled.{os.path.basename(alignment_file)}'

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Genotype File: {genotype_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Output File: {output_file}')

    # load alignment incidence matrix
    logger.info(f'Loading EMASE file: {alignment_file}')
    apm = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)
    logger.debug(f'Number Loci: {apm.num_loci}')
    logger.debug(f'Number Haplotypes: {apm.num_haplotypes}')
    logger.debug(f'Number Reads: {apm.num_reads}')

    # load genotype calls
    logger.info(f'Loading and processing genotype calls from: {genotype_file}')
    hid = dict(zip(apm.hname, np.arange(apm.num_haplotypes)))
    gid = dict(zip(apm.gname, np.arange(len(apm.gname))))
    gtmask = np.zeros((apm.num_haplotypes, apm.num_loci))
    gtcall_g = dict.fromkeys(apm.gname)

    with open(genotype_file) as fh:
        if group_file is not None:
            gtcall_t = dict.fromkeys(apm.lname)
            for line in dropwhile(utils.is_comment, fh):
                item = line.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(apm.groups[gid[g]])
                gtmask[np.meshgrid(hid2set, tid2set)] = 1.0
                for t in tid2set:
                    gtcall_t[apm.lname[t]] = gt
        else:
            for line in dropwhile(utils.is_comment, fh):
                item = line.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                gtmask[np.meshgrid(hid2set, gid[g])] = 1.0

    apm.multiply(gtmask, axis=2)
    for h in range(apm.num_haplotypes):
        apm.data[h].eliminate_zeros()

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    apm.save(h5_file=output_file)
    logger.info('Done')


def quantify(
        alignment_file: str,
        group_file: str = None,
        length_file: str = None,
        genotype_file: str = None,
        outbase: str = 'gbrs.quantified',
        multiread_model: int = 4,
        pseudocount: float = 0.0,
        max_iters: int = 999,
        tolerance: float = 0.0001,
        report_alignment_counts: bool = False,
        report_posterior: bool = False
) -> None:
    """
    Quantify gene expression using EMASE algorithm with optional genotype
    filtering.

    This function is the core quantification component of the GBRS pipeline
    that performs expression analysis on EMASE alignment data. It can operate
    in two modes:

    - Multi-way mode: Analyzes all founder haplotypes simultaneously
    - Diploid mode: Filters data by genotype calls before analysis

    The function implements the EMASE algorithm to estimate transcript
    abundances while accounting for multi-mapping reads, transcript length
    biases, and allele-specific expression differences. It can generate
    both transcript-level and gene-level estimates.

    ALGORITHM:
    - Load EMASE alignment data and optional group/length information
    - If genotype file provided:
       - Load genotype calls and create diploid filter
       - Apply genotype mask to alignment matrix
       - Set output mode to diploid
    - Initialize EMASE factory with filtered alignment data
    - Run expectation-maximization algorithm
    - Generate transcript and gene-level abundance estimates
    - Optionally generate additional reports (alignment counts, posteriors)

    MODES OF OPERATION:
    - Multi-way mode (genotype_file=None): Analyzes all founder haplotypes
      Output files: {outbase}.multiway.*
    - Diploid mode (genotype_file provided): Filters by individual genotype
      Output files: {outbase}.diploid.*

    OUTPUT FILES:
    Transcript-level TPM estimates
    - {outbase}.{multiway|diploid}.isoforms.tpm
    Transcript-level read counts
    - {outbase}.{multiway|diploid}.isoforms.expected_read_counts
    Raw alignment counts (if report_alignment_counts=True)
    - {outbase}.{multiway|diploid}.isoforms.alignment_counts
    Gene-level TPM estimates (if group_file provided)
    - {outbase}.{multiway|diploid}.genes.tpm
    Gene-level read counts (if group_file provided)
    - {outbase}.{multiway|diploid}.genes.expected_read_counts
    Gene-level alignment counts (if report_alignment_counts=True)
    - {outbase}.{multiway|diploid}.genes.alignment_counts
    Posterior probabilities (if report_posterior=True)
    - {outbase}.{multiway|diploid}.posterior.h5

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing
            alignment data.

        group_file: Path to the group file containing transcript-to-gene
            mapping. Uses default location from GBRS_DATA environment
            variable. If provided, enables gene-level analysis in addition
            to transcript-level analysis.
            Format: tab-separated, transcript ID than gene ID

        length_file: Path to the transcript length file for length bias
            correction.  Uses default location from GBRS_DATA environment
            variable. Contains transcript IDs and their lengths for
            normalization.
            Format: tab-separated, transcript ID than length

        genotype_file: Path to the genotype calls file for diploid mode
            analysis.  If None, runs in multi-way mode. If provided,
            filters alignment data by individual genotype calls before
            analysis.
            Format: tab-separated, gene ID than called haplotypes

        outbase: Base name for all output files.
        multiread_model: EMASE model to use for multi-mapping read handling.
        pseudocount: Prior pseudocount for regularization.
        max_iters: Maximum number of EM iterations.
        tolerance: Convergence tolerance for EM algorithm.
        report_alignment_counts: Whether to report raw alignment counts.
        report_posterior: Whether to report posterior probability estimates.

    Raises:
        FileNotFoundError: If any input file does not exist.
        ValueError: If the EMASE file format is invalid or parameters are
            invalid.
        RuntimeError: If the EM algorithm fails to converge or other processing
        errors.

    Notes:
        - The function automatically detects and uses default files from GBRS_DATA
          environment variable if not explicitly provided.
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
    report_group_counts = (group_file is not None)

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
    apm = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)

    # load genotype calls
    if genotype_file is not None:
        # genotype calls are at the gene level
        outbase = f'{outbase}.diploid'
        logger.debug(f'Outbase now: {outbase}')
        # haplotype as key and number as value
        hid = dict(zip(apm.hname, np.arange(apm.num_haplotypes)))
        # gene id as key and number as value
        gid = dict(zip(apm.gname, np.arange(len(apm.gname))))
        gtmask = np.zeros((apm.num_haplotypes, apm.num_loci))
        # genome genotype calls, gene_id as key
        gtcall_g = dict.fromkeys(apm.gname)
        # transcript genotype calls, transcript_id as key
        gtcall_t = dict.fromkeys(apm.lname)

        logger.info(f'Loading and processing genotype calls from: {genotype_file}')
        with open(genotype_file) as fh:
            for line in dropwhile(utils.is_comment, fh):
                item = line.rstrip().split('\t')
                g, gt = item[:2]
                gtcall_g[g] = gt
                hid2set = np.array([hid[c] for c in gt])
                tid2set = np.array(apm.groups[gid[g]])
                gtmask[tuple(np.meshgrid(hid2set, tid2set))] = 1.0
                for t in tid2set:
                    gtcall_t[apm.lname[t]] = gt

        apm.multiply(gtmask, axis=2)
        for h in range(apm.num_haplotypes):
            apm.data[h].eliminate_zeros()
    else:
        outbase = f'{outbase}.multiway'
        logger.debug(f'Outbase now: {outbase}')
        gtcall_g = None
        gtcall_t = None

    # run EMASE
    logger.info('Running EMASE')
    em_factory = EMfactory(apm)
    em_factory.prepare(pseudocount=pseudocount, lenfile=length_file)

    em_factory.run(
        model=multiread_model,
        tol=tolerance,
        max_iters=max_iters,
        verbose=True
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
        apm_counts = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)

        logger.info(f'Generating isoform Alignment Counts: {outbase}.isoforms.alignment_counts')
        apm_counts.report_alignment_counts(
            filename=f'{outbase}.isoforms.alignment_counts'
        )

        if report_group_counts:
            logger.info(f'Generating gene Alignment Counts: {outbase}.genes.alignment_counts')
            apm_counts._bundle_inline(reset=True)
            apm_counts.report_alignment_counts(
                filename=f'{outbase}.genes.alignment_counts'
            )

    logger.debug('Done')
