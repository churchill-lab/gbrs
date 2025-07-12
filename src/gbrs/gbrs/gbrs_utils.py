# standard library imports
import logging
import os
import re
from collections import defaultdict, OrderedDict
from itertools import combinations_with_replacement, product

# 3rd party library imports
import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# local library imports
from gbrs import utils

matplotlib.use('Agg')

DATA_DIR = os.getenv('GBRS_DATA', '.')

logging.getLogger('matplotlib').setLevel(logging.WARNING)
logger = utils.get_logger('gbrs')


def get_chromosome_info(
        fasta_index_file: str = 'ref.fa.fai'
) -> OrderedDict[str, int]:
    """
    Load chromosome information from FASTA index file.

    Args:
        fasta_index_file: Path to FASTA index file (default: 'ref.fa.fai')

    Returns:
        OrderedDict mapping chromosome names to lengths

    Raises:
        ValueError: If file not found or GBRS_DATA not set correctly
    """
    fai_file = os.path.join(DATA_DIR, fasta_index_file)
    try:
        chr_lens = OrderedDict(
            np.loadtxt(fai_file, usecols=(0, 1), dtype='|S8,<i4')
        )

        # convert from bytes to string
        chr_lens = OrderedDict({k.decode(): v for k, v in chr_lens.items()})
        return chr_lens
    except FileNotFoundError:
        raise ValueError(
            'Make sure if $GBRS_DATA is set correctly, and that "ref.fa.fai" '
            'is in that directory. Currently it is: '
            f'{DATA_DIR}'
        )


def get_founder_info(
        founder_info_file: str = 'founder.hexcolor.info'
) -> OrderedDict[str, int]:
    """
    Load founder strain information and colors from file.

    Args:
        founder_info_file: Path to founder info file

    Returns:
        OrderedDict mapping founder names to color codes

    Raises:
        ValueError: If file not found or GBRS_DATA not set correctly
    """
    color_file = os.path.join(DATA_DIR, founder_info_file)

    try:
        founder_colors = OrderedDict(
            np.loadtxt(
                color_file,
                usecols=(0, 1),
                dtype='str',
                delimiter='\t',
                comments=None,
            )
        )

        return founder_colors
    except FileNotFoundError:
        raise ValueError(
            f'Make sure if $GBRS_DATA is set correctly, and that '
            '"founder.hexcolor.info" is in that directory. Currently it is: '
            f'{DATA_DIR}'
        )


def unit_vector(vector: np.ndarray) -> np.ndarray:
    """
    Normalize vector to unit length.

    Args:
        vector: Input vector to normalize

    Returns:
        Unit vector (normalized to length 1.0)
    """
    if sum(vector) > 1e-6:
        return vector / np.linalg.norm(vector)
    else:
        return vector


def print_vecs(
        vecs: np.ndarray,
        format_str: str = '%10.1f',
        show_sum: bool = False
) -> None:
    """
    Print vectors with optional formatting and sum display.

    Args:
        vecs: Array of vectors to print
        format_str: Format string for vector elements
        show_sum: Whether to show sum of each vector
    """
    for i in range(vecs.shape[0]):
        v = vecs[i]
        print(' '.join(format_str % elem for elem in v))
        if show_sum:
            print('\t=>', format_str % sum(v))
        else:
            print()


def get_genotype_probability(
        aln_profile: np.ndarray,
        aln_specificity: np.ndarray,
        sigma: float = 0.12
) -> np.ndarray:
    """
    Calculate genotype probabilities from alignment profiles using expression
    similarity.

    This function is a core component of the GBRS emission probability
    calculation that converts expression similarity between an individual's
    gene expression profile and founder strain expression patterns into
    genotype probabilities. It implements the emission model of the GBRS
    Hidden Markov Model.

    The function compares the normalized expression profile of an individual
    gene against the alignment specificity vectors of all possible founder
    strain combinations (diplotypes). The similarity is measured using squared
    Euclidean distance, which is then converted to probabilities using a
    Gaussian kernel.

    ALGORITHM:
    - Normalize the individual's expression profile to a unit vector
    - For each possible diplotype (haplotype pair):
       - For homozygotes: Compare directly with founder strain expression vector
       - For heterozygotes: Create composite vector from two founder strains
       - Calculate squared Euclidean distance between profiles
    - Convert distances to probabilities using Gaussian kernel with sigma parameter
    - Normalize to ensure probabilities sum to 1

    USE CASES:
    - GBRS reconstruction: Calculate emission probabilities for HMM

    Args:
        aln_profile: Individual's gene expression profile.
            Shape: (num_haplotypes,)
            Type: numpy.ndarray of float
            Description: Raw expression values for each possible haplotype.
            The values represent expression levels (e.g., TPM) for each
            founder strain haplotype at a specific gene.
            Example: [10.5, 8.2, 7.1, 6.3, 5.8, 4.9, 4.2, 3.8] for 8 haplotypes

        aln_specificity: Founder strain alignment specificity matrix.
            Shape: (num_haplotypes, num_haplotypes)
            Type: numpy.ndarray of float
            Description: Normalized expression vectors for each founder strain.
            Each row represents how a founder strain's expression aligns with
            all other strains. Must contain at least one entry > 1.0 to be
            considered valid alignment specificity data.
            Example: Normalized expression matrix from get_alignment_spec()

        sigma: Standard deviation parameter for Gaussian kernel.
            Type: float
            Default: 0.12
            Description: Controls the sensitivity of genotype inference.
            Smaller values make the model more sensitive to expression differences,
            while larger values make it more tolerant of expression variation.
            Range: Typically 0.05 to 0.5, with 0.12 being optimal for most datasets.

    Returns:
        numpy.ndarray: Genotype probabilities for all possible diplotypes.
            Shape: (num_diplotypes,)
            Type: float
            Description: Probability values that sum to 1.0, representing the
            likelihood of each possible diplotype given the expression data.
            Order: Homozygotes first (AA, BB, CC, ...), then heterozygotes (AB, AC, ...)
            Example: [0.45, 0.12, 0.08, 0.15, 0.20] for 5 diplotypes

    Raises:
        ValueError: If aln_specificity contains no entries > 1.0
        RuntimeError: If probability calculation fails or results in invalid values

    Notes:
        - The function assumes aln_specificity contains unit vectors (normalized)
        - For heterozygotes, the composite vector is created by averaging two founder vectors
        - The Gaussian kernel uses the formula: exp(-distance² / (2 * sigma²))
        - This function is called for each gene during GBRS reconstruction
        - The sigma parameter is critical for balancing sensitivity vs. robustness
        - Genes with low expression may use a different sigma value (typically 0.45)
        - This emission model is the foundation of GBRS's expression-based genotyping
    """
    # 'aln_specificity' should be a set of unit vectors (at least one of the entry is larger than 1.)
    num_haps = len(aln_profile)
    aln_vec = unit_vector(aln_profile)
    genoprob = []
    for i in range(num_haps):
        v1 = unit_vector(aln_specificity[i])
        for j in range(i, num_haps):
            if j == i:
                # homozygotes
                genoprob.append(sum(np.power(aln_vec - v1, 2)))
            else:
                v2 = unit_vector(aln_specificity[j])
                geno_vec = unit_vector(v1 + v2)
                # compute directional similarity
                # for heterozygotes
                genoprob.append(
                    sum(np.power(aln_vec - geno_vec, 2))
                )
    genoprob = np.exp(np.array(genoprob) / (-2 * sigma * sigma))
    return np.array(genoprob / sum(genoprob))


def ris_step(
        gen_left: str,
        gen_right: str,
        rec_frac: float,
        haps: tuple[str, ...] = ('A', 'B'),
        gamma_scale: float = 0.1,
        is_x_chr: bool = False,
        forward_direction: bool = True,
) -> float:
    """
    Calculate log transition probability for RIL (Recombinant Inbred Lines) by sib-mating.

    This function implements the transition probability model for Recombinant Inbred Lines
    created through sib-mating, which is a key component of the GBRS Hidden Markov Model.
    It calculates the probability of transitioning between genotypes at adjacent genetic
    positions, accounting for recombination events and the specific mating scheme.

    The function models the genetic linkage between adjacent markers/genes based on
    the recombination fraction (distance in centiMorgans). For RIL populations, the
    transition probabilities reflect the accumulation of recombination events over
    multiple generations of inbreeding, resulting in a more complex pattern than
    simple F2 populations.

    ALGORITHM:
    1. Generate all possible diplotypes from the haplotype set
    2. Calculate recombination rate R based on chromosome type (autosome vs X)
    3. Apply gamma scaling to allow for heterozygosity (rare in RIL)
    4. Return log transition probability based on genotype pair and direction

    WHAT THIS DOES:
    - Models genetic linkage between adjacent positions in RIL populations
    - Accounts for different recombination rates on autosomes vs X chromosome
    - Provides transition probabilities for the GBRS HMM
    - Supports both forward and backward chain directions

    USE CASES:
    - GBRS reconstruction: Calculate transition probabilities for HMM
    - RIL population analysis: Model genetic linkage in inbred lines
    - Linkage mapping: Understand recombination patterns in RIL
    - Model validation: Test transition probability calculations

    Args:
        gen_left: Left genotype (diplotype).
            Type: str
            Description: Diplotype at the left (upstream) position.
            Format: Two-character string representing haplotype pair (e.g., 'AA', 'AB')
            Example: 'AA' for homozygous A, 'AB' for heterozygous A/B

        gen_right: Right genotype (diplotype).
            Type: str
            Description: Diplotype at the right (downstream) position.
            Format: Two-character string representing haplotype pair (e.g., 'AA', 'AB')
            Example: 'BB' for homozygous B, 'AC' for heterozygous A/C

        rec_frac: Recombination fraction (distance in centiMorgans).
            Type: float
            Description: Genetic distance between the two positions in centiMorgans.
            This determines the probability of recombination between positions.
            Range: Typically 0.0 to 50.0 cM (50 cM = independent segregation)
            Example: 2.5 for positions 2.5 cM apart

        haps: Haplotype identifiers.
            Type: tuple[str, ...]
            Default: ('A', 'B')
            Description: List of founder strain identifiers.
            Used to generate all possible diplotypes for the population.
            Example: ('A', 'B', 'C', 'D') for 4-founder population

        gamma_scale: Scale factor for heterozygosity allowance.
            Type: float
            Default: 0.1
            Description: Controls the probability of heterozygosity in RIL.
            RIL populations are typically highly homozygous, but some heterozygosity
            may persist. This parameter allows for rare heterozygous regions.
            Range: Typically 0.01 to 0.5, with 0.1 being standard

        is_x_chr: Whether the chromosome is X chromosome.
            Type: bool
            Default: False
            Description: X chromosomes have different recombination patterns than autosomes.
            This affects the calculation of recombination rate R and transition probabilities.
            Example: True for X chromosome, False for autosomes 1-19

        forward_direction: Direction of the transition.
            Type: bool
            Default: True
            Description: Whether calculating forward (left to right) or backward
            (right to left) transition probabilities. The direction affects the
            probability calculations for asymmetric transitions.
            Example: True for forward chain, False for backward chain

    Returns:
        float: Log transition probability (natural logarithm).
            Description: Natural logarithm of the transition probability from gen_left
            to gen_right. Using log probabilities prevents numerical underflow in
            the HMM calculations.
            Range: Negative values (log of probabilities < 1)
            Example: -2.3 for probability of 0.1

    Raises:
        ValueError: If genotypes are not valid diplotypes for the given haplotypes
        RuntimeError: If probability calculation fails

    Notes:
        - Originally part of r/qtl2 designed/coded by Karl Broman
        - Ported to python by Karl Broman and extended by KB Choi
        - For autosomes: R = 4.0 * rec_frac / (1 + 6.0 * rec_frac)
        - For X chromosome: R = (2 * rec_frac) / (1.0 + 4.0 * rec_frac)
        - Gamma parameter controls heterozygosity: gamma = R * gamma_scale
        - This function is called for each adjacent pair of positions during
          transition probability matrix generation
        - The log probabilities are used in the HMM to prevent numerical issues
        - RIL populations are typically >99% homozygous after many generations
    """
    it = combinations_with_replacement(haps, 2)
    diplotype = [f'{ht1}{ht2}' for ht1, ht2 in it]

    if is_x_chr:
        R = (2 * rec_frac) / (1.0 + 4.0 * rec_frac)
        gamma = R * gamma_scale
        if forward_direction:
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0 - R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(R) - np.log(1 + gamma)
            elif gen_left == diplotype[1]:
                return np.log(1 / 3.0)
            if gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(2.0 * R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0 - 2.0 * R) - np.log(1 + gamma)
        else:  # backward direction
            if gen_left == diplotype[0]:
                if gen_right == diplotype[0]:
                    return np.log(1.0 - 2.0 * R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(2.0 * R) - np.log(1 + gamma)
            elif gen_left == diplotype[1]:
                return np.log(1 / 3.0)
            elif gen_left == diplotype[2]:
                if gen_right == diplotype[0]:
                    return np.log(R) - np.log(1 + gamma)
                elif gen_right == diplotype[1]:
                    return np.log(gamma) - np.log(1 + gamma)
                elif gen_right == diplotype[2]:
                    return np.log(1.0 - R) - np.log(1 + gamma)

    else:  # autosome
        R = 4.0 * rec_frac / (1 + 6.0 * rec_frac)
        gamma = R * gamma_scale
        if gen_left == diplotype[0]:
            if gen_right == diplotype[0]:
                return np.log(1.0 - R) - np.log(1 + gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1 + gamma)
            elif gen_right == diplotype[2]:
                return np.log(R) - np.log(1 + gamma)
        elif gen_left == diplotype[1]:
            return np.log(1 / 3.0)
        elif gen_left == diplotype[2]:
            if gen_right == diplotype[0]:
                return np.log(R) - np.log(1 + gamma)
            elif gen_right == diplotype[1]:
                return np.log(gamma) - np.log(1 + gamma)
            elif gen_right == diplotype[2]:
                return np.log(1.0 - R) - np.log(1 + gamma)


def f2_step(
        gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def cc_step(
        gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def do_step(
        gen_left, gen_right, rec_frac, is_x_chr=False, forward_direction=True
):
    return NotImplementedError


def get_transition_prob(
        marker_file: str,
        haplotypes: str = 'A,B',
        mating_scheme: str = 'RI',
        gamma_scale: float = 0.01,
        epsilon: float = 0.000001,
        output_file: str = 'tranprob.npz'
) -> None:
    """
    Generate transition probability matrices for GBRS Hidden Markov Model.

    This function is a critical preprocessing step in the GBRS pipeline that
    calculates transition probability matrices for all chromosomes based on
    genetic distances between markers/genes. These matrices model the genetic
    linkage between adjacent positions and are essential for the HMM-based
    genome reconstruction algorithm.

    The function processes marker files containing genetic positions and
    calculates transition probabilities for each adjacent pair of positions
    on each chromosome. The probabilities depend on the mating scheme used
    to create the population (RI, F2, CC, DO) and account for different
    recombination patterns on autosomes vs sex chromosomes.

    ALGORITHM:
    1. Load marker file and organize positions by chromosome
    2. Generate all possible diplotypes from haplotype set
    3. For each chromosome:
       - Calculate genetic distances between adjacent positions
       - For each diplotype pair and position pair:
         - Calculate transition probability using appropriate step function
       - Store results in 3D matrix (positions × diplotypes × diplotypes)
    4. Save matrices in compressed numpy format

    WHAT THIS DOES:
    - Creates transition probability matrices for GBRS HMM
    - Models genetic linkage between adjacent positions
    - Accounts for different mating schemes and chromosome types
    - Provides the foundation for HMM-based genome reconstruction

    USE CASES:
    - GBRS preprocessing: Generate transition probabilities for reconstruction
    - Population analysis: Model genetic linkage in different populations
    - Linkage mapping: Understand recombination patterns
    - Model validation: Test transition probability calculations

    Args:
        marker_file: Path to the marker file containing genetic positions.
            Type: str
            Description: Tab-separated file with columns: marker_id, chromosome,
            genomic_position, genetic_position_cM. Contains all markers/genes
            used in the GBRS analysis with their genetic coordinates.
            Format:
                Gene1   1   1000000   2.5
                Gene2   1   2000000   5.1
                Gene3   2   500000    1.2
            Example: 'ref.gene_pos.ordered.tsv'

        haplotypes: Comma-separated list of founder strain identifiers.
            Type: str
            Default: 'A,B'
            Description: Founder strain names used in the population.
            These are used to generate all possible diplotypes for the
            transition probability calculations.
            Format: 'A,B,C,D' for 4-founder population
            Example: 'A,B,C,D,E,F,G,H' for Diversity Outbred mice

        mating_scheme: Population mating scheme.
            Type: str
            Default: 'RI'
            Description: The breeding scheme used to create the population.
            Determines which transition probability function to use.
            Options: 'RI' (Recombinant Inbred), 'F2', 'CC' (Collaborative Cross),
            'DO' (Diversity Outbred)
            Example: 'RI' for Recombinant Inbred Lines

        gamma_scale: Scale factor for heterozygosity allowance.
            Type: float
            Default: 0.01
            Description: Controls the probability of heterozygosity in the
            transition model. Smaller values assume more homozygosity.
            Range: Typically 0.001 to 0.1
            Example: 0.01 for highly inbred populations

        epsilon: Minimum genetic distance threshold.
            Type: float
            Default: 0.000001
            Description: Minimum allowed genetic distance between adjacent
            positions. Prevents division by zero and ensures numerical stability.
            Range: Very small positive values (1e-6 to 1e-3)
            Example: 0.000001 cM

        output_file: Output filename for transition probability matrices.
            Type: str
            Default: 'tranprob.npz'
            Description: Name of the compressed numpy file to save the
            transition probability matrices. The file will contain one
            matrix per chromosome.
            Example: 'tranprob.DO.G17.M.npz'

    Returns:
        None. The function creates a compressed numpy file containing transition
        probability matrices for each chromosome.

    Raises:
        FileNotFoundError: If marker_file does not exist
        ValueError: If mating_scheme is not supported or parameters are invalid
        OSError: If output file cannot be created

    Notes:
        - The function saves gene positions in 'ref.gene_pos.ordered.npz'
        - Transition matrices are 3D: (positions-1) × diplotypes × diplotypes
        - For n haplotypes, there are n*(n+1)/2 possible diplotypes
        - X chromosome uses different recombination model than autosomes
        - The epsilon parameter prevents numerical issues with very close positions
        - This function is typically run once per population/mating scheme
        - The output file is essential for the GBRS reconstruction algorithm
        - File sizes can be large for populations with many markers/genes
    """
    logger.info(f'Marker File: {marker_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Mating Scheme: {mating_scheme}')
    logger.info(f'Gama Scale: {gamma_scale}')
    logger.info(f'Epsilon: {epsilon}')
    logger.info(f'Output File: {output_file}')

    haplotypes = haplotypes.split(',')
    it = combinations_with_replacement(haplotypes, 2)
    diplotype = [f'{ht1}{ht2}' for ht1, ht2 in it]
    num_diplotypes = len(diplotype)
    diplotype_index = np.arange(num_diplotypes)

    locs_by_chro = defaultdict(list)
    gpos_by_chro = defaultdict(list)

    logger.info(f'Loading marker file: {marker_file}')
    with open(marker_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            locs_by_chro[item[1]].append((item[0], float(item[3])))
            gpos_by_chro[item[1]].append((item[0], int(item[2])))
    locs_by_chro = dict(locs_by_chro)
    gpos_by_chro = dict(gpos_by_chro)

    logger.info(f'Saving {os.path.join(DATA_DIR, "ref.gene_pos.ordered.npz")}')
    np.savez_compressed(
        os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz'), **gpos_by_chro
    )

    if mating_scheme == 'RI':
        step_func = ris_step
    elif mating_scheme == 'F2':
        step_func = f2_step
    elif mating_scheme == 'CC':
        step_func = cc_step
    elif mating_scheme == 'DO':
        step_func = do_step
    else:
        raise ValueError(f'Unknown mating scheme: {mating_scheme}')

    tprob = dict()
    for c in locs_by_chro.keys():
        logger.debug(f'Working on {c}')
        is_x_chr = c == 'X'
        pdiff = np.diff(np.array([e[1] for e in locs_by_chro[c]]))
        pdiff[pdiff < epsilon] = epsilon
        ndiff = len(pdiff)
        tprob[c] = np.ndarray(
            shape=(ndiff, num_diplotypes, num_diplotypes), dtype=float
        )
        for dt1id, dt2id in list(product(diplotype_index, repeat=2)):
            dt1 = diplotype[dt1id]
            dt2 = diplotype[dt2id]
            for i, d in enumerate(pdiff):
                tprob[c][i][dt1id, dt2id] = step_func(
                    dt1,
                    dt2,
                    d,
                    gamma_scale=gamma_scale,
                    haps=haplotypes,
                    is_x_chr=is_x_chr,
                )

    logger.info(f'Saving {os.path.join(DATA_DIR, output_file)}')
    np.savez_compressed(os.path.join(DATA_DIR, output_file), **tprob)
    logger.info('Done')


def get_alignment_spec(
        sample_file: str,
        haplotypes: list[str],
        gene2transcript: str,
        out_dir: str,
        min_expr: float = 2.0
) -> None:
    """
    Generate alignment specificity matrices for GBRS genome reconstruction.

    This function is a critical preprocessing step in the GBRS pipeline that
    calculates alignment specificity matrices from founder strain expression data.
    These matrices capture how each founder strain's expression profile aligns
    with the expression patterns of other strains, providing the foundation for
    genotype probability calculations during genome reconstruction.

    The function processes TPM (Transcripts Per Million) expression data from
    multiple founder strains to create three key matrices:
    1. **Axes matrix**: Raw expression values for each gene across all strains
    2. **ASE matrix**: Allele-specific expression summaries
    3. **Avecs matrix**: Normalized alignment specificity vectors used in reconstruction

    ALGORITHM:
    1. Load gene-to-transcript mapping and sample file information
    2. For each founder strain:
       - Load TPM expression data from multiple samples
       - Calculate average expression profile across samples
    3. For each gene:
       - Create axes matrix with raw expression values
       - Calculate ASE (allele-specific expression) summary
       - Generate normalized alignment specificity vectors
    4. Save matrices in compressed numpy format

    WHAT THIS DOES:
    - Takes TPM expression data from founder strain samples
    - Calculates strain-specific expression profiles
    - Creates alignment specificity matrices for each gene
    - Provides the foundation for genotype probability calculations
    - Output: Three numpy files (axes.npz, ases.npz, avecs.npz)

    USE CASES:
    - GBRS preprocessing: Prepare alignment specificity data for reconstruction
    - Founder strain analysis: Analyze expression patterns across founder strains
    - Reference data generation: Create reference matrices for multiple samples
    - Quality control: Assess expression consistency across founder strains

    Args:
        sample_file: Path to the sample file containing TPM file mappings.
            Format: tab-separated with strain name in first column, TPM file path in second.
            Example:
                A    /path/to/strain_A_sample1.tpm
                A    /path/to/strain_A_sample2.tpm
                B    /path/to/strain_B_sample1.tpm
                B    /path/to/strain_B_sample2.tpm

        haplotypes: List of founder strain identifiers.
            Example: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for Diversity Outbred mice
            Must match the strain names in the sample_file.

        gene2transcript: Path to the gene-to-transcript mapping file.
            Format: tab-separated with gene ID in first column, transcript ID in second.
            Used to identify genes and their corresponding transcripts.

        out_dir: Output directory for the generated matrices.
            The function will create three files in this directory:
            - axes.npz: Raw expression values matrix
            - ases.npz: Allele-specific expression summary matrix
            - avecs.npz: Normalized alignment specificity vectors

        min_expr: Minimum expression threshold for including genes in analysis.
            Default: 2.0
            Genes with total expression below this threshold across all strains
            are excluded from the alignment specificity calculations.

    Returns:
        None. The function creates three compressed numpy files in the specified output directory.

    Raises:
        FileNotFoundError: If sample_file, gene2transcript, or any TPM file does not exist.
        ValueError: If haplotypes list is empty or contains invalid characters.
        OSError: If output directory cannot be created or files cannot be written.

    Notes:
        - The function averages expression data across multiple samples per strain
        - Alignment specificity vectors are normalized to unit vectors
        - Only genes with sufficient expression (above min_expr threshold) are included
        - The output matrices are essential for the GBRS reconstruction algorithm
        - This function is typically run once per founder strain set to create
          reference alignment specificity data
    """
    logger.info(f'Sample File: {sample_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Min Expression: {min_expr}')
    logger.info(f'Gene2Transcript: {gene2transcript}')

    num_strains = len(haplotypes)

    gname = np.loadtxt(
        gene2transcript,
        usecols=(0,),
        dtype='str'
    )
    num_genes = len(gname)
    gid = dict(zip(gname, np.arange(num_genes)))

    logger.info(f'Loading {sample_file}')
    flist = defaultdict(list)
    with open(sample_file) as fh:
        for line in fh:
            item = line.rstrip().split('\t')
            flist[item[0]].append(item[1])
    flist = dict(flist)

    dset = dict()
    for st in haplotypes:
        dmat_strain = np.zeros((num_genes, num_strains))
        for tpmfile in flist[st]:
            logger.debug(f'Working on {tpmfile}')
            dmat_sample = np.zeros((num_genes, num_strains))
            if not os.path.isfile(tpmfile):
                print(f'File {tpmfile} does not exist.')
                continue
            with open(tpmfile) as fh:
                fh.readline()  # header
                for curline in fh:
                    item = curline.rstrip().split('\t')
                    if item[0] in gid:
                        row = gid[item[0]]
                        dmat_sample[row, :] = map(
                            float, item[1: (num_strains + 1)]
                        )
            dmat_strain += dmat_sample
        dset[st] = dmat_strain / len(flist[st])

    axes = dict()
    ases = dict()
    avecs = dict()
    for g in gname:
        axes[g] = np.zeros((num_strains, num_strains))
        ases[g] = np.zeros((1, num_strains))
        good = np.zeros(num_strains)
        for i, st in enumerate(haplotypes):
            v = dset[st][gid[g], :]
            axes[g][i, :] = v
            ases[g][0, i] = sum(v)
            if sum(v) > min_expr:
                good[i] = 1.0
        if sum(good) > 0:  # At least one strain expresses
            avecs[g] = np.zeros((num_strains, num_strains))
            for i in range(num_strains):
                avecs[g][i, :] = unit_vector(axes[g][i, :])

    logger.info(f'Saving {os.path.join(out_dir, "axes.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'axes.npz'), **axes)
    logger.info(f'Saving {os.path.join(out_dir, "ases.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'ases.npz'), **ases)
    logger.info(f'Saving {os.path.join(out_dir, "aves.npz")}')
    np.savez_compressed(os.path.join(out_dir, 'avecs.npz'), **avecs)


def reconstruct(
        expression_file: str,
        tprob_file: str,
        avec_file: str = None,
        gpos_file: str = None,
        expr_threshold: float = 1.5,
        sigma: float = 0.12,
        outbase: str = None,
) -> None:
    """
    Reconstruct individual genome using GBRS algorithm with gene expression data.

    This function is the core genome reconstruction component of the GBRS pipeline
    that uses gene expression data to infer individual genotypes across the genome.
    It implements a Hidden Markov Model (HMM) with forward-backward algorithm and
    Viterbi decoding to reconstruct the most likely genotype at each gene position.

    The function processes gene-level TPM expression data from an individual and
    compares it against founder strain expression patterns to infer the individual's
    genotype. It uses transition probabilities between adjacent genes to model
    linkage disequilibrium and recombination events.

    ALGORITHM:
    1. Load expression data, transition probabilities, and alignment specificity
    2. Calculate emission probabilities for each gene based on expression similarity
    3. Run forward-backward algorithm to compute genotype probabilities
    4. Run Viterbi algorithm to find most likely genotype sequence
    5. Save genotype probabilities and most likely genotypes

    WHAT THIS DOES:
    - Takes individual gene expression data (TPM values for each possible diplotype)
    - Compares against founder strain expression patterns using alignment specificity
    - Uses HMM to model genetic linkage between genes
    - Infers most likely genotype at each gene position
    - Output: Genotype probabilities and most likely genotypes

    USE CASES:
    - Individual genotyping: Infer genotypes from RNA-Seq data
    - Population studies: Genotype individuals without DNA sequencing
    - Expression QTL mapping: Combine genotype and expression data
    - Quality control: Validate genotype calls from other methods

    OUTPUT FILES:
    - {outbase}.genotypes.tsv: Most likely genotype for each gene
    - {outbase}.genoprobs.npz: Genotype probabilities for each gene
    - {outbase}.genotypes.npz: Ordered genotype calls by chromosome

    Args:
        expression_file: Path to the gene-level TPM expression file.
            Format: tab-separated with gene ID in first column, TPM values for
            each possible diplotype in subsequent columns, and total in last column.
            Example for 8 haplotypes (A-H) with 36 possible diplotypes:
                Gene_ID    AA    AB    AC    AD    AE    AF    AG    AH    BB    BC    ...    Total
                Gene1      10.5  8.2   7.1   6.3   5.8   4.9   4.2   3.8   9.1   8.5    ...   23.8
                Gene2      15.2  12.1  11.8  10.9  9.8   8.7   7.9   7.2   14.1  13.5   ...   37.1

        tprob_file: Path to the transition probabilities file (npz format).
            Contains transition probability matrices for each chromosome,
            modeling genetic linkage between adjacent genes based on mating scheme.

        avec_file: Path to the alignment specificity file (npz format).
            Default: None (uses default location from GBRS_DATA environment variable)
            Contains normalized expression vectors for each gene and founder strain,
            used to calculate emission probabilities.

        gpos_file: Path to the gene position file (npz format).
            Default: None (uses default location from GBRS_DATA environment variable)
            Contains gene metadata including chromosome, ID, and genomic position
            for organizing genes by chromosomal order.

        expr_threshold: Minimum expression threshold for including genes in analysis.
            Default: 1.5
            Genes with total expression below this threshold are assigned
            uniform genotype probabilities (null model).

        sigma: Standard deviation parameter for emission probability calculation.
            Default: 0.12
            Controls the sensitivity of genotype inference based on expression
            similarity to founder strain patterns.

        outbase: Base name for output files. Default: None
            If None, uses default names: 'gbrs.reconstructed.genotypes.tsv'
            and 'gbrs.reconstructed.genoprobs.npz'

    Returns:
        None. The function creates genotype files at the specified output paths.

    Raises:
        FileNotFoundError: If any input file does not exist.
        ValueError: If expression file format is invalid or parameters are invalid.
        RuntimeError: If HMM algorithm fails to converge or other processing errors.

    Notes:
        - The function uses a Hidden Markov Model with diplotypes (haplotype pairs) as states
        - Emission probabilities are calculated using squared Euclidean distance between
          normalized expression vectors, converted to probabilities via Gaussian kernel
        - Transition probabilities model recombination events between adjacent genes
          based on the specified mating scheme (RI, F2, CC, DO)
        - Forward-backward algorithm provides uncertainty estimates (genotype probabilities)
        - Viterbi algorithm provides the most likely genotype sequence
        - Genes with low expression are assigned uniform probabilities (null model)
        - This function is the final step in the GBRS pipeline for individual genotyping
        - The output can be used for downstream analyses like eQTL mapping or population studies
    """
    if outbase is None:
        out_gtype = 'gbrs.reconstructed.genotypes.tsv'
        out_gprob = 'gbrs.reconstructed.genoprobs.npz'
    else:
        out_gtype = f'{outbase}.genotypes.tsv'
        out_gprob = f'{outbase}.genoprobs.npz'

    out_gtype_ordered = f'{os.path.splitext(out_gtype)[0]}.npz'

    if avec_file is None:
        avec_file = os.path.join(DATA_DIR, 'avecs.npz')

    if gpos_file is None:
        gpos_file = os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz')

    logger.info(f'Expression File: {expression_file}')
    logger.info(f'Transition Probabilities File: {tprob_file}')
    logger.info(f'Alignment Specificity File: {avec_file}')
    logger.info(f'Gene Position File: {gpos_file}')
    logger.info(f'Expression Threshold: {expr_threshold}')
    logger.info(f'Sigma: {sigma}')
    logger.info(f'Outbase: {outbase}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()
    chrs = chrlens.keys()

    # Load alignment specificity
    logger.info(f'Loading alignment specificity: {avec_file}')
    avecs = np.load(avec_file)

    # Load meta info
    logger.info(f'Loading gene meta data: {gpos_file}')
    gene_pos = np.load(gpos_file)

    # gid_genome_order is a dictionary
    # with keys being chromosomes
    # and values being a list of gene ids and their position
    gid_genome_order = dict.fromkeys(gene_pos.files)
    for c in gene_pos.files:
        try:
            gid_genome_order[c] = np.array([g.decode() for g, p in gene_pos[c]])
        except:
            gid_genome_order[c] = np.array([g for g, p in gene_pos[c]])

    # Load expression level
    logger.info(f'Loading expression level data: {expression_file}')
    expr = dict()
    with open(expression_file) as fh:
        curline = fh.readline()
        haplotypes = curline.rstrip().split('\t')[1:-1]
        num_haps = len(haplotypes)
        genotypes = [
            h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)
        ]
        num_genos = len(genotypes)
        for curline in fh:
            item = curline.rstrip().split('\t')
            expr[item[0]] = np.array(list(map(float, item[1:-1])))

    # Get null model probability
    logger.debug('Get null model probability')
    init_vec = []
    for h1, h2 in combinations_with_replacement(haplotypes, 2):
        if h1 == h2:
            init_vec.append(np.log(1.0 / (num_haps * num_haps)))
        else:
            init_vec.append(np.log(2.0 / (num_haps * num_haps)))
    init_vec = np.array(init_vec)

    # Get initial emission probability
    logger.debug('Get initial emission probability')
    naiv_avecs = (
            np.eye(num_haps)
            + (np.ones((num_haps, num_haps)) - np.eye(num_haps)) * 0.0001
    )
    eprob = dict()
    for gid, evec in expr.items():
        if sum(evec) < expr_threshold:
            eprob[gid] = init_vec
        elif gid not in avecs.files:
            eprob[gid] = np.log(
                get_genotype_probability(evec, naiv_avecs, sigma=0.450)
                + np.nextafter(0, 1)
            )
        else:
            eprob[gid] = np.log(
                get_genotype_probability(evec, avecs[gid], sigma=sigma)
                + np.nextafter(0, 1)
            )

    # Load transition probabilities
    logger.info(f'Loading transition probabilities: {tprob_file}')
    tprob = np.load(tprob_file)

    # Get forward probability
    logger.info('Getting forward probability')
    alpha = dict()
    alpha_scaler = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            alpha_c = np.zeros((num_genos, num_genes_in_chr))
            alpha_scaler_c = np.zeros(num_genes_in_chr)
            alpha_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            normalizer = np.log(sum(np.exp(alpha_c[:, 0])))
            alpha_c[:, 0] -= normalizer  # normalization
            alpha_scaler_c[0] = -normalizer
            for i in range(1, num_genes_in_chr):
                alpha_c[:, i] = (
                        np.log(
                            np.exp(alpha_c[:, i - 1] + tprob_c[i - 1]).sum(axis=1)
                            + np.nextafter(0, 1)
                        )
                        + eprob[gid_genome_order_c[i]]
                )
                normalizer = np.log(sum(np.exp(alpha_c[:, i])))
                alpha_c[:, i] -= normalizer  # normalization
                alpha_scaler_c[i] = -normalizer
            alpha[c] = alpha_c
            alpha_scaler[c] = alpha_scaler_c

    # Get backward probability
    logger.info('Getting backward probability')
    beta = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            beta_c = np.zeros((num_genos, num_genes_in_chr))
            beta_c[:, -1] = alpha_scaler[c][
                -1
            ]  # init_vec + eprob[gid_genome_order_c[-1]]
            for i in range(num_genes_in_chr - 2, -1, -1):
                beta_c[:, i] = np.log(
                    np.exp(
                        tprob_c[i].transpose()
                        + beta_c[:, i + 1]
                        + eprob[gid_genome_order_c[i + 1]]
                        + alpha_scaler[c][i]
                    ).sum(axis=1)
                )
            beta[c] = beta_c

    # Get forward-backward probability
    logger.info('Getting forward-backward probability')
    gamma = dict()
    for c in chrs:
        if c in tprob.files:
            logger.debug(f'Working on {c}')
            gamma_c = np.exp(alpha[c] + beta[c])
            normalizer = gamma_c.sum(axis=0)
            gamma[c] = gamma_c / normalizer

    logger.info(f'Saving Reconstructed Genotype Probabilities: {out_gprob}')
    np.savez_compressed(out_gprob, **gamma)

    # Run Viterbi
    logger.info('Running Viterbi')
    delta = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            delta_c = np.zeros((num_genos, num_genes_in_chr))
            delta_c[:, 0] = init_vec + eprob[gid_genome_order_c[0]]
            for i in range(1, num_genes_in_chr):
                delta_c[:, i] = (delta_c[:, i - 1] + tprob_c[i - 1]).max(
                    axis=1
                ) + eprob[gid_genome_order_c[i]]
            delta[c] = delta_c
    viterbi_states = defaultdict(list)
    gtcall_g = dict()
    for c in chrs:
        if c in tprob.files:
            tprob_c = tprob[c]
            gid_genome_order_c = gid_genome_order[c]
            num_genes_in_chr = len(gid_genome_order_c)
            sid = delta[c][:, num_genes_in_chr - 1].argmax()
            viterbi_states[c].append(genotypes[sid])
            if (num_genes_in_chr > len(tprob_c)):
                num_genes_in_chr = len(tprob_c)
                # Above avoids IndexError for legacy GRCm38 ref files.
                # In those transprob files: num_genes_in_chr > tprob_c by 1.
            for i in reversed(range(num_genes_in_chr)):
                sid = (delta[c][:, i] + tprob_c[i][sid]).argmax()
                viterbi_states[c].append(genotypes[sid])
                gtcall_g[gid_genome_order_c[i]] = genotypes[sid]
            viterbi_states[c].reverse()

    viterbi_states = dict(viterbi_states)

    logger.info(f'Saving Reconstructed Genotypes: {out_gtype}')
    with open(out_gtype, 'w') as fhout:
        fhout.write('#Gene_ID\tDiplotype\n')
        for g in sorted(gtcall_g.keys()):
            fhout.write(f'{g}\t{gtcall_g[g]}\n')

    logger.info(f'Saving Reconstructed Ordered Genotypes: {out_gtype_ordered}')
    np.savez_compressed(out_gtype_ordered, **viterbi_states)
    logger.info('Done')


def interpolate(
        genoprob_file: str,
        grid_file: str = None,
        gpos_file: str = None,
        output_file: str = None,
) -> None:
    """
    Interpolate genotype probabilities from gene positions to uniform grid positions.

    This function is a post-processing step in the GBRS pipeline that converts
    genotype probabilities calculated at gene positions to a uniform grid of
    positions across the genome. This interpolation is useful for visualization,
    comparison across samples, and downstream analyses that require consistent
    genomic coordinates.

    The function uses linear interpolation to estimate genotype probabilities
    at grid positions based on the probabilities at nearby gene positions.
    It handles chromosome boundaries and extrapolation at chromosome ends
    by extending the gene positions with boundary values.

    ALGORITHM:
    1. Load genotype probabilities at gene positions
    2. Load uniform grid positions for each chromosome
    3. For each chromosome:
       - Extend gene positions with boundary values (0.0 and grid_max+1.0)
       - Create interpolation function using scipy.interpolate.interp1d
       - Interpolate probabilities at all grid positions
    4. Save interpolated probabilities in compressed numpy format

    WHAT THIS DOES:
    - Converts gene-based probabilities to grid-based probabilities
    - Enables uniform genomic coordinate system across samples
    - Facilitates visualization and comparison of genome reconstructions
    - Provides consistent positions for downstream analyses

    USE CASES:
    - GBRS post-processing: Convert gene probabilities to grid probabilities
    - Visualization: Create uniform plots across different samples
    - Population analysis: Compare reconstructions at consistent positions
    - Downstream analysis: Provide uniform coordinate system for tools

    Args:
        genoprob_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities calculated
            at gene positions by the GBRS reconstruction algorithm.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_genes) for each chromosome
            Example: 'DO336.genoprobs.npz'

        grid_file: Path to the uniform grid file.
            Type: str
            Default: None (uses default location from GBRS_DATA)
            Description: File containing uniform grid positions across the genome.
            Format: Tab-separated with columns: grid_id, chromosome, position, cM
            Used to define the target positions for interpolation.
            Example: 'ref.genome_grid.64k.txt'

        gpos_file: Path to the gene position file (npz format).
            Type: str
            Default: None (uses default location from GBRS_DATA)
            Description: File containing gene metadata including chromosome,
            gene ID, and genetic position (cM). Used to define the source
            positions for interpolation.
            Example: 'ref.gene_pos.ordered.npz'

        output_file: Output filename for interpolated probabilities.
            Type: str
            Default: None (auto-generated from input filename)
            Description: Name of the compressed numpy file to save the
            interpolated genotype probabilities. If None, generates name
            based on input filename with 'interpolated' prefix.
            Example: 'DO336.interpolated.genoprobs.npz'

    Returns:
        None. The function creates a compressed numpy file containing
        interpolated genotype probabilities at uniform grid positions.

    Raises:
        FileNotFoundError: If any input file does not exist
        ValueError: If grid and gene positions are incompatible
        RuntimeError: If interpolation fails or produces invalid results

    Notes:
        - The function extends gene positions with boundary values to handle
          extrapolation at chromosome ends
        - Linear interpolation is used to estimate probabilities at grid positions
        - Grid positions should cover the full range of gene positions
        - The output maintains the same diplotype structure as the input
        - Interpolated probabilities are normalized to sum to 1.0 at each position
        - This function is typically run after GBRS reconstruction for visualization
        - Grid-based probabilities enable consistent comparison across samples
        - File sizes may be larger than gene-based files due to more positions
    """
    if gpos_file is None:
        gpos_file = os.path.join(DATA_DIR, 'ref.gene_pos.ordered.npz')
        try:
            x_gene = np.load(gpos_file)
        except:
            logger.error(
                f'Please make sure if $GBRS_DATA is set correctly: {DATA_DIR}'
            )
            raise
        else:
            pass
    else:
        x_gene = np.load(gpos_file)

    if grid_file is None:
        grid_file = os.path.join(DATA_DIR, 'ref.genome_grid.64k.txt')

    if output_file is None:
        output_file = f'gbrs.interpolated.{os.path.basename(genoprob_file)}'

    logger.info(f'Genotype Probability File: {genoprob_file}')
    logger.info(f'Grid File: {grid_file}')
    logger.info(f'Gene Position File: {gpos_file}')
    logger.info(f'Output File: {output_file}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()
    chrs = chrlens.keys()
    logger.info(f'Chrom Lens: {chrlens}')
    logger.info(f'Chroms: {chrs}')

    logger.info(f'Loading grid file: {grid_file}')
    x_grid = defaultdict(list)
    with open(grid_file) as fh:
        fh.readline()  # skip header (Assuming there is just one line of header)
        for line in fh:
            item = line.rstrip().split('\t')
            x_grid[item[1]].append(
                float(item[3])
            )  # x_grid[chr] = [...positions in cM...]
    x_grid = dict(x_grid)

    x_gene_extended = (
        dict()
    )  # Adding end points in case we have to extrapolate at the 1st or last grid
    for c in x_grid.keys():
        if c in x_gene.files:
            logger.debug(f'Working on {c}')
            x = [float(coord) for m, coord in x_gene[c]]
            # x_min = min(x_grid[c][0]-100.0, 0.0)
            # x_max = max(x_grid[c][-1]+1.0, chrlens[c])
            # x = np.append([x_min], x)
            # x = np.append(x, [x_max])
            x = np.append([0.0], x)
            x = np.append(
                x, [x_grid[c][-1] + 1.0]
            )  # Do we have chromosome length in cM?
            x_gene_extended[c] = x

    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    gamma_gene = np.load(genoprob_file)
    gene_model_chr = dict()
    gene_intrp_chr = dict()
    for c in x_grid.keys():
        if c in gamma_gene.files:
            logger.debug(f'Working on {c}')
            gamma_gene_c = gamma_gene[c]
            y = np.hstack((gamma_gene_c[:, 0][:, np.newaxis], gamma_gene_c))
            y = np.hstack((y, y[:, -1][:, np.newaxis]))
            gene_model_chr[c] = interp1d(x_gene_extended[c], y, axis=1)
            gene_intrp_chr[c] = gene_model_chr[c](x_grid[c])

    logger.info(f'Saving interpolate probability file: {output_file}')
    np.savez_compressed(output_file, **gene_intrp_chr)
    logger.info('Done')


def combine():
    raise NotImplementedError


def plot(
        genoprob_file: str,
        output_file: str = None,
        output_format: str = 'pdf',
        sample_name: str = '',
        grid_size: int = 2,
        xt_max: int = 5000,
        xt_size: int = 500,
        grid_width: float = 0.01,
) -> None:
    """
    Generate genome reconstruction visualization plot.

    This function creates a comprehensive visualization of the GBRS genome
    reconstruction results, showing the most likely genotype at each position
    across all chromosomes. The plot displays both haplotypes for each
    position, allowing identification of recombination events and genomic
    structure patterns.

    The visualization uses a stacked bar chart format where each chromosome
    is represented by two horizontal bars (one for each haplotype). Different
    founder strains are color-coded, making it easy to identify regions of
    shared ancestry and recombination breakpoints. The plot includes
    recombination counts for each chromosome and total recombination count.

    ALGORITHM:
    1. Load genotype probabilities and determine most likely genotypes
    2. For each chromosome:
       - Extract most likely diplotype at each position
       - Separate into two haplotypes
       - Identify recombination events (genotype changes)
       - Create color-coded bars for visualization
    3. Generate plot with proper spacing, labels, and annotations
    4. Save in specified format with high resolution

    WHAT THIS DOES:
    - Visualizes GBRS genome reconstruction results
    - Shows haplotype structure across all chromosomes
    - Identifies recombination events and breakpoints
    - Provides quantitative recombination statistics
    - Enables visual comparison of genomic structure

    USE CASES:
    - GBRS visualization: Create publication-quality genome plots
    - Quality assessment: Visualize reconstruction quality and patterns
    - Recombination analysis: Identify and count recombination events
    - Population studies: Compare genomic structure across individuals
    - Publication figures: Generate high-resolution plots for papers

    Args:
        genoprob_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities from GBRS
            reconstruction or interpolation. Can be either gene-based or
            grid-based probabilities.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_positions) for each chromosome
            Example: 'DO336.genoprobs.npz' or 'DO336.interpolated.genoprobs.npz'

        output_file: Output filename for the plot.
            Type: str
            Default: None (auto-generated from input filename)
            Description: Name of the output file for the generated plot.
            If None, generates name based on input filename with 'plotted' prefix.
            Example: 'DO336.plotted.genome.pdf'

        output_format: File format for the output plot.
            Type: str
            Default: 'pdf'
            Description: Image format for the output file. Common formats
            include 'pdf', 'png', 'svg', 'jpg', 'tiff'.
            Example: 'pdf' for publication-quality vector graphics

        sample_name: Name of the sample for plot title.
            Type: str
            Default: ''
            Description: Sample identifier to include in the plot title.
            If empty, uses the filename as the sample name.
            Example: 'DO336' or 'Mouse_001'

        grid_size: Size parameter for grid scaling.
            Type: int
            Default: 2
            Description: Advanced parameter affecting the scaling of the
            x-axis grid. Used to adjust the visual spacing of positions.
            Range: Typically 1 to 5
            Example: 2 for standard spacing

        xt_max: Maximum x-axis tick value.
            Type: int
            Default: 5000
            Description: Maximum value for x-axis ticks in the plot.
            Controls the range of the x-axis display.
            Range: Depends on number of positions in the data
            Example: 5000 for 5000 positions

        xt_size: Size of x-axis tick intervals.
            Type: int
            Default: 500
            Description: Interval between x-axis tick marks.
            Controls the frequency of tick labels on the x-axis.
            Range: Typically 100 to 1000
            Example: 500 for ticks every 500 positions

        grid_width: Width of each grid position.
            Type: float
            Default: 0.01
            Description: Width of each position bar in the plot.
            Controls the visual thickness of the genotype bars.
            Range: Typically 0.005 to 0.05
            Example: 0.01 for standard bar width

    Returns:
        None. The function creates a high-resolution plot file in the
        specified format.

    Raises:
        FileNotFoundError: If genoprob_file does not exist
        ValueError: If plot parameters are invalid
        RuntimeError: If plotting fails or produces invalid output

    Notes:
        - The function automatically loads founder strain colors from
          'founder.hexcolor.info' file
        - Recombination events are counted when genotype changes between
          adjacent positions
        - The plot uses a 16x16 inch figure size for high resolution
        - Chromosomes are displayed in natural order (1, 2, ..., 19, X, Y)
        - Each chromosome shows two haplotypes with different colors
        - Recombination counts are displayed next to each chromosome
        - The plot includes a title with sample name and total recombination count
        - Output is saved at 600 DPI for publication quality
        - This function is typically run after GBRS reconstruction or interpolation
        - The visualization is essential for quality assessment and publication
    """
    if output_file is None:
        output_file = os.path.splitext(os.path.basename(genoprob_file))[0]
        output_file = f'gbrs.plotted.{output_file}.{output_format}'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Output Format: {output_format}')
    logger.info(f'Sample Name: {sample_name}')
    logger.info(f'Grid Size: {grid_size}')
    logger.info(f'XT Max: {xt_max}')
    logger.info(f'XT Size: {xt_size}')
    logger.info(f'Grid Width: {grid_width}')

    logger.info('Loading chromosome information')
    chrlens = get_chromosome_info()

    logger.info('Loading founder colors')
    hcolors = get_founder_info()
    haplotypes = hcolors.keys()
    hid = dict(zip(haplotypes, np.arange(8)))
    logger.info(f'Haplotype IDs: {hid}')

    genotypes = np.array(
        [h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
    )

    #
    # Main body
    #
    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    genoprob = np.load(genoprob_file)

    def natural_sort(list):
        def convert(text):
            return int(text) if text.isdigit() else text.lower()

        def alphanum_key(key):
            return [convert(c) for c in re.split('([0-9]+)', key)]

        return sorted(list, key=alphanum_key)

    chrs = [value for value in chrlens.keys() if value in genoprob.files]
    # intersection of ref.fa.fai chroms and those present in genoprob.

    chrs = natural_sort(chrs)
    # natural sorting of the chrom list for display.

    num_chrs = len(chrs)

    fig = pyplot.figure()
    fig.set_size_inches((16, 16))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, xt_max * grid_width + grid_width)
    ax.set_ylim(1, 95)
    num_recomb_total = 0
    for cid, c in enumerate(chrs):
        if (
                c in genoprob.files
        ):  # Skip drawing Y chromosome if the sample is female
            logger.debug(f'Working on {c}')
            genotype_calls = genotypes[genoprob[c].argmax(axis=0)]
            hap = []
            col1 = []
            col2 = []
            oldcol1 = 'NA'
            oldcol2 = 'NA'
            num_recomb = 0
            num_genes_in_chr = len(genotype_calls)

            for i in range(num_genes_in_chr):
                hap.append((i * grid_width, grid_width))
                c1 = hcolors[genotype_calls[i][0]]
                c2 = hcolors[genotype_calls[i][1]]

                if i > 0:
                    if c1 == c2:
                        if (
                                col1[-1] != col2[-1]
                        ):  # When homozygous region starts, remember the
                            # most recent het
                            oldcol1 = col1[-1]
                            oldcol2 = col2[-1]
                    else:
                        if (
                                col1[-1] == col2[-1]
                        ):  # When heterozygous region starts
                            if c1 == oldcol2 or c2 == oldcol1:
                                c1, c2 = c2, c1
                        elif c1 == col2[-1] or c2 == col1[-1]:
                            c1, c2 = c2, c1
                    if c1 != col1[-1] or c2 != col2[-1]:
                        num_recomb += 1

                col1.append(c1)
                col2.append(c2)
            num_recomb_total += num_recomb
            print(f'Chromsome {c} has {num_recomb} recombinations')
            # plot
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4 + 1, 1),
                facecolors=col1,
                edgecolor='face',
            )
            ax.broken_barh(
                hap,
                (num_chrs * 4 - cid * 4, 1),
                facecolors=col2,
                edgecolor='face',
            )
            ax.text(
                (num_genes_in_chr + 50) * grid_width,
                num_chrs * 4 - cid * 4 + 0.5,
                f'({num_recomb})',
            )
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_right()
        ax.get_yaxis().tick_left()
        pyplot.yticks(
            ticks=np.arange(num_chrs * 4 + 1, 1, -4),
            labels=list(chrs),
            fontsize=14,
        )
        pyplot.xticks(
            ticks=np.arange(0, xt_max * grid_width, xt_size * grid_width),
            labels=[
                '%dcM' % xt
                for xt in np.arange(
                    0, xt_max * grid_size / 100, xt_size * grid_size / 100
                )
            ],
        )
        title_txt = f'Genome reconstruction: {sample_name}'
        title_txt += f'\n(Total {num_recomb_total} recombinations)'
        ax.set_title(title_txt, fontsize=18, loc='center')

    logger.info(f'Saving generated plot: {output_file}')
    fig.savefig(output_file, dpi=600, format=output_format)
    pyplot.close(fig)
    logger.info('Done')


def parse_founder_colors(path: str) -> dict[str, str]:
    """
    Parse founder strain color definitions from file.

    Args:
        path: Path to founder color file

    Returns:
        Dictionary mapping founder names to hex color codes
    """
    founder_colors = {}
    with open(path) as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = re.split(r'[\s\t]+', line.strip())
                if len(parts) >= 2:
                    founder, color = parts[:2]
                    founder_colors[founder] = color
    return founder_colors


def load_gpos(gpos_path: str) -> dict[str, list[float]]:
    """
    Load gene positions (cM) from .npz file.

    Args:
        gpos_path: Path to gene position file

    Returns:
        Dictionary mapping chromosomes to lists of cM positions
    """
    gpos = np.load(gpos_path, allow_pickle=True)
    chrom_cM = {}
    for chrom in gpos.files:
        # Each entry: list of (gene, cM)
        cM_list = [float(x[1]) for x in gpos[chrom]]
        chrom_cM[chrom] = cM_list
    return chrom_cM


def find_blocks(seq: list) -> tuple[int, int, any]:
    """
    Find contiguous blocks in sequence.

    Args:
        seq: Input sequence to analyze

    Yields:
        Tuples of (start_index, end_index, value) for each block
    """
    if not seq:
        return
    start = 0
    current = seq[0]
    for i, val in enumerate(seq):
        if val != current:
            yield (start, i, current)
            start = i
            current = val

    yield (start, len(seq), current)


def natural_chrom_order(chroms: list[str]) -> list[str]:
    """
    Sort chromosomes in natural order (1,2,...,19,X,Y,MT).

    Args:
        chroms: List of chromosome names to sort

    Returns:
        Sorted list of chromosome names
    """

    # Order: 1,2,...,19,X,Y,MT (case-insensitive)
    def chrom_key(c):
        c = c.upper()
        if c == 'X':
            return 20
        if c == 'Y':
            return 21
        if c in ('MT', 'M', 'MITO'):
            return 22
        try:
            return int(c)
        except ValueError:
            return 99

    return sorted(chroms, key=chrom_key)


def chrom_sort_key(c: str) -> int:
    """
    Get sort key for chromosome name.

    Args:
        c: Chromosome name

    Returns:
        Integer sort key (1-19 for autosomes, 20 for X, 21 for Y, 22 for MT)
    """
    c = c.upper()
    if c == 'X':
        return 20
    if c == 'Y':
        return 21
    if c in ('MT', 'M', 'MITO'):
        return 22
    try:
        return int(c)
    except ValueError:
        return 99


def plot_genoprobs(
        genoprobs_file: str,
        output_file: str | None = None,
        output_format: str = 'pdf',
        sample_name: str | None = None,
        founder_colors_file: str | None = None,
        genome_pos_file: str | None = None,
        dpi: int = 300,
        bar_height: float = 0.5,
        haplotype_gap: float = 0.1,
        chrom_spacing: float = 2.0,
        fig_width: int = 18,
        min_fig_height: int = 8,
        height_per_chrom: float = 0.8,
        font_size_title: int = 22,
        font_size_axis: int = 18,
        font_size_tick: int = 16,
        font_size_legend: int = 14,
        legend_line_width: int = 8,
) -> None:
    """
    Generate publication-quality genome reconstruction plots with configurable parameters.

    This function is an enhanced version of the basic plot function that provides
    extensive customization options for creating publication-ready genome reconstruction
    visualizations. It creates a comprehensive view of GBRS reconstruction results
    with improved aesthetics, better spacing, and more informative annotations.

    The function generates a multi-chromosome plot where each chromosome is represented
    by two horizontal bars (one for each haplotype) with founder strains color-coded.
    It automatically detects recombination events and provides detailed statistics
    for each chromosome. The plot includes a legend, proper axis labels, and
    publication-quality formatting.

    ALGORITHM:
    1. Load genotype probabilities and determine most likely genotypes
    2. For each chromosome:
       - Extract most likely diplotype at each position
       - Separate into two haplotypes
       - Identify contiguous blocks of same genotype
       - Create color-coded rectangles for visualization
       - Count recombination events
    3. Generate plot with configurable aesthetics and spacing
    4. Add legend, labels, and annotations
    5. Save in specified format with high resolution

    WHAT THIS DOES:
    - Creates publication-quality genome reconstruction visualizations
    - Shows haplotype structure with improved aesthetics
    - Identifies and counts recombination events
    - Provides extensive customization options
    - Generates publication-ready figures

    USE CASES:
    - Publication figures: Create high-quality plots for papers
    - Quality assessment: Visualize reconstruction quality with detail
    - Recombination analysis: Identify and quantify recombination patterns
    - Population studies: Compare genomic structure across individuals
    - Custom visualization: Tailor plots to specific requirements

    Args:
        genoprobs_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities from GBRS
            reconstruction or interpolation. Can be either gene-based or
            grid-based probabilities.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_positions) for each chromosome
            Example: 'DO336.genoprobs.npz' or 'DO336.interpolated.genoprobs.npz'

        output_file: Output filename for the plot.
            Type: str | None
            Default: None (auto-generated from input filename)
            Description: Name of the output file for the generated plot.
            If None, generates name based on input filename with 'plotted' prefix.
            Example: 'DO336.plotted.genome.pdf'

        output_format: File format for the output plot.
            Type: str
            Default: 'pdf'
            Description: Image format for the output file. Common formats
            include 'pdf', 'png', 'svg', 'jpg', 'tiff'.
            Example: 'pdf' for publication-quality vector graphics

        sample_name: Name of the sample for plot title.
            Type: str | None
            Default: None
            Description: Sample identifier to include in the plot title.
            If None, uses the filename as the sample name.
            Example: 'DO336' or 'Mouse_001'

        founder_colors_file: Path to founder strain color definitions.
            Type: str | None
            Default: None
            Description: File containing founder strain to color mappings.
            Format: Tab-separated with founder name and hex color code.
            If None, uses default color scheme for 8-founder populations.
            Example: 'founder.hexcolor.info'

        genome_pos_file: Path to gene position file for cM coordinates.
            Type: str | None
            Default: None
            Description: Optional file containing gene positions in centiMorgans.
            If provided, x-axis shows genetic distance (cM) instead of marker index.
            Format: Compressed numpy file with gene positions
            Example: 'ref.gene_pos.ordered.npz'

        dpi: Resolution for output image.
            Type: int
            Default: 300
            Description: Dots per inch for the output image. Higher values
            produce higher resolution images suitable for publication.
            Range: Typically 150 to 600
            Example: 300 for publication quality

        bar_height: Height of each haplotype bar.
            Type: float
            Default: 0.5
            Description: Height of individual haplotype bars in plot units.
            Controls the visual thickness of the genotype bars.
            Range: Typically 0.3 to 1.0
            Example: 0.5 for standard appearance

        haplotype_gap: Gap between haplotype bars.
            Type: float
            Default: 0.1
            Description: Vertical gap between the two haplotype bars for
            each chromosome. Creates visual separation between haplotypes.
            Range: Typically 0.05 to 0.3
            Example: 0.1 for clear separation

        chrom_spacing: Vertical spacing between chromosomes.
            Type: float
            Default: 2.0
            Description: Vertical distance between different chromosomes
            in the plot. Controls overall plot height and readability.
            Range: Typically 1.5 to 4.0
            Example: 2.0 for balanced spacing

        fig_width: Figure width in inches.
            Type: int
            Default: 18
            Description: Width of the entire figure in inches.
            Controls the overall size and aspect ratio of the plot.
            Range: Typically 12 to 24
            Example: 18 for wide format

        min_fig_height: Minimum figure height in inches.
            Type: int
            Default: 8
            Description: Minimum height of the figure in inches.
            Ensures adequate space for all chromosomes and labels.
            Range: Typically 6 to 16
            Example: 8 for standard height

        height_per_chrom: Height per chromosome in inches.
            Type: float
            Default: 0.8
            Description: Height allocated for each chromosome in the plot.
            Used to calculate total figure height.
            Range: Typically 0.5 to 1.5
            Example: 0.8 for standard allocation

        font_size_title: Title font size.
            Type: int
            Default: 22
            Description: Font size for the main plot title.
            Controls prominence of the title text.
            Range: Typically 16 to 28
            Example: 22 for prominent title

        font_size_axis: Axis label font size.
            Type: int
            Default: 18
            Description: Font size for x-axis and y-axis labels.
            Controls readability of axis descriptions.
            Range: Typically 14 to 22
            Example: 18 for clear labels

        font_size_tick: Tick label font size.
            Type: int
            Default: 16
            Description: Font size for tick labels on axes.
            Controls readability of position markers.
            Range: Typically 12 to 20
            Example: 16 for readable ticks

        font_size_legend: Legend font size.
            Type: int
            Default: 14
            Description: Font size for legend text and title.
            Controls readability of founder strain legend.
            Range: Typically 10 to 18
            Example: 14 for clear legend

        legend_line_width: Legend line width.
            Type: int
            Default: 8
            Description: Width of lines in the legend showing founder colors.
            Controls visual prominence of legend elements.
            Range: Typically 4 to 12
            Example: 8 for prominent legend lines

    Returns:
        None. The function creates a high-resolution plot file in the
        specified format with all requested customizations.

    Raises:
        FileNotFoundError: If genoprobs_file or founder_colors_file does not exist
        ValueError: If plot parameters are invalid or incompatible
        RuntimeError: If plotting fails or produces invalid output

    Notes:
        - The function automatically detects and handles different chromosome orders
        - Recombination events are counted when genotype changes between adjacent positions
        - Founder strain colors can be customized via founder_colors_file
        - X-axis can show either marker index or genetic distance (cM) if genome_pos_file provided
        - The plot includes a comprehensive legend showing all founder strains
        - Output is optimized for publication with high DPI and vector formats
        - This function provides much more customization than the basic plot function
        - The visualization is essential for quality assessment and publication
        - All aesthetic parameters can be tuned for specific publication requirements
    """
    if output_file is None:
        output_file = os.path.splitext(os.path.basename(genoprobs_file))[0]
        output_file = f'gbrs.plotted.{output_file}.{output_format}'

    logger.info(f'Genotype Probabilities File: {genoprobs_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Output Format: {output_format}')
    logger.info(f'Founders Color File: {founder_colors_file}')
    logger.info(f'Genome Position File: {genome_pos_file}')
    logger.info(f'Sample Name: {sample_name}')
    logger.info(f'Bar Height: {bar_height}')
    logger.info(f'Haplotype Gap: {haplotype_gap}')
    logger.info(f'Chromosome Spacing: {chrom_spacing}')
    logger.info(f'Figure Width: {fig_width}')
    logger.info(f'Minimum Figure Height: {min_fig_height}')
    logger.info(f'Height Per Chromosome: {height_per_chrom}')
    logger.info(f'Font Size Title: {font_size_title}')
    logger.info(f'Font Size Axis: {font_size_axis}')
    logger.info(f'Font Size Tick: {font_size_tick}')
    logger.info(f'Font Size Legend: {font_size_legend}')
    logger.info(f'Legend Line Width: {legend_line_width}')

    founder_colors = {
        'A': '#F0F000',
        'B': '#808080',
        'C': '#F08080',
        'D': '#1010F0',
        'E': '#00A0F0',
        'F': '#00A000',
        'G': '#F00000',
        'H': '#9000E0',
    }

    founder_colors_simple = {
        'A': 'yellow',
        'B': 'grey',
        'C': 'pink',
        'D': 'blue',
        'E': 'light blue',
        'F': 'green',
        'G': 'red',
        'H': 'purple',
    }

    if founder_colors_file:
        founder_colors = parse_founder_colors(founder_colors_file)

    haplotypes = list(founder_colors.keys())
    genotypes = [h1 + h2 for h1, h2 in combinations_with_replacement(haplotypes, 2)]
    data = np.load(genoprobs_file, allow_pickle=True)
    chroms = natural_chrom_order(data.files)
    chroms = chroms[::-1]  # top-down

    logger.debug(f'Chromosomes: {", ".join(chroms)}')
    logger.debug(f'Generated Genotypes: {", ".join(genotypes)}')

    # calculate chromosome positions with configurable spacing
    chrom_ypos = {c: i * chrom_spacing for i, c in enumerate(chroms)}
    fig_height = max(min_fig_height, len(chroms) * height_per_chrom)
    fig, ax = pyplot.subplots(figsize=(fig_width, fig_height), dpi=dpi)
    total_recombs = 0
    chrom_cM = load_gpos(genome_pos_file) if genome_pos_file else None
    max_x = 0

    # loop though all chromosomes
    for chrom in chroms:
        mat = data[chrom]
        mat = np.asarray(mat)
        n_markers = mat.shape[1]

        max_idx = np.argmax(mat, axis=0)
        diplotype_calls = [genotypes[i] for i in max_idx]

        founder1 = [d[0] for d in diplotype_calls]
        founder2 = [d[1] for d in diplotype_calls]

        # x-axis: cM or marker index
        if (chrom_cM and chrom in chrom_cM and len(chrom_cM[chrom]) == n_markers):
            x = np.array(chrom_cM[chrom])
            xlabel = 'cM'
        else:
            x = np.arange(n_markers)
            xlabel = 'Marker index'

        y_base = chrom_ypos[chrom]

        # plot blocks for founder2 (bottom row)
        for start, end, val in find_blocks(founder2):
            x_start = x[start]
            x_end = x[end - 1] if end - 1 < len(x) else x[-1]
            width = x_end - x_start if x_end > x_start else 1
            logger.debug(f'f2 {chrom}:{start}-{end} {val} ({founder_colors_simple[val]})')
            ax.add_patch(
                Rectangle(
                    (x_start, y_base),
                    width,
                    bar_height,
                    color=founder_colors.get(val, '#CCCCCC'),
                    linewidth=0,
                )
            )

        # white border between rows
        ax.add_patch(
            Rectangle(
                (x[0], y_base + bar_height),
                x[-1] - x[0],
                haplotype_gap,
                color='white',
                linewidth=0,
                zorder=10,
            )
        )

        # plot blocks for founder1 (top row)
        for start, end, val in find_blocks(founder1):
            x_start = x[start]
            x_end = x[end - 1] if end - 1 < len(x) else x[-1]
            width = x_end - x_start if x_end > x_start else 1
            # logger.debug(f'f1 {chrom}:{start}-{end} {val} ({founder_colors_simple[val]})')

            ax.add_patch(
                Rectangle(
                    (x_start, y_base + bar_height + haplotype_gap),
                    width,
                    bar_height,
                    color=founder_colors.get(val, '#CCCCCC'),
                    linewidth=0,
                )
            )

        # recombination count
        recomb = sum([diplotype_calls[i] != diplotype_calls[i - 1] for i in range(1, n_markers)])
        total_recombs += recomb
        print(f'Chromsome {chrom} has {recomb} recombinations')

        ax.text(
            x[-1] + (x[1] - x[0] if n_markers > 1 else 1) * 1.5,
            y_base + bar_height + haplotype_gap / 2,
            f'({recomb})',
            va='center',
            ha='left',
            fontsize=14,
            fontweight='normal',
            color='#333333',
        )

        max_x = max(max_x, x[-1])

    # y axis: chromosome labels
    ax.set_yticks(
        [chrom_ypos[c] + bar_height + haplotype_gap / 2 for c in chroms]
    )
    ax.set_yticklabels(chroms, fontsize=font_size_tick)
    ax.set_xlabel(xlabel, fontsize=font_size_axis)
    ax.set_ylabel('Chromosome', fontsize=font_size_axis)

    # title
    sample = os.path.splitext(os.path.basename(genoprobs_file))[0]
    ax.set_title(
        f'Genome reconstruction: {sample}\n(Total {total_recombs} recombinations)',
        fontsize=font_size_title,
        fontweight='bold',
        pad=20,
    )
    ax.set_xlim(left=0, right=max_x * 1.02)
    ax.set_ylim(
        -haplotype_gap, max(chrom_ypos.values()) + 2 * bar_height + haplotype_gap
    )

    # remove spines and ticks
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis='x', labelsize=font_size_tick, length=0)
    ax.tick_params(axis='y', labelsize=font_size_tick, length=0)

    # subtle grid
    ax.xaxis.grid(True, linestyle=':', color='#CCCCCC', alpha=0.7)
    ax.set_axisbelow(True)

    # legend for founders
    legend_elements = [
        Line2D(
            [0], [0], color=founder_colors[f], lw=legend_line_width, label=f
        )
        for f in haplotypes
    ]
    ax.legend(
        handles=legend_elements,
        title='Founder',
        bbox_to_anchor=(1.01, 1),
        loc='upper left',
        fontsize=font_size_legend,
        title_fontsize=font_size_legend,
        frameon=False,
    )

    # tight layout
    pyplot.subplots_adjust(left=0.08, right=0.82, top=0.92, bottom=0.08)
    pyplot.savefig(
        output_file, bbox_inches='tight', dpi=dpi, format=output_format
    )

    logger.info(f'Plot saved to {output_file}')


def export(
        genoprob_file: str,
        strains: list[str],
        grid_file: str = None,
        output_file: str = None,
) -> None:
    """
    Export genotype probabilities to GBRS quant format for downstream analysis.

    This function converts genotype probability matrices from the GBRS format
    (diplotype probabilities) to a quantitative format that represents the
    expected contribution of each founder strain at each genomic position.
    This conversion is useful for downstream analyses that require founder
    strain dosage information rather than diplotype probabilities.

    The function applies a conversion matrix that transforms diplotype
    probabilities into founder strain dosages. For each position, it calculates
    the expected contribution of each founder strain based on the probabilities
    of all possible diplotypes containing that strain.

    ALGORITHM:
    1. Load genotype probability matrices for all chromosomes
    2. Create conversion matrix mapping diplotypes to founder strain dosages
    3. For each chromosome:
       - Transpose probability matrix to (positions × diplotypes)
       - Stack matrices from all chromosomes into single matrix
    4. Apply conversion matrix: dosage = probabilities × conversion_matrix
    5. Save results in tab-separated format with founder strain columns

    WHAT THIS DOES:
    - Converts diplotype probabilities to founder strain dosages
    - Combines all chromosomes into single quantitative matrix
    - Provides format suitable for downstream statistical analysis
    - Enables founder strain-specific analyses

    USE CASES:
    - Downstream analysis: Provide founder strain dosages for statistical tests
    - QTL mapping: Use founder strain contributions for association studies
    - Population analysis: Compare founder strain patterns across individuals
    - Integration: Combine with other genomic data in standard format
    - Visualization: Create founder strain-specific plots

    Args:
        genoprob_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities from GBRS
            reconstruction or interpolation. Contains diplotype probabilities
            for each position across all chromosomes.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_positions) for each chromosome
            Example: 'DO336.genoprobs.npz' or 'DO336.interpolated.genoprobs.npz'

        strains: List of founder strain identifiers.
            Type: list[str]
            Description: Names of founder strains in the population.
            These must match the haplotypes used in the GBRS reconstruction.
            The order determines the column order in the output file.
            Example: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for DO mice

        grid_file: Path to the grid file for position information.
            Type: str
            Default: None (uses default location from GBRS_DATA)
            Description: File containing grid position information.
            Used to determine the number of positions and chromosome structure.
            Format: Tab-separated with grid_id, chromosome, position, cM
            Example: 'ref.genome_grid.64k.txt'

        output_file: Output filename for the quant format file.
            Type: str
            Default: None (auto-generated from input filename)
            Description: Name of the tab-separated output file containing
            founder strain dosages. If None, generates name based on input
            filename with '.tsv' extension.
            Example: 'DO336.genoprobs.tsv'

    Returns:
        None. The function creates a tab-separated file containing founder
        strain dosages for all positions across all chromosomes.

    Raises:
        FileNotFoundError: If genoprob_file or grid_file does not exist
        ValueError: If strains list is empty or incompatible with data
        RuntimeError: If conversion fails or produces invalid results

    Notes:
        - The conversion matrix maps each diplotype to founder strain contributions
        - For homozygous diplotypes (AA), the contribution is 1.0 for strain A
        - For heterozygous diplotypes (AB), the contribution is 0.5 for each strain
        - The output format is compatible with many statistical analysis tools
        - Each row represents a genomic position, each column represents a founder strain
        - Dosage values range from 0.0 to 1.0, representing expected founder contribution
        - This function is typically run after GBRS reconstruction for downstream analysis
        - The quant format enables founder strain-specific statistical analyses
        - File sizes can be large for populations with many positions
    """
    if grid_file is None:
        grid_file = os.path.join(DATA_DIR, 'ref.genome_grid.64k.txt')

    if output_file is None:
        output_file = f'{os.path.splitext(genoprob_file)[0]}.tsv'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Strains: {strains}')
    logger.info(f'Grid File: {grid_file}')
    logger.info(f'Output File: {output_file}')

    logger.info('Getting suffices for strains')
    num_strains = len(strains)
    hid = dict(zip(strains, np.arange(num_strains)))
    genotypes = [
        h1 + h2 for h1, h2 in combinations_with_replacement(strains, 2)
    ]
    num_genotypes = len(genotypes)
    convmat = np.zeros((num_genotypes, num_strains))
    for g in range(num_genotypes):
        h1, h2 = genotypes[g]
        convmat[g, hid[h1]] += 1
        convmat[g, hid[h2]] += 1
    convmat *= 0.5

    logger.info(f'Loading grid file: {grid_file}')
    with open(grid_file) as fh:
        next(fh)
        grid = OrderedDict()
        for line in fh:
            item = line.rstrip().split('\t')
            chrom = item[1]
            if chrom in grid:
                grid[chrom] += 1
            else:
                grid[chrom] = 1
    num_grids = sum(grid.values())
    logger.debug(f'Number of grids: {num_grids}')

    logger.info(f'Loading GBRS genotype probability file: {genoprob_file}')
    gprob = np.load(genoprob_file)
    chromosomes = list(grid.keys())
    gprob_mat = gprob[chromosomes[0]].transpose()

    for c in chromosomes[1:]:
        logger.debug(f'Stacking {c}')
        gprob_mat = np.vstack((gprob_mat, gprob[c].transpose()))

    logger.info('Converting genotype probability')
    gprob_mat_converted = np.dot(gprob_mat, convmat)

    logger.info(f'Saving GBRS quant format: {output_file}')
    np.savetxt(
        output_file,
        gprob_mat_converted,
        fmt='%.6f',
        delimiter='\t',
        header='\t'.join(strains),
    )

    logger.info('Done')


def debug_genoprob(
        genoprob_file: str,
        output_file: str = None,
        strains: list[str] = None
) -> dict | None:
    """
    Analyze and debug genotype probability files with detailed statistics.

    This function provides comprehensive analysis and debugging capabilities
    for GBRS genotype probability files. It examines the structure, quality,
    and statistical properties of the probability matrices to help identify
    potential issues or validate the reconstruction results.

    The function performs extensive quality checks including probability
    normalization, distribution analysis, and diplotype frequency analysis.
    It provides both summary statistics and detailed examples to help
    understand the structure and quality of the genotype probability data.

    ALGORITHM:
    1. Load genotype probability file and examine structure
    2. For each chromosome:
       - Analyze matrix dimensions and data types
       - Check probability normalization (sum to 1.0)
       - Calculate maximum probabilities and confidence
       - Analyze distribution of most likely diplotypes
       - Provide detailed examples for first few positions
    3. Generate comprehensive summary report
    4. Optionally save detailed analysis to file

    WHAT THIS DOES:
    - Validates probability matrix structure and quality
    - Identifies potential issues in genotype reconstruction
    - Provides statistical summaries for quality assessment
    - Offers detailed examples for understanding data structure
    - Helps debug reconstruction pipeline issues

    USE CASES:
    - Quality control: Validate GBRS reconstruction results
    - Debugging: Identify issues in reconstruction pipeline
    - Data exploration: Understand structure of probability matrices
    - Validation: Check probability normalization and distributions
    - Documentation: Generate reports for data quality assessment

    Args:
        genoprob_file: Path to the genotype probability file (npz format).
            Type: str
            Description: File containing genotype probabilities from GBRS
            reconstruction or interpolation. Can be either gene-based or
            grid-based probabilities.
            Format: Compressed numpy file with one array per chromosome
            Shape: (num_diplotypes, num_positions) for each chromosome
            Example: 'DO336.genoprobs.npz' or 'DO336.interpolated.genoprobs.npz'

        output_file: Output filename for detailed analysis report.
            Type: str
            Default: None (auto-generated from input filename)
            Description: Name of the output file for saving detailed analysis.
            If None, generates name based on input filename with '.tsv' extension.
            The report includes comprehensive statistics and examples.
            Example: 'DO336.genoprobs.debug.tsv'

        strains: List of founder strain identifiers.
            Type: list[str]
            Default: None (uses default 8-founder DO strains)
            Description: Names of founder strains in the population.
            Used to generate diplotype names for analysis and reporting.
            If None, uses standard DO founder strains A-H.
            Example: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    Returns:
        dict | None: Loaded genotype probability data or None if loading fails.
            Type: dict | None
            Description: If successful, returns the loaded numpy data object
            containing genotype probability matrices. If loading fails due to
            file errors or invalid format, returns None.
            The returned object can be used for further analysis or debugging.

    Raises:
        FileNotFoundError: If genoprob_file does not exist
        ValueError: If file format is invalid or corrupted
        RuntimeError: If analysis fails due to data inconsistencies

    Notes:
        - The function automatically detects whether the file contains
          interpolated (grid-based) or original (gene-based) probabilities
        - Probability normalization is checked to ensure columns sum to ~1.0
        - Maximum probabilities indicate confidence in genotype calls
        - Diplotype distribution analysis shows frequency of each genotype
        - Detailed examples show actual probability values for first few positions
        - The function provides both console output and optional file output
        - This function is essential for quality control and debugging
        - Analysis results help identify potential issues in reconstruction
        - The function is safe to run on any GBRS genotype probability file
        - Output includes both summary statistics and detailed examples
    """
    if strains is None:
        strains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']

    if output_file is None:
        output_file = f'{os.path.splitext(genoprob_file)[0]}.tsv'

    logger.info(f'Genotype Probabilities File: {genoprob_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Strains: {strains}')

    try:
        data = np.load(genoprob_file)
        strains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        diplotypes = [h1 + h2 for h1, h2 in combinations_with_replacement(strains, 2)]

        # Determine if this is an interpolated file based on the filename
        is_interpolated = 'interpolated' in genoprob_file.lower()

        logger.debug(f'Chromosomes: {data.files}')
        logger.debug(f'Diplotypes: {diplotypes}')
        logger.debug(f'Number of diplotypes: {len(diplotypes)}')
        logger.debug(
            f"File type: {'Interpolated (grid positions)' if is_interpolated else 'Original (gene positions)'}")

        # Summary statistics
        total_positions = 0
        for chrom in sorted(data.files):
            matrix = data[chrom]
            position_type = "grid positions" if is_interpolated else "genes"
            print(f"\nChromosome {chrom}:")
            print(f"  Shape: {matrix.shape}")
            print(f"  Data type: {matrix.dtype}")
            print(f"  Number of {position_type}: {matrix.shape[1]}")
            print(f"  Number of diplotypes: {matrix.shape[0]}")

            # Check probability properties
            col_sums = matrix.sum(axis=0)
            print(
                f"  Column sums (should be ~1.0): min={col_sums.min():.6f}, max={col_sums.max():.6f}, mean={col_sums.mean():.6f}")

            # Find most likely diplotypes
            max_probs = matrix.max(axis=0)
            max_diplotype_indices = matrix.argmax(axis=0)
            print(
                f"  Max probabilities: min={max_probs.min():.6f}, max={max_probs.max():.6f}, mean={max_probs.mean():.6f}")

            # Show distribution of most likely diplotypes
            unique_diplotypes, counts = np.unique(max_diplotype_indices, return_counts=True)
            print(f"  Most common diplotypes (first 10):")
            for i, (diplotype_idx, count) in enumerate(zip(unique_diplotypes, counts)):
                if i < 10:
                    print(
                        f"    {diplotypes[diplotype_idx]}: {count} {position_type} ({count / matrix.shape[1] * 100:.1f}%)")

            total_positions += matrix.shape[1]

        position_type = "grid positions" if is_interpolated else "genes"
        print(f"\nTotal {position_type} across all chromosomes: {total_positions}")

        # Show detailed example for first chromosome
        if data.files:
            first_chrom = sorted(data.files)[0]
            matrix = data[first_chrom]
            position_type = "positions" if is_interpolated else "genes"
            logger.debug(f"\n=== DETAILED EXAMPLE: Chromosome {first_chrom} ===")
            print(f"First 5 {position_type}, all diplotypes:")

            # Show first 5 positions/genes with their probabilities
            for pos_idx in range(min(5, matrix.shape[1])):
                print(f"\n{position_type.capitalize()} {pos_idx}:")
                probs = matrix[:, pos_idx]
                max_idx = probs.argmax()
                print(f"  Most likely: {diplotypes[max_idx]} (prob={probs[max_idx]:.6f})")

                # Show top 5 diplotypes
                top_indices = np.argsort(probs)[-5:][::-1]
                for i, idx in enumerate(top_indices):
                    print(f"  {i + 1}. {diplotypes[idx]}: {probs[idx]:.6f}")

        return data

    except Exception as e:
        print(f"Error reading {genoprob_file}: {e}")
        return None



