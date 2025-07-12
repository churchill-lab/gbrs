# standard library imports
import gzip
import logging
import os
import string
import subprocess
from collections import defaultdict
from itertools import dropwhile

# 3rd party library imports
import numpy as np
import tables
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.sparse import csc_matrix

# local library imports
from gbrs import utils
from gbrs.emase.AlignmentMatrixFactory import AlignmentMatrixFactory
from gbrs.emase.AlignmentPropertyMatrix import AlignmentPropertyMatrix
from gbrs.emase.EMfactory import EMfactory

logging.getLogger('numexpr').setLevel(logging.WARNING)
logger = utils.get_logger('gbrs')


def bam2emase(
        alignment_files: list[str],
        haplotypes: list[str],
        locusid_file: str,
        output_file: str = 'alignments.transcriptome.h5',
        delim: str = '_',
        index_dtype: str = 'uint32',
        data_dtype: str = 'uint8'
) -> None:
    """
    Convert BAM file to EMASE format (HDF5) for allele-specific expression analysis.

    This function is a core component of the GBRS pipeline that transforms standard
    RNA-Seq alignment data (BAM format) into the specialized EMASE format required
    for allele-specific expression analysis in multiparent populations.

    The function processes BAM files where reads have been aligned to a hybrid
    transcriptome containing transcripts from all founder strains. Reference names
    in the BAM files must follow the pattern: {transcript_id}{delim}{haplotype}
    (e.g., "ENSMUST00000000001_A" for transcript ENSMUST00000000001 from haplotype A).

    The output is an HDF5 file containing sparse matrices that efficiently represent
    the alignment relationships between paired reads, loci (transcripts), and haplotypes.

    Args:
        alignment_files: List of paths to input BAM files containing paired-end RNA-Seq alignments.
            Typically contains two files: one for R1 reads and one for R2 reads.
            Example: ['sample_R1.bam', 'sample_R2.bam']
            All BAM files must be aligned to the same hybrid transcriptome.

        haplotypes: List of haplotype identifiers representing founder strains.
            For example, Diversity Outbred mice: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
            These identifiers must match those used in the BAM reference names.

        locusid_file: Path to the transcript/locus information file (TSV format).
            Expected format: tab-separated values with transcript IDs in the first column.
            Example:
                ENSMUST00000000001    0
                ENSMUST00000000002    0
                ENSMUST00000000003    0
            The second column is typically 0 and used for additional metadata.

        output_file: Path for the output EMASE file (HDF5 format).
            Default: 'alignments.transcriptome.h5'
            The file will contain sparse matrices organized by haplotype.

        delim: Delimiter string between transcript ID and haplotype in BAM reference names.
            Default: '_'
            Example: If BAM references are "ENSMUST00000000001_A", use delim='_'

        index_dtype: Data type for matrix indices in the output HDF5 file.
            Default: 'uint32'
            Options: 'uint16', 'uint32', 'uint64'
            Choose based on the scale of your dataset:
            - uint16: Up to 65,535 loci/reads (small datasets)
            - uint32: Up to 4.3 billion loci/reads (most datasets)
            - uint64: For extremely large datasets

        data_dtype: Data type for matrix values in the output HDF5 file.
            Default: 'uint8'
            Options: 'uint8', 'uint16', 'uint32', 'float32'
            - uint8: Binary presence/absence (most memory efficient)
            - uint16/uint32: For read counts or alignment scores
            - float32: For alignment probabilities or weights

    Raises:
        FileNotFoundError: If any alignment_file or locusid_file does not exist.
        ValueError: If haplotypes list is empty or contains invalid characters.
        OSError: If output_file cannot be created due to permission issues.

    Output Format:
        The function generates an HDF5 file with the following structure:

        / (root)
        ├── @incidence_only: True/False (based on data_dtype)
        ├── @mtype: 'csc_matrix'
        ├── @shape: (num_loci, num_haplotypes, num_reads)
        ├── @hname: ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
        ├── /lname: Array of locus/transcript names
        ├── /rname: Array of read names (combined from all input files)
        ├── /h0: Sparse matrix for haplotype 0 (A)
        │   ├── /indptr: CSC matrix indptr array
        │   ├── /indices: CSC matrix indices array
        │   └── /data: CSC matrix data array (if incidence_only=False)
        ├── /h1: Sparse matrix for haplotype 1 (B)
        │   └── ...
        └── /hN: Sparse matrix for haplotype N
            └── ...

    Implementation Details:
        1. Parses the locus ID file to extract transcript names
        2. Creates a PairedAlignmentMatrixFactory to handle multiple BAM files
        3. Extracts all unique read names from fits BAM file (assumes all reads are the same).
        4. Creates temporary binary files for each haplotype and BAM file.
        5. Processes alignments from all files and builds sparse matrices
        6. Saves data in compressed HDF5 format
        7. Cleans up temporary files

    Notes:
        - Output files are typically smaller than the combined input BAM files due to compression
    """
    logger.info(f'BAM Files: {alignment_files}')
    logger.info(f'Locus ID File: {locusid_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Haplotypes: {haplotypes}')
    logger.info(f'Delimiter: {delim}')
    logger.info(f'Index dtype: {index_dtype}')
    logger.info(f'Data dtype: {data_dtype}')

    logger.info(f'Parsing Locus ID File: {locusid_file}')
    loci = utils.get_names(locusid_file)

    logger.info(f'Parsing BAM Files: {alignment_files}')
    amf = AlignmentMatrixFactory(alignment_files)
    amf.prepare(haplotypes, loci, delim=delim, out_dir=os.path.dirname(output_file))

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    amf.produce(output_file, index_dtype=index_dtype, data_dtype=data_dtype)
    amf.cleanup()
    logger.info('Done')


def combine(
        emase_files: list[str],
        output_file: str,
        comp_lib: str = 'zlib'
) -> None:
    """
    Combine multiple EMASE files into a single file by concatenating read data.

    This function merges multiple EMASE files by concatenating their read data
    along the read dimension. It performs a UNION operation on read data,
    combining all reads from all input files into a single EMASE file. This is
    different from `get_common_alignments()` which performs an INTERSECTION
    operation.

    ALGORITHM:
    1. Load the first EMASE file as the base matrix
    2. For each subsequent file:
       - Load the EMASE file
       - Verify compatibility (same loci, haplotypes, structure)
       - Concatenate read data: combined_reads = [reads1, reads2, reads3, ...]
    3. Save the combined result

    WHAT THIS DOES:
    - Takes multiple EMASE files (can have different reads)
    - Combines all read data into a single file
    - Preserves all alignment information from all files
    - Result: Single EMASE file with all reads from all input files
    - Output: EMASE file with same structure but more reads

    USE CASES:
    - Sample pooling: Combine multiple biological samples
    - Replicate merging: Combine technical or biological replicates
    - Batch processing: Merge files from different sequencing runs
    - Data consolidation: Create single file for downstream analysis

    DIFFERENCE FROM OTHER FUNCTIONS:
    - combine(): UNION operation - combines all reads from all files
    - get_common_alignments(): INTERSECTION operation - only keeps reads that align
      consistently across all files
    - compress(): Groups identical alignment patterns (equivalence classes)

    REQUIREMENTS:
    - All input files must have identical loci (lname) and haplotypes (hname)
    - All files must have the same matrix structure and dimensions
    - Files can have different reads (unlike get_common_alignments)

    Args:
        emase_files: List of EMASE file paths to combine. Files must have
            compatible structure (same loci and haplotypes).

        output_file: Path for the output combined EMASE file. The file will
            contain all reads from all input files.

        comp_lib: Compression library to use for the output HDF5 file.
            Default: 'zlib' (good balance of compression and speed)
            Options: 'zlib', 'lzo', 'bzip2', 'blosc'

    Returns:
        None. The function creates a combined EMASE file at the specified output_file path.

    Raises:
        FileNotFoundError: If any EMASE file does not exist.
        ValueError: If files have incompatible structure (different loci/haplotypes).
        RuntimeError: If matrices are not finalized or other processing errors.

    Implementation Details:
        1. Loads each EMASE file using AlignmentPropertyMatrix
        2. Validates compatibility between files (loci, haplotypes)
        3. Uses the combine() method to concatenate matrices along read dimension
        4. Preserves metadata from the first file (loci, haplotypes)
        5. Concatenates read names from all files
        6. Saves the combined result with specified compression
    """
    for f in emase_files:
        logger.info(f'EMASE file: {f}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    amf_list = list()
    for f in emase_files:
        logger.info(f'Loading EMASE file: {f}')
        amf = AlignmentPropertyMatrix(h5_file=f)
        logger.debug(f'Number Loci: {amf.num_loci}')
        logger.debug(f'Number Haplotypes: {amf.num_haplotypes}')
        logger.debug(f'Number Reads: {amf.num_reads}')
        amf_list.append(amf)

    amf_final = amf_list[0].copy()
    for i, amf in enumerate(amf_list[1:]):
        logger.info(f'Combining EMASE file: {emase_files[i]}')
        amf_final = amf_final.combine(amf)
        logger.debug(f'Combined Number Loci: {amf_final.num_loci}')
        logger.debug(f'CombinedNumber Haplotypes: {amf_final.num_haplotypes}')
        logger.debug(f'Combined Number Reads: {amf_final.num_reads}')

    logger.info(f'Saving EMASE file {output_file}')
    amf_final.save(h5_file=output_file, complib=comp_lib)
    logger.info('Dome')


def count_alignments(
        alignment_file: str,
        group_file: str,
        outbase: str = 'emase'
) -> None:
    """
    Count the number of alignments for each locus and gene in EMASE format data.

    This function analyzes EMASE alignment data to generate alignment count reports
    at both the transcript (isoform) and gene levels. It is part of the GBRS pipeline
    for quantifying allele-specific expression in multiparent populations.

    The function processes an EMASE file containing sparse matrices of read alignments
    to transcripts from different founder haplotypes. It generates two output files:
    1. Isoform-level alignment counts: counts per transcript
    2. Gene-level alignment counts: aggregated counts per gene

    ALGORITHM:
    1. Load the EMASE file with group information
    2. Generate isoform-level alignment counts
    3. Bundle data to save memory
    4. Generate gene-level alignment counts

    OUTPUT FILES:
    - {outbase}.isoforms.alignment_counts: Alignment counts per transcript
    - {outbase}.genes.alignment_counts: Alignment counts per gene

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing alignment data.
            The file should contain sparse matrices representing read alignments to
            transcripts from different founder haplotypes.

        group_file: Path to the group file containing transcript-to-gene mapping.
            This file enables aggregation of transcript-level counts to gene-level counts.
            Format: tab-separated with transcript ID in first column, gene ID in second.

        outbase: Base name for output files. Default: 'emase'
            Output files will be named:
            - {outbase}.isoforms.alignment_counts
            - {outbase}.genes.alignment_counts

    Returns:
        None. The function creates alignment count files at the specified outbase path.

    Raises:
        FileNotFoundError: If alignment_file or group_file does not exist.
        ValueError: If the EMASE file format is invalid or incompatible.
        RuntimeError: If the AlignmentPropertyMatrix cannot be loaded or processed.

    Notes:
        - The function uses memory optimization by bundling inline data between
          isoform and gene-level processing.
        - Gene-level counts are aggregated from transcript-level counts using
          the group file mapping.
        - This function is typically used after bam2emase() to analyze alignment
          patterns in the processed data.
    """
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Outbase: {outbase}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    amf = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)
    logger.debug(f'Number Loci: {amf.num_loci}')
    logger.debug(f'Number Haplotypes: {amf.num_haplotypes}')
    logger.debug(f'Number Reads: {amf.num_reads}')

    tmp_filename = f'{outbase}.isoforms.alignment_counts'
    logger.info(f'Generating isoform Alignment Counts: {tmp_filename}')
    amf.report_alignment_counts(filename=tmp_filename)

    # bundle the inline data to save memory
    amf._bundle_inline(reset=True)

    tmp_filename = f'{outbase}.genes.alignment_counts'
    logger.info(f'Generating gene Alignment Counts: {tmp_filename}')
    amf.report_alignment_counts(filename=tmp_filename)

    logger.info('Done')


def get_num_shared_multireads(amf: AlignmentPropertyMatrix) -> np.ndarray:
    """
    Calculate the number of shared multireads between all pairs of loci.

    This function computes a pairwise similarity matrix that counts how many reads
    are shared between each pair of loci (transcripts/genes) across all haplotypes.
    It is used for analyzing read distribution patterns and identifying loci with
    similar alignment profiles.

    ALGORITHM:
    1. Sum alignment matrix across haplotype dimension to get total alignments per locus
    2. Convert to binary matrix (presence/absence of alignments)
    3. Compute pairwise similarity: cnt_mat = haplotype_sum.T * haplotype_sum
    4. Result: symmetric matrix where cnt_mat[i,j] = number of reads shared between loci i and j

    WHAT THIS DOES:
    - Takes an AlignmentPropertyMatrix with dimensions (loci, haplotypes, reads)
    - Collapses haplotype dimension to get locus-read relationships
    - Computes pairwise similarity between all loci based on shared reads
    - Returns symmetric matrix of shared read counts

    Args:
        amf: AlignmentPropertyMatrix object containing alignment data.
            Should have dimensions (num_loci, num_haplotypes, num_reads).
            The matrix represents read alignments to transcripts from different haplotypes.

    Returns:
        np.ndarray: Symmetric matrix of shape (num_loci, num_loci) where each element
            [i,j] represents the number of reads shared between loci i and j.
            Diagonal elements represent the total number of reads aligned to each locus.

    Notes:
        - The result is a symmetric matrix: cnt_mat[i,j] = cnt_mat[j,i]
        - Diagonal elements give the total number of reads per locus
        - Off-diagonal elements give the number of reads shared between locus pairs
        - This function is typically used for analyzing read distribution patterns
          and identifying loci with similar alignment profiles.
    """
    haplotype_sum = amf.sum(axis=AlignmentPropertyMatrix.Axis.HAPLOTYPE)
    haplotype_sum.data = np.ones(haplotype_sum.nnz)
    cnt_mat = haplotype_sum.transpose() * haplotype_sum
    return cnt_mat


def count_shared_multireads_pairwise(
        alignment_file: str,
        group_file: str,
        outbase: str = 'emase'
) -> None:
    """
    Count shared multireads between all pairs of loci at both transcript and gene levels.

    This function analyzes EMASE alignment data to identify and quantify reads that
    align to multiple loci (transcripts or genes). It generates pairwise similarity
    matrices that show how many reads are shared between each pair of loci, which is
    useful for understanding read distribution patterns and identifying loci with
    similar expression profiles.

    The function processes alignment data at two levels:
    1. Transcript (isoform) level: analyzes shared reads between individual transcripts
    2. Gene level: analyzes shared reads between genes (aggregated from transcripts)

    ALGORITHM:
    1. Load EMASE file with group information
    2. Generate isoform-level shared read count matrix
    3. Bundle data to save memory
    4. Generate gene-level shared read count matrix

    OUTPUT FILES:
    - {outbase}.isoforms.shared_read_counts.npz: Compressed numpy array containing
      pairwise shared read counts between transcripts
    - {outbase}.genes.shared_read_counts.npz: Compressed numpy array containing
      pairwise shared read counts between genes

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing alignment data.
            The file should contain sparse matrices representing read alignments to
            transcripts from different founder haplotypes.

        group_file: Path to the group file containing transcript-to-gene mapping.
            This file enables aggregation of transcript-level analysis to gene-level analysis.
            Format: tab-separated with transcript ID in first column, gene ID in second.

        outbase: Base name for output files. Default: 'emase'
            Output files will be named:
            - {outbase}.isoforms.shared_read_counts.npz
            - {outbase}.genes.shared_read_counts.npz

    Returns:
        None. The function creates compressed numpy files containing shared read count matrices.

    Raises:
        FileNotFoundError: If alignment_file or group_file does not exist.
        ValueError: If the EMASE file format is invalid or incompatible.
        RuntimeError: If the AlignmentPropertyMatrix cannot be loaded or processed.

    Notes:
        - The function uses memory optimization by bundling inline data between
          isoform and gene-level processing.
        - Output files are compressed numpy arrays (.npz format) for efficient storage.
        - Each output file contains a symmetric matrix where element [i,j] represents
          the number of reads shared between loci i and j.
        - This function is useful for analyzing read distribution patterns and
          identifying loci with similar alignment profiles.
    """
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Outbase: {outbase}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    aln_mat = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)
    logger.debug(f'Number Loci: {aln_mat.num_loci}')
    logger.debug(f'Number Haplotypes: {aln_mat.num_haplotypes}')
    logger.debug(f'Number Reads: {aln_mat.num_reads}')

    outfile1 = f'{outbase}.isoforms.shared_read_counts'
    logger.info(f'Generating isoform Shared Read Counts: {outfile1}')
    cnt_mat = get_num_shared_multireads(aln_mat)
    np.savez_compressed(outfile1, counts=cnt_mat)

    aln_mat._bundle_inline(reset=True)

    outfile2 = f'{outbase}.isoforms.shared_read_counts'
    logger.info(f'Generating genes Shared Read Counts: {outfile2}')
    cnt_mat = get_num_shared_multireads(aln_mat)
    np.savez_compressed(outfile2, counts=cnt_mat)

    logger.info('Done')


def create_hybrid(
        fasta_list: list[str],
        haplotypes: list[str],
        output_file: str = 'emase.pooled.targets.fa',
        build_bowtie_index: bool = False
) -> None:
    """
    Create a hybrid transcriptome by combining FASTA files from multiple founder haplotypes.

    This function is a core component of the GBRS pipeline that creates a pooled
    transcriptome for RNA-Seq alignment. It combines transcript sequences from
    multiple founder strains (haplotypes) into a single FASTA file, with each
    transcript ID suffixed by its haplotype identifier.

    The function processes multiple FASTA files representing transcriptomes from
    different founder strains (e.g., Diversity Outbred mouse founder strains A-H).
    Each transcript in the output file is uniquely identified by appending the
    haplotype suffix to the original transcript ID.

    ALGORITHM:
    1. Create output directory if it doesn't exist
    2. For each haplotype FASTA file:
       - Read transcript sequences
       - Append haplotype suffix to transcript IDs
       - Write to pooled FASTA file
       - Record transcript lengths
    3. Optionally build Bowtie index for alignment

    OUTPUT FILES:
    - {output_file}: Pooled FASTA file containing all transcripts from all haplotypes
    - {outbase}.info: Transcript length information file
    - {outbase}.bowtie1.*: Bowtie index files (if build_bowtie_index=True)

    Args:
        fasta_list: List of paths to input FASTA files, one per founder haplotype.
            Each file should contain transcript sequences from a single founder strain.
            Example: ['strain_A.fa', 'strain_B.fa', 'strain_C.fa']

        haplotypes: List of haplotype identifiers corresponding to each FASTA file.
            These identifiers will be appended to transcript IDs in the output.
            Example: ['A', 'B', 'C'] for Diversity Outbred founder strains
            Must have same length as fasta_list.

        output_file: Path for the output pooled FASTA file. Default: 'emase.pooled.targets.fa'
            The file will contain all transcripts from all haplotypes with unique IDs.

        build_bowtie_index: Whether to build a Bowtie index for the pooled transcriptome.
            Default: False
            If True, creates Bowtie index files for RNA-Seq read alignment.

    Returns:
        None. The function creates the pooled FASTA file and optional index files.

    Raises:
        FileNotFoundError: If any FASTA file does not exist.
        ValueError: If fasta_list and haplotypes have different lengths.
        OSError: If output directory cannot be created or files cannot be written.
        subprocess.CalledProcessError: If Bowtie index building fails.

    Notes:
        - Transcript IDs in the output follow the pattern: {original_id}_{haplotype}
        - The function creates a transcript length file for downstream processing
        - Bowtie index building can be time-consuming for large transcriptomes
        - This function is typically used as the first step in the GBRS pipeline
          to prepare the reference transcriptome for RNA-Seq alignment.
    """
    out_dir = os.path.dirname(output_file)
    if out_dir != '' and not os.path.exists(out_dir):
        os.mkdir(out_dir)

    for x in fasta_list:
        logger.info(f'Fasta File: {x}')
    logger.info(f'Haplotype List: {haplotypes}')
    logger.info(f'Output File: {output_file}')

    # Get pooled transcriptome
    outbase = os.path.splitext(output_file)[0]
    num_haps = len(fasta_list)
    lenfile = f'{outbase}.info'
    seqout = open(output_file, 'w')
    lenout = open(lenfile, 'w')

    logger.debug('Looping through haplotypes')
    for hid in range(num_haps):
        fasta = fasta_list[hid]
        hapname = haplotypes[hid]
        logger.info(
            f'Adding suffix "_{hapname}" to the sequence ID\'s of {fasta}'
        )
        fh = open(fasta)
        line = fh.readline()  # the first fasta header
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
        logger.info('Building bowtie1 index (could take some time)')
        subprocess.call(f'bowtie-build {output_file} {out_index}', shell=True)

    logger.info('Done')


def get_common_alignments(
        emase_files: list[str],
        output_file: str = None,
        comp_lib: str = 'zlib'
) -> None:
    """
    Find reads that align to the same loci across multiple EMASE files.

    This function performs ELEMENT-WISE MULTIPLICATION of sparse matrices to
    find reads that have identical alignment patterns across all input files.

    The function performs an INTERSECTION operation on alignment data - it
    only keeps reads that align consistently across all input files. This is
    different from `combine()` which performs a UNION operation.

    ALGORITHM:
    1. Load the first EMASE file as the base matrix
    2. For each subsequent file:
       - Load the EMASE file
       - Verify read IDs are identical across files
       - Perform element-wise multiplication: aln_mat = aln_mat * aln_mat_next
    3. Save the result

    WHAT THIS DOES:
    - Takes multiple EMASE files with the SAME reads
    - For each read, only keeps alignments that exist in ALL files
    - Result: Reads that align to the same loci in all input files
    - Output: EMASE file with same structure but fewer non-zero elements

    DIFFERENCE FROM OTHER FUNCTIONS:
    - combine(): UNION operation - combines all reads from all files
    - compress(): Groups identical alignment patterns (equivalence classes)

    REQUIREMENTS:
    - All input files must have identical read IDs (rname)
    - All files must have identical loci (lname) and haplotypes (hname)
    - All files must have the same matrix structure and dimensions
    - Files should represent same sample

    Args:
        emase_files: List of EMASE file paths to process. All files must have
            identical read IDs and compatible structure (same loci and
            haplotypes).

        output_file: Path for the output common alignments EMASE file. If None,
            auto-generates filename based on first input file.

        comp_lib: Compression library to use for the output HDF5 file.
            Options: 'zlib', 'lzo', 'bzip2', 'blosc'

    Raises:
        FileNotFoundError: If any EMASE file does not exist.
        ValueError: If read IDs are not identical across files.
        RuntimeError: If matrices are not finalized or other processing errors.

    Notes:
        Only alignments that exist in ALL files are retained in the output.
        For combining different samples or replicates, use combine() instead.
        For grouping identical alignment patterns, use compress() instead.
    """
    if output_file is None:
        output_file = f'alignments.common.{os.path.basename(emase_files[0])}'

    for f in emase_files:
        logger.info(f'EMASE file: {f}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')

    logger.info(f'Loading EMASE file: {emase_files[0]}')
    apm = AlignmentPropertyMatrix(h5_file=emase_files[0])
    logger.debug(f'Number Loci: {apm.num_loci}')
    logger.debug(f'Number Haplotypes: {apm.num_haplotypes}')
    logger.debug(f'Number Reads: {apm.num_reads}')

    for f in emase_files[1:]:
        logger.info(f'Loading EMASE file: {f}')
        apm_next = AlignmentPropertyMatrix(h5_file=f)
        logger.debug(f'Number Loci: {apm_next.num_loci}')
        logger.debug(f'Number Haplotypes: {apm_next.num_haplotypes}')
        logger.debug(f'Number Reads: {apm_next.num_reads}')

        if np.all(apm.rname == apm_next.rname):
            apm = apm * apm_next
        else:
            logger.error('The read ID\'s are not compatible.')
            raise ValueError('The read ID\'s are not compatible.')

        logger.debug(f'Total Number Loci: {apm.num_loci}')
        logger.debug(f'Total Number Haplotypes: {apm.num_haplotypes}')
        logger.debug(f'Total Number Reads: {apm.num_reads}')

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    apm.save(h5_file=output_file, complib=comp_lib)
    logger.info('Done')


def process_haplotype_optimized(
        haplotype_id: int,
        emase_files: list[str],
        num_loci: int,
        num_reads: int,
) -> (int, np.ndarray, np.ndarray, np.ndarray):
    """
    Process a single haplotype across all EMASE files.

    Args:
        haplotype_id: Index of the haplotype to process
        emase_files: List of EMASE files
        num_loci: Number of loci
        num_reads: Number of reads

    Returns:
        Tuple of (haplotype_id, indices, indptr, data) for the processed haplotype
    """
    logger.debug(f'Processing haplotype {haplotype_id}')

    # load first file's haplotype data
    with tables.open_file(emase_files[0], 'r') as f:
        hap_node = f.get_node(f'/h{haplotype_id}')
        indices = hap_node.indices.read()
        indptr = hap_node.indptr.read()
        # for incidence matrices, data is all ones
        data = np.ones(len(indices), dtype=np.float64)

    # create CSC matrix for first file
    current_matrix = csc_matrix((data, indices, indptr), shape=(num_reads, num_loci))

    # process remaining files
    for file_idx, emase_file in enumerate(emase_files[1:], 1):
        logger.debug(f'Processing file {file_idx + 1}/{len(emase_files)} for haplotype {haplotype_id}')

        with tables.open_file(emase_file, 'r') as f:
            hap_node = f.get_node(f'/h{haplotype_id}')
            next_indices = hap_node.indices.read()
            next_indptr = hap_node.indptr.read()
            next_data = np.ones(len(next_indices), dtype=np.float64)

        next_matrix = csc_matrix((next_data, next_indices, next_indptr), shape=(num_reads, num_loci))

        # wlement-wise multiplication (common alignments)
        current_matrix = current_matrix.multiply(next_matrix)

    # extract final result
    final_indices = current_matrix.indices
    final_indptr = current_matrix.indptr
    final_data = current_matrix.data

    return haplotype_id, final_indices, final_indptr, final_data


def save_optimized_result(
        haplotype_results: list,
        output_file: str,
        lname: list,
        rname: list,
        num_loci: int,
        num_haplotypes: int,
        num_reads: int,
        comp_lib: str
) -> None:
    """
    Save the optimized result to HDF5 file.

    Args:
        haplotype_results: List of (haplotype_id, indices, indptr, data) tuples
        output_file: Output file path
        lname: List of locus names
        rname: List of read names
        num_loci: Number of loci
        num_haplotypes: Number of haplotypes
        num_reads: Number of reads
        comp_lib: Compression library
    """
    # sort results by haplotype_id
    haplotype_results.sort(key=lambda x: x[0])

    with tables.open_file(output_file, 'w') as f:
        # set root attributes (matching original format exactly)
        f.set_node_attr('/', 'shape', (num_loci, num_haplotypes, num_reads))
        f.set_node_attr('/', 'mtype', 'csc_matrix')  # String, not bytes
        f.set_node_attr('/', 'incidence_only', True)

        # set hname attribute
        hname = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'][:num_haplotypes]
        f.set_node_attr('/', 'hname', hname)

        # save locus names using create_carray (matching original)
        lname_array = np.array([name.encode() for name in lname])
        fil = tables.Filters(complevel=1, complib=comp_lib)
        f.create_carray('/', 'lname', obj=lname_array, title='Locus Names', filters=fil)

        # save read names using create_carray (matching original)
        rname_array = np.array([name.encode() for name in rname])
        f.create_carray('/', 'rname', obj=rname_array, title='Read Names', filters=fil)

        # save haplotype data
        for haplotype_id, indices, indptr, data in haplotype_results:
            logger.debug(f'Saving haplotype {haplotype_id}')
            hap_group = f.create_group('/', f'h{haplotype_id}')

            # ensure proper data types
            indices = np.asarray(indices, dtype=np.uint32)
            indptr = np.asarray(indptr, dtype=np.uint32)

            # save indices with compression
            f.create_carray(
                hap_group, 'indices',
                tables.UInt32Atom(),
                indices.shape,
                obj=indices,
                filters=tables.Filters(complevel=6, complib=comp_lib)
            )

            # save indptr with compression
            f.create_carray(
                hap_group, 'indptr',
                tables.UInt32Atom(),
                indptr.shape,
                obj=indptr,
                filters=tables.Filters(complevel=6, complib=comp_lib)
            )


def get_common_alignments_optimized(
        emase_files: list[str],
        output_file: str = None,
        comp_lib: str = 'zlib',
        validate: bool = True
) -> None:
    """
    Optimized version of get_common_alignments that finds reads with
    identical alignment patterns across multiple EMASE files.

    This function performs the same operation as get_common_alignments
    but with significant memory and performance optimizations. It finds
    reads that align to the same loci across all input files using
    element-wise multiplication of sparse matrices, but processes data
    haplotype-by-haplotype to reduce memory usage.

    ALGORITHM:
    1. Load metadata (shape, locus names, read names) from the first file
    2. Validate all files have identical dimensions and read IDs
    3. Process each haplotype separately:
       - Load haplotype data from all files using direct HDF5 access
       - Perform element-wise multiplication of sparse matrices
       - Store results for each haplotype
    4. Combine all haplotype results and save to output file

    DIFFERENCES FROM ORIGINAL get_common_alignments:

    MEMORY OPTIMIZATION:
       - Original: Loads full AlignmentPropertyMatrix objects for all files simultaneously
       - Optimized: Processes one haplotype at a time, using direct HDF5 access
       - Result: Significantly lower peak memory usage, especially for large datasets

    DATA ACCESS PATTERN:
       - Original: Uses AlignmentPropertyMatrix class with full object overhead
       - Optimized: Direct HDF5 file access using PyTables for minimal memory footprint
       - Result: Faster data loading and reduced memory allocations

    PROCESSING:
       - Original: Loads all files into memory, then performs matrix multiplication
       - Optimized: Processes haplotype-by-haplotype, allowing better memory management
       - Result: Can handle larger datasets that would exceed available memory

    4. IMPLEMENTATION:
       - Original: Uses high-level AlignmentPropertyMatrix operations
       - Optimized: Uses low-level scipy.sparse operations for better performance
       - Result: More efficient sparse matrix operations

    5. ERROR HANDLING:
       - Original: Validates read IDs during AlignmentPropertyMatrix loading
       - Optimized: Validates dimensions and read IDs before processing (added validation)
       - Result: Better error detection while maintaining performance benefits

    Args:
        emase_files: List of EMASE file paths to process. All files must have
            identical read IDs and compatible structure (same loci and
            haplotypes).

        output_file: Path for the output common alignments EMASE file. If None,
            auto-generates filename based on first input file.

        comp_lib: Compression library to use for the output HDF5 file.
            Options: 'zlib', 'lzo', 'bzip2', 'blosc'

        validate: Whether to validate file shapes and read IDs
    """
    if output_file is None:
        output_file = f'alignments.common.optimized.{os.path.basename(emase_files[0])}'

    for x in emase_files:
        logger.info(f'EMASE file: {x}')

    logger.info(f'Output File: {output_file}')
    logger.info(f'Compression Library: {comp_lib}')
    logger.info('Loading metadata from first file')

    # get information from first file
    with tables.open_file(emase_files[0], 'r') as f:
        shape = f.get_node_attr('/', 'shape')
        num_loci, num_haplotypes, num_reads = shape
        logger.info(f'Matrix shape: {shape}')
        logger.info(f'Number Loci: {num_loci}')
        logger.info(f'Number Haplotypes: {num_haplotypes}')
        logger.info(f'Number Reads: {num_reads}')

        # load metadata (locus names, read names) from first file
        lname = f.get_node('/', 'lname').read()
        rname = f.get_node('/', 'rname').read()
        # Convert from bytes to string
        lname = [x.decode() for x in lname]
        rname = [x.decode() for x in rname]

    # VALIDATION: Check all files have identical dimensions and read IDs
    if validate:
        logger.info('Validating file compatibility...')
        for file_idx, emase_file in enumerate(emase_files[1:], 1):
            logger.debug(f'Validating file {file_idx + 1}/{len(emase_files)}: {emase_file}')

            with tables.open_file(emase_file, 'r') as f:
                # Check shape/dimensions
                file_shape = f.get_node_attr('/', 'shape')
                if file_shape != shape:
                    error_msg = (
                        f'File {emase_file} has incompatible dimensions. '
                        f'Expected shape {shape}, got {file_shape}. '
                        f'All files must have identical number of loci, haplotypes, and reads.'
                    )
                    logger.error(error_msg)
                    raise ValueError(error_msg)

                # Check read IDs
                file_rname = f.get_node('/', 'rname').read()
                file_rname = [x.decode() for x in file_rname]
                if file_rname != rname:
                    error_msg = (
                        f'File {emase_file} has incompatible read IDs. '
                        f'Read IDs must be identical across all files for common alignments analysis. '
                        f'First file has {len(rname)} reads, this file has {len(file_rname)} reads.'
                    )
                    logger.error(error_msg)
                    raise ValueError(error_msg)

                # Check haplotype names (optional but good for consistency)
                try:
                    file_hname = f.get_node_attr('/', 'hname')
                    if file_idx == 1:  # Get hname from first file for comparison
                        with tables.open_file(emase_files[0], 'r') as f0:
                            hname = f0.get_node_attr('/', 'hname')
                    if file_hname != hname:
                        logger.warning(
                            f'File {emase_file} has different haplotype names than first file. '
                            f'This may indicate different haplotype ordering.'
                        )
                except (AttributeError, tables.exceptions.NoSuchNodeError):
                    logger.debug(f'File {emase_file} does not have haplotype names attribute')

        logger.info('All files validated successfully')
    else:
        logger.warning('WARNING: Not performing file validation for compatibility check')

    # process each haplotype
    haplotype_results = []
    for haplotype_id in range(num_haplotypes):
        result = process_haplotype_optimized(haplotype_id, emase_files, num_loci, num_reads)
        haplotype_results.append(result)

    # combine results and save
    logger.debug('Combining haplotype results and saving')
    logger.info(f'Saving EMASE Formatted File: {output_file}')

    save_optimized_result(
        haplotype_results, output_file, lname, rname,
        num_loci, num_haplotypes, num_reads, comp_lib
    )

    logger.info('Done')


def pull_out_unique_reads(
        alignment_file: str,
        output_file: str,
        group_file: str = None,
        shallow: bool = False,
        ignore_alleles: bool = False
) -> None:
    """
    Extract unique reads from EMASE alignment data, optionally aggregating by gene groups.

    This function filters EMASE alignment data to retain only unique reads, which can
    be useful for reducing redundancy and focusing analysis on distinct read sequences.
    It can operate at both transcript and gene levels, depending on whether a group
    file is provided.

    The function processes EMASE alignment data to identify and extract reads that
    have unique alignment patterns. When a group file is provided, it first aggregates
    alignments by gene groups before identifying unique reads.

    ALGORITHM:
    1. Load EMASE file with optional group information
    2. If group file provided:
       - Bundle alignments by gene groups
       - Identify unique reads at gene level
       - Extract corresponding transcript-level alignments
    3. If no group file:
       - Identify unique reads directly at transcript level
    4. Save filtered alignment data

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing alignment data.
            The file should contain sparse matrices representing read alignments to
            transcripts from different founder haplotypes.

        output_file: Path for the output EMASE file containing only unique reads.
            The output will have the same structure as the input but with fewer reads.

        group_file: Path to the group file containing transcript-to-gene mapping.
            If provided, unique reads are identified at the gene level and then
            mapped back to transcript-level alignments.
            Format: tab-separated with transcript ID in first column, gene ID in second.

        shallow: Whether to use shallow processing mode. Default: False
            If True, uses memory-efficient processing that may be slower.
            If False, loads full data into memory for faster processing.

        ignore_alleles: Whether to ignore haplotype information when identifying unique reads.
            If True, reads are considered unique based only on locus alignment,
            ignoring which haplotype they align to.
            If False, haplotype information is considered in uniqueness determination.

    Returns:
        None. The function creates an EMASE file containing only unique reads.

    Raises:
        FileNotFoundError: If alignment_file or group_file does not exist.
        ValueError: If the EMASE file format is invalid or incompatible.
        RuntimeError: If the AlignmentPropertyMatrix cannot be loaded or processed.

    Notes:
        - The function can significantly reduce the number of reads in the output file
        - When ignore_alleles=True, reads aligning to the same locus across different
          haplotypes are considered duplicates
        - Shallow mode is useful for large datasets that don't fit in memory
        - This function is useful for removing redundant reads before downstream analysis
        - The output maintains the same structure as the input but with fewer reads
    """
    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Output File: {output_file}')
    logger.info(f'Shallow: {shallow}')
    logger.info(f'Ignore Alleles: {ignore_alleles}')

    logger.info(f'Loading EMASE file: {alignment_file}')
    apm = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)
    logger.debug(f'Number Loci: {apm.num_loci}')
    logger.debug(f'Number Haplotypes: {apm.num_haplotypes}')
    logger.debug(f'Number Reads: {apm.num_reads}')

    logger.info('Getting unique reads')
    if group_file:
        logger.debug('Using group file')
        aln_mat_g = apm.bundle(reset=True, shallow=shallow)
        aln_mat_g_uniq = aln_mat_g.get_unique_reads(
            ignore_haplotype=ignore_alleles, shallow=shallow
        )
        num_alns_per_read = aln_mat_g_uniq.sum(axis=AlignmentPropertyMatrix.Axis.LOCUS).sum(
            axis=AlignmentPropertyMatrix.Axis.HAPLOTYPE
        )
        aln_mat_uniq = apm.pull_alignments_from(
            (num_alns_per_read > 0), shallow=shallow
        )
    else:
        logger.debug('Not using group file')
        aln_mat_uniq = apm.get_unique_reads(
            ignore_haplotype=ignore_alleles, shallow=shallow
        )

    logger.info(f'Saving EMASE Formatted File: {output_file}')
    aln_mat_uniq.save(h5_file=output_file, shallow=shallow)
    logger.info('Done')


def parse_gtf(gtf_fh) -> tuple[dict, dict]:
    """
    Parse a GTF (Gene Transfer Format) file to extract gene and transcript annotations.

    This function parses a GTF file to build comprehensive gene and transcript databases
    that contain genomic coordinates, exon structures, and other annotation features.
    It is used in the GBRS pipeline to extract transcript sequences from genome files
    based on gene annotations.

    The function processes GTF entries for different feature types:
    - gene: Gene-level annotations with genomic coordinates
    - transcript: Transcript-level annotations and structure
    - exon: Exon coordinates and numbering
    - UTR: Untranslated region coordinates
    - start_codon/stop_codon: Translation start/stop sites

    ALGORITHM:
    1. Skip comment lines at the beginning of the file
    2. Parse each GTF line into tab-separated fields
    3. Extract attributes from the 9th field (semicolon-separated key-value pairs)
    4. Process different feature types:
       - gene: Create gene database entry
       - transcript: Create transcript database entry
       - exon: Add exon coordinates to transcript
       - UTR/start_codon/stop_codon: Add feature coordinates
    5. Build gene-transcript relationships

    Args:
        gtf_fh: File handle or iterable containing GTF file lines.
            The file should be in standard GTF format with tab-separated fields:
            chromosome, source, feature, start, end, score, strand, frame, attributes

    Returns:
        tuple[dict, dict]: A tuple containing two dictionaries:
            - gdb (gene database): Dictionary mapping gene IDs to gene information
                Keys: gene IDs (str)
                Values: dict containing:
                    - 'chr': chromosome name (str)
                    - 'strand': strand (+ or -) (str)
                    - 'start': gene start position (int)
                    - 'end': gene end position (int)
                    - 'isoform': set of transcript IDs (set)
                    - Additional attributes from GTF file

            - tdb (transcript database): Dictionary mapping transcript IDs to transcript information
                Keys: transcript IDs (str)
                Values: dict containing:
                    - 'chr': chromosome name (str)
                    - 'strand': strand (+ or -) (str)
                    - 'start': transcript start position (int)
                    - 'end': transcript end position (int)
                    - 'exon': list of (start, end) exon coordinates (list)
                    - 'exon_number': list of exon numbers (list)
                    - 'UTR': list of UTR coordinates (list)
                    - 'five_prime_utr': list of 5' UTR coordinates (list)
                    - 'three_prime_utr': list of 3' UTR coordinates (list)
                    - 'start_codon': list of start codon coordinates (list)
                    - 'stop_codon': list of stop codon coordinates (list)
                    - Additional attributes from GTF file

    Raises:
        ValueError: If GTF file format is invalid or required attributes are missing.
        IndexError: If GTF line has insufficient fields.

    Notes:
        - The function handles both standard and non-standard GTF files
        - Duplicate gene or transcript entries are reported as errors
        - Exon coordinates are stored as (start, end) tuples
        - The function builds gene-transcript relationships automatically
        - CDS features are parsed but not stored (only coordinates are used)
        - Unknown feature types are logged as errors but processing continues
    """
    gdb = dict()
    tdb = dict()
    for line in dropwhile(utils.is_comment, gtf_fh):
        item = line.rstrip().split('\t')
        attribute = defaultdict(list)
        for e in item[8].split('; '):
            e1, e2 = e.replace(';', '').split('"')[:2]
            attribute[e1.strip()].append(e2.strip())
        attribute = dict(attribute)
        # leave the chromosome names as is in gtf file
        chrom = item[0]
        feature = item[2]
        s = int(item[3])
        e = int(item[4])
        strand = item[6]
        if feature == 'gene':
            # Seqnature currently has some issue processing this entry
            gid = attribute.pop('gene_id')[0]
            if gid in gdb:
                print(f'[Error] Duplicate entry: {gid}')
            else:
                gdb[gid] = {
                    'chr': chrom,
                    'strand': strand,
                    'start': s,
                    'end': e,
                    'isoform': set(),
                }
                for k, v in attribute.items():
                    if len(v) == 1:
                        v = v.pop()
                    gdb[gid][k] = v
        elif feature == 'transcript':
            gid = attribute['gene_id'][0]
            tid = attribute.pop('transcript_id')[0]
            if tid in tdb:
                print(f'[Error] Duplicate entry: {tid}')
            else:
                tdb[tid] = {
                    'chr': chrom,
                    'strand': strand,
                    'start': s,
                    'end': e,
                    'eid': [],
                    'exon': [],
                    'exon_number': [],  # Do we need this?
                    'UTR': [],
                    'five_prime_utr': [],
                    'three_prime_utr': [],
                    'Selenocysteine': [],
                    'start_codon': [],
                    'stop_codon': [],
                    'start_codon_frame': [],
                    'stop_codon_frame': [],
                }
                for k, v in attribute.items():
                    if len(v) == 1:
                        v = v.pop()
                    tdb[tid][k] = v
                gid = tdb[tid]['gene_id']
                if gid in gdb:
                    gdb[gid]['isoform'].add(tid)
                else:  # This is a non-standard case where the input gtf does not have detailed gene-level annotation
                    gdb[gid] = dict()
                    gdb[gid]['chr'] = chrom
                    gdb[gid]['isoform'] = set(tid)
        else:
            gid = attribute['gene_id'][0]
            tid = attribute['transcript_id'][0]
            if feature == 'exon':
                if tid in tdb:
                    try:
                        enu = int(attribute.pop('exon_number')[0])
                        tdb[tid]['exon_number'].append(enu)
                    except:
                        pass
                    try:
                        eid = attribute.pop('exon_id')[0]
                        tdb[tid]['eid'].append(eid)
                    except:
                        pass
                    tdb[tid]['exon'].append((s, e))
                else:  # This is a non-standard case where input gtf does not have detailed transcript-level annotation
                    tdb[tid] = dict()
                    tdb[tid]['chr'] = chrom
                    tdb[tid]['strand'] = strand
                    tdb[tid]['exon'] = [(s, e)]
            elif feature in (
                    'UTR',
                    'five_prime_utr',
                    'three_prime_utr',
                    'Selenocysteine',
            ):
                tdb[tid][feature].append((s, e))
            elif feature in ('start_codon', 'stop_codon'):
                tdb[tid][feature].append((s, e))
                tdb[tid][feature + '_frame'].append(int(item[7]))
            elif feature == 'CDS':
                pass
            else:
                logger.error(f'[Unknown feature: {feature}]/n{line}')
    return gdb, tdb


def get_fragment(start: int, end: int, chro: str, strand: str, genome: dict) -> Seq:
    """
    Extract a genomic fragment and return the appropriate DNA sequence.

    This function extracts a DNA sequence from a specified genomic region and
    returns the appropriate strand based on the gene orientation. It is used
    in the GBRS pipeline to extract transcript sequences from genome files
    based on gene annotations.

    The function handles both positive and negative strand genes:
    - Positive strand (+): Returns the sequence as-is from the genome
    - Negative strand (-): Returns the reverse complement of the sequence

    Args:
        start: Start position of the genomic region (1-based coordinates).
            Should be within the bounds of the specified chromosome.

        end: End position of the genomic region (1-based coordinates).
            Should be greater than or equal to start.

        chro: Chromosome name to extract the sequence from.
            Must be a key in the genome dictionary.

        strand: Strand orientation of the gene ('+' or '-').
            Determines whether to return the sequence as-is or its reverse complement.

        genome: Dictionary containing chromosome sequences.
            Keys: chromosome names (str)
            Values: Bio.SeqRecord.SeqRecord objects containing chromosome sequences

    Returns:
        Bio.Seq.Seq: DNA sequence of the specified genomic region.
            For positive strand genes: sequence as extracted from genome
            For negative strand genes: reverse complement of the extracted sequence

    Raises:
        KeyError: If the specified chromosome is not found in the genome dictionary.
        IndexError: If start or end positions are outside the chromosome bounds.
        ValueError: If start is greater than end or if strand is not '+' or '-'.

    Notes:
        - Uses 1-based coordinates (GTF standard) but converts to 0-based for Python indexing
        - The function assumes the genome dictionary contains Bio.SeqRecord objects
        - For negative strand genes, the reverse complement is computed automatically
        - This function is typically used in conjunction with parse_gtf() to extract
          transcript sequences from genome files
    """
    fragment = genome[chro].seq[(start - 1): end]
    if strand == '-':
        fragment = fragment.reverse_complement()
    return fragment


def prepare(
        genome_files: list[str],
        haplotypes: list[str] = None,
        gtf_files: list[str] = None,
        out_dir: str = None,
        save_g2tmap: bool = False,
        save_dbs: bool = False,
        no_bowtie_index: bool = False
) -> None:
    """
    Prepare transcriptome files for EMASE analysis from genome and annotation files.

    This function is a comprehensive preparation step in the GBRS pipeline that
    extracts transcript sequences from genome files using gene annotations and
    creates the necessary files for downstream RNA-Seq analysis. It processes
    multiple founder haplotypes to create a pooled transcriptome suitable for
    allele-specific expression analysis.

    The function performs several key operations:
    1. Parses GTF annotation files to extract gene and transcript information
    2. Extracts transcript sequences from genome files based on exon coordinates
    3. Creates a pooled transcriptome with haplotype-specific transcript IDs
    4. Optionally builds Bowtie indices for RNA-Seq read alignment
    5. Generates transcript length files and gene-transcript mappings

    ALGORITHM:
    1. Validate input files and create output directory
    2. For each haplotype:
       - Load genome sequence from FASTA file
       - Parse GTF annotation file
       - Extract transcript sequences using exon coordinates
       - Append haplotype suffix to transcript IDs
       - Write to pooled transcriptome file
    3. Generate transcript length information
    4. Optionally create gene-to-transcript mapping
    5. Optionally build Bowtie index

    OUTPUT FILES:
    - emase.transcripts.fa: Pooled transcriptome FASTA file
    - emase.transcripts.info: Transcript length information
    - emase.gene2transcripts.tsv: Gene-to-transcript mapping (if save_g2tmap=True)
    - emase.gdb.pkl: Gene database pickle file (if save_dbs=True)
    - emase.tdb.pkl: Transcript database pickle file (if save_dbs=True)
    - bowtie.transcripts.*: Bowtie index files (if no_bowtie_index=False)

    Args:
        genome_files: List of paths to genome FASTA files, one per founder haplotype.
            Each file should contain chromosome sequences from a single founder strain.
            Example: ['strain_A.fa', 'strain_B.fa', 'strain_C.fa']

        haplotypes: List of haplotype identifiers corresponding to each genome file.
            If None, uses uppercase letters starting from 'A' for each genome file.
            Example: ['A', 'B', 'C'] for Diversity Outbred founder strains

        gtf_files: List of paths to GTF annotation files, one per genome file.
            Each file should contain gene annotations for the corresponding genome.
            Must have same length as genome_files if provided.

        out_dir: Output directory for generated files. Default: None (current directory)
            If None, files are created in the current working directory.

        save_g2tmap: Whether to save gene-to-transcript mapping file.
            If True, creates emase.gene2transcripts.tsv with gene-transcript relationships.

        save_dbs: Whether to save gene and transcript databases as pickle files.
            If True, creates emase.gdb.pkl and emase.tdb.pkl for later use.

        no_bowtie_index: Whether to skip Bowtie index creation.
            If True, does not create Bowtie index files for RNA-Seq alignment.

    Returns:
        None. The function creates transcriptome files and optional index files.

    Raises:
        FileNotFoundError: If any genome or GTF file does not exist.
        ValueError: If genome_files and gtf_files have different lengths.
        OSError: If output directory cannot be created or files cannot be written.
        subprocess.CalledProcessError: If Bowtie index building fails.

    Notes:
        - Transcript IDs in the output follow the pattern: {original_id}_{haplotype}
        - For single genome analysis, no haplotype suffix is added
        - The function handles both compressed (.gz) and uncompressed files
        - Bowtie index building can be time-consuming for large transcriptomes
        - This function is typically used as the first step in the GBRS pipeline
          to prepare reference transcriptomes for RNA-Seq alignment
    """
    for x in genome_files:
        logger.info(f'Genome File: {x}')
    logger.info(f'Haplotypes: {haplotypes}')
    for x in gtf_files:
        logger.info(f'GTF File: {x}')
    logger.info(f'Output Dir: {out_dir}')
    logger.info(f'Save Gene 2 Transcript Map: {save_g2tmap}')
    logger.info(f'Save DBs: {save_dbs}')
    logger.info(f'No bowtie Index: {no_bowtie_index}')

    num_haps = len(genome_files)

    if haplotypes is None:
        haplotypes = list(string.ascii_uppercase[:num_haps])
        if num_haps == 1:
            logger.info('Assuming single genome analysis. No suffix will be added to ID\'s')
        else:
            logger.info(f'Default haplotype names will be used: {", ".join(haplotypes)}')

    if gtf_files is None:
        gtf_files = []
        for genome_file in genome_files:
            gtf_files.append(f'{os.path.splitext(genome_file)[0]}.gtf')
        gtf_files_str = '\n'.join(gtf_files)
        logger.info(f'Assuming there exist the following GTF files:\n{gtf_files_str}')

    if len(haplotypes) != num_haps or len(gtf_files) != num_haps:
        logger.warning('The number of gtf files or specified haplotypes is not matching to the number of genomes.')

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # Get pooled transcriptome
    if num_haps > 1:
        transcriptome_file = os.path.join(out_dir, 'emase.pooled.transcripts.fa')
        len_file = os.path.join(out_dir, 'emase.pooled.transcripts.info')
    else:
        transcriptome_file = os.path.join(out_dir, 'emase.transcripts.fa')
        len_file = os.path.join(out_dir, 'emase.transcripts.info')

    seq_out = open(transcriptome_file, 'w')
    len_out = open(len_file, 'w')

    for hid in range(num_haps):
        genome_file = genome_files[hid]
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
        hap_name = haplotypes[hid]
        gtf_file = gtf_files[hid]
        logger.info(f'Loading {genome_name} genome...')

        if os.path.splitext(genome_file)[1] == '.gz':
            genome_fh = gzip.open(genome_file, 'rb')
        else:
            genome_fh = open(genome_file)

        genome = SeqIO.to_dict(SeqIO.parse(genome_fh, 'fasta'))
        genome_fh.close()

        if os.path.splitext(gtf_file)[1] == '.gz':
            anno_fh = gzip.open(gtf_file, 'rb')
        else:
            anno_fh = open(gtf_file)

        logger.info(f'Parsing {os.path.basename(gtf_file)}...')
        gdb, tdb = parse_gtf(anno_fh)
        anno_fh.close()

        if num_haps == 1:
            print(
                f'Building {genome_name} transcriptome (Note: No suffix added to ID\'s)...'
            )
            if save_dbs:
                import _pickle as cPickle

                cPickle.dump(
                    gdb, open(os.path.join(out_dir, 'emase.gdb.pkl'), 'wb')
                )
                cPickle.dump(
                    tdb, open(os.path.join(out_dir, 'emase.tdb.pkl'), 'wb')
                )
        elif num_haps > 1:
            logger.info(
                f'Building {genome_name} transcriptome using suffix "_{hap_name}"...',
            )

        for tid in sorted(list(tdb.keys())):
            tinfo = tdb[tid]
            if tinfo['chr'] in genome:
                # Filter out transcripts from chromosome that the input genome does not contain
                if num_haps > 1:
                    # No need to add suffix if we deal with a single genome
                    tid = f'{tid}_{hap_name}'
                fragment = Seq('')
                for exon in tinfo['exon']:
                    fragment += get_fragment(
                        exon[0],
                        exon[1],
                        tinfo['chr'],
                        tinfo['strand'],
                        genome,
                    )
                if len(fragment) > 0:
                    SeqIO.write(
                        SeqRecord(fragment, tid, '', ''), seq_out, 'fasta'
                    )
                len_out.write(f'{tid}\t{len(fragment)}\n')
            else:
                print(
                    f'Skipping Transcript {tid} of Chromosome {tinfo["chr"]}...',
                )
    seq_out.close()
    len_out.close()

    if save_g2tmap:
        with open(
                os.path.join(out_dir, 'emase.gene2transcripts.tsv'), 'w'
        ) as fhout:
            logger.info('Recording mapping of gene id to transcript id\'s...', )
            for gid in sorted(list(gdb.keys())):
                if gdb[gid]['chr'] in genome:
                    item = [gid]
                    item = item + list(gdb[gid]['isoform'])
                    fhout.write('\t'.join(item) + '\n')

    #
    # Build bowtie index for the pooled transcriptome
    if not no_bowtie_index:
        out_index = os.path.join(
            os.path.dirname(transcriptome_file), 'bowtie.transcripts'
        )
        logger.info('Building bowtie index...')
        subprocess.call(
            f'bowtie-build {transcriptome_file} {out_index}', shell=True
        )

    logger.info('Done')


def run(
        alignment_file: str,
        group_file: str = None,
        length_file: str = None,
        outbase: str = 'emase',
        multiread_model: int = 4,
        read_length: int = 100,
        pseudocount: float = 0.0,
        max_iters: int = 999,
        tolerance: float = 0.0001,
        report_alignment_counts: bool = False,
        report_posterior: bool = False
) -> None:
    """
    Run EMASE (Expectation-Maximization for Allele-Specific Expression) analysis.

    This function is the core analysis component of the GBRS pipeline that performs
    allele-specific expression quantification using an expectation-maximization
    algorithm. It processes EMASE alignment data to estimate transcript abundances
    and allele-specific expression levels in multiparent populations.

    The function implements the EMASE algorithm described in the scientific literature
    for quantifying gene expression from RNA-Seq data while accounting for:
    - Multi-mapping reads that align to multiple transcripts
    - Allele-specific expression differences between founder haplotypes
    - Transcript length biases in RNA-Seq quantification
    - Sequencing depth normalization

    ALGORITHM:
    1. Load EMASE alignment data and optional group information
    2. Initialize EMASE factory with alignment matrix
    3. Prepare EMASE with transcript lengths and read length information
    4. Run expectation-maximization algorithm with specified model
    5. Generate transcript abundance estimates (TPM and read counts)
    6. Optionally generate gene-level estimates and additional reports

    OUTPUT FILES:
    - {outbase}.isoforms.tpm: Transcript-level TPM (Transcripts Per Million) estimates
    - {outbase}.isoforms.expected_read_counts: Transcript-level read count estimates
    - {outbase}.genes.tpm: Gene-level TPM estimates (if group_file provided)
    - {outbase}.genes.expected_read_counts: Gene-level read count estimates (if group_file provided)
    - {outbase}.isoforms.alignment_counts: Raw alignment counts per transcript (if report_alignment_counts=True)
    - {outbase}.genes.alignment_counts: Raw alignment counts per gene (if report_alignment_counts=True)
    - {outbase}.posterior.h5: Posterior probability estimates (if report_posterior=True)

    Args:
        alignment_file: Path to the EMASE file (HDF5 format) containing alignment data.
            The file should contain sparse matrices representing read alignments to
            transcripts from different founder haplotypes.

        group_file: Path to the group file containing transcript-to-gene mapping.
            If provided, enables gene-level analysis in addition to transcript-level analysis.
            Format: tab-separated with transcript ID in first column, gene ID in second.

        length_file: Path to the transcript length file.
            Contains transcript IDs and their lengths for length bias correction.
            Format: tab-separated with transcript ID in first column, length in second.

        outbase: Base name for all output files. Default: 'emase'
            All output files will be prefixed with this base name.

        multiread_model: EMASE model to use for multi-mapping read handling.
            Options: 1-4, where higher numbers represent more sophisticated models
            for handling reads that align to multiple transcripts.

        read_length: Average read length for length bias correction.
            Used in the EMASE algorithm to account for transcript length biases.

        pseudocount: Prior pseudocount for regularization.
            Adds a small prior to prevent zero probabilities in the EM algorithm.

        max_iters: Maximum number of EM iterations.
            The algorithm will stop if convergence is reached before this limit.

        tolerance: Convergence tolerance for EM algorithm.
            The algorithm stops when the change in TPM estimates is below this threshold.

        report_alignment_counts: Whether to report raw alignment counts.
            If True, generates files with raw alignment counts per transcript/gene.

        report_posterior: Whether to report posterior probability estimates.
            If True, generates HDF5 file with posterior probabilities for each read-transcript pair.

    Raises:
        FileNotFoundError: If alignment_file, group_file, or length_file does not exist.
        ValueError: If the EMASE file format is invalid or parameters are invalid.
        RuntimeError: If the EM algorithm fails to converge or other processing errors.

    Notes:
        - The function uses the EMASE algorithm for robust expression quantification
        - TPM values are normalized for transcript length and sequencing depth
        - Gene-level estimates are aggregated from transcript-level estimates
        - The algorithm handles multi-mapping reads using the specified model
        - Posterior probabilities provide uncertainty estimates for expression levels
        - This function is the final step in the GBRS pipeline for expression quantification
    """
    report_group_counts = (
            group_file is not None
    )

    logger.info(f'Alignment File: {alignment_file}')
    logger.info(f'Group File: {group_file}')
    logger.info(f'Read Length File: {length_file}')
    logger.info(f'Read Length: {read_length}')
    logger.info(f'Outbase: {outbase}')
    logger.info(f'Multiread Model: {multiread_model}')
    logger.info(f'Pseudocount: {pseudocount}')
    logger.info(f'Tolerance: {tolerance}')
    logger.info(f'Report Alignment Counts: {report_alignment_counts}')
    logger.info(f'Report Posterior: {report_posterior}')

    # load alignment incidence matrix ('alignment_file' is assumed to be in multiway transcriptome)
    logger.info(f'Loading EMASE file: {alignment_file}')
    apm = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)
    logger.debug(f'Number Loci: {apm.num_loci}')
    logger.debug(f'Number Haplotypes: {apm.num_haplotypes}')
    logger.debug(f'Number Reads: {apm.num_reads}')

    logger.info('Running EMASE')
    em_factory = EMfactory(apm)
    em_factory.prepare(pseudocount=pseudocount, lenfile=length_file, read_length=read_length)
    em_factory.run(
        model=multiread_model, tol=tolerance, max_iters=max_iters, verbose=True
    )

    logger.info(f'Generating isoform TPMs: {outbase}.isoforms.tpm')
    em_factory.report_depths(
        filename=f'{outbase}.isoforms.tpm', tpm=True
    )

    logger.info(f'Generating isoform Read Counts: {outbase}.isoforms.expected_read_counts')
    em_factory.report_read_counts(
        filename=f'{outbase}.isoforms.expected_read_counts'
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
            grp_wise=True
        )

        logger.info(f'Generating gene Read Counts: {outbase}.genes.expected_read_counts')
        em_factory.report_read_counts(
            filename=f'{outbase}.genes.expected_read_counts',
            grp_wise=True
        )

    if report_alignment_counts:
        apm = AlignmentPropertyMatrix(h5_file=alignment_file, grp_file=group_file)

        logger.info(f'Generating isoform Alignment Counts: {outbase}.isoforms.alignment_counts')
        apm.report_alignment_counts(filename=f'{outbase}.isoforms.alignment_counts')

        if report_group_counts:
            logger.info(f'Generating gene Alignment Counts: {outbase}.genes.alignment_counts')
            apm._bundle_inline(reset=True)
            apm.report_alignment_counts(filename=f'{outbase}.genes.alignment_counts')
    logger.debug('Done')
