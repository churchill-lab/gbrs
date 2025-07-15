# Standard library imports
import os
import struct
from collections import Counter

# 3rd party library imports
import numpy as np
import pysam
import tables
from scipy.sparse import coo_matrix, csc_matrix

# Local library imports
from gbrs import h5_utils
from gbrs import utils

logger = utils.get_logger('gbrs')


class AlignmentMatrixFactory:
    """
    Factory class for creating alignment matrices from BAM files.

    This class handles the conversion of paired-end BAM alignment data to the EMASE format used by
    the GBRS pipeline. It processes multiple BAM files (typically R1 and R2 from paired-end
    sequencing) and creates HDF5 files containing sparse alignment matrices organized by haplotypes.

    The factory works in two phases:

    **Preparation**: Parse multiple BAM files, extract read and locus information, create temporary
        files for each haplotype and read end
    **Production**: Convert temporary files to HDF5 format with sparse matrices in CSC (Compressed
        Sparse Column) format

    Key Features:
    - Supports multiple BAM files (typically R1 and R2)
    - Handles haplotype-specific reference sequences
    - Creates memory-efficient sparse matrices
    - Combines data from multiple read ends using element-wise multiplication
    - Supports compression and custom data types
    - Automatic cleanup of temporary files

    The resulting HDF5 files can be loaded directly by AlignmentPropertyMatrix for analysis.

    Attributes:
        alignment_files: List of paths to input BAM files
        hname: List of haplotype names
        lname: List of locus (transcript) names
        rname: Array of read names (sorted)
        tmp_files: Nested mapping of haplotype names and file indices to temporary file paths
    """

    def __init__(
        self,
        alignment_files: list[str]
    ) -> None:
        """
        Initialize the AlignmentMatrixFactory, validate the input BAM files and prepare the factory
        for processing RNA-seq alignment data.

        Args:
            alignment_files: List of paths to input BAM files containing RNA-seq alignments

        Raises:
            FileNotFoundError: If any of the specified BAM files do not exist
            ValueError: If alnfile is not a list or is empty

        Note:
            The BAM files should contain alignments to a haplotype-specific reference genome where
            reference names follow the format: 'locus_haplotype' (e.g., 'ENSMUST00000123456_A').

            Typically, this would be used with two BAM files:
            - R1 alignments (first read of each pair)
            - R2 alignments (second read of each pair)

            The factory assumes that both BAM files contain the same read names and will use the
            read names from the first BAM file.
        """
        if not isinstance(alignment_files, list) or len(alignment_files) == 0:
            raise ValueError('alignment_files must be a non-empty list of BAM file paths')

        if not all(os.path.exists(bam) for bam in alignment_files):
            raise FileNotFoundError('One or more BAM files do not exist')

        self.alignment_files = alignment_files
        self.hname: list[str] | None = None
        self.lname: list[str] | None = None
        self.rname: np.ndarray | None = None
        self.tmp_files: dict | None = None


    def prepare(
        self,
        haplotypes: list[str],
        loci: list[str],
        delim: str = '_',
        out_dir: str | None = None
    ) -> None:
        """
        Prepare the factory for matrix production by processing BAM files.

        This method performs the preparation phase of the factory workflow:
        - Extracts read names from the first BAM file
        - Creates temporary binary files for each haplotype and read end
        - Processes alignments and writes read-locus pairs to temporary files

        The preparation phase is memory-efficient as it processes BAM files sequentially and stores
        intermediate results in temporary binary files rather than keeping all data in memory.

        Args:
            haplotypes: List of haplotype names to process.
            loci: List of locus (transcript) names
            delim: Delimiter used in reference names to separate locus and haplotype.
            out_dir: Directory to write temporary files. If None, uses directory of first BAM file.

        Raises:
            RuntimeError: If BAM files cannot be read or processed
            ValueError: If loci list is empty or invalid

        Note:
            This method creates temporary binary files with names like:
                '{haplotype}_{file_index}_{process_id}.bin'.
            Each temporary file contains binary data with pairs of uint32 values:
                (read_index, locus_index) for each alignment.
        """
        if len(loci) == 0:
            raise ValueError('loci list cannot be empty')

        # set haplotype names
        if len(haplotypes) > 0:
            self.hname = haplotypes
        else:
            self.hname = ['h0']

        self.lname = loci

        # extract read names from first BAM file
        logger.debug('Gathering all read names...')
        save = pysam.set_verbosity(0)
        fh = pysam.AlignmentFile(self.alignment_files[0], 'rb')
        pysam.set_verbosity(save)

        # NOTE: Single BAM file is used, as the assumption is that both have the same reads
        # not looking at flags, just taking all reads
        self.rname = {aln.query_name for aln in fh.fetch(until_eof=True)}
        logger.debug(f'Retrieved {len(self.rname):,} read names')
        fh.close()

        # sort read names for consistent ordering
        logger.debug('Sorting read names...')
        self.rname = np.array(list(self.rname), dtype='S')
        self.rname.sort()
        logger.debug(f'Sorted {len(self.rname):,} read names')

        # create ID mappings
        lid = {self.lname[i]: i for i in range(len(self.lname))}
        rid = {self.rname[i].decode(): i for i in range(len(self.rname))}

        # set output directory
        if out_dir is None:
            out_dir = os.path.dirname(self.alignment_files[0])

        # initialize temporary file structure
        fhout = dict.fromkeys(self.hname)
        for hap in self.hname:
            fhout[hap] = {idx: [] for idx in range(len(self.alignment_files))}
        self.tmp_files = fhout

        # create temporary files
        for hap in self.hname:
            for idx, bam in enumerate(self.alignment_files):
                outfile = os.path.join(out_dir, f'{hap}_{idx}_{os.getpid()}.bin')
                logger.debug(f'Initializing temp file: {outfile}')
                self.tmp_files[hap][idx] = outfile
                fhout[hap][idx] = open(outfile, 'wb')

        # pre-compile struct.pack for better performance
        pack_uint32 = struct.Struct('>I').pack

        alignment_count_total = Counter()
        alignment_count_skipped = Counter()

        # process each BAM file
        for idx, bam in enumerate(self.alignment_files):
            logger.debug(f'Processing BAM file: {idx + 1}/{len(self.alignment_files)}: {bam}')
            save = pysam.set_verbosity(0)
            fh = pysam.AlignmentFile(bam, 'rb')
            pysam.set_verbosity(save)

            has_haplotypes = len(haplotypes) > 0

            for alignment in fh.fetch(until_eof=True):
                alignment_count_total[bam] += 1

                # Skip reads that should never be processed:
                # - Unmapped (0x4), no alignment to reference
                # - Secondary (0x100), not the primary alignment (e.g., multimapping)
                # - Supplementary (0x800), split or chimeric alignments
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    alignment_count_skipped[bam] += 1
                    continue

                # For paired-end reads, apply stricter filters to ensure clean, concordant pairs
                if alignment.is_paired:
                    if (
                            # 0x80: skip second read to avoid double-counting
                            alignment.is_read2 or
                            # !0x2: improperly paired (bad orientation/insert)
                            not alignment.is_proper_pair or
                            # mapped to different chromosomes
                            alignment.reference_id != alignment.next_reference_id or
                            # 0x8: mate not mapped
                            alignment.mate_is_unmapped
                    ):
                        alignment_count_skipped[bam] += 1

                # extract reference name and parse locus/haplotype if needed
                reference_name = fh.get_reference_name(alignment.reference_id)
                if has_haplotypes:
                    # assumes reference_name = "locus(delim)hap"
                    locus, hap = reference_name.split(delim)
                else:
                    # use default haplotype name
                    hap = self.hname[0]
                    locus = reference_name

                # write encoded query name and locus ID to appropriate haplotype stream
                fhout[hap][idx].write(pack_uint32(rid[alignment.query_name]))
                fhout[hap][idx].write(pack_uint32(lid[locus]))

            # close temporary files for this BAM file
            for hap in self.hname:
                fhout[hap][idx].close()

        # output file statistics
        alignments_total = 0
        alignments_skipped = 0

        logger.info(f'{"TOTAL":>12} {"VALID":>12}  FILE')
        logger.info('-' * 80)

        for idx, bam in enumerate(self.alignment_files):
            total = alignment_count_total[bam]
            valid = alignment_count_total[bam] - alignment_count_skipped[bam]
            if len(self.alignment_files) > 1:
                logger.info(f'{total:>12,} {valid:>12,}  {os.path.basename(bam)}')

            alignments_total += alignment_count_total[bam]
            alignments_skipped += alignment_count_skipped[bam]

        logger.info('-' * 80)
        logger.info(f'{alignments_total:>12,} {(alignments_total - alignments_skipped):>12,}  TOTAL')

        logger.debug('Finished temp file creation')

    def produce(
        self,
        h5_file: str,
        title: str = 'Alignments',
        index_dtype: str = 'uint32',
        data_dtype: type = float,
        complib: str = 'zlib',
        incidence_only: bool = True
    ) -> None:
        """
        Produce an HDF5 file containing the alignment matrix.

        This method performs the production phase of the factory workflow:
        - Reads temporary binary files created during preparation
        - Converts binary data to sparse matrices
        - Combines matrices from multiple read ends using element-wise multiplication
        - Saves final sparse matrices in HDF5 format

        Args:
            h5_file: Path to the output HDF5 file
            title: Title for the HDF5 file. Defaults to 'Alignments'
            index_dtype: Data type for sparse matrix indices (indptr, indices).
            data_dtype: Data type for sparse matrix data. Defaults to float
            complib: Compression library to use. Defaults to 'zlib'
            incidence_only: If True, store only binary incidence (1.0 for alignments). If False,
                store actual alignment counts.

        Raises:
            RuntimeError: If temporary files are missing or cannot be read
            ValueError: If factory has not been prepared (prepare() not called)

        Note:
            The method combines data from multiple BAM files using element-wise multiplication of
            sparse matrices. This ensures that only reads present in ALL input files are included
            in the final matrix.
        """
        if any(v is None for v in (self.tmp_files, self.hname, self.lname, self.rname)):
            raise ValueError(
                'Factory must be prepared before producing matrix. Call prepare() first.'
            )

        logger.debug(f'Constructing: {h5_file}')
        h5fh = tables.open_file(h5_file, 'w', title=title)
        fil = tables.Filters(complevel=1, complib=complib)

        # set root attributes
        logger.debug(f'Creating attribute /incidence_only: {incidence_only}')
        h5fh.set_node_attr(h5fh.root, 'incidence_only', incidence_only)

        logger.debug('Creating attribute /mtype: csc_matrix')
        h5fh.set_node_attr(h5fh.root, 'mtype', 'csc_matrix')

        logger.debug(
            f'Creating attribute /shape: {(len(self.lname), len(self.hname), len(self.rname))}'
        )
        h5fh.set_node_attr(
            h5fh.root,
            'shape',
            (len(self.lname), len(self.hname), len(self.rname)),
        )

        logger.debug(f'Creating attribute /hname: {self.hname}')
        h5fh.set_node_attr(h5fh.root, 'hname', self.hname)

        # create metadata arrays
        logger.debug('Creating array /lname...')
        h5fh.create_carray(
            h5fh.root,
            'lname',
            obj=self.lname,
            title='Locus Names',
            filters=fil,
        )

        logger.debug('Creating array /rname...')
        h5fh.create_carray(
            h5fh.root,
            'rname',
            obj=self.rname,
            title='Read Names',
            filters=fil
        )

        spmat = dict()

        # process each haplotype
        logger.debug('Looping through haplotypes')
        for hid in range(len(self.hname)):
            logger.debug(f'Processing haplotype {hid}: {self.hname[hid]}')

            # process each BAM file for this haplotype
            for idx, bam in enumerate(self.alignment_files):
                hap = self.hname[hid]
                infile = self.tmp_files[hap][idx]
                logger.debug(f'Reading file: {infile.name}')

                # Read binary data and reshape
                dmat = np.fromfile(open(infile.name, 'rb'), dtype='>I')
                dmat = dmat.reshape((int(len(dmat) / 2), 2)).T

                # set data values (1.0 for incidence, or counts if available)
                if dmat.shape[0] > 2:
                    dvec = dmat[2]
                else:
                    dvec = np.ones(dmat.shape[1])

                # create sparse matrix
                #
                # spmat contains the read name and locus name of each read in the alignment file,
                # in the form: ((read name, locus name), read count))
                #
                # spmat is a sparse matrix with the read names as the row indices and the locus
                # names as the column indices.
                #
                # The data in the matrix is the number of reads that align to a specific locus.
                spmat[idx] = coo_matrix(
                    (dvec, dmat[:2]), shape=(len(self.rname), len(self.lname))
                )

                spmat[idx] = spmat[idx].tocsc()

            # combine matrices from multiple BAM files using element-wise multiplication
            spmat_mul = spmat[0].multiply(spmat[len(self.alignment_files) - 1])

            # create haplotype group
            logger.debug(f'Creating group /h{hid}')
            hgroup = h5fh.create_group(
                h5fh.root,
                f'h{hid}',
                f'Sparse matrix components for Haplotype {hid}',
            )

            # save sparse matrix components
            logger.debug(f'Creating array /h{hid}/indptr')
            h5fh.create_carray(
                hgroup,
                'indptr',
                obj=spmat_mul.indptr.astype(index_dtype),
                filters=fil,
            )

            logger.debug(f'Creating array /h{hid}/indices')
            h5fh.create_carray(
                hgroup,
                'indices',
                obj=spmat_mul.indices.astype(index_dtype),
                filters=fil,
            )

            if not incidence_only:
                logger.debug(f'Creating array /h{hid}/data')
                h5fh.create_carray(
                    hgroup,
                    'data',
                    obj=spmat_mul.data.astype(data_dtype),
                    filters=fil,
                )

        h5fh.flush()
        h5fh.close()

        logger.debug('File created')


    def cleanup(self) -> None:
        """
        Clean up temporary binary files created by the prepare() method. It should be called after
        produce() to free up disk space.
        """
        if self.tmp_files is None:
            logger.debug('No temporary files to clean up')
            return

        try:
            for hap, file_dict in self.tmp_files.items():
                for idx, tmp_file in file_dict.items():
                    if os.path.exists(tmp_file.name):
                        os.remove(tmp_file.name)
                    logger.debug(f'Removing temporary file: {tmp_file.name}')
        except Exception as e:
            logger.warning(f'Error during cleanup: {e}')

