# standard library imports
import copy
from enum import IntEnum
from typing import Self

# 3rd party library imports
import numpy as np
import tables
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, lil_matrix, spmatrix

# local library imports
from gbrs import utils
from gbrs.emase.Sparse3DMatrix import Sparse3DMatrix

logger = utils.get_logger('gbrs')


class AxisEnum(IntEnum):
    LOCUS = 0
    HAPLOTYPE = 1
    READ = 2
    GROUP = 3
    HAPLOGROUP = 4


class AlignmentPropertyMatrix(Sparse3DMatrix):
    """
    A specialized 3D sparse matrix for RNA-seq alignment data with metadata support. It is used to
    store and manipulate the alignment probabilities for each locus, haplotype, and read.

    This class extends Sparse3DMatrix to provide additional functionality specific to RNA-seq
    alignment analysis, including:
    - Named dimensions (loci, haplotypes, reads)
    - Group-based operations for transcript isoforms
    - Normalization methods for read-level analysis
    - Specialized summarization and filtering operations

    The matrix represents RNA-seq alignment data where:
    - Axis 0 (LOCUS): Transcripts or genes
    - Axis 1 (HAPLOTYPE): Founder strains or haplotypes
    - Axis 2 (READ): Sequencing reads
    """
    Axis = AxisEnum

    def __init__(
        self,
        other: Self | None = None,
        h5_file: str | None = None,
        h5_object: tables.File | None = None,
        datanode: str = '/',
        metanode: str = '/',
        shallow: bool = False,
        shape: tuple[int, int, int] | None = None,
        dtype: type = float,
        haplotype_names: list[str] | None = None,
        locus_names: list[str] | None = None,
        read_names: list[str] | None = None,
        grp_file: str | None = None
    ) -> None:
        """
        Initialize an AlignmentPropertyMatrix object. It extends Sparse3DMatrix initialization with
        metadata support for RNA-seq alignment analysis.

        It supports multiple initialization modes:

        - Copy from existing matrix: Provide 'other' parameter
        - Load from HDF5 file: Provide 'h5file' parameter
        - Create new matrix: Provide 'shape' and optional name parameters

        Args:
            other: Existing matrix to copy from
            h5_file: Path to HDF5 file containing matrix data
            h5_object: Open HDF5 file object
            datanode: HDF5 node path for matrix data. Defaults to '/'
            metanode: HDF5 node path for metadata. Defaults to '/'
            shallow: If True, copy only data without metadata
            shape: 3-tuple (num_loci, num_haplotypes, num_reads)
            dtype: Data type for matrix elements. Defaults to float
            haplotype_names: Names for haplotypes
            locus_names: Names for loci (transcripts/genes)
            read_names: Names for reads
            grp_file: Path to transcript group definition file

        Raises:
            RuntimeError: If copying from non-finalized matrix or invalid parameters

        Note:
            The grp_file should contain tab-separated values with group name in first column and
            transcript names in subsequent columns.

        """
        Sparse3DMatrix.__init__(
            self,
            other=other,
            h5_file=h5_file,
            datanode=datanode,
            shape=shape,
            dtype=dtype
        )

        # number of loci, haplotypes, and reads
        self.num_loci, self.num_haplotypes, self.num_reads = self.shape
        # number of groups
        self.num_groups = 0
        # count matrix
        self.count = None
        # haplotype names
        self.hname = None
        # locus names
        self.lname = None
        # read names
        self.rname = None
        # locus IDs
        self.lid = None
        # read IDs
        self.rid = None
        # group names
        self.gname = None
        # groups in terms of locus IDs
        self.groups = None

        if other is not None:
            # use for copying from other existing AlignmentPropertyMatrix object
            if other.count is not None:
                self.count = copy.copy(other.count)

            if not shallow:
                self.__copy_names(other)
                self.__copy_group_info(other)

        elif h5_file is not None:
            # use for loading from a pytables file
            h5fh = tables.open_file(h5_file, 'r')

            if h5fh.__contains__(f'{datanode}/count'):
                self.count = h5fh.get_node(datanode, 'count').read()

            if not shallow:
                self.hname = h5fh.get_node_attr(datanode, 'hname')
                self.lname = h5fh.get_node(metanode, 'lname').read()

                # convert from bytes to string
                self.lname = [x.decode() for x in self.lname]
                self.lid = dict(zip(self.lname, np.arange(self.num_loci)))

                if h5fh.__contains__(f'{metanode}/rname'):
                    self.rname = h5fh.get_node(metanode, 'rname').read()
                    self.rid = dict(zip(self.rname, np.arange(self.num_reads)))

            h5fh.close()

        elif shape is not None:
            # use for initializing an empty matrix
            if haplotype_names is not None:
                if len(haplotype_names) == self.num_haplotypes:
                    self.hname = haplotype_names
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')

            if locus_names is not None:
                if len(locus_names) == self.num_loci:
                    self.lname = np.array(locus_names)
                    self.lid = dict(zip(self.lname, np.arange(self.num_loci)))
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')

            if read_names is not None:
                if len(read_names) == self.num_reads:
                    self.rname = np.array(read_names)
                    self.rid = dict(zip(self.rname, np.arange(self.num_reads)))
                else:
                    raise RuntimeError('The number of names does not match to the matrix shape.')

        if grp_file is not None:
            self.__load_groups(grp_file)


    def __load_groups(self, grp_file: str) -> None:
        """
        Load transcript group definitions from a file from a tab-separated file and create mappings
        for group-based operations.

        Args:
            grp_file: Path to the transcript group definition file

        Raises:
            RuntimeError: If locus IDs are not available (lname/lid not set)

        Note:
            The group file format should be:
            ```
            Group1    Transcript1    Transcript2    Transcript3
            Group2    Transcript4    Transcript5
            Group3    Transcript6
            ```
        """
        if self.lid is not None:
            self.gname = list()
            self.groups = list()

            with open(grp_file) as fh:
                for line in fh:
                    item = line.rstrip().split('\t')
                    self.gname.append(item[0])
                    tid_list = [self.lid[t] for t in item[1:]]
                    self.groups.append(tid_list)

            self.gname = np.array(self.gname)
            self.num_groups = len(self.gname)
        else:
            raise RuntimeError('Locus IDs are not available.')

    load_groups = __load_groups


    def __copy_names(self, other: Self) -> None:
        """
        Copy all the naming information from another AlignmentPropertyMatrix object in place.

        Args:
            other: The source matrix to copy names from
        """
        self.hname = other.hname
        self.lname = copy.copy(other.lname)
        self.rname = copy.copy(other.rname)
        self.lid = copy.copy(other.lid)
        self.rid = copy.copy(other.rid)


    def __copy_group_info(self, other: Self) -> None:
        """
        Copy all the group information from another AlignmentPropertyMatrix object in place.

        Args:
            other: The source matrix to copy group info from
        """
        if other.groups is not None and other.gname is not None:
            self.groups = copy.deepcopy(other.groups)
            self.gname = copy.copy(other.gname)
            self.num_groups = other.num_groups


    def copy(self, shallow: bool = False) -> Self:
        """
        Create a copy of the AlignmentPropertyMatrix object including all data and metadata. The
        copy operation preserves the matrix structure, names, group information, and count data.

        Args:
            shallow: If True, copy only the sparse matrix data without metadata (names, groups,
                counts). If False, copy everything.

        Returns:
            A new matrix object with copied data

        Raises:
            RuntimeError: If the source matrix is not finalized

        Note:
            Shallow copies are useful for temporary operations where metadata is not needed, saving
            memory and computation time.
        """
        dmat = Sparse3DMatrix.copy(self)
        dmat.count = copy.copy(self.count)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape

        if not shallow:
            dmat.__copy_names(self)
            dmat.__copy_group_info(self)

        return dmat

    def _bundle_inline(self, reset: bool = False) -> None:
        """
        Perform inline bundling of transcript groups in place.

        - Changes shape from:
            (num_transcripts, num_haplotypes, num_reads) to (num_groups, num_haplotypes, num_reads)
        - Updates locus names to group names
        - Clears group information after bundling

        Inline bundling is memory-efficient but destructive - the original transcript-level data is
        lost. Use bundle() for non-destructive bundling.

        Args:
            reset: If True, reset all values to 1.0 after bundling.
                If False, preserve the original values.

        Returns:
            None: Modifies the current object in-place

        Raises:
            RuntimeError: If matrix is not finalized or no group information available
        """
        if self.finalized:
            if self.num_groups > 0 and self.groups is not None and self.gname is not None:
                grp_conv_mat = lil_matrix((self.num_loci, self.num_groups))

                for i in range(self.num_groups):
                    grp_conv_mat[self.groups[i], i] = 1.0

                grp_conv_mat = grp_conv_mat.tocsc()

                for hid in range(self.num_haplotypes):
                    # TODO: Is there any better way to save memory?
                    self.data[hid] = (self.data[hid] * grp_conv_mat)

                self.num_loci = self.num_groups
                self.shape = (self.num_groups, self.num_haplotypes, self.num_reads)
                self.lname = copy.copy(self.gname)
                self.lid = dict(zip(self.gname, np.arange(self.num_groups)))
                self.num_groups = 0
                self.groups = None
                self.gname = None

                if reset:
                    self.reset()
            else:
                raise RuntimeError('No group information is available for bundling.')
        else:
            raise RuntimeError('The matrix is not finalized.')


    def bundle(self, reset: bool = False, shallow: bool = False) -> Self:
        """
        Create a bundled AlignmentPropertyMatrix with transcript groups. This method creates a new
        matrix where transcript-level alignments are aggregated into group-level alignments (e.g.,
        gene-level from transcript-level). The bundling operation uses the group conversion matrix
        to combine transcript alignments within each group.

        The bundled matrix has shape (num_groups, num_haplotypes, num_reads) instead of
        (num_transcripts, num_haplotypes, num_reads).

        For memory-efficient bundling, consider using _bundle_inline() if the original data can be
        discarded.

        Args:
            reset: If True, set all bundled values to 1.0 (binary incidence). If False, preserve
                the original alignment values.
            shallow: If True, copy only the sparse matrix data without metadata. If False, copy all
                metadata including names.

        Returns:
            A new matrix with group-level structure

        Raises:
            RuntimeError: If matrix is not finalized or no group information available

        Note:
            Group information must be loaded via load_groups() before bundling.
        """
        if self.finalized:
            if self.groups is not None and self.gname is not None:
                grp_conv_mat = lil_matrix((self.num_loci, self.num_groups))

                for i in range(self.num_groups):
                    grp_conv_mat[self.groups[i], i] = 1.0

                grp_align = Sparse3DMatrix.__mul__(self, grp_conv_mat)
                grp_align.num_loci = self.num_groups
                grp_align.num_haplotypes = self.num_haplotypes
                grp_align.num_reads = self.num_reads

                grp_align.shape = (
                    grp_align.num_loci,
                    grp_align.num_haplotypes,
                    grp_align.num_reads,
                )

                if not shallow:
                    grp_align.lname = copy.copy(self.gname)
                    grp_align.hname = self.hname
                    grp_align.rname = copy.copy(self.rname)
                    grp_align.lid = dict(zip(grp_align.lname, np.arange(grp_align.num_loci)))
                    grp_align.rid = copy.copy(self.rid)

                if reset:
                    grp_align.reset()

                return grp_align
            else:
                raise RuntimeError(
                    'No group information is available for bundling.'
                )
        else:
            raise RuntimeError('The matrix is not finalized.')

    #
    # Binary Operators
    #

    def __add__(self, other: Self) -> Self:
        """
        Add 2 AlignmentPropertyMatrix objects element-wise. Both matrices must have the same shape.

        Args:
            other: Matrix to add to self

        Returns:
            New matrix with sum of elements
        """
        dmat = Sparse3DMatrix.__add__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = self.shape
        dmat.__copy_names(self)
        dmat.__copy_group_info(self)
        return dmat

    def __sub__(self, other: Self) -> Self:
        """
        Subtract another AlignmentPropertyMatrix from self element-wise. Both matrices must have
        the same shape.

        Args:
            other: Matrix to subtract from self

        Returns:
            New matrix with difference of elements
        """
        dmat = Sparse3DMatrix.__sub__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = self.shape
        dmat.__copy_names(self)
        dmat.__copy_group_info(self)
        return dmat


    def __mul__(self, other: Self | np.ndarray | spmatrix | float) -> Self:
        """
        Multiply AlignmentPropertyMatrix with another matrix or scalar.

        Args:
            other: Matrix or scalar to multiply with. Can be:
                - AlignmentPropertyMatrix: Element-wise multiplication
                - numpy.ndarray: Matrix multiplication
                - scipy.sparse matrix: Matrix multiplication
                - scalar: Element-wise scaling

        Returns:
            Result of multiplication
        """
        dmat = Sparse3DMatrix.__mul__(self, other)
        dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape

        if isinstance(other, (np.ndarray, csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
            dmat.hname = self.hname
            dmat.rname = copy.copy(self.rname)
            dmat.rid = copy.copy(self.rid)
            dmat.num_groups = 0
        else:
            dmat.__copy_names(self)
            dmat.__copy_group_info(self)

        return dmat

    #
    # Helper functions
    #

    def sum(self, axis: AxisEnum) -> np.ndarray | spmatrix:
        """
        Sum the AlignmentPropertyMatrix along a specified axis. The result provides aggregated
        alignment data for downstream analysis.

        Args:
            axis: Axis along which to sum:
                - Axis.LOCUS: Sum across loci for each haplotype-read pair
                - Axis.HAPLOTYPE: Sum across haplotypes for each locus-read pair
                - Axis.READ: Sum across reads for each locus-haplotype pair

        Returns:
            Summed data with shape depending on the axis:
            - LOCUS: (num_reads, num_haplotypes) dense array
            - HAPLOTYPE: (num_loci, num_reads) sparse matrix
            - READ: (num_haplotypes, num_loci) dense array

        Raises:
            RuntimeError: If matrix is not finalized or invalid axis specified

        Note:
            For READ axis summation, if count data is available, it is used to weight the summation.
            Otherwise, binary incidence is assumed.

            The HAPLOTYPE axis returns a sparse matrix to preserve sparsity of the original data,
            while other axes return dense arrays.
        """
        if self.finalized:
            if axis == self.Axis.LOCUS:
                # sum along loci
                sum_mat = []
                for hid in range(self.num_haplotypes):
                    sum_mat.append(self.data[hid].sum(axis=1).A)

                sum_mat = np.hstack(sum_mat)
            elif axis == self.Axis.HAPLOTYPE:
                # sum along haplotypes
                sum_mat = self.data[0]
                for hid in range(1, self.num_haplotypes):
                    # unlike others, this sum_mat is still sparse matrix
                    sum_mat = (sum_mat + self.data[hid])
            elif axis == self.Axis.READ:
                # sum along reads
                sum_mat = []
                for hid in range(self.num_haplotypes):
                    if self.count is None:
                        sum_hap = self.data[hid].sum(axis=0).A
                    else:
                        hap_mat = self.data[hid].copy()
                        hap_mat.data *= self.count[hap_mat.indices]
                        sum_hap = hap_mat.sum(axis=0).A
                    sum_mat.append(sum_hap)
                sum_mat = np.vstack(sum_mat)
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
            return sum_mat
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def normalize_reads(
        self,
        axis: AxisEnum,
        grouping_mat: spmatrix | None = None
    ) -> None:
        """
        Normalize read-level alignment probabilities along specified axis.

        This method performs read-wise normalization to convert raw alignment counts to
        probabilities. The normalization ensures that for each read, the sum of probabilities
        across the specified dimension equals 1.0.

        Args:
            axis: Axis along which to normalize:
                - Axis.LOCUS: Normalize across loci for each read
                - Axis.HAPLOTYPE: Normalize across haplotypes for each read
                - Axis.READ: Normalize each read as a whole
                - Axis.GROUP: Normalize across transcript groups for each read
                - Axis.HAPLOGROUP: Normalize across haplotype-groups for each read
            grouping_mat: Incidence matrix specifying which transcripts belong to the same
            gene/group. Required for GROUP and HAPLOGROUP normalization.

        Raises:
            RuntimeError: If matrix is not finalized, invalid axis, or missing grouping matrix when
                required

        Note:
            This method modifies the matrix in-place. The normalization is essential for the EMASE
            algorithm to work correctly, as it ensures that alignment probabilities sum to 1.0 for
            each read.
        """
        if self.finalized:
            if axis == self.Axis.LOCUS:
                # locus-wise normalization on each read
                # sparse matrix of |reads| x |loci|
                normalizer = self.sum(axis=self.Axis.HAPLOTYPE)
                normalizer.eliminate_zeros()

                for hid in range(self.num_haplotypes):
                    # trying to avoid numerical problem (inf or nan)
                    self.data[hid].eliminate_zeros()

                    # element-wise division
                    self.data[hid] = np.divide(self.data[hid], normalizer)
            elif axis == self.Axis.HAPLOTYPE:
                # haplotype-wise normalization on each read
                for hid in range(self.num_haplotypes):
                    # 1-dim Sparse matrix of |reads| x 1
                    normalizer = self.data[hid].sum(axis=self.Axis.HAPLOTYPE)
                    normalizer = normalizer.A.flatten()
                    self.data[hid].data /= normalizer[self.data[hid].indices]
            elif axis == self.Axis.READ:
                # normalization each read as a whole
                sum_mat = self.sum(axis=self.Axis.LOCUS)
                normalizer = sum_mat.sum(axis=self.Axis.HAPLOTYPE)
                normalizer = normalizer.ravel()

                for hid in range(self.num_haplotypes):
                    self.data[hid].data /= normalizer[self.data[hid].indices]
            elif axis == self.Axis.GROUP:
                # group-wise normalization on each read
                if grouping_mat is None:
                    raise RuntimeError('Group information matrix is missing.')

                normalizer = self.sum(axis=self.Axis.HAPLOTYPE) * grouping_mat

                for hid in range(self.num_haplotypes):
                    # trying to avoid numerical problem (inf or nan)
                    self.data[hid].eliminate_zeros()
                    self.data[hid] = np.divide(self.data[hid], normalizer)
            elif axis == self.Axis.HAPLOGROUP:
                # haplotype-wise & group-wise normalization on each read
                if grouping_mat is None:
                    raise RuntimeError('Group information matrix is missing.')

                for hid in range(self.num_haplotypes):
                    # normalizer is different hap-by-hap

                    # Sparse matrix of |reads| x |loci|
                    normalizer = (self.data[hid] * grouping_mat)

                    # Trying to avoid numerical problem (inf or nan)
                    self.data[hid].eliminate_zeros()
                    self.data[hid] = np.divide(self.data[hid], normalizer)
            else:
                raise RuntimeError('The axis should be 0, 1, 2, or 3.')
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def pull_alignments_from(self, reads_to_use: np.ndarray, shallow: bool = False) -> Self:
        """
        Extract alignments for a subset of reads.

        A new matrix containing only the alignments for the specified reads is created. It
        filters the matrix along the read dimension while preserving the locus and haplotype
        structure.

        Args:
            reads_to_use: Boolean array of length num_reads specifying which reads to include
                (True) or exclude (False).
            shallow: If True, copy only the sparse matrix data without any metadata. If False, copy
                all metadata and update read names and IDs.

        Returns:
            New matrix with filtered read data

        Raises:
            RuntimeError: If reads_to_use length doesn't match num_reads

        Note:
            This method is commonly used in conjunction with get_unique_reads() to analyze specific
            subsets of reads, such as uniquely-aligning reads or reads with specific properties.
        """
        new_alnmat = self.copy(shallow=shallow)

        for hid in range(self.num_haplotypes):
            hdata = new_alnmat.data[hid]
            hdata.data *= reads_to_use[hdata.indices]
            hdata.eliminate_zeros()

        if new_alnmat.count is not None:
            new_alnmat.count[np.logical_not(reads_to_use)] = 0

        return new_alnmat


    def get_unique_reads(
        self,
        ignore_haplotype: bool = False,
        shallow: bool = False
    ) -> Self:
        """
        Extract reads that align to exactly one locus-haplotype combination, which are crucial for
        accurate expression quantification in RNA-seq analysis.

        Args:
            ignore_haplotype: If True, consider reads as unique if they align to only one locus
                regardless of haplotype.  If False, reads must align to exactly one locus-haplotype
                combination to be considered unique.
            shallow: If True, copy only the sparse matrix data without metadata. If False, copy all
                metadata.

        Returns:
            New matrix containing only unique reads

        Raises:
            RuntimeError: If matrix is not finalized
        """
        if self.finalized:
            if ignore_haplotype:
                summat = self.sum(axis=self.Axis.HAPLOTYPE)
                nnz_per_read = np.diff(summat.tocsr().indptr)
                unique_reads = np.logical_and(nnz_per_read > 0, nnz_per_read < 2)
            else:
                # allelic multireads should be removed
                alncnt_per_read = self.sum(axis=self.Axis.LOCUS).sum(axis=self.Axis.HAPLOTYPE)
                unique_reads = np.logical_and(alncnt_per_read > 0, alncnt_per_read < 2)
            return self.pull_alignments_from(unique_reads, shallow=shallow)
        else:
            raise RuntimeError('The matrix is not finalized.')


    def count_unique_reads(
        self,
        ignore_haplotype: bool = False
    ) -> np.ndarray:
        """
        Count uniquely-aligning reads per locus or locus-haplotype combination.

        Args:
            ignore_haplotype: If True, count reads unique to each locus regardless of haplotype.
                If False, count reads unique to each locus-haplotype combination.

        Returns:
            Count array with shape.

        Raises:
            RuntimeError: If matrix is not finalized

        Note:
            When count data is available, it is used to weight the counts. Otherwise, binary
            incidence is assumed (each alignment counts as 1).
        """
        if self.finalized:
            unique_reads = self.get_unique_reads(ignore_haplotype=ignore_haplotype, shallow=True)

            if ignore_haplotype:
                numaln_per_read = unique_reads.sum(axis=self.Axis.HAPLOTYPE)

                if self.count is None:
                    numaln_per_read.data = np.ones(numaln_per_read.nnz)
                else:
                    numaln_per_read.data = self.count[numaln_per_read.indices]

                # an array of size |num_loci|
                return numaln_per_read.sum(axis=0).A.ravel()
            else:
                # an array of size |num_haplotypes|x|num_loci|
                return unique_reads.sum(axis=self.Axis.READ)
        else:
            raise RuntimeError('The matrix is not finalized.')


    def count_alignments(self) -> np.ndarray:
        """
        Count all alignments (both unique and multireads) for each locus-haplotype combination,
        providing the raw alignment counts before uniqueness filtering.

        Returns:
            Alignment count array with shape (num_haplotypes, num_loci) containing total alignment
            counts for each locus-haplotype combination.

        Raises:
            RuntimeError: If matrix is not finalized

        Note:
            This method returns the same result as sum(axis=Axis.READ) but is provided for clarity
            and convenience.
        """
        if self.finalized:
            return self.sum(axis=self.Axis.READ)
        else:
            raise RuntimeError('The matrix is not finalized.')


    def report_alignment_counts(
        self,
        filename: str
    ) -> None:
        """
        Generate a comprehensive alignment count report. Creates a tab-separated text file
        containing detailed alignment statistics for each locus, including total alignments,
        unique alignments per haplotype, and locus-level unique counts.

        Args:
            filename: Path to the output file where the report will be saved

        Raises:
            RuntimeError: If matrix is not finalized or names are not available

        Output Format:
            locus    aln_A    aln_B    aln_C    uniq_A    uniq_B    uniq_C    locus_uniq
            Gene1    100      95       105      80        75        85        240
            Gene2    50       45       55       40        35        45        120
        """
        alignment_counts = self.count_alignments()
        allelic_unique_counts = self.count_unique_reads(ignore_haplotype=False)
        locus_unique_counts = self.count_unique_reads(ignore_haplotype=True)
        cntdata = np.vstack((alignment_counts, allelic_unique_counts))
        cntdata = np.vstack((cntdata, locus_unique_counts))

        fhout = open(filename, 'w')
        fhout.write('locus\t' + '\t'.join([f'aln_{h}' for h in self.hname]) + '\t')
        fhout.write('\t'.join([f'uniq_{h}' for h in self.hname]) + '\t')
        fhout.write('locus_uniq' + '\n')

        for locus_id in range(self.num_loci):
            lout = [self.lname[locus_id]]
            lout.extend(list(map(str, cntdata[:, locus_id].ravel())))
            fhout.write('\t'.join(lout))
            fhout.write('\n')

        fhout.close()


    def combine(
        self,
        other: Self,
        shallow: bool = False
    ) -> Self:
        """
        Combine two AlignmentPropertyMatrix objects along the read dimension. This method
        concatenates two matrices along the read axis, effectively merging their read data while
        preserving the locus and haplotype structure. This is useful for combining data from
        multiple samples or sequencing runs.

        Args:
            other: Matrix to combine with self.
            shallow: If True, copy only the sparse matrix data without metadata. If False, copy all
                metadata and update read names and IDs.

        Returns:
            Combined matrix with concatenated read data

        Raises:
            RuntimeError: If either matrix is not finalized or shapes are incompatible

        Note:
            The matrices must have compatible shapes:
            - Same number of loci (num_loci)
            - Same number of haplotypes (num_haplotypes)
            - Read dimensions are concatenated (num_reads = self.num_reads + other.num_reads)
        """
        if self.finalized and other.finalized:
            dmat = Sparse3DMatrix.combine(self, other)
            dmat.num_loci, dmat.num_haplotypes, dmat.num_reads = dmat.shape

            if self.count is not None and other.count is not None:
                dmat.count = np.concatenate((self.count, other.count))

            if not shallow:
                dmat.hname = self.hname
                dmat.lname = copy.copy(self.lname)
                dmat.rname = np.concatenate((self.rname, other.rname))
                dmat.lid = copy.copy(self.lid)
                dmat.rid = dict(zip(dmat.rname, np.arange(dmat.num_reads)))
                dmat.__copy_group_info(self)

            return dmat
        else:
            raise RuntimeError('Both matrices must be finalized.')


    def save(
        self,
        h5_file: str,
        title: str | None = None,
        index_dtype: str = 'uint32',
        data_dtype: type = float,
        incidence_only: bool = True,
        complib: str = 'zlib',
        shallow: bool = False
    ) -> None:
        """
        Save the AlignmentPropertyMatrix to an HDF5 file. The file includes the sparse matrix data,
        count information, and metadata about loci, haplotypes, and reads.

        Args:
            h5_file: Path to the output HDF5 file
            title: Title for the HDF5 file.
            index_dtype: Data type for matrix indices.
            data_dtype: Data type for matrix values.
            incidence_only: If True, store only binary incidence (0/1). If False, store actual
                alignment values.
            complib: Compression library to use ('zlib', 'lzo', 'bzip2', 'blosc').
            shallow: True to store lname and rname.

        Raises:
            RuntimeError: If matrix is not finalized

        Note:
            The HDF5 file structure includes:
            - /data: Sparse matrix data (inherited from Sparse3DMatrix)
            - /count: Equivalence class counts (if available)
            - /lname: Locus names (if shallow=False)
            - /rname: Read names (if shallow=False and available)
            - hname: Haplotype names as file attribute (if shallow=False)
        """
        Sparse3DMatrix.save(
            self,
            h5_file=h5_file,
            title=title,
            index_dtype=index_dtype,
            data_dtype=data_dtype,
            incidence_only=incidence_only,
            complib=complib
        )

        h5fh = tables.open_file(h5_file, 'a')
        fil = tables.Filters(complevel=1, complib=complib)

        if self.count is not None:
            h5fh.create_carray(
                h5fh.root,
                'count',
                obj=self.count,
                title='Equivalence Class Counts',
                filters=fil,
            )

        if not shallow:
            h5fh.set_node_attr(h5fh.root, 'hname', self.hname)
            h5fh.create_carray(
                h5fh.root,
                'lname',
                obj=self.lname,
                title='Locus Names',
                filters=fil,
            )

            if self.rname is not None:
                h5fh.create_carray(
                    h5fh.root,
                    'rname',
                    obj=self.rname,
                    title='Read Names',
                    filters=fil,
                )

        h5fh.flush()
        h5fh.close()


    def get_read_data(self, rid: int) -> spmatrix:
        """
        Extract alignment data for a specific read. Returns a 2D sparse matrix containing all
        alignment information for the specified read across all loci and haplotypes.

        Args:
            rid: Read ID (index) to extract data for

        Returns:
            2D sparse matrix with shape (num_haplotypes, num_loci) containing alignment data for
                the specified read.

        Raises:
            IndexError: If rid is out of range
            RuntimeError: If matrix is not finalized
        """
        return self.get_cross_section(index=rid, axis=self.Axis.READ)


    def print_read(self, rid: int) -> None:
        """
        Print detailed alignment information for a specific read.  Prints a human-readable summary
        of all alignments for the specified read, showing which loci it aligns to and the
        corresponding alignment values for each haplotype.

        Args:
            rid: Read ID (index) to print information for

        Raises:
            IndexError: If rid is out of range
            RuntimeError: If matrix is not finalized
        """
        if self.rname is not None:
            print(self.rname[rid])
            print('--')

        r = self.get_read_data(rid)
        aligned_loci = np.unique(r.nonzero()[1])

        for locus in aligned_loci:
            nzvec = r[:, locus].todense().transpose()[0].A.flatten()

            if self.lname is not None:
                print(self.lname[locus], end=' ')
            else:
                print(locus, end=' ')

            print(nzvec)

    #
    # For future use
    #
    def get_reads_aligned_to_locus(self, lid: int, hid: int | None = None) -> list[int]:
        """
        Get read IDs that align to a specific locus or locus-haplotype combination. It can filter
        by haplotype or return reads aligned to the locus across all haplotypes.

        Args:
            lid: Locus ID (index) to query
            hid: Haplotype ID (index) to filter by. If None, returns reads aligned to the locus
                across all haplotypes.

        Returns:
            Sorted list of read IDs that align to the specified locus

        Raises:
            IndexError: If lid or hid is out of range
            RuntimeError: If matrix is not finalized
        """
        ridset = set()
        if hid is None:
            for hid in range(self.num_haplotypes):
                curset = set(np.nonzero(self.data[hid][:, lid])[0])
                ridset = ridset.union(curset)
            return sorted(list(ridset))
        else:
            return sorted(np.nonzero(self.data[hid][:, lid])[0])
