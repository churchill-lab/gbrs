# Standard library imports
import copy
from numbers import Number
from typing import Self

# Third-party library imports
import numpy as np
import tables
from scipy.sparse import coo_matrix, csc_matrix, csr_matrix, hstack, lil_matrix, vstack

# Local library imports
from gbrs import utils

logger = utils.get_logger('gbrs')


class Sparse3DMatrix:
    """
    A 3-dimensional sparse matrix designed for efficient storage and manipulation of RNA-seq
    alignment data.
    
    This class implements a sparse 3D matrix structure optimized for RNA-seq alignments, where the
    three dimensions represent:
    - Axis 0: Loci (transcripts/genes)
    - Axis 1: Haplotypes (founder strains)
    - Axis 2: Reads (sequencing reads)
    
    The matrix uses scipy sparse matrices internally for memory efficiency, storing only non-zero
    elements. This is particularly important for RNA-seq data where most reads align to only a small
    subset of possible loci-haplotype combinations.
    
    Key Features:
    - Memory-efficient sparse storage with CSC (Compressed Sparse Column) format
    - Support for binary operations (addition, subtraction, multiplication)
    - Cross-sectional slicing along any axis
    - HDF5 file I/O for persistence
    - Two-phase workflow: population then finalization
    
    The matrix follows a two-phase workflow:
    **Population phase**: Add data using set_value() or add_value()
    **Finalization phase**: Call finalize() to convert to optimized CSC format
    
    Attributes:
        shape: 3-tuple of (num_loci, num_haplotypes, num_reads)
        ndim: Always 3 for 3D matrix
        data: List of scipy.sparse.csc_matrix objects, one per haplotype
        finalized: Whether the matrix has been converted to CSC format
    """

    def __init__(
        self,
        other: Self | None = None,
        h5_file: str | None = None,
        datanode: str = '/',
        shape: tuple[int, int, int] | None = None,
        dtype: type = float
    ) -> None:
        """
        Initialize a Sparse3DMatrix object.
        
        This constructor supports multiple initialization modes:
        
        1. **Copy from existing matrix**: Provide 'other' parameter
        2. **Load from HDF5 file**: Provide 'h5file' parameter  
        3. **Create new empty matrix**: Provide 'shape' parameter
        
        Args:
            other: Existing matrix to copy from. The source matrix must be finalized.
            h5_file: Path to HDF5 file containing matrix data in EMASE format.
            datanode: HDF5 node path containing matrix data. Defaults to '/'
            shape: 3-tuple specifying (num_loci, num_haplotypes, num_reads).
                Required for creating new empty matrix
            dtype: Data type for matrix elements.
                
        Raises:
            RuntimeError: If copying from non-finalized matrix or invalid shape provided
            ValueError: If shape is not a 3-tuple of positive integers
        """
        self.shape = (0, 0, 0)
        self.ndim = 3
        #self.data: list[csc_matrix] = []
        self.data: list = list()
        self.finalized = False
        
        # copy from existing matrix
        if other is not None:
            if other.finalized:
                self.shape = other.shape
                self.data = copy.deepcopy(other.data)
                self.finalized = True
            else:
                raise RuntimeError('The original matrix must be finalized.')
        
        # load from HDF5 file
        elif h5_file is not None:
            h5fh = tables.open_file(h5_file, 'r')
            self.shape = h5fh.get_node_attr(datanode, 'shape')
            for hid in range(self.shape[1]):
                self.data.append(self._reconstruct_spmat(h5fh, hid, datanode, dtype))
            h5fh.close()
            self.finalize()  # Convert to CSC format
        
        # create new empty matrix
        elif shape is not None:
            if len(shape) != 3:
                raise ValueError('The shape must be a tuple of three positive integers.')
            if any(dim < 1 for dim in shape):
                raise ValueError('The shape must be a tuple of three positive integers.')
            
            self.shape = shape
            for hid in range(self.shape[1]):
                self.data.append(
                    lil_matrix((self.shape[2], self.shape[0]), dtype=dtype)
                )

    def _reconstruct_spmat(
        self,
        h5_fh: tables.File,
        hid: int,
        datanode: str,
        dtype: type
    ) -> csc_matrix | coo_matrix:
        """
        Reconstruct a sparse matrix from HDF5 file data for a specific haplotype.
        
        This private method reads sparse matrix data from an HDF5 file and reconstructs the
        corresponding scipy sparse matrix object. It supports both CSC and COO formats.
        
        Args:
            h5_fh: Open HDF5 file handle.
            hid: Haplotype index to reconstruct.
            datanode: HDF5 node path containing matrix data.
            dtype: Data type for matrix elements.
            
        Returns:
            Reconstructed sparse matrix in CSC or COO format.
            
        Raises:
            RuntimeError: If unsupported matrix type is encountered
            
        Note:
            For CSC format, it reads indptr, indices, and data arrays.
            For COO format, it reads coordinate and data arrays.
            If 'incidence_only' is True, data values are set to 1.0.
        """
        try:
            mtype = h5_fh.get_node_attr(datanode, 'mtype')
            try:
                mtype = mtype.decode()
            except Exception:
                pass
            incidence_only = h5_fh.get_node_attr(datanode, 'incidence_only')
        except AttributeError:
            mtype = 'coo_matrix'
            incidence_only = False
        
        hapnode = h5_fh.get_node(f'{datanode}/h{hid}')
        
        if mtype == 'csc_matrix':
            indptr = h5_fh.get_node(hapnode, 'indptr').read().astype(int)
            indices = h5_fh.get_node(hapnode, 'indices').read().astype(int)
            
            if not incidence_only:
                data = h5_fh.get_node(hapnode, 'data').read().astype(dtype)
            else:
                data = np.ones(len(indices), dtype=dtype)
            
            spmat = csc_matrix((data, indices, indptr), shape=(self.shape[2], self.shape[0]))
        
        elif mtype == 'coo_matrix':
            coor = h5_fh.get_node(hapnode, 'coor').read()
            data = h5_fh.get_node(hapnode, 'data').read().astype(dtype)
            spmat = coo_matrix((data, coor), shape=(self.shape[2], self.shape[0]))
        
        else:
            raise RuntimeError('Only csc or coo matrices are supported.')
        
        return spmat


    def copy(self) -> Self:
        """
        Create a deep copy of the Sparse3DMatrix. The copy inherits the finalized state of the
        original matrix.
        
        Returns:
            A new matrix with identical data and properties
            
        Raises:
            RuntimeError: If the original matrix is not finalized
        """
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            dmat.data = copy.deepcopy(self.data)
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def __add__(
        self,
        other: Self | csc_matrix | csr_matrix | coo_matrix | lil_matrix
    ) -> Self:
        """
        Add another matrix or sparse matrix to this Sparse3DMatrix.
        
        This method supports element-wise addition between:
        - Two Sparse3DMatrix objects (both must be finalized)
        - Sparse3DMatrix and scipy sparse matrix (csc, csr, coo, lil)
        
        Args:
            other: Matrix to add. Can be Sparse3DMatrix or scipy sparse matrix
            
        Returns:
            New matrix containing the sum.
            
        Raises:
            RuntimeError: If matrices are not finalized or incompatible
            TypeError: If unsupported type is provided
        """
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            
            if isinstance(other, Sparse3DMatrix):
                if other.finalized:
                    for hid in range(self.shape[1]):
                        dmat.data.append(self.data[hid] + other.data[hid])
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            
            elif isinstance(other, (csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
                other_csc = other.tocsc()

                for hid in range(self.shape[1]):
                    dmat.data.append(self.data[hid] + other_csc)
            
            else:
                raise TypeError('This operator is not supported between the given types.')
            
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def __sub__(
        self,
        other: Self | csc_matrix | csr_matrix | coo_matrix | lil_matrix
    ) -> Self:
        """
        Subtract another matrix or sparse matrix from this Sparse3DMatrix.
        
        This method supports element-wise subtraction between:
        - Two Sparse3DMatrix objects (both must be finalized)
        - Sparse3DMatrix and scipy sparse matrix (csc, csr, coo, lil)
        
        Args:
            other: Matrix to subtract. Can be Sparse3DMatrix or scipy sparse
                matrix
            
        Returns:
            New matrix containing the difference.
            
        Raises:
            RuntimeError: If matrices are not finalized or incompatible
            TypeError: If unsupported type is provided
        """
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            
            if isinstance(other, Sparse3DMatrix):
                if other.finalized:
                    for hid in range(self.shape[1]):
                        dmat.data.append(self.data[hid] - other.data[hid])
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            
            elif isinstance(other, (csc_matrix, csr_matrix, coo_matrix, lil_matrix)):
                other_csc = other.tocsc()
                for hid in range(self.shape[1]):
                    dmat.data.append(self.data[hid] - other_csc)
            
            else:
                raise TypeError('This operator is not supported between the given types.')
            
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def __mul__(
        self,
        other: Self | np.ndarray | csc_matrix | csr_matrix | coo_matrix | lil_matrix | Number
    ) -> Self:
        """
        Multiply this Sparse3DMatrix by another matrix, sparse matrix, or scalar.
        
        This method supports multiple multiplication modes:
        - Element-wise multiplication with another Sparse3DMatrix
        - Matrix multiplication with numpy array or scipy sparse matrix
        - Scalar multiplication (rescaling)
        
        Args:
            other: Multiplier. Can be:
                - Sparse3DMatrix: Element-wise multiplication
                - numpy.ndarray or scipy sparse matrix: Matrix multiplication
                - Number: Scalar multiplication
                
        Returns:
            Result of multiplication.
            
        Raises:
            RuntimeError: If matrices are not finalized or incompatible
            TypeError: If unsupported type is provided
            
        Note:
            Matrix multiplication changes the shape of the result matrix.
            For matrix multiplication, the result shape becomes:
            (other.shape[1], self.shape[1], self.shape[2])
        """
        if self.finalized:
            dmat = self.__class__()
            dmat.shape = self.shape
            
            if isinstance(other, Sparse3DMatrix):
                # element-wise multiplication between same kind
                if other.finalized:
                    for hid in range(self.shape[1]):
                        dmat.data.append(self.data[hid].multiply(other.data[hid]))
                else:
                    raise RuntimeError('Both matrices must be finalized.')
            
            elif isinstance(other, (np.ndarray, csc_matrix, csr_matrix)):
                # matrix-matrix multiplication
                for hid in range(self.shape[1]):
                    dmat.data.append(self.data[hid] * other)
                dmat.shape = (other.shape[1], self.shape[1], self.shape[2])
            
            elif isinstance(other, (coo_matrix, lil_matrix)):
                # matrix-matrix multiplication
                other_csc = other.tocsc()
                for hid in range(self.shape[1]):
                    dmat.data.append(self.data[hid] * other_csc)
                dmat.shape = (other_csc.shape[1], self.shape[1], self.shape[2])
            
            elif isinstance(other, Number):
                # rescaling of matrix
                for hid in range(self.shape[1]):
                    dmat.data.append(self.data[hid] * other)
            
            else:
                raise TypeError('This operator is not supported between the given types.')
            
            dmat.finalized = True
            return dmat
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def set_value(self, lid: int, hid: int, rid: int, value: float) -> None:
        """
        Set a specific value in the matrix at the given coordinates.
        
        This method sets the value at position (lid, hid, rid) to the specified value.  Use this
        during the population phase before calling finalize().
        
        Args:
            lid: Locus index (axis 0).
            hid: Haplotype index (axis 1).
            rid: Read index (axis 2).
            value: Value to set at the specified position.
            
        Raises:
            RuntimeError: If the matrix has not been properly initialized with shape
        """
        if np.all(self.shape > 0):
            self.data[hid][rid, lid] = value
        else:
            raise RuntimeError('The Sparse3DMatrix has only been declared.')

    def add_value(self, lid: int, hid: int, rid: int, value: float) -> None:
        """
        Add a value to the existing value at the given coordinates.
        
        This method adds the specified value to the current value at position (lid, hid, rid). Use
        this during the population phase before calling finalize().
        
        Args:
            lid: Locus index (axis 0)
            hid: Haplotype index (axis 1)
            rid: Read index (axis 2) 
            value: Value to add to the current value at the specified position
            
        Raises:
            RuntimeError: If the matrix has not been properly initialized with shape
        """
        if np.all(self.shape > 0):
            self.data[hid][rid, lid] += value
        else:
            raise RuntimeError('The Sparse3DMatrix has only been declared.')


    def reset(self) -> None:
        """
        Reset all non-zero values in the matrix to 1.0.
        
        This method sets all non-zero elements in the matrix to 1.0, effectively converting the
        matrix to an incidence matrix (binary matrix).
        
        Raises:
            RuntimeError: If the matrix is not finalized
            
        Note:
            This is useful for converting alignment count matrices to binary incidence matrices,
            where only the presence/absence of alignments is important, not their counts.
        """
        if self.finalized:
            for hid in range(self.shape[1]):
                # TODO: inherit the dtype from orig
                self.data[hid].data = np.ones(self.data[hid].nnz, dtype=self.data[hid].dtype)
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def finalize(self) -> None:
        """
        Convert the matrix to optimized CSC format.
        
        This method converts all internal sparse matrices from LIL (List of Lists) format to
        CSC (Compressed Sparse Column) format for optimal performance.

        After finalization, the matrix becomes read-only for element-wise operations.
        
        Raises:
            RuntimeError: If the matrix has not been properly initialized
            
        Note:
            Finalization is a required step before performing matrix operations like addition,
            subtraction, or multiplication. It optimizes memory usage and computational performance
            for large matrices.
            
            After finalization, you can no longer use set_value() or add_value(). Use reset() to
            modify values if needed.
        """
        if not self.finalized:
            for hid in range(self.shape[1]):
                self.data[hid] = self.data[hid].tocsc()
            self.finalized = True

    def sum(
        self,
        axis: int = 2
    ) -> np.ndarray | csc_matrix:
        """
        Sum the matrix along a specified axis.
        
        This method performs summation along one of the three dimensions:
            loci (axis 0), haplotypes (axis 1), or reads (axis 2).
        
        Args:
            axis: Axis along which to sum:
                - 0: Sum across loci for each haplotype-read pair
                - 1: Sum across haplotypes for each locus-read pair  
                - 2: Sum across reads for each locus-haplotype pair (default)
                
        Returns:
            Summed data with shape depending on the axis:
            - Axis 0: (num_reads, num_haplotypes) dense array
            - Axis 1: (num_loci, num_reads) sparse matrix
            - Axis 2: (num_haplotypes, num_loci) dense array
                
        Raises:
            RuntimeError: If matrix is not finalized or invalid axis specified
        """
        if self.finalized:
            if axis == 0:
                # sum along loci
                sum_mat = []
                for hid in range(self.shape[1]):
                    sum_mat.append(self.data[hid].sum(axis=1).A)
                sum_mat = np.hstack(sum_mat)
            elif axis == 1:
                # sum along haplotypes
                sum_mat = self.data[0]
                for hid in range(1, self.shape[1]):
                    # still sparse matrix
                    sum_mat = sum_mat + self.data[hid]
            elif axis == 2:
                # sum along reads
                sum_mat = []
                for hid in range(self.shape[1]):
                    sum_mat.append(self.data[hid].sum(axis=0).A)
                sum_mat = np.vstack(sum_mat)
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
            return sum_mat
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def get_cross_section(
        self,
        index: int,
        axis: int = 0
    ) -> csc_matrix:
        """
        Get a 2D cross-section of the 3D matrix along a specified axis.
        
        This method extracts a 2D slice from the 3D matrix, returning a sparse matrix representing
        the cross-section at the specified index along the given axis.
        
        Args:
            index: Index along the specified axis for the cross-section
            axis: Axis along which to take the cross-section:
                - 0: Cross-section across loci (returns haplotype × read matrix)
                - 1: Cross-section across haplotypes (returns locus × read matrix)
                - 2: Cross-section across reads (returns locus × haplotype matrix)
                
        Returns:
            2D sparse matrix in CSC format representing the cross-section
            
        Raises:
            RuntimeError: If matrix is not finalized or invalid axis specified
        """
        if self.finalized:
            if axis == 0:
                cols = []
                for hid in range(len(self.data)):
                    cols.append(self.data[hid][:, index])
                return hstack(cols).tocsc()
            elif axis == 1:
                return self.data[axis]
            elif axis == 2:
                rows = []
                for hid in range(len(self.data)):
                    rows.append(self.data[hid][index, :])
                return vstack(rows).tocsc()
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def add(
        self,
        addend_mat: Self,
        axis: int = 1
    ) -> Self:
        """
        Add another matrix along a specified axis.
        
        This method performs in-place addition of another Sparse3DMatrix along the specified axis.
        The operation modifies the current matrix directly without creating a new matrix.
        
        Args:
            addend_mat: Matrix to add along the specified axis
            axis: Axis along which to add:
                - 0: Add along loci dimension (not implemented)
                - 1: Add along haplotypes dimension (default)
                - 2: Add along reads dimension (not implemented)
                
        Returns:
            Self (for method chaining)
            
        Raises:
            RuntimeError: If matrix is not finalized or invalid axis specified
            NotImplementedError: If axis 0 or 2 is specified (not yet implemented)
            
        Note:
            This method modifies the matrix in-place. Only axis 1 (haplotypes) is currently
            implemented. The addend_mat should be compatible with the current matrix structure.
        """
        if self.finalized:
            if axis == 0:
                raise NotImplementedError('The method is not yet implemented for the axis.')
            elif axis == 1:
                for hid in range(self.shape[1]):
                    self.data[hid] = self.data[hid] + addend_mat
            elif axis == 2:
                raise NotImplementedError('The method is not yet implemented for the axis.')
            else:
                raise RuntimeError('The axis should be 0, 1, or 2.')
        else:
            raise RuntimeError('The original matrix must be finalized.')

    def multiply(
        self,
        multiplier: np.ndarray | csc_matrix | csr_matrix | coo_matrix | lil_matrix,
        axis: int | None = None
    ) -> Self:
        """
        Multiply the matrix by a multiplier along a specified axis.
        
        This method performs in-place multiplication of the matrix by different types of multipliers
        along specified axes. The operation modifies the current matrix directly without creating
        a new matrix. This is a key method used in the EMASE algorithm for updating alignment
        probabilities.
        
        Args:
            multiplier: Multiplier matrix or array. Can be:
                - 1D numpy array: Vector multiplication along specified axis
                - 2D numpy array: Matrix multiplication along specified axis
                - scipy sparse matrix: Sparse matrix multiplication
                - Sparse3DMatrix: Element-wise multiplication with another 3D matrix
            axis: Axis along which to multiply:
                - None: Not used (for Sparse3DMatrix multipliers)
                - 0: Multiply along loci dimension (not implemented for 1D)
                - 1: Multiply along haplotypes dimension
                - 2: Multiply along reads dimension
                
        Returns:
            Self (for method chaining)
            
        Raises:
            RuntimeError: If matrix is not finalized, invalid axis, or incompatible multiplier
            NotImplementedError: If axis 0 is specified for 1D multipliers
            
        Note:
            This method modifies the matrix in-place and is critical for the EMASE algorithm.
            
            For 1D multipliers:
            - Axis 1: Multiplier length must match number of loci
            - Axis 2: Multiplier length must match number of reads
            
            For 2D multipliers:
            - Axis 0: Shape should be (num_reads, num_haplotypes)
            - Axis 1: Shape should be (num_reads, num_loci)  
            - Axis 2: Shape should be (num_haplotypes, num_loci) - used in EM algorithm
            
            For Sparse3DMatrix: Performs element-wise multiplication with compatible 3D matrix
        """
        if self.finalized:
            if multiplier.ndim == 1:
                if axis == 0:
                    # multiplier is np.array of length |haplotypes|
                    raise NotImplementedError('The method is not yet implemented for the axis.')
                elif axis == 1:
                    # multiplier is np.array of length |loci|
                    sz = len(multiplier)
                    multiplier_mat = lil_matrix((sz, sz))
                    multiplier_mat.setdiag(multiplier)
                    for hid in range(self.shape[1]):
                        self.data[hid] = self.data[hid] * multiplier_mat
                elif axis == 2:
                    # multiplier is np.array of length |reads|
                    for hid in range(self.shape[1]):
                        self.data[hid].data *= multiplier[self.data[hid].indices]
                else:
                    raise RuntimeError('The axis should be 0, 1, or 2.')
            elif multiplier.ndim == 2:
                if axis == 0:
                    # multiplier is sp.sparse matrix of shape |reads| x |haplotypes|
                    for hid in range(self.shape[1]):
                        self.data[hid].data *= multiplier[self.data[hid].indices, hid]
                elif axis == 1:
                    # multiplier is sp.sparse matrix of shape |reads| x |loci|
                    for hid in range(self.shape[1]):
                        self.data[hid] = self.data[hid].multiply(multiplier)
                elif axis == 2:
                    # multiplier is np.matrix of shape |haplotypes| x |loci|
                    for hid in range(self.shape[1]):
                        multiplier_vec = multiplier[hid, :]
                        multiplier_vec = multiplier_vec.ravel()
                        self.data[hid].data *= multiplier_vec.repeat(np.diff(self.data[hid].indptr))
                else:
                    raise RuntimeError('The axis should be 0, 1, or 2.')
            elif isinstance(multiplier, Sparse3DMatrix):
                # multiplier is Sparse3DMatrix object
                for hid in range(self.shape[1]):
                    self.data[hid] = self.data[hid].multiply(multiplier.data[hid])
            else:
                raise RuntimeError(
                    'The multiplier should be 1, 2 dimensional numpy array or '
                    'a Sparse3DMatrix object.'
                )
        else:
            raise RuntimeError('The original matrix must be finalized.')


    def combine(
        self,
        other: Self
    ) -> Self:
        """
        Combine this matrix with another matrix along the read dimension.
        
        This method combines two Sparse3DMatrix objects along the read dimension, effectively
        concatenating their read data. This is useful for merging data from different samples
        or experiments.
        
        Args:
            other: Matrix to combine with along the read dimension
            
        Returns:
            New matrix with combined read data
            
        Raises:
            RuntimeError: If matrices are not finalized or incompatible

        Note:
            The matrices must have the same locus and haplotype dimensions.
            The resulting matrix will have the sum of the read dimensions.
        """
        if self.finalized and other.finalized:
            dmat = self.__class__()
            dmat.shape = (self.shape[0], self.shape[1], self.shape[2] + other.shape[2])
            
            for hid in range(self.shape[1]):
                dmat.data.append(vstack([self.data[hid], other.data[hid]]))
            
            dmat.finalized = True
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
        complib: str = 'zlib'
    ) -> None:
        """
        Save the matrix to an HDF5 file in EMASE format.
        
        This method saves the Sparse3DMatrix to an HDF5 file using the EMASE format, which is
        compatible with the GBRS pipeline and other EMASE-based tools.
        
        Args:
            h5_file: Path to the output HDF5 file.
            title: Title for the HDF5 file. If None, uses default.
            index_dtype: Data type for sparse matrix indices (indptr, indices).
            data_dtype: Data type for sparse matrix data. Defaults to float
            incidence_only: If True, store only binary incidence (1.0 for alignments). If False,
                store actual values.
            complib: Compression library to use.
                
        Raises:
            RuntimeError: If matrix is not finalized
            
        Note:
            The output HDF5 file contains:
            - Root attributes: shape, mtype, incidence_only
            - Groups: h0, h1, ... (one per haplotype)
            - Each haplotype group contains:
                indptr, indices,  data (if not incidence_only)
            
            The file format is compatible with AlignmentPropertyMatrix.load().
        """
        if not self.finalized:
            raise RuntimeError('The matrix must be finalized before saving.')
        
        logger.debug(f'Saving matrix to: {h5_file}')
        h5fh = tables.open_file(h5_file, 'w', title=title or 'Sparse3DMatrix')
        fil = tables.Filters(complevel=1, complib=complib)
        
        # set root attributes
        h5fh.set_node_attr(h5fh.root, 'shape', self.shape)
        h5fh.set_node_attr(h5fh.root, 'mtype', 'csc_matrix')
        h5fh.set_node_attr(h5fh.root, 'incidence_only', incidence_only)
        
        # save each haplotype matrix
        for hid in range(self.shape[1]):
            logger.debug(f'Saving haplotype {hid}')
            hgroup = h5fh.create_group(
                h5fh.root,
                f'h{hid}',
                f'Sparse matrix components for Haplotype {hid}',
            )
            
            # save sparse matrix components
            h5fh.create_carray(
                hgroup,
                'indptr',
                obj=self.data[hid].indptr.astype(index_dtype),
                filters=fil,
            )
            
            h5fh.create_carray(
                hgroup,
                'indices',
                obj=self.data[hid].indices.astype(index_dtype),
                filters=fil,
            )
            
            if not incidence_only:
                h5fh.create_carray(
                    hgroup,
                    'data',
                    obj=self.data[hid].data.astype(data_dtype),
                    filters=fil,
                )
        
        h5fh.flush()
        h5fh.close()
        logger.debug('Matrix saved successfully')
