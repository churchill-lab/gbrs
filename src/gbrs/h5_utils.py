"""
GBRS HDF5 Utilities

This module provides utilities for working with HDF5 files in the GBRS pipeline,
including file comparison, inspection, and validation tools.

Key Functions:
    compare_h5_files: Compare two HDF5 files for structural and content differences
    compare_h5_structure: Compare only the structure (nodes, attributes) of HDF5 files
    compare_h5_content: Compare only the content (array data) of HDF5 files
    h5_inspect: Inspect and display contents of EMASE HDF5 files
    get_file_info: Get detailed information about HDF5 file structure and contents

File Comparison:
    The module provides both strict and tolerant comparison modes:
    - Strict: Exact equality checking for all data
    - Tolerant: Allows small numerical differences for floating-point data
    
    Comparisons include:
    - File structure (nodes, groups, attributes)
    - Array shapes, dtypes, and values
    - Node attributes and metadata

HDF5 Inspection:
    Tools for examining EMASE format HDF5 files:
    - Display file structure and metadata
    - Show sparse matrix contents
    - Validate file format and integrity
    - Debug alignment data issues

EMASE Format Support:
    Specialized functions for working with EMASE format HDF5 files:
    - Haplotype-specific matrix inspection
    - Sparse matrix visualization
    - Alignment data validation
    
    EMASE files contain:
    - Root attributes: shape, mtype, incidence_only
    - Haplotype groups: /h0, /h1, ... with sparse matrices
    - Metadata arrays: lname, rname, hname
"""

import os
import numpy as np
import tables
import argparse
from typing import Dict, Any, List, Tuple, Optional
import logging

# Set up logging
from gbrs import utils

DATA_DIR = os.getenv('GBRS_DATA', '.')
logger = utils.get_logger('gbrs')



def get_array_size(node) -> int | None:
    """Get the size of an array node."""
    if hasattr(node, 'size'):
        return node.size
    elif hasattr(node, 'nrows'):
        if len(node.shape) == 1:
            return node.nrows
        else:
            return int(np.prod(node.shape))
    else:
        return None


def compare_attributes(node1, node2, path: str) -> list[str]:
    """Compare attributes between two nodes."""
    differences = []
    
    try:
        # Get attributes safely
        attrs1 = {}
        attrs2 = {}
        
        try:
            # Get attributes more safely by iterating through them
            for attr_name in node1._v_attrs._f_list():
                try:
                    attrs1[attr_name] = node1._v_attrs[attr_name]
                except Exception:
                    # Skip problematic attributes
                    pass
        except Exception as e:
            # some nodes might have weird attributes, skip them for now
            pass
            
        try:
            # Get attributes more safely by iterating through them
            for attr_name in node2._v_attrs._f_list():
                try:
                    attrs2[attr_name] = node2._v_attrs[attr_name]
                except Exception:
                    # Skip problematic attributes
                    pass
        except Exception as e:
            # some nodes might have weird attributes, skip them for now
            pass
        
        # compare attribute keys
        keys1 = set(attrs1.keys())
        keys2 = set(attrs2.keys())
        
        if keys1 != keys2:
            missing_in_2 = keys1 - keys2
            missing_in_1 = keys2 - keys1
            if missing_in_2:
                differences.append(f'{path}: Missing attributes in file2: {missing_in_2}')
            if missing_in_1:
                differences.append(f'{path}: Missing attributes in file1: {missing_in_1}')
        
        # compare common attributes
        common_keys = keys1 & keys2
        for key in common_keys:
            val1 = attrs1[key]
            val2 = attrs2[key]
            
            if isinstance(val1, np.ndarray) and isinstance(val2, np.ndarray):
                if not np.array_equal(val1, val2):
                    differences.append(f"{path}.{key}: Array values differ")
            elif val1 != val2:
                differences.append(f"{path}.{key}: Values differ - {val1} vs {val2}")
                
    except Exception as e:
        # skip problematic nodes instead of failing
        differences.append(f'{path}: Error comparing attributes: {e}')

    return differences


def compare_arrays(array1, array2, path: str, tolerant: bool = False, tolerance: float = 1e-6) -> list[str]:
    """Compare two arrays for equivalence."""
    differences = []
    
    try:
        if array1.shape != array2.shape:
            differences.append(f'{path}: Shapes differ - {array1.shape} vs {array2.shape}')
            return differences
        
        if array1.dtype != array2.dtype:
            differences.append(f'{path}: Dtypes differ - {array1.dtype} vs {array2.dtype}')
            return differences
        
        if array1.size != array2.size:
            differences.append(f'{path}: Sizes differ - {array1.size} vs {array2.size}')
            return differences
        
        # compare values
        if tolerant:
            # tolerant comparison for floating point????
            if np.issubdtype(array1.dtype, np.floating):
                if not np.allclose(array1, array2, rtol=tolerance, atol=tolerance):
                    max_diff = np.max(np.abs(array1 - array2))
                    differences.append(f'{path}: Values differ (tolerant comparison, max diff: {max_diff})')
            else:
                if not np.array_equal(array1, array2):
                    differences.append(f'{path}: Values differ (non-floating point)')
        else:
            if not np.array_equal(array1, array2):
                differences.append(f'{path}: Values differ (strict comparison)')
                    
    except Exception as e:
        differences.append(f'{path}: Error comparing arrays: {e}')
    
    return differences


def compare_h5_structure(file1: str, file2: str) -> (bool, list[str]):
    """Compare the structure of two H5 files."""
    differences = []
    
    try:
        with tables.open_file(file1, 'r') as f1, tables.open_file(file2, 'r') as f2:
            # compare root attributes
            root_diffs = compare_attributes(f1.root, f2.root, '/')
            differences.extend(root_diffs)
            
            # get all nodes from both files
            nodes1 = {node._v_pathname: node for node in f1.walk_nodes()}
            nodes2 = {node._v_pathname: node for node in f2.walk_nodes()}
            
            # check for missing nodes
            missing_in_2 = set(nodes1.keys()) - set(nodes2.keys())
            missing_in_1 = set(nodes2.keys()) - set(nodes1.keys())
            
            if missing_in_2:
                differences.append(f'Missing nodes in file2: {missing_in_2}')
            if missing_in_1:
                differences.append(f'Missing nodes in file1: {missing_in_1}')
            
            # compare common nodes
            common_nodes = set(nodes1.keys()) & set(nodes2.keys())
            
            for node_path in sorted(common_nodes):
                node1 = nodes1[node_path]
                node2 = nodes2[node_path]
                
                # compare node types
                if type(node1) != type(node2):
                    differences.append(f'{node_path}: Node types differ - {type(node1)} vs {type(node2)}')
                    continue
                
                # compare attributes
                attr_diffs = compare_attributes(node1, node2, node_path)
                differences.extend(attr_diffs)
                
                # compare array data if applicable
                if hasattr(node1, 'read') and hasattr(node2, 'read'):
                    try:
                        data1 = node1.read()
                        data2 = node2.read()
                        
                        # for now, we'll do strict comparison here, the main comparison will 
                        # be done in compare_h5_content
                        if not np.array_equal(data1, data2):
                            differences.append(f'{node_path}: Array data differs')
                    except Exception as e:
                        differences.append(f'{node_path}: Error reading array data: {e}')
    
    except Exception as e:
        differences.append(f"Error comparing file structure: {e}")
    
    return len(differences) == 0, differences


def compare_h5_content(file1: str, file2: str, tolerant: bool = False, tolerance: float = 1e-6) -> (bool, list[str]):
    """Compare the content of two H5 files."""
    differences = []
    
    try:
        with tables.open_file(file1, 'r') as f1, tables.open_file(file2, 'r') as f2:
            # cet all array nodes from both files
            arrays1 = {node._v_pathname: node for node in f1.walk_nodes() 
                      if hasattr(node, 'read') and hasattr(node, 'shape')}
            arrays2 = {node._v_pathname: node for node in f2.walk_nodes() 
                      if hasattr(node, 'read') and hasattr(node, 'shape')}
            
            # check for missing arrays
            missing_in_2 = set(arrays1.keys()) - set(arrays2.keys())
            missing_in_1 = set(arrays2.keys()) - set(arrays1.keys())
            
            if missing_in_2:
                differences.append(f'Missing arrays in file2: {missing_in_2}')
            if missing_in_1:
                differences.append(f'Missing arrays in file1: {missing_in_1}')
            
            # compare common arrays
            common_arrays = set(arrays1.keys()) & set(arrays2.keys())
            
            for array_path in sorted(common_arrays):
                try:
                    array1 = arrays1[array_path].read()
                    array2 = arrays2[array_path].read()
                    
                    array_diffs = compare_arrays(array1, array2, array_path, tolerant, tolerance)
                    differences.extend(array_diffs)
                    
                except Exception as e:
                    differences.append(f'{array_path}: Error comparing array: {e}')
    
    except Exception as e:
        differences.append(f'Error comparing file content: {e}')
    
    return len(differences) == 0, differences


def get_file_info(filename: str) -> Dict[str, Any]:
    """Get basic information about an H5 file."""
    info = {
        'filename': filename,
        'size': os.path.getsize(filename) / (1024**3),  # GB
        'nodes': [],
        'arrays': [],
        'total_array_size': 0
    }
    
    try:
        with tables.open_file(filename, 'r') as f:
            # get root attributes
            try:
                shape = f.get_node_attr('/', 'shape')
                info['shape'] = shape
                info['num_loci'], info['num_haplotypes'], info['num_reads'] = shape
            except:
                info['shape'] = None
            
            # list all nodes
            for node in f.walk_nodes():
                node_info = {
                    'path': node._v_pathname,
                    'type': type(node).__name__,
                    'size': get_array_size(node)
                }
                
                if hasattr(node, 'shape'):
                    node_info['shape'] = node.shape
                if hasattr(node, 'dtype'):
                    node_info['dtype'] = node.dtype
                
                info['nodes'].append(node_info)
                
                # track arrays specifically
                if hasattr(node, 'read') and hasattr(node, 'shape'):
                    info['arrays'].append(node_info)
                    if node_info['size']:
                        info['total_array_size'] += node_info['size']
    
    except Exception as e:
        info['error'] = str(e)
    
    return info


def compare_h5_files(file1: str, file2: str, tolerant: bool = True, tolerance: float = 1e-6) -> None:
    logger.warning(f'Comparing files:')
    logger.warning(f'  File 1: {file1}')
    logger.warning(f'  File 2: {file2}')
    logger.warning(f'  Mode: {"Tolerant" if tolerant else "Strict"}')
    if tolerant:
        logger.warning(f"  Tolerance: {tolerance}")
    
    
    # get file information
    logger.warning('=== File Information ===')
    info1 = get_file_info(file1)
    info2 = get_file_info(file2)
    
    logger.warning(f'File 1: {info1["filename"]}')
    logger.warning(f'  Size: {info1["size"]:.2f} GB')
    if info1['shape']:
        logger.warning(f'  Shape: {info1["shape"]}')
    logger.warning(f'  Arrays: {len(info1["arrays"])}')
    logger.warning(f'  Total array size: {info1["total_array_size"]:,}')
    
    logger.warning(f'File 2: {info2["filename"]}')
    logger.warning(f'  Size: {info2["size"]:.2f} GB')
    if info2['shape']:
        logger.warning(f'  Shape: {info2["shape"]}')
    logger.warning(f'  Arrays: {len(info2["arrays"])}')
    logger.warning(f'  Total array size: {info2["total_array_size"]:,}')
    
    # compare structure
    logger.warning('=== Structure Comparison ===')
    structure_match, structure_diffs = compare_h5_structure(file1, file2)
    
    if structure_match:
        logger.warning('Structure: Files have identical structure')
    else:
        logger.error('Structure: Files have different structure')
        for diff in structure_diffs:
            logger.error(f'  {diff}')
    
    # compare content
    logger.warning('=== Content Comparison ==')
    content_match, content_diffs = compare_h5_content(file1, file2, tolerant, tolerance)
    
    if content_match:
        logger.warning('Content: Files have identical content')
    else:
        logger.error('Content: Files have different content')
        for diff in content_diffs:
            logger.error(f'  {diff}')
    
    # Summary
    logger.warning('=== Summary ===')
    if structure_match and content_match:
        logger.warning('Files are EQUIVALENT')
    else:
        logger.error('Files are NOT EQUIVALENT')
        if not structure_match:
            logger.error(f'Structure differences: {len(structure_diffs)}')
        if not content_match:
            logger.error(f'Content differences: {len(content_diffs)}')



def h5_inspect(h5_file: str, haplotypes: list[str] = None, show_matrix: bool = False, 
               max_loci: int = 10,  max_reads: int = 10, dense: bool = False) -> None:
    """
    Inspect a bam2emase HDF5 file with both debug information and matrix display.
    
    Args:
        h5_file: Path to the HDF5 file
        haplotypes: List of haplotype names to display (e.g., ['A', 'B']). If None, shows all.
        show_matrix: If True, show the matrix for each haplotype (default: False).
        max_loci: Maximum number of loci to display (default: 10)
        max_reads: Maximum number of reads to display (default: 10)
        dense: If True, convert sparse matrices to dense format for display (default: False)
    """
    logger.warning(f'Inspecting h5 file: {h5_file}')
    logger.debug(f'Haplotypes: {haplotypes}')
    logger.debug(f'Show matrix: {show_matrix}')
    logger.debug(f'Max loci: {max_loci}')
    logger.debug(f'Max reads: {max_reads}')
    logger.debug(f'Dense: {dense}')
    
    # Check if file exists
    if not os.path.exists(h5_file):
        logger.error(f'File does not exist: {h5_file}')
        return
    
    # Detect file type (compressed vs per-read)
    file_type = None
    matrix_row_label = 'rows'
    try:
        with tables.open_file(h5_file, 'r') as f:
            has_count = '/count' in f
            has_rname = '/rname' in f
            if has_count and not has_rname:
                file_type = 'Compressed (EC)'
                count_len = len(f.root.count)
                logger.warning(f'File type: {file_type}')
                logger.warning(f'ECs: {count_len}')
                matrix_row_label = 'ECs'
            elif has_rname and not has_count:
                file_type = 'Per-read alignment (bam2emase)'
                rname_len = len(f.root.rname)
                logger.warning(f'File type: {file_type}')
                logger.warning(f'Reads: {rname_len}')
                matrix_row_label = 'reads'
            elif has_count and has_rname:
                logger.warning('File type: WARNING: Both /count and /rname present (hybrid or unexpected file)')
                count_len = len(f.root.count)
                rname_len = len(f.root.rname)
                logger.warning(f'ECs: {count_len}')
                logger.warning(f'Reads: {rname_len}')
                matrix_row_label = 'rows'  # ambiguous
            else:
                logger.warning('File type: Unknown (neither /count nor /rname present)')
                matrix_row_label = 'rows'
    except Exception as e:
        logger.error(f'Error detecting file type: {e}')
        matrix_row_label = 'rows'
    
    # get file information
    logger.warning('=== File Information ===')
    info = get_file_info(h5_file)
    
    logger.warning(f'File: {info["filename"]}')
    logger.warning(f'  Size: {info["size"]:.2f} GB')
    if info['shape']:
        logger.warning(f'  Shape: {info["shape"]}')
    logger.warning(f'  Arrays: {len(info["arrays"])}')
    logger.warning(f'  Total array size: {info["total_array_size"]:,}')
    logger.warning('')
    
    try:
        with tables.open_file(h5_file, 'r') as f:
            # Get root attributes
            root = f.root
            incidence_only = getattr(root._v_attrs, 'incidence_only', False)
            mtype = getattr(root._v_attrs, 'mtype', 'unknown')
            shape = getattr(root._v_attrs, 'shape', None)
            hname = getattr(root._v_attrs, 'hname', [])
            
            logger.warning('=== Matrix Information ===')
            logger.warning(f'Matrix type: {mtype}')
            logger.warning(f'Shape: {shape}')
            logger.warning(f'Haplotypes: {hname}')
            logger.warning(f'Incidence only: {incidence_only}')
            logger.warning('')
            
            # Filter haplotypes if specified
            if haplotypes is not None:
                # Convert to uppercase for case-insensitive matching
                haplotypes_upper = [h.upper() for h in haplotypes]
                hname_upper = [h.upper() for h in hname]
                
                # Find matching haplotypes
                matching_indices = []
                for i, h in enumerate(hname_upper):
                    if h in haplotypes_upper:
                        matching_indices.append(i)
                
                if not matching_indices:
                    logger.error(f'No matching haplotypes found. Available: {hname}')
                    return
                
                logger.warning(f'Displaying haplotypes: {[hname[i] for i in matching_indices]}')
                logger.warning('')
            else:
                matching_indices = list(range(len(hname)))
            
            # Get locus and read names if available
            lname = None
            rname = None
            try:
                if '/lname' in f:
                    lname = f.root.lname.read()
                if '/rname' in f:
                    rname = f.root.rname.read()
            except Exception as e:
                logger.warning(f'Could not read locus/read names: {e}')
            
            # Display locus and read name information
            logger.warning('=== Dataset Information ===')
            if lname is not None:
                logger.warning(f'Locus names (lname):')
                logger.warning(f'  Shape: {lname.shape}')
                logger.warning(f'  Dtype: {lname.dtype}')
                logger.warning(f'  Elements: {len(lname):,}')
                logger.warning(f'  Memory: {len(lname) * lname.dtype.itemsize / (1024**2):.2f} MB')
                # Show sample locus names based on max_loci parameter
                if len(lname) > 0:
                    sample_loci = lname[:min(max_loci, len(lname))]
                    logger.warning(f'  Sample loci ({len(sample_loci)} of {len(lname):,}): {[x.decode() if isinstance(x, bytes) else x for x in sample_loci]}')
                logger.warning('')
            
            if rname is not None:
                logger.warning(f'Read names (rname):')
                logger.warning(f'  Shape: {rname.shape}')
                logger.warning(f'  Dtype: {rname.dtype}')
                logger.warning(f'  Elements: {len(rname):,}')
                logger.warning(f'  Memory: {len(rname) * rname.dtype.itemsize / (1024**2):.2f} MB')
                # Show sample read names based on max_reads parameter
                if len(rname) > 0:
                    sample_reads = rname[:min(max_reads, len(rname))]
                    logger.warning(f'  Sample reads ({len(sample_reads)} of {len(rname):,}): {[x.decode() if isinstance(x, bytes) else x for x in sample_reads]}')
                logger.warning('')
            
            # Process each selected haplotype
            for i in matching_indices:
                haplotype = hname[i]
                haplotype_path = f'/h{i}'
                
                if haplotype_path not in f:
                    logger.warning(f'Haplotype {haplotype} ({haplotype_path}) not found, skipping')
                    continue
                
                logger.warning(f'=== Haplotype {haplotype} (h{i}) ===')
                
                # Debug information for this haplotype
                try:
                    h_node = f.get_node(haplotype_path)
                    logger.warning(f'  Node type: {type(h_node).__name__}')
                    logger.warning(f'  Node title: {h_node._v_title}')
                    logger.warning(f'  Children: {len(h_node._v_children)}')
                    
                    # Show child nodes
                    for child_name, child_node in h_node._v_children.items():
                        if hasattr(child_node, 'shape'):
                            shape = child_node.shape
                            dtype = child_node.dtype if hasattr(child_node, 'dtype') else 'Unknown'
                            total_elements = np.prod(shape) if shape else 0
                            memory_usage = total_elements * child_node.dtype.itemsize if hasattr(child_node, 'dtype') else 0
                            
                            logger.warning(f'    {child_name}:')
                            logger.warning(f'      Shape: {shape}')
                            logger.warning(f'      Dtype: {dtype}')
                            logger.warning(f'      Elements: {total_elements:,}')
                            logger.warning(f'      Memory: {memory_usage / (1024**2):.2f} MB')
                    
                except Exception as e:
                    logger.error(f'  Error getting debug info for haplotype {haplotype}: {e}')

                if show_matrix:
                    
                    logger.warning('')
                    
                    # Matrix display for this haplotype
                    try:
                        # Read CSC matrix components
                        h_node = f.get_node(haplotype_path)
                        indptr = h_node.indptr.read()
                        indices = h_node.indices.read()
                        
                        # Check if data exists
                        has_data = False
                        data = None
                        try:
                            data = h_node.data.read()
                            has_data = True
                        except Exception:
                            pass
                        
                        if not has_data and not incidence_only:
                            logger.warning('  No data array found and not incidence-only mode')
                            continue
                        
                        # Convert to scipy sparse matrix
                        from scipy.sparse import csc_matrix
                        
                        # Get the expected shape from the file metadata
                        expected_shape = None
                        try:
                            shape_attr = f.root._v_attrs.shape
                            if isinstance(shape_attr, (tuple, list)) and len(shape_attr) == 3:
                                num_loci_file, num_haplotypes_file, num_reads_file = shape_attr
                                expected_shape = (num_reads_file, num_loci_file)
                                logger.debug(f'  Using expected shape from tuple metadata: {expected_shape}')
                            else:
                                logger.debug(f'  Shape attribute is not in expected format: {type(shape_attr)} = {shape_attr}')
                        except Exception as e:
                            logger.debug(f'  Could not get shape from file metadata: {e}')
                        
                        # Always use the expected shape if available, otherwise fall back to data inference
                        if expected_shape is None:
                            logger.debug(f'  Falling back to shape inference from data')
                        
                        if has_data:
                            if expected_shape:
                                sparse_matrix = csc_matrix((data, indices, indptr), shape=expected_shape)
                            else:
                                sparse_matrix = csc_matrix((data, indices, indptr))
                        else:
                            # For incidence-only, create matrix with ones
                            if expected_shape:
                                sparse_matrix = csc_matrix((np.ones_like(indices), indices, indptr), shape=expected_shape)
                            else:
                                sparse_matrix = csc_matrix((np.ones_like(indices), indices, indptr))
                        
                        # Get actual dimensions
                        num_reads, num_loci = sparse_matrix.shape
                        logger.warning(f'  Matrix shape: {num_reads} {matrix_row_label} × {num_loci} loci')
                        logger.warning(f'  Total non-zero elements: {sparse_matrix.nnz:,}')
                        
                        # Check if we should use dense format
                        max_dense_size = 20  # Maximum size for dense display
                        use_dense = dense and num_loci <= max_dense_size and num_reads <= max_dense_size
                        
                        if use_dense or dense:
                            # Always display the dense matrix for the specified region
                            dense_cap = 1000
                            display_reads = min(max_reads, num_reads, dense_cap)
                            display_loci = min(max_loci, num_loci, dense_cap)
                            if max_reads > dense_cap or max_loci > dense_cap:
                                logger.warning(f'  Requested dense display region ({max_reads}×{max_loci}) exceeds {dense_cap}×{dense_cap}. Showing only first {dense_cap}×{dense_cap}.')
                            dense_matrix = sparse_matrix.toarray()
                            if display_reads == 0 or display_loci == 0:
                                logger.warning(f'  No data to display in the specified region (reads: {display_reads}, loci: {display_loci})')
                            else:
                                logger.warning(f'  Dense matrix ({display_reads}×{display_loci} of {num_reads}×{num_loci}):')
                                for row in range(display_reads):
                                    row_str = '  '
                                    for col in range(display_loci):
                                        val = dense_matrix[row, col]
                                        if incidence_only or not has_data:
                                            row_str += '1 ' if val > 0 else '0 '
                                        else:
                                            row_str += f'{val} '
                                    logger.warning(row_str)
                                if display_reads < num_reads or display_loci < num_loci:
                                    logger.warning(f'  ... (showing {display_reads} of {num_reads} reads, {display_loci} of {num_loci} loci)')
                        else:
                            # Use sparse format (original logic)
                            if dense and (num_loci > max_dense_size or num_reads > max_dense_size):
                                logger.warning(f'  Matrix too large for dense display (>{max_dense_size}×{max_dense_size}), showing sparse format')
                            
                            # Find some actual data to display
                            coo_matrix = sparse_matrix.tocoo()
                            
                            if coo_matrix.nnz == 0:
                                logger.warning('  (completely empty matrix)')
                            else:
                                # Filter elements to only show those within max_reads and max_loci range
                                valid_mask = (coo_matrix.row < max_reads) & (coo_matrix.col < max_loci)
                                filtered_rows = coo_matrix.row[valid_mask]
                                filtered_cols = coo_matrix.col[valid_mask]
                                filtered_vals = coo_matrix.data[valid_mask]
                                
                                if len(filtered_rows) == 0:
                                    logger.warning(f'  No non-zero elements in first {max_reads} reads × {max_loci} loci region')
                                    logger.warning(f'  Total non-zero elements in full matrix: {coo_matrix.nnz:,}')
                                    
                                    # Show a sample from the full matrix for context
                                    sample_size = min(10, coo_matrix.nnz)
                                    if sample_size < coo_matrix.nnz:
                                        # Take a random sample
                                        import random
                                        sample_indices = random.sample(range(coo_matrix.nnz), sample_size)
                                        sample_rows = coo_matrix.row[sample_indices]
                                        sample_cols = coo_matrix.col[sample_indices]
                                        sample_vals = coo_matrix.data[sample_indices]
                                    else:
                                        sample_rows = coo_matrix.row
                                        sample_cols = coo_matrix.col
                                        sample_vals = coo_matrix.data
                                    
                                    # Sort by row, then column for consistent display
                                    sorted_indices = np.lexsort((sample_cols, sample_rows))
                                    rows = sample_rows[sorted_indices]
                                    cols = sample_cols[sorted_indices]
                                    vals = sample_vals[sorted_indices]
                                    
                                    logger.warning(f'  Sample of {len(rows)} non-zero elements from full matrix:')
                                    for idx in range(len(rows)):
                                        row, col, val = rows[idx], cols[idx], vals[idx]
                                        if incidence_only or not has_data:
                                            logger.warning(f'  [{row},{col}]: 1')
                                        else:
                                            logger.warning(f'  [{row},{col}]: {val}')
                                    
                                    if len(rows) < coo_matrix.nnz:
                                        logger.warning(f'  ... and {coo_matrix.nnz - len(rows):,} more elements')
                                else:
                                    # Sort by row, then column for consistent display
                                    sorted_indices = np.lexsort((filtered_cols, filtered_rows))
                                    rows = filtered_rows[sorted_indices]
                                    cols = filtered_cols[sorted_indices]
                                    vals = filtered_vals[sorted_indices]
                                    
                                    logger.warning(f'  Non-zero elements in first {max_reads} reads × {max_loci} loci:')
                                    
                                    # Show all elements in the specified range
                                    for idx in range(len(rows)):
                                        row, col, val = rows[idx], cols[idx], vals[idx]
                                        if incidence_only or not has_data:
                                            logger.warning(f'  [{row},{col}]: 1')
                                        else:
                                            logger.warning(f'  [{row},{col}]: {val}')
                                    
                                    # Show statistics about the data distribution
                                    logger.warning(f'  Found {len(rows)} elements in specified range')
                                    if coo_matrix.nnz > len(rows):
                                        logger.warning(f'  ... and {coo_matrix.nnz - len(rows):,} more elements outside this range')
                                
                                # Also show some statistics about the data distribution
                                logger.warning(f'  Full matrix index ranges: rows 0-{coo_matrix.shape[0]-1}, cols 0-{coo_matrix.shape[1]-1}')
                                
                                # Check if there's data in the first few rows/cols
                                early_data = (coo_matrix.row < 100) & (coo_matrix.col < 100)
                                if early_data.any():
                                    logger.warning(f'  Found {early_data.sum()} elements in first 100x100 region')
                                else:
                                    logger.warning('  No data in first 100x100 region (very sparse matrix)')
                    
                    except Exception as e:
                        logger.error(f'  Error processing matrix for haplotype {haplotype}: {e}')
                    
                    logger.warning('')
    
    except Exception as e:
        logger.error(f'Error opening file: {e}')
        return

