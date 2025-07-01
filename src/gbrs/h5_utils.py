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
            attrs1 = dict(node1._v_attrs)
        except Exception as e:
            # some nodes might have weird attributes, skip them for now
            pass
            
        try:
            attrs2 = dict(node2._v_attrs)
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




def debug_h5_file(h5_file: str) -> None:
    logger.warning(f'Comparing files:')
    logger.warning(f'Debugging h5 file: {h5_file}')
    
    
    # get file information
    logger.warning('=== File Information ===')
    info1 = get_file_info(h5_file)
    
    logger.warning(f'File: {info1["filename"]}')
    logger.warning(f'  Size: {info1["size"]:.2f} GB')
    if info1['shape']:
        logger.warning(f'  Shape: {info1["shape"]}')
    logger.warning(f'  Arrays: {len(info1["arrays"])}')
    logger.warning(f'  Total array size: {info1["total_array_size"]:,}')
    


