"""
Molecular Processor Module

This module provides functionality for batch processing of molecular data in the
Molecular Data Management and CRO Integration Platform. It includes functions for
validating, filtering, sorting, and enriching molecular data, as well as calculating
molecular properties.

The module is optimized for processing large datasets with parallel execution and
batch processing to ensure performance requirements are met.
"""

from typing import List, Dict, Optional, Union, Any, Tuple, Callable
import concurrent.futures
import logging
from tqdm import tqdm  # version 4.64+

# RDKit version 2023.03+
from rdkit import Chem

from ..exceptions import MolecularProcessingException, ValidationException
from ..constants import DEFAULT_MOLECULE_PROPERTIES, BATCH_SIZE
from .validator import (
    validate_smiles_structure, validate_molecule, validate_molecules_batch, PROPERTY_RANGES
)
from .molecule_converter import (
    smiles_to_mol, smiles_to_mol_with_exception
)
from .property_calculator import (
    calculate_properties, calculate_properties_from_smiles, 
    batch_calculate_properties, PROPERTY_CALCULATORS
)

# Set up module logger
logger = logging.getLogger(__name__)

# Default number of workers for parallel processing
DEFAULT_NUM_WORKERS = 4


def process_molecule(smiles: str, properties: List[str] = None) -> Dict[str, Any]:
    """
    Processes a single molecule with validation and property calculation.
    
    Args:
        smiles: SMILES string representation of the molecule
        properties: List of property names to calculate (default: None)
        
    Returns:
        Dict[str, Any]: Dictionary with SMILES, validation results, and calculated properties
    """
    # Initialize result with SMILES
    result = {'smiles': smiles}
    
    # Validate molecule structure and properties
    validation_result = validate_molecule(smiles)
    result['validation'] = validation_result
    
    # If structure is valid and properties are requested, calculate them
    if validation_result['structure_valid'] and properties:
        try:
            calculated_properties = calculate_properties_from_smiles(smiles, properties)
            result['properties'] = calculated_properties
        except MolecularProcessingException as e:
            logger.warning(f"Failed to calculate properties for molecule: {e.message}")
            result['properties'] = {}
            result['property_error'] = e.message
    
    return result


def process_molecules_batch(smiles_list: List[str], properties: List[str] = None, num_workers: int = None) -> List[Dict[str, Any]]:
    """
    Processes a batch of molecules in parallel.
    
    Args:
        smiles_list: List of SMILES strings
        properties: List of property names to calculate (default: None)
        num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
        
    Returns:
        List[Dict[str, Any]]: List of processed molecule results
    """
    if properties is None:
        properties = DEFAULT_MOLECULE_PROPERTIES
    
    if num_workers is None:
        num_workers = DEFAULT_NUM_WORKERS
    
    if not smiles_list:
        return []
    
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit tasks for each SMILES string
        future_to_index = {}
        for i, smiles in enumerate(smiles_list):
            future = executor.submit(process_molecule, smiles, properties)
            future_to_index[future] = i
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            try:
                result = future.result()
                result['index'] = index
                results.append(result)
            except Exception as e:
                logger.error(f"Error processing molecule at index {index}: {str(e)}")
                results.append({
                    'index': index,
                    'smiles': smiles_list[index],
                    'error': str(e),
                    'validation': {'structure_valid': False}
                })
    
    # Sort results by original index
    results.sort(key=lambda x: x['index'])
    return results


def process_molecules_stream(smiles_list: List[str], properties: List[str] = None, batch_size: int = None, num_workers: int = None) -> List[Dict[str, Any]]:
    """
    Processes large datasets of molecules efficiently in batches with progress tracking.
    
    Args:
        smiles_list: List of SMILES strings
        properties: List of property names to calculate (default: None)
        batch_size: Size of each processing batch (default: BATCH_SIZE from constants)
        num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
        
    Returns:
        List[Dict[str, Any]]: List of processed molecule results
    """
    if properties is None:
        properties = DEFAULT_MOLECULE_PROPERTIES
    
    if batch_size is None:
        batch_size = BATCH_SIZE
    
    if num_workers is None:
        num_workers = DEFAULT_NUM_WORKERS
    
    if not smiles_list:
        return []
    
    results = []
    total_molecules = len(smiles_list)
    
    # Split into batches
    batches = [smiles_list[i:i + batch_size] for i in range(0, total_molecules, batch_size)]
    
    logger.info(f"Processing {total_molecules} molecules in {len(batches)} batches")
    
    for i, batch in enumerate(tqdm(batches, desc="Processing molecules")):
        batch_results = process_molecules_batch(batch, properties, num_workers)
        results.extend(batch_results)
        logger.info(f"Processed batch {i+1}/{len(batches)} ({len(batch)} molecules)")
    
    logger.info(f"Completed processing {total_molecules} molecules")
    return results


def filter_molecules(molecules: List[Dict[str, Any]], property_filters: Dict[str, Dict[str, float]] = None, substructure: str = None) -> List[Dict[str, Any]]:
    """
    Filters molecules based on properties and/or substructure.
    
    Args:
        molecules: List of molecule dictionaries with properties
        property_filters: Dictionary mapping property names to min/max ranges
                          e.g., {'molecular_weight': {'min': 100, 'max': 500}}
        substructure: SMILES or SMARTS pattern to match as a substructure
        
    Returns:
        List[Dict[str, Any]]: Filtered list of molecules
    """
    if not molecules:
        return []
    
    filtered_molecules = molecules
    
    # Apply property filters if provided
    if property_filters:
        for property_name, value_range in property_filters.items():
            min_value = value_range.get('min', float('-inf'))
            max_value = value_range.get('max', float('inf'))
            
            filtered_molecules = [
                mol for mol in filtered_molecules
                if 'properties' in mol
                and property_name in mol['properties']
                and min_value <= mol['properties'][property_name] <= max_value
            ]
    
    # Apply substructure filter if provided
    if substructure:
        try:
            # Try to use RDKit directly for substructure matching
            pattern_mol = Chem.MolFromSmarts(substructure)
            if not pattern_mol:
                pattern_mol = Chem.MolFromSmiles(substructure)
            
            if pattern_mol:
                filtered_molecules = [
                    mol for mol in filtered_molecules
                    if 'smiles' in mol and smiles_to_mol(mol['smiles']) and 
                    smiles_to_mol(mol['smiles']).HasSubstructMatch(pattern_mol)
                ]
            else:
                logger.warning(f"Invalid substructure pattern: {substructure}")
        except Exception as e:
            logger.warning(f"Substructure filtering failed: {str(e)}")
    
    return filtered_molecules


def sort_molecules(molecules: List[Dict[str, Any]], property_name: str, ascending: bool = True) -> List[Dict[str, Any]]:
    """
    Sorts molecules based on a property value.
    
    Args:
        molecules: List of molecule dictionaries with properties
        property_name: Name of the property to sort by
        ascending: Sort in ascending order if True, descending if False
        
    Returns:
        List[Dict[str, Any]]: Sorted list of molecules
    """
    if not molecules:
        return []
    
    # Check if property exists in at least one molecule
    property_exists = any(
        'properties' in mol and property_name in mol['properties']
        for mol in molecules
    )
    
    if not property_exists:
        logger.warning(f"Property '{property_name}' not found in any molecules")
        return molecules
    
    # Define a key function that handles missing properties
    def sort_key(mol):
        if 'properties' in mol and property_name in mol['properties']:
            return mol['properties'][property_name]
        return float('-inf') if ascending else float('inf')
    
    return sorted(molecules, key=sort_key, reverse=not ascending)


def calculate_missing_properties(molecules: List[Dict[str, Any]], required_properties: List[str]) -> List[Dict[str, Any]]:
    """
    Calculates missing properties for molecules.
    
    Args:
        molecules: List of molecule dictionaries
        required_properties: List of property names that should be present
        
    Returns:
        List[Dict[str, Any]]: Molecules with calculated missing properties
    """
    if not molecules or not required_properties:
        return molecules
    
    # Identify molecules with valid structures that need property calculation
    molecules_needing_calculation = []
    indices_needing_calculation = []
    missing_properties_per_molecule = []
    
    for i, molecule in enumerate(molecules):
        if 'smiles' not in molecule:
            continue
        
        # Skip if molecule has invalid structure
        if 'validation' in molecule and not molecule['validation'].get('structure_valid', False):
            continue
        
        # Initialize properties dict if it doesn't exist
        if 'properties' not in molecule:
            molecule['properties'] = {}
        
        # Identify missing properties
        existing_properties = set(molecule['properties'].keys())
        missing_props = [prop for prop in required_properties if prop not in existing_properties]
        
        if missing_props:
            molecules_needing_calculation.append(molecule['smiles'])
            indices_needing_calculation.append(i)
            missing_properties_per_molecule.append(missing_props)
    
    # Calculate missing properties in batch if needed
    if molecules_needing_calculation:
        try:
            # Calculate all required properties for efficiency
            batch_results = batch_calculate_properties_from_smiles(
                molecules_needing_calculation, required_properties
            )
            
            # Update original molecules with calculated properties
            for i, props in enumerate(batch_results):
                orig_index = indices_needing_calculation[i]
                missing_props = missing_properties_per_molecule[i]
                
                # Only update the missing properties
                for prop in missing_props:
                    if prop in props:
                        molecules[orig_index]['properties'][prop] = props[prop]
        
        except Exception as e:
            logger.error(f"Batch property calculation failed: {str(e)}")
            # Fall back to individual calculations
            for i, smiles in enumerate(molecules_needing_calculation):
                try:
                    props = calculate_properties_from_smiles(smiles, missing_properties_per_molecule[i])
                    orig_index = indices_needing_calculation[i]
                    molecules[orig_index]['properties'].update(props)
                except MolecularProcessingException as e:
                    logger.warning(f"Failed to calculate properties for {smiles}: {e.message}")
    
    return molecules


def deduplicate_molecules(molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """
    Removes duplicate molecules based on SMILES strings.
    
    Args:
        molecules: List of molecule dictionaries
        
    Returns:
        List[Dict[str, Any]]: Deduplicated list of molecules
    """
    if not molecules:
        return []
    
    unique_molecules = {}
    
    for molecule in molecules:
        if 'smiles' not in molecule:
            continue
        
        smiles = molecule['smiles']
        
        # If SMILES not yet in unique molecules, add it
        if smiles not in unique_molecules:
            unique_molecules[smiles] = molecule
        else:
            # If duplicate found, keep the one with more properties or information
            existing_props = unique_molecules[smiles].get('properties', {})
            current_props = molecule.get('properties', {})
            
            if len(current_props) > len(existing_props):
                unique_molecules[smiles] = molecule
    
    return list(unique_molecules.values())


def enrich_molecules(molecules: List[Dict[str, Any]], additional_properties: List[str]) -> List[Dict[str, Any]]:
    """
    Enriches molecules with additional calculated properties.
    
    Args:
        molecules: List of molecule dictionaries
        additional_properties: List of additional property names to calculate
        
    Returns:
        List[Dict[str, Any]]: Enriched molecules with additional properties
    """
    if not molecules or not additional_properties:
        return molecules
    
    enriched_molecules = []
    
    for molecule in molecules:
        if 'smiles' not in molecule:
            enriched_molecules.append(molecule)
            continue
        
        # Skip if molecule has invalid structure
        if 'validation' in molecule and not molecule['validation'].get('structure_valid', False):
            enriched_molecules.append(molecule)
            continue
        
        # Initialize properties dict if it doesn't exist
        if 'properties' not in molecule:
            molecule['properties'] = {}
        
        # Calculate additional properties
        try:
            calculated_properties = calculate_properties_from_smiles(
                molecule['smiles'], additional_properties
            )
            molecule['properties'].update(calculated_properties)
        except MolecularProcessingException as e:
            logger.warning(f"Failed to enrich molecule {molecule['smiles']}: {e.message}")
        
        enriched_molecules.append(molecule)
    
    return enriched_molecules


def validate_and_standardize_molecules(molecules: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Validates and standardizes a list of molecules, separating valid and invalid molecules.
    
    Args:
        molecules: List of molecule dictionaries
        
    Returns:
        Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]: Tuple of valid and invalid molecules
    """
    if not molecules:
        return [], []
    
    valid_molecules = []
    invalid_molecules = []
    
    for molecule in molecules:
        if 'smiles' not in molecule:
            invalid_molecules.append({
                **molecule,
                'error': 'Missing SMILES string'
            })
            continue
        
        smiles = molecule['smiles']
        
        # Validate structure
        is_valid = validate_smiles_structure(smiles)
        
        if is_valid:
            # Attempt to standardize (canonicalize) the SMILES string
            try:
                from .molecule_converter import canonicalize_smiles
                canonical_smiles = canonicalize_smiles(smiles)
                
                # Update molecule with canonical SMILES and validation result
                validated_molecule = {
                    **molecule,
                    'smiles': canonical_smiles,
                    'validation': {'structure_valid': True}
                }
                
                valid_molecules.append(validated_molecule)
            except MolecularProcessingException as e:
                # If canonicalization fails, consider the molecule invalid
                invalid_molecules.append({
                    **molecule,
                    'validation': {'structure_valid': False},
                    'error': f"Standardization failed: {e.message}"
                })
        else:
            invalid_molecules.append({
                **molecule,
                'validation': {'structure_valid': False},
                'error': 'Invalid SMILES structure'
            })
    
    return valid_molecules, invalid_molecules


def generate_processing_summary(processed_molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Generates a summary of molecule processing results.
    
    Args:
        processed_molecules: List of processed molecule dictionaries
        
    Returns:
        Dict[str, Any]: Summary statistics of processing results
    """
    if not processed_molecules:
        return {
            'total_molecules': 0,
            'valid_molecules': 0,
            'invalid_molecules': 0,
            'success_rate': 0,
            'errors': {}
        }
    
    total_molecules = len(processed_molecules)
    
    # Count valid and invalid molecules
    valid_molecules = sum(
        1 for mol in processed_molecules
        if 'validation' in mol and mol['validation'].get('structure_valid', False)
    )
    
    invalid_molecules = total_molecules - valid_molecules
    
    # Calculate success rate
    success_rate = (valid_molecules / total_molecules) * 100 if total_molecules > 0 else 0
    
    # Collect error information
    errors = {}
    for mol in processed_molecules:
        if 'error' in mol:
            error_message = mol['error']
            errors[error_message] = errors.get(error_message, 0) + 1
    
    # Calculate property statistics for valid molecules
    property_stats = {}
    if valid_molecules > 0:
        # Find all available properties
        all_properties = set()
        for mol in processed_molecules:
            if 'properties' in mol:
                all_properties.update(mol['properties'].keys())
        
        # Calculate statistics for each property
        for prop in all_properties:
            values = [mol['properties'][prop] for mol in processed_molecules
                     if 'properties' in mol and prop in mol['properties']]
            
            if values:
                property_stats[prop] = {
                    'min': min(values),
                    'max': max(values),
                    'average': sum(values) / len(values),
                    'count': len(values)
                }
    
    return {
        'total_molecules': total_molecules,
        'valid_molecules': valid_molecules,
        'invalid_molecules': invalid_molecules,
        'success_rate': success_rate,
        'errors': errors,
        'property_statistics': property_stats
    }


class MoleculeProcessor:
    """
    Class for processing and managing molecular data with caching capability.
    """
    
    def __init__(self, num_workers: int = None, batch_size: int = None):
        """
        Initializes the MoleculeProcessor with configuration parameters.
        
        Args:
            num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
            batch_size: Size of each processing batch (default: BATCH_SIZE from constants)
        """
        self._num_workers = num_workers if num_workers is not None else DEFAULT_NUM_WORKERS
        self._batch_size = batch_size if batch_size is not None else BATCH_SIZE
        self._cache = {}  # Cache for processed results
        logger.info(f"MoleculeProcessor initialized with {self._num_workers} workers and batch size {self._batch_size}")
    
    def process(self, smiles: str, properties: List[str] = None) -> Dict[str, Any]:
        """
        Processes a single molecule with caching.
        
        Args:
            smiles: SMILES string representation of the molecule
            properties: List of property names to calculate (default: None)
            
        Returns:
            Dict[str, Any]: Processed molecule data
        """
        # Check if result is already cached
        cache_key = (smiles, tuple(properties) if properties else None)
        if cache_key in self._cache:
            return self._cache[cache_key]
        
        # Process molecule
        result = process_molecule(smiles, properties)
        
        # Cache result
        self._cache[cache_key] = result
        
        return result
    
    def process_batch(self, smiles_list: List[str], properties: List[str] = None) -> List[Dict[str, Any]]:
        """
        Processes a batch of molecules with caching.
        
        Args:
            smiles_list: List of SMILES strings
            properties: List of property names to calculate (default: None)
            
        Returns:
            List[Dict[str, Any]]: List of processed molecule data
        """
        # Find which molecules are already cached
        uncached_smiles = []
        uncached_indices = []
        results = [None] * len(smiles_list)  # Pre-allocate results list
        
        for i, smiles in enumerate(smiles_list):
            cache_key = (smiles, tuple(properties) if properties else None)
            if cache_key in self._cache:
                results[i] = self._cache[cache_key]
            else:
                uncached_smiles.append(smiles)
                uncached_indices.append(i)
        
        # Process uncached molecules in batch
        if uncached_smiles:
            batch_results = process_molecules_batch(uncached_smiles, properties, self._num_workers)
            
            # Update cache and results
            for i, result in enumerate(batch_results):
                orig_idx = uncached_indices[i]
                smiles = uncached_smiles[i]
                cache_key = (smiles, tuple(properties) if properties else None)
                self._cache[cache_key] = result
                results[orig_idx] = result
        
        return results
    
    def process_stream(self, smiles_list: List[str], properties: List[str] = None) -> List[Dict[str, Any]]:
        """
        Processes a large stream of molecules in batches with caching.
        
        Args:
            smiles_list: List of SMILES strings
            properties: List of property names to calculate (default: None)
            
        Returns:
            List[Dict[str, Any]]: List of processed molecule data
        """
        # Process in batches
        batch_size = self._batch_size
        num_batches = (len(smiles_list) + batch_size - 1) // batch_size
        
        all_results = []
        for i in range(num_batches):
            start_idx = i * batch_size
            end_idx = min((i + 1) * batch_size, len(smiles_list))
            batch = smiles_list[start_idx:end_idx]
            
            logger.info(f"Processing batch {i+1}/{num_batches} ({len(batch)} molecules)")
            batch_results = self.process_batch(batch, properties)
            all_results.extend(batch_results)
        
        return all_results
    
    def filter(self, molecules: List[Dict[str, Any]], property_filters: Dict[str, Dict[str, float]] = None, substructure: str = None) -> List[Dict[str, Any]]:
        """
        Filters molecules based on properties and/or substructure.
        
        Args:
            molecules: List of molecule dictionaries with properties
            property_filters: Dictionary mapping property names to min/max ranges
            substructure: SMILES or SMARTS pattern to match as a substructure
            
        Returns:
            List[Dict[str, Any]]: Filtered list of molecules
        """
        return filter_molecules(molecules, property_filters, substructure)
    
    def sort(self, molecules: List[Dict[str, Any]], property_name: str, ascending: bool = True) -> List[Dict[str, Any]]:
        """
        Sorts molecules based on a property value.
        
        Args:
            molecules: List of molecule dictionaries with properties
            property_name: Name of the property to sort by
            ascending: Sort in ascending order if True, descending if False
            
        Returns:
            List[Dict[str, Any]]: Sorted list of molecules
        """
        return sort_molecules(molecules, property_name, ascending)
    
    def deduplicate(self, molecules: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Removes duplicate molecules based on SMILES strings.
        
        Args:
            molecules: List of molecule dictionaries
            
        Returns:
            List[Dict[str, Any]]: Deduplicated list of molecules
        """
        return deduplicate_molecules(molecules)
    
    def enrich(self, molecules: List[Dict[str, Any]], additional_properties: List[str]) -> List[Dict[str, Any]]:
        """
        Enriches molecules with additional calculated properties.
        
        Args:
            molecules: List of molecule dictionaries
            additional_properties: List of additional property names to calculate
            
        Returns:
            List[Dict[str, Any]]: Enriched molecules with additional properties
        """
        return enrich_molecules(molecules, additional_properties)
    
    def validate_and_standardize(self, molecules: List[Dict[str, Any]]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
        """
        Validates and standardizes a list of molecules, separating valid and invalid molecules.
        
        Args:
            molecules: List of molecule dictionaries
            
        Returns:
            Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]: Tuple of valid and invalid molecules
        """
        return validate_and_standardize_molecules(molecules)
    
    def generate_summary(self, processed_molecules: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Generates a summary of molecule processing results.
        
        Args:
            processed_molecules: List of processed molecule dictionaries
            
        Returns:
            Dict[str, Any]: Summary statistics of processing results
        """
        return generate_processing_summary(processed_molecules)
    
    def clear_cache(self) -> None:
        """
        Clears the processor's result cache.
        """
        cache_size = len(self._cache)
        self._cache = {}
        logger.info(f"Cleared cache containing {cache_size} entries")