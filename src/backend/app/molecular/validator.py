"""
Molecular Validator Module

This module provides utilities for validating molecular structures and properties 
in the Molecular Data Management and CRO Integration Platform. It includes functions 
for validating SMILES strings, checking property values against expected ranges, 
and batch validation of multiple molecules.

It also provides functions for checking drug-likeness and bioavailability using 
established rules such as Lipinski's Rule of Five and Veber's rules.
"""

from typing import List, Dict, Optional, Union, Any, Tuple
import concurrent.futures
import logging

# RDKit version 2023.03+
from rdkit import Chem
from rdkit.Chem import Descriptors

from ..exceptions import MolecularProcessingException, ValidationException
from ..constants import DEFAULT_MOLECULE_PROPERTIES
from .molecule_converter import smiles_to_mol, smiles_to_mol_with_exception

# Set up module logger
logger = logging.getLogger(__name__)

# Default number of workers for parallel processing
DEFAULT_NUM_WORKERS = 4

# Define valid ranges for molecular properties
PROPERTY_RANGES = {
    'molecular_weight': {'min': 0, 'max': 1000},
    'logp': {'min': -10, 'max': 10},
    'h_bond_donors': {'min': 0, 'max': 10},
    'h_bond_acceptors': {'min': 0, 'max': 20},
    'rotatable_bonds': {'min': 0, 'max': 15},
    'polar_surface_area': {'min': 0, 'max': 200},
    'heavy_atom_count': {'min': 0, 'max': 100},
    'ring_count': {'min': 0, 'max': 10},
    'aromatic_rings': {'min': 0, 'max': 7},
    'solubility': {'min': -10, 'max': 10}
}


def validate_smiles_structure(smiles: str) -> bool:
    """
    Validates the structure of a SMILES string by attempting to convert it to an RDKit molecule.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        bool: True if the SMILES string is valid, False otherwise
    """
    if not smiles:
        return False
    
    mol = smiles_to_mol(smiles)
    return mol is not None


def validate_smiles_structure_with_exception(smiles: str) -> bool:
    """
    Validates the structure of a SMILES string and raises an exception if invalid.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        bool: True if the SMILES string is valid
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    if not smiles:
        raise ValidationException(
            "SMILES string cannot be empty", 
            details={"smiles": smiles}
        )
    
    try:
        mol = smiles_to_mol_with_exception(smiles)
        return True
    except MolecularProcessingException as e:
        raise ValidationException(
            f"Invalid SMILES structure: {e.message}",
            details=e.details
        )


def validate_property_value(property_name: str, value: float) -> bool:
    """
    Validates if a property value is within the expected range.
    
    Args:
        property_name: Name of the property to validate
        value: Value to validate
        
    Returns:
        bool: True if value is within range, False otherwise
    """
    if property_name not in PROPERTY_RANGES:
        # If property is not in our range dictionary, assume it's valid
        return True
    
    prop_range = PROPERTY_RANGES[property_name]
    return prop_range['min'] <= value <= prop_range['max']


def validate_molecule_properties(properties: Dict[str, float]) -> Dict[str, bool]:
    """
    Validates multiple properties of a molecule against expected ranges.
    
    Args:
        properties: Dictionary of property names and their values
        
    Returns:
        Dict[str, bool]: Dictionary of property names and their validation results
    """
    results = {}
    
    for prop_name, value in properties.items():
        results[prop_name] = validate_property_value(prop_name, value)
    
    return results


def validate_molecule(smiles: str, properties: Dict[str, float] = None) -> Dict[str, Any]:
    """
    Validates both the structure and properties of a molecule.
    
    Args:
        smiles: SMILES string to validate
        properties: Dictionary of property values to validate (optional)
        
    Returns:
        Dict[str, Any]: Dictionary with structure validation result and property validation results
    """
    # Initialize results with structure validation
    results = {
        'structure_valid': validate_smiles_structure(smiles)
    }
    
    # If structure is valid and properties are provided, validate properties
    if results['structure_valid'] and properties:
        property_results = validate_molecule_properties(properties)
        results['property_validation'] = property_results
        results['all_properties_valid'] = all(property_results.values())
    
    return results


def validate_molecules_batch(molecules: List[Dict[str, Any]], num_workers: int = None) -> List[Dict[str, Any]]:
    """
    Validates a batch of molecules in parallel.
    
    Args:
        molecules: List of dictionaries, each containing at least a 'smiles' key and optionally a 'properties' key
        num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
        
    Returns:
        List[Dict[str, Any]]: List of validation results for each molecule
    """
    if num_workers is None:
        num_workers = DEFAULT_NUM_WORKERS
    
    if not molecules:
        return []
    
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit validation tasks
        future_to_index = {}
        for i, molecule in enumerate(molecules):
            smiles = molecule.get('smiles', '')
            properties = molecule.get('properties', {})
            future = executor.submit(validate_molecule, smiles, properties)
            future_to_index[future] = i
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            try:
                validation_result = future.result()
                # Add original molecule data along with validation results
                result = {
                    'index': index,
                    'smiles': molecules[index].get('smiles', ''),
                    'validation': validation_result
                }
                if 'properties' in molecules[index]:
                    result['properties'] = molecules[index]['properties']
                results.append(result)
            except Exception as e:
                logger.error(f"Error validating molecule at index {index}: {str(e)}")
                results.append({
                    'index': index,
                    'smiles': molecules[index].get('smiles', ''),
                    'validation': {
                        'structure_valid': False,
                        'error': str(e)
                    }
                })
    
    # Sort results by original index
    results.sort(key=lambda x: x['index'])
    return results


def check_lipinski_rule_of_five(mol: Chem.rdchem.Mol) -> Dict[str, Any]:
    """
    Checks if a molecule follows Lipinski's Rule of Five for drug-likeness.
    
    Lipinski's rules:
    - Molecular weight < 500 Da
    - LogP < 5
    - H-bond donors < 5
    - H-bond acceptors < 10
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dict[str, Any]: Dictionary with Lipinski properties and rule compliance
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot check Lipinski rules with None molecule",
            details={"mol": "None"}
        )
    
    # Calculate properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    
    # Check rules
    mw_pass = mw <= 500
    logp_pass = logp <= 5
    donors_pass = h_donors <= 5
    acceptors_pass = h_acceptors <= 10
    
    # Count violations
    violations = sum(not x for x in [mw_pass, logp_pass, donors_pass, acceptors_pass])
    
    # Generally, up to 1 violation is still considered drug-like
    is_drug_like = violations <= 1
    
    return {
        'molecular_weight': mw,
        'logp': logp,
        'h_bond_donors': h_donors,
        'h_bond_acceptors': h_acceptors,
        'mw_pass': mw_pass,
        'logp_pass': logp_pass,
        'donors_pass': donors_pass,
        'acceptors_pass': acceptors_pass,
        'violations': violations,
        'is_drug_like': is_drug_like
    }


def check_veber_rules(mol: Chem.rdchem.Mol) -> Dict[str, Any]:
    """
    Checks if a molecule follows Veber's rules for oral bioavailability.
    
    Veber's rules:
    - Rotatable bonds ≤ 10
    - Polar surface area ≤ 140 Å²
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dict[str, Any]: Dictionary with Veber properties and rule compliance
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot check Veber rules with None molecule",
            details={"mol": "None"}
        )
    
    # Calculate properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    polar_surface_area = Descriptors.TPSA(mol)
    
    # Check rules
    rotatable_pass = rotatable_bonds <= 10
    psa_pass = polar_surface_area <= 140
    
    # Both rules must be satisfied
    passes_veber_rules = rotatable_pass and psa_pass
    
    return {
        'rotatable_bonds': rotatable_bonds,
        'polar_surface_area': polar_surface_area,
        'rotatable_pass': rotatable_pass,
        'psa_pass': psa_pass,
        'passes_veber_rules': passes_veber_rules
    }


def get_property_range(property_name: str) -> Dict[str, float]:
    """
    Gets the expected range for a molecular property.
    
    Args:
        property_name: Name of the property
        
    Returns:
        Dict[str, float]: Dictionary with min and max values for the property
    """
    if property_name in PROPERTY_RANGES:
        return PROPERTY_RANGES[property_name]
    else:
        # Return a default wide range if property not found
        return {'min': float('-inf'), 'max': float('inf')}