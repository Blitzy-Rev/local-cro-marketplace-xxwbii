"""
Utility module providing helper functions for molecular data processing in the
Molecular Data Management and CRO Integration Platform.

This module contains functions for molecular structure manipulation, property calculation,
similarity comparison, and drug-likeness assessment. These functions are used throughout
the application to support processing and validation of molecular data during CSV import,
filtering, and organization of molecules.
"""

import logging
from typing import List, Dict, Optional, Union, Any, Tuple

# RDKit version 2023.03+
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import DataStructs

from ..exceptions import MolecularProcessingException, ValidationException
from ..constants import DEFAULT_MOLECULE_PROPERTIES
from ..molecular.molecule_converter import (
    smiles_to_mol,
    smiles_to_mol_with_exception,
    mol_to_smiles
)
from ..molecular.property_calculator import calculate_properties

# Set up module logger
logger = logging.getLogger(__name__)

# Define valid ranges for molecular properties
MOLECULAR_WEIGHT_RANGE = {"min": 0, "max": 2000}
LOGP_RANGE = {"min": -10, "max": 10}

# Lipinski's Rule of Five thresholds
LIPINSKI_RULES = {
    "molecular_weight": 500,
    "logp": 5,
    "h_bond_donors": 5,
    "h_bond_acceptors": 10
}

# Veber rules thresholds
VEBER_RULES = {
    "rotatable_bonds": 10,
    "polar_surface_area": 140
}


def is_valid_smiles(smiles: str) -> bool:
    """
    Checks if a SMILES string is valid by attempting to create an RDKit molecule.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        True if SMILES is valid, False otherwise
    """
    if not smiles or not isinstance(smiles, str):
        return False
    
    mol = smiles_to_mol(smiles)
    return mol is not None


def validate_smiles(smiles: str) -> bool:
    """
    Validates a SMILES string and raises an exception if invalid.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        True if SMILES is valid
        
    Raises:
        ValidationException: If the SMILES string is invalid or empty
    """
    if not smiles or not isinstance(smiles, str):
        raise ValidationException(
            "SMILES string cannot be empty or must be a string",
            details={"smiles": smiles}
        )
    
    try:
        smiles_to_mol_with_exception(smiles)
        return True
    except MolecularProcessingException as e:
        raise ValidationException(
            f"Invalid SMILES string: {e.message}",
            details=e.details
        )


def canonicalize_smiles(smiles: str) -> str:
    """
    Converts a SMILES string to its canonical form.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Canonical SMILES string
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_smiles(mol, canonical=True)


def calculate_molecular_properties(smiles: str, properties: List[str] = None) -> Dict[str, float]:
    """
    Calculates molecular properties for a given SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        properties: List of property names to calculate (uses DEFAULT_MOLECULE_PROPERTIES if None)
        
    Returns:
        Dictionary of property names and their calculated values
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    
    if properties is None:
        properties = DEFAULT_MOLECULE_PROPERTIES
        
    return calculate_properties(mol, properties)


def is_property_in_range(property_name: str, value: float, ranges: Dict[str, Dict[str, float]] = None) -> bool:
    """
    Checks if a molecular property value is within a specified range.
    
    Args:
        property_name: Name of the property
        value: Value to check
        ranges: Dictionary of property ranges (uses predefined ranges if None)
        
    Returns:
        True if property value is within range, False otherwise
    """
    if ranges is None:
        # Use predefined ranges for common properties
        if property_name == "molecular_weight":
            min_val, max_val = MOLECULAR_WEIGHT_RANGE["min"], MOLECULAR_WEIGHT_RANGE["max"]
        elif property_name == "logp":
            min_val, max_val = LOGP_RANGE["min"], LOGP_RANGE["max"]
        else:
            # No validation for other properties
            return True
    else:
        # Use provided ranges
        if property_name not in ranges:
            # No validation for properties not in ranges
            return True
        min_val, max_val = ranges[property_name]["min"], ranges[property_name]["max"]
    
    return min_val <= value <= max_val


def filter_molecules_by_property(molecules: List[Dict[str, Any]], property_name: str, 
                                min_value: float = None, max_value: float = None) -> List[Dict[str, Any]]:
    """
    Filters a list of molecules based on property values.
    
    Args:
        molecules: List of molecule dictionaries, each containing properties
        property_name: Name of the property to filter by
        min_value: Minimum property value (inclusive)
        max_value: Maximum property value (inclusive)
        
    Returns:
        Filtered list of molecules
    """
    result = []
    
    for molecule in molecules:
        # Check if the molecule has the property
        if "properties" in molecule and property_name in molecule["properties"]:
            value = molecule["properties"][property_name]
            
            # Check if value is within range
            in_range = True
            if min_value is not None and value < min_value:
                in_range = False
            if max_value is not None and value > max_value:
                in_range = False
                
            if in_range:
                result.append(molecule)
    
    return result


def calculate_similarity(smiles1: str, smiles2: str, fingerprint_type: str = 'morgan') -> float:
    """
    Calculates the similarity between two molecules specified by SMILES strings.
    
    Args:
        smiles1: SMILES string of the first molecule
        smiles2: SMILES string of the second molecule
        fingerprint_type: Type of fingerprint to use ('morgan', 'maccs', 'topological', 'pattern')
        
    Returns:
        Similarity score between 0 and 1
        
    Raises:
        ValidationException: If either SMILES string is invalid
    """
    # Validate both SMILES strings
    validate_smiles(smiles1)
    validate_smiles(smiles2)
    
    # Convert SMILES to molecules
    mol1 = smiles_to_mol_with_exception(smiles1)
    mol2 = smiles_to_mol_with_exception(smiles2)
    
    # Generate fingerprints based on the specified type
    if fingerprint_type == 'morgan':
        fp1 = Chem.AllChem.GetMorganFingerprintAsBitVect(mol1, 2, 2048)
        fp2 = Chem.AllChem.GetMorganFingerprintAsBitVect(mol2, 2, 2048)
    elif fingerprint_type == 'maccs':
        fp1 = Chem.MACCSkeys.GenMACCSKeys(mol1)
        fp2 = Chem.MACCSkeys.GenMACCSKeys(mol2)
    elif fingerprint_type == 'topological':
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
    elif fingerprint_type == 'pattern':
        fp1 = Chem.PatternFingerprint(mol1)
        fp2 = Chem.PatternFingerprint(mol2)
    else:
        raise ValidationException(
            f"Invalid fingerprint type: {fingerprint_type}",
            details={"valid_types": ['morgan', 'maccs', 'topological', 'pattern']}
        )
    
    # Calculate Tanimoto similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


def check_lipinski_rule_of_five(smiles: str) -> Dict[str, Union[bool, int]]:
    """
    Checks if a molecule follows Lipinski's Rule of Five for drug-likeness.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Dictionary with rule compliance results and violation count
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    
    # Calculate relevant properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    
    # Check each rule
    mw_pass = mw <= LIPINSKI_RULES["molecular_weight"]
    logp_pass = logp <= LIPINSKI_RULES["logp"]
    donors_pass = h_donors <= LIPINSKI_RULES["h_bond_donors"]
    acceptors_pass = h_acceptors <= LIPINSKI_RULES["h_bond_acceptors"]
    
    # Count violations
    violations = sum(not x for x in [mw_pass, logp_pass, donors_pass, acceptors_pass])
    
    return {
        "molecular_weight": mw,
        "logp": logp,
        "h_bond_donors": h_donors,
        "h_bond_acceptors": h_acceptors,
        "mw_pass": mw_pass,
        "logp_pass": logp_pass,
        "donors_pass": donors_pass,
        "acceptors_pass": acceptors_pass,
        "violations": violations,
        "is_drug_like": violations <= 1  # Up to 1 violation is often still considered drug-like
    }


def check_veber_rules(smiles: str) -> Dict[str, Union[bool, int]]:
    """
    Checks if a molecule follows Veber's rules for oral bioavailability.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Dictionary with rule compliance results and violation count
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    
    # Calculate relevant properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    polar_surface_area = Descriptors.TPSA(mol)
    
    # Check each rule
    rotatable_pass = rotatable_bonds <= VEBER_RULES["rotatable_bonds"]
    psa_pass = polar_surface_area <= VEBER_RULES["polar_surface_area"]
    
    # Count violations
    violations = sum(not x for x in [rotatable_pass, psa_pass])
    
    return {
        "rotatable_bonds": rotatable_bonds,
        "polar_surface_area": polar_surface_area,
        "rotatable_pass": rotatable_pass,
        "psa_pass": psa_pass,
        "violations": violations,
        "passes_veber_rules": violations == 0  # Both rules must be satisfied
    }


def check_drug_likeness(smiles: str) -> Dict[str, Any]:
    """
    Checks if a molecule follows both Lipinski's Rule of Five and Veber's rules.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Dictionary with combined drug-likeness assessment
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    
    # Get Lipinski and Veber results
    lipinski_results = check_lipinski_rule_of_five(smiles)
    veber_results = check_veber_rules(smiles)
    
    # Combine results
    combined_results = {
        "lipinski": lipinski_results,
        "veber": veber_results,
        "total_violations": lipinski_results["violations"] + veber_results["violations"],
        "is_drug_like": lipinski_results["is_drug_like"] and veber_results["passes_veber_rules"]
    }
    
    return combined_results


def has_substructure(smiles: str, smarts_pattern: str) -> bool:
    """
    Checks if a molecule contains a specific substructure defined by SMARTS pattern.
    
    Args:
        smiles: SMILES string representation of a molecule
        smarts_pattern: SMARTS pattern defining the substructure
        
    Returns:
        True if substructure is found, False otherwise
        
    Raises:
        ValidationException: If the SMILES string is invalid
        MolecularProcessingException: If the SMARTS pattern is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    
    # Parse SMARTS pattern
    substructure = Chem.MolFromSmarts(smarts_pattern)
    if substructure is None:
        raise MolecularProcessingException(
            "Invalid SMARTS pattern",
            details={"smarts_pattern": smarts_pattern}
        )
    
    # Check if molecule contains the substructure
    return mol.HasSubstructMatch(substructure)


def get_molecular_formula(smiles: str) -> str:
    """
    Calculates the molecular formula for a given SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Molecular formula (e.g., C6H12O6)
        
    Raises:
        ValidationException: If the SMILES string is invalid
    """
    validate_smiles(smiles)
    mol = smiles_to_mol_with_exception(smiles)
    
    # Get atom dictionary
    atom_dict = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol in atom_dict:
            atom_dict[symbol] += 1
        else:
            atom_dict[symbol] = 1
    
    # Order atoms: C, H, then alphabetically
    ordered_atoms = []
    if 'C' in atom_dict:
        ordered_atoms.append(('C', atom_dict['C']))
        del atom_dict['C']
    
    if 'H' in atom_dict:
        ordered_atoms.append(('H', atom_dict['H']))
        del atom_dict['H']
    
    # Add remaining atoms in alphabetical order
    for symbol in sorted(atom_dict.keys()):
        ordered_atoms.append((symbol, atom_dict[symbol]))
    
    # Build formula string
    formula = ''
    for symbol, count in ordered_atoms:
        formula += symbol
        if count > 1:
            formula += str(count)
    
    return formula