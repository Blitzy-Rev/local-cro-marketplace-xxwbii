"""
Unit tests for molecular_utils module.

This module contains tests for the molecular utility functions used in the 
Molecular Data Management and CRO Integration Platform. It tests SMILES validation,
property calculation, filtering, similarity comparison, and drug-likeness assessment.
"""

import pytest
from typing import List, Dict
from rdkit import Chem  # version 2023.03+
from rdkit.Chem import Descriptors  # version 2023.03+

from ../../app/utils/molecular_utils import is_valid_smiles
from ../../app/utils/molecular_utils import validate_smiles
from ../../app/utils/molecular_utils import canonicalize_smiles
from ../../app/utils/molecular_utils import calculate_molecular_properties
from ../../app/utils/molecular_utils import is_property_in_range
from ../../app/utils/molecular_utils import filter_molecules_by_property
from ../../app/utils/molecular_utils import calculate_similarity
from ../../app/utils/molecular_utils import check_lipinski_rule_of_five
from ../../app/utils/molecular_utils import check_veber_rules
from ../../app/utils/molecular_utils import check_drug_likeness
from ../../app/utils/molecular_utils import has_substructure
from ../../app/utils/molecular_utils import get_molecular_formula
from ../../app/utils/molecular_utils import MOLECULAR_WEIGHT_RANGE
from ../../app/utils/molecular_utils import LOGP_RANGE
from ../../app/utils/molecular_utils import LIPINSKI_RULES
from ../../app/utils/molecular_utils import VEBER_RULES
from ../../app/exceptions import ValidationException
from ../../app/exceptions import MolecularProcessingException

# Test data - valid SMILES strings
VALID_SMILES = [
    'CCO',              # Ethanol 
    'c1ccccc1',         # Benzene
    'CC(=O)O',          # Acetic acid
    'C1CCCCC1',         # Cyclohexane
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'  # Caffeine
]

# Test data - invalid SMILES strings
INVALID_SMILES = [
    '',                # Empty string
    'X',               # Invalid atom
    'C1CC',            # Unclosed ring
    'c1ccccc',         # Unclosed aromatic ring
    'malformedsmiles'  # Not a SMILES string
]

# Test data - molecular properties
TEST_PROPERTIES = {
    "molecular_weight": 46.07,
    "logp": -0.14,
    "h_bond_donors": 1,
    "h_bond_acceptors": 1,
    "rotatable_bonds": 1,
    "polar_surface_area": 20.23
}

# Test data - molecules with properties
TEST_MOLECULES = [
    {"smiles": "CCO", "properties": {"molecular_weight": 46.07, "logp": -0.14}},
    {"smiles": "c1ccccc1", "properties": {"molecular_weight": 78.11, "logp": 1.90}},
    {"smiles": "CC(=O)O", "properties": {"molecular_weight": 60.05, "logp": -0.17}}
]

# Test data - property ranges
TEST_RANGES = {
    "molecular_weight": {"min": 0, "max": 500},
    "logp": {"min": -3, "max": 5}
}

# Test data - SMARTS patterns for substructure searching
TEST_SMARTS_PATTERNS = {
    "alcohol": "[OH]C",
    "carboxylic_acid": "C(=O)[OH]",
    "benzene": "c1ccccc1"
}

def test_is_valid_smiles_valid():
    """Tests that valid SMILES strings are correctly identified."""
    for smiles in VALID_SMILES:
        assert is_valid_smiles(smiles) is True

def test_is_valid_smiles_invalid():
    """Tests that invalid SMILES strings are correctly identified."""
    for smiles in INVALID_SMILES:
        assert is_valid_smiles(smiles) is False

def test_validate_smiles_valid():
    """Tests that valid SMILES strings are correctly validated with exception handling."""
    for smiles in VALID_SMILES:
        assert validate_smiles(smiles) is True

def test_validate_smiles_invalid():
    """Tests that invalid SMILES strings raise ValidationException."""
    for smiles in INVALID_SMILES:
        with pytest.raises(ValidationException):
            validate_smiles(smiles)

def test_canonicalize_smiles():
    """Tests conversion of SMILES strings to canonical form."""
    test_cases = [
        ('C=C', 'C=C'),          # Ethylene
        ('C-C-O', 'CCO'),        # Ethanol with explicit bonds
        ('C(C)O', 'CCO'),        # Ethanol with branch
        ('OCC', 'CCO'),          # Ethanol from different atom ordering
        ('c1ccccc1', 'c1ccccc1') # Benzene
    ]
    
    for input_smiles, expected_canonical in test_cases:
        assert canonicalize_smiles(input_smiles) == expected_canonical

def test_canonicalize_smiles_invalid():
    """Tests that invalid SMILES strings raise ValidationException during canonicalization."""
    for smiles in INVALID_SMILES:
        with pytest.raises(ValidationException):
            canonicalize_smiles(smiles)

def test_calculate_molecular_properties():
    """Tests calculation of molecular properties from SMILES."""
    # Test with a simple molecule (ethanol)
    properties = calculate_molecular_properties('CCO')
    
    # Check that result is a dictionary
    assert isinstance(properties, dict)
    
    # Check that common properties are included
    assert 'molecular_weight' in properties
    assert 'logp' in properties
    assert 'h_bond_donors' in properties
    assert 'h_bond_acceptors' in properties
    
    # Check that property values are of the correct type
    assert isinstance(properties['molecular_weight'], float)
    assert isinstance(properties['logp'], float)

def test_calculate_molecular_properties_specific():
    """Tests calculation of specific molecular properties."""
    # Test with a specific property list
    properties = calculate_molecular_properties('CCO', ['molecular_weight', 'logp'])
    
    # Check that only requested properties are returned
    assert set(properties.keys()) == {'molecular_weight', 'logp'}
    
    # Check approximate values (may vary slightly by RDKit version)
    assert abs(properties['molecular_weight'] - 46.07) < 0.01
    assert abs(properties['logp'] - (-0.14)) < 0.1

def test_calculate_molecular_properties_invalid():
    """Tests that invalid SMILES strings raise ValidationException during property calculation."""
    for smiles in INVALID_SMILES:
        with pytest.raises(ValidationException):
            calculate_molecular_properties(smiles)

def test_is_property_in_range():
    """Tests checking if property values are within specified ranges."""
    # Test with built-in ranges
    assert is_property_in_range('molecular_weight', 100) is True
    assert is_property_in_range('molecular_weight', 3000) is False
    
    assert is_property_in_range('logp', 2) is True
    assert is_property_in_range('logp', 15) is False
    
    # Test with non-existent property (should return True as there's no validation)
    assert is_property_in_range('non_existent_property', 100) is True

def test_is_property_in_range_custom():
    """Tests checking property values against custom ranges."""
    custom_ranges = {
        'molecular_weight': {'min': 0, 'max': 300},
        'logp': {'min': -2, 'max': 3}
    }
    
    # Test within custom ranges
    assert is_property_in_range('molecular_weight', 200, custom_ranges) is True
    assert is_property_in_range('logp', 1, custom_ranges) is True
    
    # Test outside custom ranges
    assert is_property_in_range('molecular_weight', 400, custom_ranges) is False
    assert is_property_in_range('logp', 4, custom_ranges) is False

def test_filter_molecules_by_property():
    """Tests filtering molecules based on property values."""
    molecules = TEST_MOLECULES
    
    # Test filtering by molecular weight with a range
    filtered = filter_molecules_by_property(molecules, 'molecular_weight', 50, 100)
    assert len(filtered) == 1
    assert filtered[0]['smiles'] == 'c1ccccc1'
    
    # Test filtering by logp with only min value
    filtered = filter_molecules_by_property(molecules, 'logp', 0)
    assert len(filtered) == 1
    assert filtered[0]['smiles'] == 'c1ccccc1'
    
    # Test filtering by logp with only max value
    filtered = filter_molecules_by_property(molecules, 'logp', max_value=0)
    assert len(filtered) == 2
    assert filtered[0]['smiles'] == 'CCO'
    assert filtered[1]['smiles'] == 'CC(=O)O'

def test_filter_molecules_by_property_no_range():
    """Tests filtering molecules with default range (no min/max specified)."""
    molecules = TEST_MOLECULES
    
    # Should return all molecules that have the property
    filtered = filter_molecules_by_property(molecules, 'molecular_weight')
    assert len(filtered) == 3
    
    # Test with property that doesn't exist in the data
    filtered = filter_molecules_by_property(molecules, 'non_existent_property')
    assert len(filtered) == 0

def test_filter_molecules_by_property_empty():
    """Tests filtering an empty list of molecules."""
    filtered = filter_molecules_by_property([], 'molecular_weight')
    assert filtered == []

def test_calculate_similarity():
    """Tests calculation of similarity between molecules."""
    # Test same molecule (should have similarity 1.0)
    similarity = calculate_similarity('CCO', 'CCO')
    assert similarity == 1.0
    
    # Test similar molecules
    similarity = calculate_similarity('CCO', 'CCCO')  # Ethanol vs. Propanol
    assert 0.5 <= similarity <= 0.9  # Similarity should be high but less than 1
    
    # Test dissimilar molecules
    similarity = calculate_similarity('CCO', 'c1ccccc1')  # Ethanol vs. Benzene
    assert 0.0 <= similarity <= 0.3  # Similarity should be low
    
    # Test with caffeine and theophylline (structurally related)
    caffeine = 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'
    theophylline = 'CN1C=NC2=C1C(=O)NC(=O)N2'
    similarity = calculate_similarity(caffeine, theophylline)
    assert 0.7 <= similarity <= 1.0  # Should be quite similar

def test_calculate_similarity_different_fingerprints():
    """Tests similarity calculation with different fingerprint types."""
    # Define test SMILES
    smiles1 = 'CCO'
    smiles2 = 'CCCO'
    
    # Test with different fingerprint types
    morgan_sim = calculate_similarity(smiles1, smiles2, 'morgan')
    maccs_sim = calculate_similarity(smiles1, smiles2, 'maccs')
    topo_sim = calculate_similarity(smiles1, smiles2, 'topological')
    
    # All similarities should be between 0 and 1
    assert 0 <= morgan_sim <= 1
    assert 0 <= maccs_sim <= 1
    assert 0 <= topo_sim <= 1
    
    # Different fingerprints may give different similarity values
    # but we just check they're valid similarity scores

def test_calculate_similarity_invalid():
    """Tests that invalid SMILES strings raise ValidationException during similarity calculation."""
    valid_smiles = 'CCO'
    
    # Test with first SMILES invalid
    with pytest.raises(ValidationException):
        calculate_similarity('invalid', valid_smiles)
    
    # Test with second SMILES invalid
    with pytest.raises(ValidationException):
        calculate_similarity(valid_smiles, 'invalid')
    
    # Test with both SMILES invalid
    with pytest.raises(ValidationException):
        calculate_similarity('invalid', 'also_invalid')

def test_check_lipinski_rule_of_five():
    """Tests checking Lipinski's Rule of Five for drug-likeness."""
    # Ethanol - should pass all Lipinski rules
    result = check_lipinski_rule_of_five('CCO')
    assert result['is_drug_like'] is True
    assert result['violations'] == 0
    assert result['mw_pass'] is True
    assert result['logp_pass'] is True
    assert result['donors_pass'] is True
    assert result['acceptors_pass'] is True
    
    # Test a larger molecule that might violate some rules
    # Example: Cyclosporine A - large cyclic peptide, violates multiple Lipinski rules
    large_molecule = 'CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1)C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C'
    result = check_lipinski_rule_of_five(large_molecule)
    assert result['violations'] > 0  # Should violate at least one rule

def test_check_veber_rules():
    """Tests checking Veber's rules for oral bioavailability."""
    # Ethanol - should pass all Veber rules
    result = check_veber_rules('CCO')
    assert result['passes_veber_rules'] is True
    assert result['rotatable_pass'] is True
    assert result['psa_pass'] is True
    
    # Test a molecule that might violate Veber rules
    # Example: A molecule with many rotatable bonds and high PSA
    complex_molecule = 'O=C(O)CCCCCCC(=O)NCCOCCOCCNC(=O)c1ccc(OCCOCCOCCOC)cc1'
    result = check_veber_rules(complex_molecule)
    assert result['passes_veber_rules'] is False  # Should violate at least one rule

def test_check_drug_likeness():
    """Tests checking combined drug-likeness criteria."""
    # Ethanol - should be drug-like
    result = check_drug_likeness('CCO')
    assert result['is_drug_like'] is True
    assert result['lipinski']['violations'] == 0
    assert result['veber']['passes_veber_rules'] is True
    
    # Aspirin - should be drug-like
    result = check_drug_likeness('CC(=O)Oc1ccccc1C(=O)O')
    assert result['is_drug_like'] is True
    
    # Test a non-drug-like molecule
    # Example: A large, complex molecule that violates many rules
    non_drug_like = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
    result = check_drug_likeness(non_drug_like)
    assert result['is_drug_like'] is False

def test_has_substructure():
    """Tests checking if molecules contain specific substructures."""
    # Test alcohol substructure
    assert has_substructure('CCO', '[OH]') is True
    assert has_substructure('CCC', '[OH]') is False
    
    # Test benzene substructure
    assert has_substructure('c1ccccc1', 'c1ccccc1') is True
    assert has_substructure('CCO', 'c1ccccc1') is False
    
    # Test carboxylic acid substructure
    assert has_substructure('CC(=O)O', 'C(=O)[OH]') is True
    assert has_substructure('CCO', 'C(=O)[OH]') is False

def test_has_substructure_invalid():
    """Tests that invalid SMILES or SMARTS patterns are handled correctly."""
    # Test with invalid SMILES
    with pytest.raises(ValidationException):
        has_substructure('invalid', '[OH]')
    
    # Test with invalid SMARTS pattern
    with pytest.raises(MolecularProcessingException):
        has_substructure('CCO', 'invalid_smarts')

def test_get_molecular_formula():
    """Tests generation of molecular formula from SMILES."""
    # Test simple molecules
    assert get_molecular_formula('CCO') == 'C2H6O'
    assert get_molecular_formula('c1ccccc1') == 'C6H6'
    assert get_molecular_formula('CC(=O)O') == 'C2H4O2'
    
    # Test with multiple elements
    assert get_molecular_formula('CCN') == 'C2H7N'
    assert get_molecular_formula('c1ccccc1Cl') == 'C6H5Cl'
    
    # Test ordering (C, H, then alphabetically)
    assert get_molecular_formula('BrCCl') == 'C1H2BrCl'

def test_get_molecular_formula_invalid():
    """Tests that invalid SMILES strings raise ValidationException during formula generation."""
    for smiles in INVALID_SMILES:
        with pytest.raises(ValidationException):
            get_molecular_formula(smiles)