"""
Unit tests for the molecular validator module.

This module contains tests for validating SMILES structures and molecular properties,
including structure validation, property range validation, and batch processing.
"""

import pytest
from typing import List, Dict

from rdkit import Chem
from rdkit.Chem import Descriptors

from ../../app/molecular.validator import (
    validate_smiles_structure,
    validate_smiles_structure_with_exception,
    validate_property_value,
    validate_molecule_properties,
    validate_molecule,
    validate_molecules_batch,
    check_lipinski_rule_of_five,
    check_veber_rules,
    get_property_range,
    PROPERTY_RANGES
)
from ../../app/molecular.molecule_converter import smiles_to_mol
from ../../app/exceptions import ValidationException, MolecularProcessingException

# Test data
VALID_SMILES = [
    'CCO',  # Ethanol
    'c1ccccc1',  # Benzene
    'CC(=O)O',  # Acetic acid
    'C1CCCCC1',  # Cyclohexane
    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C'  # Caffeine
]

INVALID_SMILES = [
    '',  # Empty string
    'X',  # Invalid atom
    'C1CC',  # Unclosed ring
    'c1ccccc',  # Unclosed aromatic ring
    'malformedsmiles'  # Not a SMILES string
]

TEST_PROPERTIES = {
    "molecular_weight": 46.07,
    "logp": -0.14,
    "h_bond_donors": 1,
    "h_bond_acceptors": 1,
    "rotatable_bonds": 1,
    "polar_surface_area": 20.23,
    "heavy_atom_count": 3,
    "ring_count": 0
}

OUT_OF_RANGE_PROPERTIES = {
    "molecular_weight": 1500,
    "logp": 15,
    "h_bond_donors": 20,
    "h_bond_acceptors": 30,
    "rotatable_bonds": 25,
    "polar_surface_area": 300,
    "heavy_atom_count": 150,
    "ring_count": 15
}


def test_validate_smiles_structure_valid():
    """Tests that valid SMILES strings are correctly validated."""
    for smiles in VALID_SMILES:
        result = validate_smiles_structure(smiles)
        assert result is True, f"Expected '{smiles}' to be valid, but it was not."


def test_validate_smiles_structure_invalid():
    """Tests that invalid SMILES strings are correctly identified."""
    for smiles in INVALID_SMILES:
        result = validate_smiles_structure(smiles)
        assert result is False, f"Expected '{smiles}' to be invalid, but it was not."


def test_validate_smiles_structure_with_exception_valid():
    """Tests that valid SMILES strings are correctly validated with exception handling."""
    for smiles in VALID_SMILES:
        result = validate_smiles_structure_with_exception(smiles)
        assert result is True, f"Expected '{smiles}' to be valid, but it was not."


def test_validate_smiles_structure_with_exception_invalid():
    """Tests that invalid SMILES strings raise ValidationException."""
    for smiles in INVALID_SMILES:
        with pytest.raises(ValidationException) as excinfo:
            validate_smiles_structure_with_exception(smiles)
        assert "Invalid SMILES structure" in str(excinfo.value) or "SMILES string cannot be empty" in str(excinfo.value)


def test_validate_property_value_in_range():
    """Tests that property values within range are correctly validated."""
    for prop_name, value in TEST_PROPERTIES.items():
        if prop_name in PROPERTY_RANGES:  # Only test properties with defined ranges
            result = validate_property_value(prop_name, value)
            assert result is True, f"Expected property '{prop_name}' with value {value} to be in range, but it was not."


def test_validate_property_value_out_of_range():
    """Tests that property values outside range are correctly identified."""
    for prop_name, value in OUT_OF_RANGE_PROPERTIES.items():
        if prop_name in PROPERTY_RANGES:  # Only test properties with defined ranges
            result = validate_property_value(prop_name, value)
            assert result is False, f"Expected property '{prop_name}' with value {value} to be out of range, but it was not."


def test_validate_property_value_unknown_property():
    """Tests that unknown properties are handled correctly."""
    result = validate_property_value("unknown_property", 42)
    assert result is True, "Expected unknown property to be considered valid, but it was not."


def test_validate_molecule_properties():
    """Tests validation of multiple molecular properties."""
    result = validate_molecule_properties(TEST_PROPERTIES)
    assert isinstance(result, dict), "Expected result to be a dictionary."
    
    # Check that all properties are valid
    for prop_name, is_valid in result.items():
        assert is_valid is True, f"Expected property '{prop_name}' to be valid, but it was not."


def test_validate_molecule_properties_mixed():
    """Tests validation of properties with mixed validity."""
    # Create a dictionary with both in-range and out-of-range properties
    mixed_properties = {
        "molecular_weight": 46.07,  # in range
        "logp": 15,  # out of range
        "h_bond_donors": 1,  # in range
        "h_bond_acceptors": 30  # out of range
    }
    
    result = validate_molecule_properties(mixed_properties)
    assert isinstance(result, dict), "Expected result to be a dictionary."
    
    # Check specific properties
    assert result["molecular_weight"] is True, "Expected 'molecular_weight' to be valid."
    assert result["h_bond_donors"] is True, "Expected 'h_bond_donors' to be valid."
    assert result["logp"] is False, "Expected 'logp' to be invalid."
    assert result["h_bond_acceptors"] is False, "Expected 'h_bond_acceptors' to be invalid."


def test_validate_molecule_valid():
    """Tests validation of a valid molecule with valid properties."""
    smiles = VALID_SMILES[0]  # Use first valid SMILES
    result = validate_molecule(smiles, TEST_PROPERTIES)
    
    assert "structure_valid" in result, "Expected 'structure_valid' key in result."
    assert result["structure_valid"] is True, "Expected structure to be valid."
    assert "property_validation" in result, "Expected 'property_validation' key in result."
    
    # Check that all properties are valid
    for prop_name, is_valid in result["property_validation"].items():
        assert is_valid is True, f"Expected property '{prop_name}' to be valid, but it was not."
    
    assert result["all_properties_valid"] is True, "Expected all properties to be valid."


def test_validate_molecule_invalid_structure():
    """Tests validation of a molecule with invalid structure."""
    smiles = INVALID_SMILES[1]  # Use an invalid SMILES
    result = validate_molecule(smiles, TEST_PROPERTIES)
    
    assert "structure_valid" in result, "Expected 'structure_valid' key in result."
    assert result["structure_valid"] is False, "Expected structure to be invalid."
    assert "property_validation" not in result, "Expected no property validation for invalid structure."
    assert "all_properties_valid" not in result, "Expected no all_properties_valid flag for invalid structure."


def test_validate_molecule_no_properties():
    """Tests validation of a molecule without properties."""
    smiles = VALID_SMILES[0]  # Use first valid SMILES
    result = validate_molecule(smiles)
    
    assert "structure_valid" in result, "Expected 'structure_valid' key in result."
    assert result["structure_valid"] is True, "Expected structure to be valid."
    assert "property_validation" not in result, "Expected no property validation when no properties provided."
    assert "all_properties_valid" not in result, "Expected no all_properties_valid flag when no properties provided."


def test_validate_molecules_batch():
    """Tests batch validation of multiple molecules."""
    # Create a list of molecules with SMILES and properties
    molecules = [
        {"smiles": VALID_SMILES[0], "properties": TEST_PROPERTIES},
        {"smiles": INVALID_SMILES[1], "properties": TEST_PROPERTIES},
        {"smiles": VALID_SMILES[1], "properties": OUT_OF_RANGE_PROPERTIES}
    ]
    
    results = validate_molecules_batch(molecules)
    
    assert len(results) == len(molecules), "Expected same number of results as input molecules."
    
    # First molecule should have valid structure and properties
    assert results[0]["validation"]["structure_valid"] is True, "Expected first molecule structure to be valid."
    assert "property_validation" in results[0]["validation"], "Expected property validation for first molecule."
    
    # Second molecule should have invalid structure
    assert results[1]["validation"]["structure_valid"] is False, "Expected second molecule structure to be invalid."
    
    # Third molecule should have valid structure but some invalid properties
    assert results[2]["validation"]["structure_valid"] is True, "Expected third molecule structure to be valid."
    assert "property_validation" in results[2]["validation"], "Expected property validation for third molecule."
    assert results[2]["validation"]["all_properties_valid"] is False, "Expected some invalid properties for third molecule."


def test_validate_molecules_batch_empty():
    """Tests batch validation with an empty list."""
    results = validate_molecules_batch([])
    assert results == [], "Expected empty result list for empty input."


def test_check_lipinski_rule_of_five():
    """Tests Lipinski's Rule of Five checking."""
    # Ethanol should pass Lipinski's rules
    mol = smiles_to_mol('CCO')
    result = check_lipinski_rule_of_five(mol)
    
    assert isinstance(result, dict), "Expected result to be a dictionary."
    assert "molecular_weight" in result, "Expected 'molecular_weight' key in result."
    assert "logp" in result, "Expected 'logp' key in result."
    assert "h_bond_donors" in result, "Expected 'h_bond_donors' key in result."
    assert "h_bond_acceptors" in result, "Expected 'h_bond_acceptors' key in result."
    assert "mw_pass" in result, "Expected 'mw_pass' key in result."
    assert "logp_pass" in result, "Expected 'logp_pass' key in result."
    assert "donors_pass" in result, "Expected 'donors_pass' key in result."
    assert "acceptors_pass" in result, "Expected 'acceptors_pass' key in result."
    assert "violations" in result, "Expected 'violations' key in result."
    assert "is_drug_like" in result, "Expected 'is_drug_like' key in result."
    
    assert result["violations"] == 0, "Expected 0 violations for ethanol."
    assert result["is_drug_like"] is True, "Expected ethanol to be drug-like."


def test_check_lipinski_rule_of_five_violations():
    """Tests Lipinski's Rule of Five with violations."""
    # Create a molecule that violates Lipinski's rules
    # This is a large molecule that should violate multiple rules
    large_molecule = smiles_to_mol('CCCCCCCCCCCCCCCCCCCCC(=O)NCCCCCCCCCCCCCCCCCCCC')
    
    result = check_lipinski_rule_of_five(large_molecule)
    
    assert result["molecular_weight"] > 500, "Expected molecular weight > 500."
    assert result["logp"] > 5, "Expected logP > 5."
    assert result["violations"] > 0, "Expected at least one Lipinski violation."
    assert result["is_drug_like"] is False, "Expected molecule to not be drug-like."


def test_check_veber_rules():
    """Tests Veber's rules checking."""
    # Ethanol should pass Veber's rules
    mol = smiles_to_mol('CCO')
    result = check_veber_rules(mol)
    
    assert isinstance(result, dict), "Expected result to be a dictionary."
    assert "rotatable_bonds" in result, "Expected 'rotatable_bonds' key in result."
    assert "polar_surface_area" in result, "Expected 'polar_surface_area' key in result."
    assert "rotatable_pass" in result, "Expected 'rotatable_pass' key in result."
    assert "psa_pass" in result, "Expected 'psa_pass' key in result."
    assert "passes_veber_rules" in result, "Expected 'passes_veber_rules' key in result."
    
    assert result["rotatable_pass"] is True, "Expected rotatable bonds to pass."
    assert result["psa_pass"] is True, "Expected polar surface area to pass."
    assert result["passes_veber_rules"] is True, "Expected ethanol to pass Veber's rules."


def test_check_veber_rules_violations():
    """Tests Veber's rules with violations."""
    # Create a molecule that violates Veber's rules
    # This molecule has many rotatable bonds and high polar surface area
    violation_molecule = smiles_to_mol('CNC(=O)C(CNC(=O)C(N)CCCCN)NC(=O)C(CC(C)C)NC(=O)C(CC(=O)O)NC(=O)C(CCSC)NC(=O)C(CCCCN)NC(=O)CCCCCCCCCNC(=O)CCCC')
    
    result = check_veber_rules(violation_molecule)
    
    assert result["rotatable_bonds"] > 10, "Expected rotatable bonds > 10."
    assert result["polar_surface_area"] > 140, "Expected polar surface area > 140."
    assert result["rotatable_pass"] is False, "Expected rotatable bonds to fail."
    assert result["psa_pass"] is False, "Expected polar surface area to fail."
    assert result["passes_veber_rules"] is False, "Expected molecule to fail Veber's rules."


def test_get_property_range_known():
    """Tests getting range for a known property."""
    property_name = "molecular_weight"
    result = get_property_range(property_name)
    
    assert isinstance(result, dict), "Expected result to be a dictionary."
    assert "min" in result, "Expected 'min' key in result."
    assert "max" in result, "Expected 'max' key in result."
    assert result["min"] == PROPERTY_RANGES[property_name]["min"], "Expected min value to match PROPERTY_RANGES."
    assert result["max"] == PROPERTY_RANGES[property_name]["max"], "Expected max value to match PROPERTY_RANGES."


def test_get_property_range_unknown():
    """Tests getting range for an unknown property."""
    property_name = "unknown_property"
    result = get_property_range(property_name)
    
    assert isinstance(result, dict), "Expected result to be a dictionary."
    assert "min" in result, "Expected 'min' key in result."
    assert "max" in result, "Expected 'max' key in result."
    assert result["min"] == float('-inf'), "Expected min value to be negative infinity for unknown property."
    assert result["max"] == float('inf'), "Expected max value to be infinity for unknown property."