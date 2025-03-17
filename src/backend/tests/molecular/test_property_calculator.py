import pytest
from typing import List, Dict, Any

from rdkit import Chem  # version 2023.03+
from rdkit.Chem import Descriptors  # version 2023.03+

from ../../app/molecular/property_calculator import (
    calculate_property,
    calculate_properties,
    calculate_properties_from_smiles,
    batch_calculate_properties,
    batch_calculate_properties_from_smiles,
    check_lipinski_rule_of_five,
    check_veber_rules,
    register_property_calculator,
    get_available_properties,
    PropertyCalculator,
    PROPERTY_CALCULATORS
)
from ../../app/molecular/molecule_converter import smiles_to_mol
from ../../app/exceptions import MolecularProcessingException

# Valid SMILES strings for testing
VALID_SMILES = ['CCO', 'c1ccccc1', 'CC(=O)O', 'C1CCCCC1', 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C']  # ethanol, benzene, acetic acid, cyclohexane, caffeine

# Invalid SMILES strings for testing
INVALID_SMILES = ['', 'X', 'C1CC', 'c1ccccc', 'malformedsmiles']

# Properties to test
TEST_PROPERTIES = ["molecular_weight", "logp", "h_bond_donors", "h_bond_acceptors", "rotatable_bonds"]

# Custom property for testing
CUSTOM_PROPERTY_NAME = "custom_property"
CUSTOM_PROPERTY_FUNCTION = lambda mol: 42.0  # Simple function that always returns 42.0


def test_calculate_property_valid():
    """Tests calculation of a single property for a valid molecule"""
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test calculating molecular weight
    result = calculate_property(mol, "molecular_weight")
    
    # Assert the result is a float
    assert isinstance(result, float)
    
    # Assert the result is approximately the expected value for ethanol (46.07)
    assert 45.0 < result < 47.0


def test_calculate_property_invalid_molecule():
    """Tests that calculate_property raises exception for invalid molecule"""
    # Test that it raises an exception when molecule is None
    with pytest.raises(MolecularProcessingException) as excinfo:
        calculate_property(None, "molecular_weight")
    
    # Assert that the exception message contains expected text
    assert "Cannot calculate property for None molecule" in str(excinfo.value)


def test_calculate_property_invalid_property():
    """Tests that calculate_property raises exception for invalid property name"""
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test that it raises an exception when property name is invalid
    with pytest.raises(MolecularProcessingException) as excinfo:
        calculate_property(mol, "invalid_property")
    
    # Assert that the exception message contains expected text
    assert "Unknown property" in str(excinfo.value)
    
    # Assert that the exception details contain available properties
    assert "available_properties" in excinfo.value.details


def test_calculate_properties_valid():
    """Tests calculation of multiple properties for a valid molecule"""
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test calculating multiple properties
    results = calculate_properties(mol, TEST_PROPERTIES)
    
    # Assert the result is a dictionary
    assert isinstance(results, dict)
    
    # Assert all requested properties are in the result
    for prop in TEST_PROPERTIES:
        assert prop in results
        
    # Assert all property values are floats
    for value in results.values():
        assert isinstance(value, float)
    
    # Assert some specific property values for ethanol
    assert 45.0 < results["molecular_weight"] < 47.0  # ~46.07
    assert -0.5 < results["logp"] < 0.5  # ~-0.14
    assert 0.9 < results["h_bond_donors"] < 1.1  # 1
    assert 0.9 < results["h_bond_acceptors"] < 1.1  # 1


def test_calculate_properties_all():
    """Tests calculation of all available properties when none specified"""
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test calculating all properties (None = all)
    results = calculate_properties(mol)
    
    # Assert the result is a dictionary
    assert isinstance(results, dict)
    
    # Assert all available properties are in the result
    for prop in PROPERTY_CALCULATORS.keys():
        assert prop in results
        
    # Assert all property values are floats
    for value in results.values():
        assert isinstance(value, float)


def test_calculate_properties_invalid_molecule():
    """Tests that calculate_properties raises exception for invalid molecule"""
    # Test that it raises an exception when molecule is None
    with pytest.raises(MolecularProcessingException) as excinfo:
        calculate_properties(None, TEST_PROPERTIES)
    
    # Assert that the exception message contains expected text
    assert "Cannot calculate properties for None molecule" in str(excinfo.value)


def test_calculate_properties_from_smiles_valid():
    """Tests calculation of properties from a valid SMILES string"""
    # Test calculating properties from a valid SMILES
    results = calculate_properties_from_smiles(VALID_SMILES[0], TEST_PROPERTIES)
    
    # Assert the result is a dictionary
    assert isinstance(results, dict)
    
    # Assert all requested properties are in the result
    for prop in TEST_PROPERTIES:
        assert prop in results
        
    # Assert all property values are floats
    for value in results.values():
        assert isinstance(value, float)
    
    # Assert some specific property values for ethanol
    assert 45.0 < results["molecular_weight"] < 47.0  # ~46.07


def test_calculate_properties_from_smiles_invalid():
    """Tests that calculate_properties_from_smiles raises exception for invalid SMILES"""
    # Test that it raises an exception for invalid SMILES
    with pytest.raises(MolecularProcessingException) as excinfo:
        calculate_properties_from_smiles(INVALID_SMILES[1], TEST_PROPERTIES)
    
    # Assert that the exception message contains expected text
    assert "Invalid SMILES" in str(excinfo.value)


def test_batch_calculate_properties():
    """Tests batch calculation of properties for multiple molecules"""
    # Create a list of molecules from valid SMILES
    mols = [smiles_to_mol(smiles) for smiles in VALID_SMILES]
    
    # Test batch calculation of properties
    results = batch_calculate_properties(mols, TEST_PROPERTIES)
    
    # Assert the result is a list of the same length as input
    assert isinstance(results, list)
    assert len(results) == len(mols)
    
    # Assert each result is a dictionary with all requested properties
    for result in results:
        assert isinstance(result, dict)
        for prop in TEST_PROPERTIES:
            assert prop in result
            assert isinstance(result[prop], float)


def test_batch_calculate_properties_empty():
    """Tests batch calculation with an empty list"""
    # Test batch calculation with an empty list
    results = batch_calculate_properties([], TEST_PROPERTIES)
    
    # Assert the result is an empty list
    assert isinstance(results, list)
    assert len(results) == 0


def test_batch_calculate_properties_from_smiles():
    """Tests batch calculation of properties from multiple SMILES strings"""
    # Test batch calculation from valid SMILES strings
    results = batch_calculate_properties_from_smiles(VALID_SMILES, TEST_PROPERTIES)
    
    # Assert the result is a list of the same length as input
    assert isinstance(results, list)
    assert len(results) == len(VALID_SMILES)
    
    # Assert each result is a dictionary with all requested properties
    for result in results:
        assert isinstance(result, dict)
        for prop in TEST_PROPERTIES:
            assert prop in result
            assert isinstance(result[prop], float)


def test_batch_calculate_properties_from_smiles_with_invalid():
    """Tests batch calculation with some invalid SMILES strings"""
    # Create a mixed list with valid and invalid SMILES
    mixed_smiles = VALID_SMILES[:2] + INVALID_SMILES[:2] + VALID_SMILES[2:]
    
    # Test batch calculation with mixed valid/invalid SMILES
    results = batch_calculate_properties_from_smiles(mixed_smiles, TEST_PROPERTIES)
    
    # Assert the result is a list of the same length as input
    assert isinstance(results, list)
    assert len(results) == len(mixed_smiles)
    
    # Check valid SMILES have properties calculated
    for i, smiles in enumerate(mixed_smiles):
        if smiles in VALID_SMILES:
            assert isinstance(results[i], dict)
            assert len(results[i]) > 0
        else:
            # Invalid SMILES should have empty dictionaries
            assert results[i] == {} or results[i] is None


def test_check_lipinski_rule_of_five():
    """Tests Lipinski's Rule of Five checking"""
    # Create a molecule for testing (ethanol, which passes Lipinski's rules)
    mol = smiles_to_mol(VALID_SMILES[0])
    
    # Test Lipinski rule checking
    result = check_lipinski_rule_of_five(mol)
    
    # Assert the result is a dictionary with expected keys
    assert isinstance(result, dict)
    assert 'molecular_weight' in result
    assert 'logp' in result
    assert 'h_bond_donors' in result
    assert 'h_bond_acceptors' in result
    assert 'mw_pass' in result
    assert 'logp_pass' in result
    assert 'donors_pass' in result
    assert 'acceptors_pass' in result
    assert 'violations' in result
    assert 'is_drug_like' in result
    
    # Assert ethanol passes Lipinski's rules
    assert result['mw_pass'] is True
    assert result['logp_pass'] is True
    assert result['donors_pass'] is True
    assert result['acceptors_pass'] is True
    assert result['violations'] == 0
    assert result['is_drug_like'] is True


def test_check_veber_rules():
    """Tests Veber's rules checking"""
    # Create a molecule that passes Veber's rules
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test Veber rule checking
    result = check_veber_rules(mol)
    
    # Assert the result is a dictionary with expected keys
    assert isinstance(result, dict)
    assert 'rotatable_bonds' in result
    assert 'polar_surface_area' in result
    assert 'rotatable_pass' in result
    assert 'psa_pass' in result
    assert 'passes_veber_rules' in result
    
    # Assert ethanol passes Veber's rules
    assert result['rotatable_pass'] is True
    assert result['psa_pass'] is True
    assert result['passes_veber_rules'] is True


def test_register_property_calculator():
    """Tests registration of a custom property calculator"""
    # Define a custom property calculator function
    def custom_calculator(mol):
        return 42.0
    
    # Register the custom calculator
    register_property_calculator(CUSTOM_PROPERTY_NAME, custom_calculator)
    
    # Create a molecule for testing
    mol = smiles_to_mol(VALID_SMILES[0])
    
    # Test the custom property calculator
    result = calculate_property(mol, CUSTOM_PROPERTY_NAME)
    
    # Assert the result is the expected value from the custom function
    assert result == 42.0


def test_get_available_properties():
    """Tests getting the list of available property calculators"""
    # Get available properties
    properties = get_available_properties()
    
    # Assert the result is a list
    assert isinstance(properties, list)
    
    # Assert it contains expected property names
    for prop in TEST_PROPERTIES:
        assert prop in properties
    
    # Assert it matches the keys in PROPERTY_CALCULATORS
    assert set(properties) == set(PROPERTY_CALCULATORS.keys())


def test_property_calculator_class_init():
    """Tests initialization of PropertyCalculator class"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Assert the instance is created successfully
    assert calculator is not None
    
    # Assert it has the expected attributes
    assert hasattr(calculator, '_calculators')
    assert hasattr(calculator, '_property_ranges')
    assert hasattr(calculator, '_cache')


def test_property_calculator_register_calculator():
    """Tests registration of a custom calculator in PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Define a custom property calculator function
    def custom_calculator(mol):
        return 42.0
    
    # Register the custom calculator
    calculator.register_calculator(CUSTOM_PROPERTY_NAME, custom_calculator)
    
    # Create a molecule for testing
    mol = smiles_to_mol(VALID_SMILES[0])
    
    # Test the custom property calculator
    result = calculator.calculate(mol, CUSTOM_PROPERTY_NAME)
    
    # Assert the result is the expected value from the custom function
    assert result == 42.0


def test_property_calculator_calculate():
    """Tests calculation of a single property using PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test calculating molecular weight
    result = calculator.calculate(mol, "molecular_weight")
    
    # Assert the result is a float
    assert isinstance(result, float)
    
    # Assert the result is approximately the expected value for ethanol (46.07)
    assert 45.0 < result < 47.0


def test_property_calculator_calculate_all():
    """Tests calculation of multiple properties using PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Test calculating multiple properties
    results = calculator.calculate_all(mol, TEST_PROPERTIES)
    
    # Assert the result is a dictionary
    assert isinstance(results, dict)
    
    # Assert all requested properties are in the result
    for prop in TEST_PROPERTIES:
        assert prop in results
        
    # Assert all property values are floats
    for value in results.values():
        assert isinstance(value, float)
    
    # Assert some specific property values for ethanol
    assert 45.0 < results["molecular_weight"] < 47.0  # ~46.07


def test_property_calculator_calculate_from_smiles():
    """Tests calculation of properties from SMILES using PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Test calculating properties from a valid SMILES
    results = calculator.calculate_from_smiles(VALID_SMILES[0], TEST_PROPERTIES)
    
    # Assert the result is a dictionary
    assert isinstance(results, dict)
    
    # Assert all requested properties are in the result
    for prop in TEST_PROPERTIES:
        assert prop in results
        
    # Assert all property values are floats
    for value in results.values():
        assert isinstance(value, float)
    
    # Assert some specific property values for ethanol
    assert 45.0 < results["molecular_weight"] < 47.0  # ~46.07


def test_property_calculator_batch_calculate():
    """Tests batch calculation using PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Create a list of molecules from VALID_SMILES
    mols = [smiles_to_mol(smiles) for smiles in VALID_SMILES]
    
    # Test batch calculation of properties
    results = calculator.batch_calculate(mols, TEST_PROPERTIES)
    
    # Assert the result is a list of the same length as input
    assert isinstance(results, list)
    assert len(results) == len(mols)
    
    # Assert each result is a dictionary with all requested properties
    for result in results:
        assert isinstance(result, dict)
        for prop in TEST_PROPERTIES:
            assert prop in result
            assert isinstance(result[prop], float)


def test_property_calculator_set_property_range():
    """Tests setting custom property range in PropertyCalculator"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Set a custom property range for molecular_weight
    calculator.set_property_range("molecular_weight", 10.0, 100.0)
    
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Calculate molecular weight
    mw = calculator.calculate(mol, "molecular_weight")
    
    # Test validation within the custom range
    assert calculator.validate_property("molecular_weight", mw) is True
    
    # Test validation outside the custom range
    assert calculator.validate_property("molecular_weight", 5.0) is False
    assert calculator.validate_property("molecular_weight", 150.0) is False


def test_property_calculator_clear_cache():
    """Tests clearing the property calculation cache"""
    # Create a PropertyCalculator instance
    calculator = PropertyCalculator()
    
    # Create a molecule from a valid SMILES
    mol = smiles_to_mol(VALID_SMILES[0])  # ethanol
    
    # Calculate a property to populate cache
    calculator.calculate(mol, "molecular_weight")
    
    # Assert the cache is not empty
    assert len(calculator._cache) > 0
    
    # Clear the cache
    calculator.clear_cache()
    
    # Assert the cache is empty
    assert len(calculator._cache) == 0