import pytest
import random
from typing import List, Dict

from app.crud.crud_molecule import molecule
from app.models.molecule import Molecule
from app.models.molecule_property import MoleculeProperty
from app.schemas.molecule import MoleculeCreate


def test_get_by_smiles(db_session, create_test_molecule):
    """Test retrieving a molecule by its SMILES string"""
    # Create a test molecule with a specific SMILES string
    test_smiles = "CCO"  # Ethanol
    test_molecule = create_test_molecule(test_smiles, 1)
    
    # Retrieve the molecule using get_by_smiles function
    db_molecule = molecule.get_by_smiles(db_session, test_smiles)
    
    # Assert that the retrieved molecule has the correct SMILES string
    assert db_molecule is not None
    assert db_molecule.smiles == test_smiles
    
    # Test with a non-existent SMILES string
    non_existent_smiles = "CC(=O)C"  # Acetone (not in the database)
    result = molecule.get_by_smiles(db_session, non_existent_smiles)
    
    # Assert that None is returned for non-existent SMILES
    assert result is None


def test_get_with_properties(db_session, create_test_molecule):
    """Test retrieving a molecule with its properties"""
    # Create a test molecule with properties
    test_smiles = "CCO"  # Ethanol
    properties = {
        "MW": 46.07,
        "LogP": -0.14,
        "TPSA": 20.23
    }
    test_molecule = create_test_molecule(test_smiles, 1, properties)
    
    # Retrieve the molecule with properties using get_with_properties function
    db_molecule = molecule.get_with_properties(db_session, test_molecule.id)
    
    # Assert that the molecule has the correct number of properties
    assert db_molecule is not None
    assert len(db_molecule.properties) == len(properties)
    
    # Assert that the property values match the expected values
    property_dict = {prop.property_name: prop.property_value for prop in db_molecule.properties}
    for name, value in properties.items():
        assert name in property_dict
        assert float(property_dict[name]) == value
    
    # Test with a non-existent molecule ID
    non_existent_id = 9999
    result = molecule.get_with_properties(db_session, non_existent_id)
    
    # Assert that None is returned for non-existent ID
    assert result is None


def test_get_by_property_range(db_session, create_test_molecule):
    """Test filtering molecules by property range"""
    # Create multiple test molecules with different property values
    test_molecules = [
        create_test_molecule("CCO", 1, {"MW": 46.07, "LogP": -0.14}),  # Ethanol
        create_test_molecule("CCCO", 1, {"MW": 60.10, "LogP": 0.25}),  # Propanol
        create_test_molecule("CCCCO", 1, {"MW": 74.12, "LogP": 0.88}),  # Butanol
        create_test_molecule("CCCCCO", 1, {"MW": 88.15, "LogP": 1.51}),  # Pentanol
    ]
    
    # Test filtering with min value only
    results = molecule.get_by_property_range(db_session, "MW", min_value=60.0)
    assert len(results) == 3  # Should return Propanol, Butanol, Pentanol
    
    # Test filtering with max value only
    results = molecule.get_by_property_range(db_session, "MW", max_value=60.10)
    assert len(results) == 2  # Should return Ethanol, Propanol
    
    # Test filtering with both min and max values
    results = molecule.get_by_property_range(db_session, "LogP", min_value=0.0, max_value=1.0)
    assert len(results) == 2  # Should return Propanol, Butanol
    
    # Test with sorting (ascending)
    results = molecule.get_by_property_range(
        db_session, "MW", sort_by="MW", sort_desc=False
    )
    assert results[0].smiles == "CCO"  # Ethanol (lowest MW)
    assert results[-1].smiles == "CCCCCO"  # Pentanol (highest MW)
    
    # Test with sorting (descending)
    results = molecule.get_by_property_range(
        db_session, "MW", sort_by="MW", sort_desc=True
    )
    assert results[0].smiles == "CCCCCO"  # Pentanol (highest MW)
    assert results[-1].smiles == "CCO"  # Ethanol (lowest MW)
    
    # Test with pagination
    results = molecule.get_by_property_range(
        db_session, "MW", skip=1, limit=2
    )
    assert len(results) == 2  # Should return 2 molecules
    assert results[0].smiles == "CCCO"  # Propanol (second lowest MW)


def test_count_by_property_range(db_session, create_test_molecule):
    """Test counting molecules by property range"""
    # Create multiple test molecules with different property values
    test_molecules = [
        create_test_molecule("CCO", 1, {"MW": 46.07, "LogP": -0.14}),  # Ethanol
        create_test_molecule("CCCO", 1, {"MW": 60.10, "LogP": 0.25}),  # Propanol
        create_test_molecule("CCCCO", 1, {"MW": 74.12, "LogP": 0.88}),  # Butanol
        create_test_molecule("CCCCCO", 1, {"MW": 88.15, "LogP": 1.51}),  # Pentanol
    ]
    
    # Test counting with min value only
    count = molecule.count_by_property_range(db_session, "MW", min_value=60.0)
    assert count == 3  # Should count Propanol, Butanol, Pentanol
    
    # Test counting with max value only
    count = molecule.count_by_property_range(db_session, "MW", max_value=60.10)
    assert count == 2  # Should count Ethanol, Propanol
    
    # Test counting with both min and max values
    count = molecule.count_by_property_range(db_session, "LogP", min_value=0.0, max_value=1.0)
    assert count == 2  # Should count Propanol, Butanol
    
    # Test with a property that doesn't exist
    count = molecule.count_by_property_range(db_session, "NonExistentProperty")
    assert count == 0


def test_get_by_flag_status(db_session, create_test_molecule):
    """Test retrieving molecules by flag status"""
    # Create multiple test molecules with different flag statuses
    flagged_molecules = [
        create_test_molecule("CCO", 1, flag_status="HIGH"),  # Ethanol - High priority
        create_test_molecule("CCCO", 1, flag_status="HIGH"),  # Propanol - High priority
        create_test_molecule("CCCCO", 1, flag_status="MEDIUM"),  # Butanol - Medium priority
        create_test_molecule("CCCCCO", 1, flag_status="LOW"),  # Pentanol - Low priority
    ]
    
    # Test retrieving molecules with a specific flag status
    results = molecule.get_by_flag_status(db_session, "HIGH")
    assert len(results) == 2  # Should return Ethanol, Propanol
    for mol in results:
        assert mol.flag_status == "HIGH"
    
    # Test with pagination
    results = molecule.get_by_flag_status(db_session, "HIGH", skip=1, limit=1)
    assert len(results) == 1  # Should return only 1 molecule
    
    # Test with a flag status that no molecules have
    results = molecule.get_by_flag_status(db_session, "CRITICAL")
    assert len(results) == 0  # Should return empty list


def test_update_flag_status(db_session, create_test_molecule):
    """Test updating the flag status of a molecule"""
    # Create a test molecule with an initial flag status
    test_molecule = create_test_molecule("CCO", 1, flag_status="LOW")  # Ethanol - Low priority
    
    # Update the flag status using update_flag_status function
    updated_molecule = molecule.update_flag_status(db_session, test_molecule.id, "HIGH")
    
    # Assert that the flag status was updated correctly
    assert updated_molecule is not None
    assert updated_molecule.flag_status == "HIGH"
    
    # Verify in database
    db_molecule = molecule.get(db_session, test_molecule.id)
    assert db_molecule.flag_status == "HIGH"
    
    # Test with a non-existent molecule ID
    non_existent_id = 9999
    result = molecule.update_flag_status(db_session, non_existent_id, "HIGH")
    
    # Assert that None is returned for non-existent ID
    assert result is None


def test_get_by_user(db_session, create_test_molecule, test_pharma_user):
    """Test retrieving molecules created by a specific user"""
    # Create multiple test molecules with different user IDs
    user1_molecules = [
        create_test_molecule("CCO", test_pharma_user.id),  # Ethanol - User 1
        create_test_molecule("CCCO", test_pharma_user.id),  # Propanol - User 1
    ]
    
    user2_id = test_pharma_user.id + 1  # Different user ID
    user2_molecules = [
        create_test_molecule("CCCCO", user2_id),  # Butanol - User 2
    ]
    
    # Test retrieving molecules for a specific user
    results = molecule.get_by_user(db_session, test_pharma_user.id)
    assert len(results) == 2  # Should return User 1's molecules
    for mol in results:
        assert mol.created_by == test_pharma_user.id
    
    # Test with pagination
    results = molecule.get_by_user(db_session, test_pharma_user.id, skip=1, limit=1)
    assert len(results) == 1  # Should return only 1 molecule
    
    # Test with a user ID that has no molecules
    non_existent_user_id = 9999
    results = molecule.get_by_user(db_session, non_existent_user_id)
    assert len(results) == 0  # Should return empty list


def test_get_property_ranges(db_session, create_test_molecule):
    """Test retrieving the min and max values for a property"""
    # Create multiple test molecules with different property values
    test_molecules = [
        create_test_molecule("CCO", 1, {"MW": 46.07, "LogP": -0.14}),  # Ethanol
        create_test_molecule("CCCO", 1, {"MW": 60.10, "LogP": 0.25}),  # Propanol
        create_test_molecule("CCCCO", 1, {"MW": 74.12, "LogP": 0.88}),  # Butanol
        create_test_molecule("CCCCCO", 1, {"MW": 88.15, "LogP": 1.51}),  # Pentanol
    ]
    
    # Test retrieving the property range for MW
    min_val, max_val = molecule.get_property_ranges(db_session, "MW")
    assert min_val == 46.07  # Ethanol MW
    assert max_val == 88.15  # Pentanol MW
    
    # Test retrieving the property range for LogP
    min_val, max_val = molecule.get_property_ranges(db_session, "LogP")
    assert min_val == -0.14  # Ethanol LogP
    assert max_val == 1.51  # Pentanol LogP
    
    # Test with a property that doesn't exist
    min_val, max_val = molecule.get_property_ranges(db_session, "NonExistentProperty")
    assert min_val == 0.0  # Default values
    assert max_val == 0.0  # Default values


def test_get_available_properties(db_session, create_test_molecule):
    """Test retrieving all available property names"""
    # Create multiple test molecules with different properties
    test_molecules = [
        create_test_molecule("CCO", 1, {"MW": 46.07, "LogP": -0.14, "TPSA": 20.23}),  # Ethanol
        create_test_molecule("CCCO", 1, {"MW": 60.10, "LogP": 0.25, "RotBonds": 1}),  # Propanol
        create_test_molecule("CCCCO", 1, {"MW": 74.12, "HBA": 1, "HBD": 1}),  # Butanol
    ]
    
    # Test retrieving all available property names
    properties = molecule.get_available_properties(db_session)
    
    # Assert that all expected property names are returned
    expected_properties = ["MW", "LogP", "TPSA", "RotBonds", "HBA", "HBD"]
    for prop_name in expected_properties:
        assert prop_name in properties
    
    # Assert that duplicate property names are not returned
    assert len(properties) == len(set(properties))  # Should have no duplicates


def test_bulk_create(db_session, test_pharma_user):
    """Test creating multiple molecules in a single transaction"""
    # Prepare a list of molecule data for bulk creation
    molecules_data = [
        {
            "smiles": "CCO",  # Ethanol
            "properties": [
                {"property_name": "MW", "property_value": 46.07},
                {"property_name": "LogP", "property_value": -0.14}
            ]
        },
        {
            "smiles": "CCCO",  # Propanol
            "properties": [
                {"property_name": "MW", "property_value": 60.10},
                {"property_name": "LogP", "property_value": 0.25}
            ]
        },
        {
            "smiles": "CCCCO",  # Butanol
            "properties": [
                {"property_name": "MW", "property_value": 74.12},
                {"property_name": "LogP", "property_value": 0.88}
            ]
        }
    ]
    
    # Call the bulk_create function with the list
    result = molecule.bulk_create(db_session, molecules_data, test_pharma_user.id)
    
    # Assert that the correct number of molecules were created
    assert result["total_created"] == 3
    assert len(result["created_molecules"]) == 3
    
    # Assert that the molecules have the correct properties
    for i, mol_data in enumerate(molecules_data):
        # Get the molecule from database
        db_mol = molecule.get_by_smiles(db_session, mol_data["smiles"])
        assert db_mol is not None
        
        # Check properties
        db_mol_with_props = molecule.get_with_properties(db_session, db_mol.id)
        assert db_mol_with_props is not None
        
        # Create a dict of property name -> value for easier comparison
        db_props = {prop.property_name: prop.property_value for prop in db_mol_with_props.properties}
        for expected_prop in mol_data["properties"]:
            prop_name = expected_prop["property_name"]
            assert prop_name in db_props
            assert float(db_props[prop_name]) == expected_prop["property_value"]
    
    # Test with duplicate SMILES strings
    duplicate_molecules = [
        {
            "smiles": "CCO",  # Ethanol (already exists)
            "properties": [
                {"property_name": "TPSA", "property_value": 20.23}
            ]
        },
        {
            "smiles": "c1ccccc1",  # Benzene (new molecule)
            "properties": [
                {"property_name": "MW", "property_value": 78.11}
            ]
        }
    ]
    
    result = molecule.bulk_create(db_session, duplicate_molecules, test_pharma_user.id)
    
    # Assert that only the new molecule was created
    assert result["total_created"] == 1
    assert result["total_duplicates"] == 1
    assert result["total_processed"] == 2
    
    # Verify that the new molecule was created
    benzene = molecule.get_by_smiles(db_session, "c1ccccc1")
    assert benzene is not None
    assert benzene.smiles == "c1ccccc1"
    
    # Verify property was added to benzene
    benzene_with_props = molecule.get_with_properties(db_session, benzene.id)
    prop_names = [prop.property_name for prop in benzene_with_props.properties]
    assert "MW" in prop_names