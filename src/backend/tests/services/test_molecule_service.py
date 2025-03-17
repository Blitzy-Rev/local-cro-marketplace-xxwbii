# src/backend/tests/services/test_molecule_service.py
import pytest
from unittest.mock import MagicMock, patch
import uuid
import io
from io import BytesIO
from rdkit import Chem

from src.backend.app.services.molecule_service import (
    get_molecule,
    get_molecule_by_smiles,
    create_molecule,
    update_molecule,
    delete_molecule,
    get_molecules,
    update_molecule_flag,
    get_property_ranges,
    get_available_properties,
    bulk_create_molecules,
    generate_molecule_image,
    check_molecule_drug_likeness
)
from src.backend.app.crud.crud_molecule import molecule
from src.backend.app.molecular.processor import process_molecule, process_molecules_batch
from src.backend.app.models.molecule import Molecule
from src.backend.app.models.molecule_property import MoleculeProperty
from src.backend.app.schemas.molecule import MoleculeCreate, MoleculeUpdate, MoleculeFilter
from src.backend.app.exceptions import MoleculeServiceException, ValidationException, ResourceNotFoundException

# Mock UUID for testing purposes
TEST_USER_ID = uuid.uuid4()

def test_get_molecule_success(Session):
    """Test successful retrieval of a molecule by ID"""
    # Create a test molecule with properties
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule.properties = [MoleculeProperty(property_name="MW", property_value=46.07)]
    Session.add(test_molecule)
    Session.commit()
    Session.refresh(test_molecule)

    # Call get_molecule with the molecule ID
    retrieved_molecule = get_molecule(test_molecule.id)

    # Assert that the returned molecule matches the created molecule
    assert retrieved_molecule is not None
    assert retrieved_molecule["id"] == test_molecule.id
    assert retrieved_molecule["smiles"] == test_molecule.smiles

    # Assert that properties are included in the result
    assert len(retrieved_molecule["properties"]) == 1
    assert retrieved_molecule["properties"][0]["property_name"] == "MW"
    assert retrieved_molecule["properties"][0]["property_value"] == 46.07

def test_get_molecule_not_found(Session):
    """Test molecule retrieval with non-existent ID"""
    # Call get_molecule with a non-existent ID
    retrieved_molecule = get_molecule(uuid.uuid4())

    # Assert that None is returned
    assert retrieved_molecule is None

def test_get_molecule_by_smiles_success(Session):
    """Test successful retrieval of a molecule by SMILES"""
    # Create a test molecule with a specific SMILES
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    Session.add(test_molecule)
    Session.commit()
    Session.refresh(test_molecule)

    # Call get_molecule_by_smiles with the SMILES string
    retrieved_molecule = get_molecule_by_smiles("CCO")

    # Assert that the returned molecule matches the created molecule
    assert retrieved_molecule is not None
    assert retrieved_molecule["id"] == test_molecule.id
    assert retrieved_molecule["smiles"] == test_molecule.smiles

def test_get_molecule_by_smiles_not_found(Session):
    """Test molecule retrieval with non-existent SMILES"""
    # Call get_molecule_by_smiles with a non-existent SMILES
    retrieved_molecule = get_molecule_by_smiles("NonExistentSMILES")

    # Assert that None is returned
    assert retrieved_molecule is None

@patch('src.backend.app.services.molecule_service.process_molecule')
def test_create_molecule_success(mock_process_molecule, Session):
    """Test successful creation of a molecule"""
    # Create MoleculeCreate object with valid SMILES and properties
    molecule_data = MoleculeCreate(smiles="CCO", properties=[{"property_name": "MW", "property_value": 46.07}])

    # Mock process_molecule to return processed data
    mock_process_molecule.return_value = {"validation": {"structure_valid": True}, "properties": {"MW": 46.07}}

    # Call create_molecule with the molecule data and user ID
    created_molecule = create_molecule(molecule_data, TEST_USER_ID)

    # Assert that the returned molecule has correct SMILES and properties
    assert created_molecule["smiles"] == "CCO"
    assert len(created_molecule["properties"]) == 1
    assert created_molecule["properties"][0]["property_name"] == "MW"
    assert created_molecule["properties"][0]["property_value"] == 46.07

    # Assert that molecule is stored in the database
    db_molecule = Session.query(Molecule).filter(Molecule.id == created_molecule["id"]).first()
    assert db_molecule is not None
    assert db_molecule.smiles == "CCO"

def test_create_molecule_duplicate(Session):
    """Test molecule creation failure with duplicate SMILES"""
    # Create a test molecule with a specific SMILES
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    Session.add(test_molecule)
    Session.commit()

    # Create MoleculeCreate object with the same SMILES
    molecule_data = MoleculeCreate(smiles="CCO", properties=[{"property_name": "MW", "property_value": 46.07}])

    # Call create_molecule with the molecule data
    with pytest.raises(MoleculeServiceException) as exc_info:
        create_molecule(molecule_data, TEST_USER_ID)

    # Assert that MoleculeServiceException is raised with 'already exists' message
    assert "already exists" in str(exc_info.value)

@patch('src.backend.app.services.molecule_service.process_molecule')
def test_create_molecule_invalid(mock_process_molecule, Session):
    """Test molecule creation failure with invalid SMILES"""
    # Create MoleculeCreate object with invalid SMILES
    molecule_data = MoleculeCreate(smiles="InvalidSMILES", properties=[{"property_name": "MW", "property_value": 46.07}])

    # Mock process_molecule to raise ValidationException
    mock_process_molecule.side_effect = ValidationException("Invalid SMILES")

    # Call create_molecule with the molecule data
    with pytest.raises(ValidationException) as exc_info:
        create_molecule(molecule_data, TEST_USER_ID)

    # Assert that ValidationException is raised
    assert "Invalid SMILES" in str(exc_info.value)

def test_update_molecule_success(Session):
    """Test successful update of a molecule"""
    # Create a test molecule
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID, flag_status="LOW")
    Session.add(test_molecule)
    Session.commit()
    Session.refresh(test_molecule)

    # Create MoleculeUpdate object with new flag_status and properties
    update_data = MoleculeUpdate(flag_status="HIGH", properties=[{"property_name": "LogP", "property_value": 0.5}])

    # Call update_molecule with molecule ID and update data
    updated_molecule = update_molecule(test_molecule.id, update_data)

    # Assert that the returned molecule has updated flag_status and properties
    assert updated_molecule["id"] == test_molecule.id
    assert updated_molecule["flag_status"] == "HIGH"
    assert len(updated_molecule["properties"]) == 1
    assert updated_molecule["properties"][0]["property_name"] == "LogP"
    assert updated_molecule["properties"][0]["property_value"] == 0.5

    # Assert that changes are stored in the database
    db_molecule = Session.query(Molecule).filter(Molecule.id == test_molecule.id).first()
    assert db_molecule.flag_status == "HIGH"
    assert len(db_molecule.properties) == 1
    assert db_molecule.properties[0].property_name == "LogP"
    assert db_molecule.properties[0].property_value == 0.5

def test_update_molecule_not_found(Session):
    """Test molecule update failure with non-existent ID"""
    # Create MoleculeUpdate object with flag_status and properties
    update_data = MoleculeUpdate(flag_status="HIGH", properties=[{"property_name": "LogP", "property_value": 0.5}])

    # Call update_molecule with non-existent ID and update data
    with pytest.raises(MoleculeServiceException) as exc_info:
        update_molecule(uuid.uuid4(), update_data)

    # Assert that MoleculeServiceException is raised with 'not found' message
    assert "not found" in str(exc_info.value)

def test_delete_molecule_success(Session):
    """Test successful deletion of a molecule"""
    # Create a test molecule
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    Session.add(test_molecule)
    Session.commit()
    Session.refresh(test_molecule)

    # Call delete_molecule with the molecule ID
    deleted = delete_molecule(test_molecule.id)

    # Assert that True is returned
    assert deleted is True

    # Verify that molecule is no longer in the database
    db_molecule = Session.query(Molecule).filter(Molecule.id == test_molecule.id).first()
    assert db_molecule is None

def test_delete_molecule_not_found(Session):
    """Test molecule deletion with non-existent ID"""
    # Call delete_molecule with a non-existent ID
    deleted = delete_molecule(uuid.uuid4())

    # Assert that False is returned
    assert deleted is False

def test_get_molecules_no_filters(Session):
    """Test retrieval of molecules without filters"""
    # Create multiple test molecules
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Call get_molecules without filters
    molecules_data = get_molecules()

    # Assert that all molecules are returned
    assert len(molecules_data["items"]) == 2
    assert molecules_data["items"][0]["smiles"] == "CCO"
    assert molecules_data["items"][1]["smiles"] == "CCC"

    # Assert that total count matches number of molecules
    assert molecules_data["total"] == 2

    # Assert that pagination info is correct
    assert molecules_data["skip"] == 0
    assert molecules_data["limit"] == 100

def test_get_molecules_with_smiles_filter(Session):
    """Test retrieval of molecules filtered by SMILES"""
    # Create multiple test molecules with different SMILES
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Create MoleculeFilter with specific SMILES
    molecule_filter = MoleculeFilter(smiles="CCO")

    # Call get_molecules with the filter
    molecules_data = get_molecules(filters=molecule_filter)

    # Assert that only matching molecules are returned
    assert len(molecules_data["items"]) == 1
    assert molecules_data["items"][0]["smiles"] == "CCO"

    # Assert that total count matches number of matching molecules
    assert molecules_data["total"] == 1

def test_get_molecules_with_user_filter(Session):
    """Test retrieval of molecules filtered by user"""
    # Create test molecules with different user IDs
    user_id1 = uuid.uuid4()
    user_id2 = uuid.uuid4()
    test_molecule1 = Molecule(smiles="CCO", created_by=user_id1)
    test_molecule2 = Molecule(smiles="CCC", created_by=user_id2)
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Create MoleculeFilter with specific user ID
    molecule_filter = MoleculeFilter(created_by=user_id1)

    # Call get_molecules with the filter
    molecules_data = get_molecules(filters=molecule_filter)

    # Assert that only molecules created by the user are returned
    assert len(molecules_data["items"]) == 1
    assert molecules_data["items"][0]["smiles"] == "CCO"

    # Assert that total count matches number of user's molecules
    assert molecules_data["total"] == 1

def test_get_molecules_with_flag_filter(Session):
    """Test retrieval of molecules filtered by flag status"""
    # Create test molecules with different flag statuses
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID, flag_status="HIGH")
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID, flag_status="LOW")
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Create MoleculeFilter with specific flag status
    molecule_filter = MoleculeFilter(flag_status="HIGH")

    # Call get_molecules with the filter
    molecules_data = get_molecules(filters=molecule_filter)

    # Assert that only molecules with matching flag status are returned
    assert len(molecules_data["items"]) == 1
    assert molecules_data["items"][0]["smiles"] == "CCO"

    # Assert that total count matches number of flagged molecules
    assert molecules_data["total"] == 1

def test_get_molecules_with_property_filter(Session):
    """Test retrieval of molecules filtered by property values"""
    # Create test molecules with different property values
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule1.properties = [MoleculeProperty(property_name="MW", property_value=46.07)]
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    test_molecule2.properties = [MoleculeProperty(property_name="MW", property_value=60.0)]
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Create MoleculeFilter with property range filter
    molecule_filter = MoleculeFilter(properties=[{"name": "MW", "min_value": 40, "max_value": 50}])

    # Call get_molecules with the filter
    molecules_data = get_molecules(filters=molecule_filter)

    # Assert that only molecules with properties in range are returned
    assert len(molecules_data["items"]) == 1
    assert molecules_data["items"][0]["smiles"] == "CCO"

    # Assert that total count matches number of matching molecules
    assert molecules_data["total"] == 1

def test_get_molecules_with_sorting(Session):
    """Test retrieval of molecules with sorting"""
    # Create test molecules with different property values
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule1.properties = [MoleculeProperty(property_name="MW", property_value=60.0)]
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    test_molecule2.properties = [MoleculeProperty(property_name="MW", property_value=46.07)]
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Create MoleculeFilter with sort_by and sort_desc parameters
    molecule_filter_asc = MoleculeFilter(sort_by="MW", sort_desc=False, properties=[{"name": "MW"}])
    molecule_filter_desc = MoleculeFilter(sort_by="MW", sort_desc=True, properties=[{"name": "MW"}])

    # Call get_molecules with the filter
    molecules_data_asc = get_molecules(filters=molecule_filter_asc)
    molecules_data_desc = get_molecules(filters=molecule_filter_desc)

    # Assert that molecules are returned in correct order
    assert molecules_data_asc["items"][0]["smiles"] == "CCC"
    assert molecules_data_desc["items"][0]["smiles"] == "CCO"

    # Test both ascending and descending order
    assert molecules_data_asc["items"][1]["smiles"] == "CCO"
    assert molecules_data_desc["items"][1]["smiles"] == "CCC"

def test_get_molecules_with_pagination(Session):
    """Test retrieval of molecules with pagination"""
    # Create multiple test molecules (more than default limit)
    for i in range(150):
        test_molecule = Molecule(smiles=f"C{i}", created_by=TEST_USER_ID)
        Session.add(test_molecule)
    Session.commit()

    # Call get_molecules with skip and limit parameters
    molecules_data = get_molecules(skip=50, limit=25)

    # Assert that correct subset of molecules is returned
    assert len(molecules_data["items"]) == 25
    assert molecules_data["items"][0]["smiles"] == "C50"
    assert molecules_data["items"][-1]["smiles"] == "C74"

    # Assert that total count includes all molecules
    assert molecules_data["total"] == 150

    # Test multiple pages to ensure correct pagination
    molecules_data_page2 = get_molecules(skip=75, limit=25)
    assert molecules_data_page2["items"][0]["smiles"] == "C75"

def test_update_molecule_flag_success(Session):
    """Test successful update of molecule flag status"""
    # Create a test molecule with initial flag status
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID, flag_status="LOW")
    Session.add(test_molecule)
    Session.commit()
    Session.refresh(test_molecule)

    # Call update_molecule_flag with molecule ID and new flag status
    updated_molecule = update_molecule_flag(test_molecule.id, "HIGH")

    # Assert that the returned molecule has updated flag status
    assert updated_molecule["id"] == test_molecule.id
    assert updated_molecule["flag_status"] == "HIGH"

    # Assert that change is stored in the database
    db_molecule = Session.query(Molecule).filter(Molecule.id == test_molecule.id).first()
    assert db_molecule.flag_status == "HIGH"

def test_update_molecule_flag_not_found(Session):
    """Test flag update failure with non-existent molecule"""
    # Call update_molecule_flag with non-existent ID and flag status
    with pytest.raises(MoleculeServiceException) as exc_info:
        update_molecule_flag(uuid.uuid4(), "HIGH")

    # Assert that MoleculeServiceException is raised with 'not found' message
    assert "not found" in str(exc_info.value)

def test_get_property_ranges_success(Session):
    """Test successful retrieval of property ranges"""
    # Create test molecules with varying property values
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule1.properties = [MoleculeProperty(property_name="MW", property_value=46.07)]
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    test_molecule2.properties = [MoleculeProperty(property_name="MW", property_value=60.0)]
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Call get_property_ranges with property name
    property_range = get_property_ranges("MW")

    # Assert that returned min/max values match expected range
    assert property_range["min"] == 46.07
    assert property_range["max"] == 60.0

    # Test with multiple different properties
    property_range = get_property_ranges("NonExistentProperty")
    assert property_range["min"] == 0.0
    assert property_range["max"] == 0.0

def test_get_property_ranges_no_data(Session):
    """Test property range retrieval with no data"""
    # Call get_property_ranges with property name that doesn't exist
    property_range = get_property_ranges("NonExistentProperty")

    # Assert that default min/max values are returned (0, 0)
    assert property_range["min"] == 0.0
    assert property_range["max"] == 0.0

def test_get_available_properties_success(Session):
    """Test successful retrieval of available properties"""
    # Create test molecules with various properties
    test_molecule1 = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    test_molecule1.properties = [MoleculeProperty(property_name="MW", property_value=46.07),
                                 MoleculeProperty(property_name="LogP", property_value=0.5)]
    test_molecule2 = Molecule(smiles="CCC", created_by=TEST_USER_ID)
    test_molecule2.properties = [MoleculeProperty(property_name="MW", property_value=60.0),
                                 MoleculeProperty(property_name="TPSA", property_value=20.0)]
    Session.add_all([test_molecule1, test_molecule2])
    Session.commit()

    # Call get_available_properties
    properties = get_available_properties()

    # Assert that returned list contains all expected property names
    assert "MW" in properties
    assert "LogP" in properties
    assert "TPSA" in properties

    # Assert that each property name appears only once (distinct)
    assert len(properties) == 3

def test_get_available_properties_no_data(Session):
    """Test property retrieval with no data"""
    # Call get_available_properties with empty database
    properties = get_available_properties()

    # Assert that empty list is returned
    assert properties == []

@patch('src.backend.app.services.molecule_service.process_molecules_batch')
def test_bulk_create_molecules_success(mock_process_molecules_batch, Session):
    """Test successful bulk creation of molecules"""
    # Create list of molecule data dictionaries with SMILES and properties
    molecules_data = [
        {"smiles": "CCO", "properties": [{"property_name": "MW", "property_value": 46.07}]},
        {"smiles": "CCC", "properties": [{"property_name": "MW", "property_value": 60.0}]}
    ]

    # Mock process_molecules_stream to return processed data
    mock_process_molecules_batch.return_value = [
        {"smiles": "CCO", "validation": {"structure_valid": True}, "properties": {"MW": 46.07}},
        {"smiles": "CCC", "validation": {"structure_valid": True}, "properties": {"MW": 60.0}}
    ]

    # Call bulk_create_molecules with the data and user ID
    result = bulk_create_molecules(molecules_data, TEST_USER_ID)

    # Assert that correct number of molecules are created
    assert result["total_created"] == 2

    # Assert that returned statistics match expected values
    assert result["total_duplicates"] == 0
    assert result["total_processed"] == 2
    assert result["total_failed"] == 0

    # Verify that molecules are stored in the database
    db_molecules = Session.query(Molecule).all()
    assert len(db_molecules) == 2
    assert db_molecules[0].smiles == "CCO"
    assert db_molecules[1].smiles == "CCC"

def test_bulk_create_molecules_with_duplicates(Session):
    """Test bulk creation with duplicate molecules"""
    # Create a test molecule in the database
    test_molecule = Molecule(smiles="CCO", created_by=TEST_USER_ID)
    Session.add(test_molecule)
    Session.commit()

    # Create list of molecule data including the existing SMILES
    molecules_data = [
        {"smiles": "CCO", "properties": [{"property_name": "MW", "property_value": 46.07}]},
        {"smiles": "CCC", "properties": [{"property_name": "MW", "property_value": 60.0}]}
    ]

    # Call bulk_create_molecules with the data
    result = bulk_create_molecules(molecules_data, TEST_USER_ID)

    # Assert that only new molecules are created
    assert result["total_created"] == 1

    # Assert that statistics show correct count of duplicates
    assert result["total_duplicates"] == 1
    assert result["total_processed"] == 2
    assert result["total_failed"] == 0

    # Verify that existing molecule is unchanged
    db_molecule = Session.query(Molecule).filter(Molecule.smiles == "CCO").first()
    assert db_molecule is not None
    assert db_molecule.created_by == TEST_USER_ID

@patch('src.backend.app.services.molecule_service.process_molecules_batch')
def test_bulk_create_molecules_without_processing(mock_process_molecules_batch, Session):
    """Test bulk creation without property processing"""
    # Create list of molecule data dictionaries
    molecules_data = [
        {"smiles": "CCO", "properties": [{"property_name": "MW", "property_value": 46.07}]},
        {"smiles": "CCC", "properties": [{"property_name": "MW", "property_value": 60.0}]}
    ]

    # Call bulk_create_molecules with process_properties=False
    result = bulk_create_molecules(molecules_data, TEST_USER_ID, process_properties=False)

    # Assert that process_molecules_stream is not called
    mock_process_molecules_batch.assert_not_called()

    # Assert that molecules are created with original data
    assert result["total_created"] == 2
    assert result["total_duplicates"] == 0
    assert result["total_processed"] == 2
    assert result["total_failed"] == 0

    # Verify that molecules are stored in the database
    db_molecules = Session.query(Molecule).all()
    assert len(db_molecules) == 2
    assert db_molecules[0].smiles == "CCO"
    assert db_molecules[1].smiles == "CCC"

@patch('src.backend.app.services.molecule_service.Chem.MolFromSmiles')
@patch('src.backend.app.services.molecule_service.Draw.MolToImage')
@patch('src.backend.app.services.molecule_service.upload_molecule_image')
@patch('src.backend.app.services.molecule_service.get_molecule_image_url')
def test_generate_molecule_image_success(mock_get_molecule_image_url, mock_upload_molecule_image, mock_mol_to_image, mock_mol_from_smiles):
    """Test successful generation of molecule image"""
    # Mock Chem.MolFromSmiles to return a valid molecule
    mock_mol = MagicMock()
    mock_mol_from_smiles.return_value = mock_mol

    # Mock Draw.MolToImage to return a mock image
    mock_image = MagicMock()
    mock_mol_to_image.return_value = mock_image

    # Mock upload_molecule_image to return a file path
    mock_upload_molecule_image.return_value = "test_image.png"

    # Mock get_molecule_image_url to return a URL
    mock_get_molecule_image_url.return_value = "http://example.com/test_image.png"

    # Call generate_molecule_image with valid SMILES
    image_data = generate_molecule_image("CCO")

    # Assert that image URL is returned
    assert image_data["image_url"] == "http://example.com/test_image.png"

    # Verify that all mocked functions were called with correct arguments
    mock_mol_from_smiles.assert_called_once_with("CCO")
    mock_mol_to_image.assert_called_once()
    mock_upload_molecule_image.assert_called_once()
    mock_get_molecule_image_url.assert_called_once_with("test_image.png")

@patch('src.backend.app.services.molecule_service.Chem.MolFromSmiles')
def test_generate_molecule_image_invalid_smiles(mock_mol_from_smiles):
    """Test image generation failure with invalid SMILES"""
    # Mock Chem.MolFromSmiles to return None (invalid molecule)
    mock_mol_from_smiles.return_value = None

    # Call generate_molecule_image with invalid SMILES
    with pytest.raises(MoleculeServiceException) as exc_info:
        generate_molecule_image("InvalidSMILES")

    # Assert that MoleculeServiceException is raised with 'Invalid SMILES' message
    assert "Invalid SMILES" in str(exc_info.value)

@patch('src.backend.app.services.molecule_service.Chem.MolFromSmiles')
@patch('src.backend.app.services.molecule_service.Draw.MolToImage')
@patch('src.backend.app.services.molecule_service.upload_molecule_image')
def test_generate_molecule_image_without_storage(mock_upload_molecule_image, mock_mol_to_image, mock_mol_from_smiles):
    """Test image generation without storage"""
    # Mock Chem.MolFromSmiles to return a valid molecule
    mock_mol = MagicMock()
    mock_mol_from_smiles.return_value = mock_mol

    # Mock Draw.MolToImage to return a mock image
    mock_image = MagicMock()
    mock_mol_to_image.return_value = mock_image

    # Call generate_molecule_image with store_image=False
    image_data = generate_molecule_image("CCO", store_image=False)

    # Assert that image data is returned instead of URL
    assert "image_data" in image_data

    # Verify that upload_molecule_image is not called
    mock_upload_molecule_image.assert_not_called()

@patch('src.backend.app.services.molecule_service.process_molecule')
@patch('src.backend.app.services.molecule_service.check_drug_likeness')
def test_check_molecule_drug_likeness_success(mock_check_drug_likeness, mock_process_molecule):
    """Test successful drug-likeness check"""
    # Mock process_molecule to return processed data
    mock_process_molecule.return_value = {"validation": {"structure_valid": True}}

    # Mock check_drug_likeness to return drug-likeness assessment
    mock_check_drug_likeness.return_value = {"is_drug_like": True}

    # Call check_molecule_drug_likeness with valid SMILES
    drug_likeness = check_molecule_drug_likeness("CCO")

    # Assert that drug-likeness assessment is returned
    assert drug_likeness["is_drug_like"] is True

    # Verify that mocked functions were called with correct arguments
    mock_process_molecule.assert_called_once_with("CCO")
    mock_check_drug_likeness.assert_called_once()

@patch('src.backend.app.services.molecule_service.process_molecule')
def test_check_molecule_drug_likeness_invalid_smiles(mock_process_molecule):
    """Test drug-likeness check failure with invalid SMILES"""
    # Mock process_molecule to raise ValidationException
    mock_process_molecule.side_effect = ValidationException("Invalid SMILES")

    # Call check_molecule_drug_likeness with invalid SMILES
    with pytest.raises(ValidationException) as exc_info:
        check_molecule_drug_likeness("InvalidSMILES")

    # Assert that ValidationException is raised
    assert "Invalid SMILES" in str(exc_info.value)