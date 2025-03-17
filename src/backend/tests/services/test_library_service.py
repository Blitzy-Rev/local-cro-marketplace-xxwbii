import pytest
import uuid
from unittest.mock import MagicMock, patch
from datetime import datetime

from ...app.services.library_service import (
    get_library,
    get_library_by_name,
    create_library,
    update_library,
    delete_library,
    get_libraries,
    get_user_libraries,
    get_library_with_molecules,
    add_molecules_to_library,
    remove_molecules_from_library,
    check_molecule_in_library,
    get_library_molecules
)
from ...app.crud.crud_library import library
from ...app.crud.crud_molecule import molecule
from ...app.models.library import Library
from ...app.schemas.library import LibraryCreate, LibraryUpdate, LibraryFilter
from ...app.exceptions import LibraryServiceException, ResourceNotFoundException


def test_get_library_success(db_session):
    """Test successful retrieval of a library by ID."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Mock the library.get function to return our test library
    with patch.object(library, 'get', return_value=test_library):
        # Call the function being tested
        result = get_library(library_id)
        
        # Assert that the returned library matches the created library
        assert result is not None
        assert result["id"] == test_library.id
        assert result["name"] == test_library.name
        assert result["description"] == test_library.description
        assert result["created_by"] == test_library.created_by
        assert "created_at" in result
        assert "updated_at" in result


def test_get_library_not_found(db_session):
    """Test library retrieval with non-existent ID."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Mock the library.get function to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested
        result = get_library(non_existent_id)
        
        # Assert that None is returned
        assert result is None


def test_get_library_by_name_success(db_session):
    """Test successful retrieval of a library by name and user ID."""
    # Create a test library with a specific name and user ID
    user_id = uuid.uuid4()
    library_name = "Test Library"
    test_library = Library(
        id=uuid.uuid4(),
        name=library_name,
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Mock the library.get_by_name function to return our test library
    with patch.object(library, 'get_by_name', return_value=test_library):
        # Call the function being tested
        result = get_library_by_name(library_name, user_id)
        
        # Assert that the returned library matches the created library
        assert result is not None
        assert result["id"] == test_library.id
        assert result["name"] == library_name
        assert result["description"] == test_library.description
        assert result["created_by"] == user_id
        assert "created_at" in result
        assert "updated_at" in result


def test_get_library_by_name_not_found(db_session):
    """Test library retrieval with non-existent name."""
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    non_existent_name = "Non-existent Library"
    
    # Mock the library.get_by_name function to return None (not found)
    with patch.object(library, 'get_by_name', return_value=None):
        # Call the function being tested
        result = get_library_by_name(non_existent_name, user_id)
        
        # Assert that None is returned
        assert result is None


def test_create_library_success(db_session):
    """Test successful creation of a library."""
    # Create LibraryCreate object with name and description
    library_data = LibraryCreate(
        name="New Library",
        description="New Library Description"
    )
    
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    
    # Create a mock library that will be returned by library.create
    created_library = Library(
        id=uuid.uuid4(),
        name=library_data.name,
        description=library_data.description,
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Mock library.get_by_name to return None (no existing library) and library.create to return our mock
    with patch.object(library, 'get_by_name', return_value=None), \
         patch.object(library, 'create', return_value=created_library):
        
        # Call the function being tested
        result = create_library(library_data, user_id)
        
        # Assert that the returned library has the correct name and description
        assert result is not None
        assert result["name"] == library_data.name
        assert result["description"] == library_data.description
        assert result["created_by"] == user_id
        assert "created_at" in result
        assert "updated_at" in result
        
        # Assert that library was created with correct data
        library.create.assert_called_once()


def test_create_library_duplicate_name(db_session):
    """Test library creation failure with duplicate name for same user."""
    # Create a test library with a specific name and user ID
    user_id = uuid.uuid4()
    existing_library = Library(
        id=uuid.uuid4(),
        name="Duplicate Library",
        description="Existing Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create LibraryCreate object with the same name
    library_data = LibraryCreate(
        name="Duplicate Library",
        description="New Description"
    )
    
    # Mock library.get_by_name to return the existing library
    with patch.object(library, 'get_by_name', return_value=existing_library):
        # Call the function being tested and assert that LibraryServiceException is raised
        with pytest.raises(LibraryServiceException) as excinfo:
            create_library(library_data, user_id)
        
        # Assert that the exception message contains "already exists"
        assert "already exists" in str(excinfo.value)


def test_create_library_same_name_different_user(db_session):
    """Test library creation with same name but different user."""
    # Create a test library with a specific name and user ID
    user1_id = uuid.uuid4()
    existing_library = Library(
        id=uuid.uuid4(),
        name="Shared Library Name",
        description="First User's Library",
        created_by=user1_id,
        created_at=datetime.utcnow()
    )
    
    # Create LibraryCreate object with the same name
    library_data = LibraryCreate(
        name="Shared Library Name",
        description="Second User's Library"
    )
    
    # Generate a different user ID
    user2_id = uuid.uuid4()
    
    # Create a mock library that will be returned by library.create
    created_library = Library(
        id=uuid.uuid4(),
        name=library_data.name,
        description=library_data.description,
        created_by=user2_id,
        created_at=datetime.utcnow()
    )
    
    # Mock library.get_by_name to return None (no existing library for second user)
    # and library.create to return our created_library
    with patch.object(library, 'get_by_name', return_value=None), \
         patch.object(library, 'create', return_value=created_library):
        
        # Call the function being tested
        result = create_library(library_data, user2_id)
        
        # Assert that the library is created successfully
        assert result is not None
        assert result["name"] == library_data.name
        assert result["description"] == library_data.description
        assert result["created_by"] == user2_id
        
        # Assert that library was created with correct data
        library.create.assert_called_once()


def test_update_library_success(db_session):
    """Test successful update of a library."""
    # Create a test library in the database
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    existing_library = Library(
        id=library_id,
        name="Original Library",
        description="Original Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create LibraryUpdate object with new name and description
    update_data = LibraryUpdate(
        name="Updated Library",
        description="Updated Description"
    )
    
    # Create the updated library that will be returned by library.update
    updated_library = Library(
        id=library_id,
        name=update_data.name,
        description=update_data.description,
        created_by=user_id,
        created_at=existing_library.created_at,
        updated_at=datetime.utcnow()
    )
    
    # Mock library.get to return the existing library
    # Mock library.get_by_name to return None (no duplicate name conflict)
    # Mock library.update to return the updated library
    with patch.object(library, 'get', return_value=existing_library), \
         patch.object(library, 'get_by_name', return_value=None), \
         patch.object(library, 'update', return_value=updated_library):
        
        # Call the function being tested
        result = update_library(library_id, update_data, user_id)
        
        # Assert that the returned library has updated name and description
        assert result["id"] == library_id
        assert result["name"] == update_data.name
        assert result["description"] == update_data.description
        assert result["created_by"] == user_id
        assert "updated_at" in result
        
        # Assert that library.update was called with the right parameters
        library.update.assert_called_once()


def test_update_library_not_found(db_session):
    """Test library update failure with non-existent ID."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Create LibraryUpdate object with name and description
    update_data = LibraryUpdate(
        name="Updated Library",
        description="Updated Description"
    )
    
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    
    # Mock library.get to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            update_library(non_existent_id, update_data, user_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)


def test_update_library_unauthorized(db_session):
    """Test library update failure with unauthorized user."""
    # Create a test library with a specific user ID
    library_id = uuid.uuid4()
    creator_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=creator_id,
        created_at=datetime.utcnow()
    )
    
    # Create LibraryUpdate object with new name and description
    update_data = LibraryUpdate(
        name="Updated Library",
        description="Updated Description"
    )
    
    # Generate a different user ID
    different_user_id = uuid.uuid4()
    
    # Mock library.get to return the test library
    with patch.object(library, 'get', return_value=test_library):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            update_library(library_id, update_data, different_user_id)
        
        # Assert that the exception message contains "unauthorized"
        assert "unauthorized" in str(excinfo.value)


def test_update_library_duplicate_name(db_session):
    """Test library update failure with duplicate name."""
    # Create two test libraries with different names but same user ID
    user_id = uuid.uuid4()
    
    library1_id = uuid.uuid4()
    library1 = Library(
        id=library1_id,
        name="First Library",
        description="First Library Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    library2 = Library(
        id=uuid.uuid4(),
        name="Second Library",
        description="Second Library Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create LibraryUpdate object with the name of the second library
    update_data = LibraryUpdate(
        name="Second Library",
        description="Updated Description"
    )
    
    # Mock library.get to return the first library
    # Mock library.get_by_name to return the second library (duplicate name)
    with patch.object(library, 'get', return_value=library1), \
         patch.object(library, 'get_by_name', return_value=library2):
        
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            update_library(library1_id, update_data, user_id)
        
        # Assert that the exception message contains "already exists"
        assert "already exists" in str(excinfo.value)


def test_delete_library_success(db_session):
    """Test successful deletion of a library."""
    # Create a test library in the database
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Mock library.get to return the test library
    # Mock library.remove to succeed
    with patch.object(library, 'get', return_value=test_library), \
         patch.object(library, 'remove', return_value=test_library):
        
        # Call the function being tested
        result = delete_library(library_id, user_id)
        
        # Assert that True is returned
        assert result is True
        
        # Verify that library is no longer in the database
        library.remove.assert_called_once_with(db_session, library_id)


def test_delete_library_not_found(db_session):
    """Test library deletion with non-existent ID."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    
    # Mock library.get to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            delete_library(non_existent_id, user_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)


def test_delete_library_unauthorized(db_session):
    """Test library deletion failure with unauthorized user."""
    # Create a test library with a specific user ID
    library_id = uuid.uuid4()
    creator_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=creator_id,
        created_at=datetime.utcnow()
    )
    
    # Generate a different user ID
    different_user_id = uuid.uuid4()
    
    # Mock library.get to return the test library
    with patch.object(library, 'get', return_value=test_library):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            delete_library(library_id, different_user_id)
        
        # Assert that the exception message contains "unauthorized"
        assert "unauthorized" in str(excinfo.value)
        
        # Verify that library still exists in the database
        library.remove.assert_not_called()


def test_get_libraries_no_filters(db_session):
    """Test retrieval of libraries without filters."""
    # Create multiple test libraries
    libraries = [
        Library(
            id=uuid.uuid4(),
            name=f"Library {i}",
            description=f"Description {i}",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ) for i in range(3)
    ]
    
    total_count = len(libraries)
    
    # Mock library.filter_libraries to return our libraries and total count
    # Mock library.count_molecules for each library
    with patch.object(library, 'filter_libraries', return_value=(libraries, total_count)), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries()
        
        # Assert that all libraries are returned
        assert len(result["items"]) == len(libraries)
        assert result["total"] == total_count
        
        # Assert that pagination info is correct
        assert result["skip"] == 0
        assert result["limit"] == 100
        
        # Check each library has expected fields
        for i, lib in enumerate(result["items"]):
            assert lib["id"] == libraries[i].id
            assert lib["name"] == libraries[i].name
            assert lib["description"] == libraries[i].description
            assert lib["molecule_count"] == 5


def test_get_libraries_with_name_filter(db_session):
    """Test retrieval of libraries filtered by name."""
    # Create multiple test libraries with different names
    libraries = [
        Library(
            id=uuid.uuid4(),
            name="Test Library",
            description="Description 1",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ),
        Library(
            id=uuid.uuid4(),
            name="Another Library",
            description="Description 2",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        )
    ]
    
    # Create a filter with name parameter
    name_filter = "Test"
    filter_obj = LibraryFilter(name=name_filter)
    
    # Only the first library should match the filter
    matching_libraries = [libraries[0]]
    total_matching = len(matching_libraries)
    
    # Mock library.filter_libraries to return our matching libraries
    # Mock library.count_molecules for each library
    with patch.object(library, 'filter_libraries', return_value=(matching_libraries, total_matching)), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries(filter_obj)
        
        # Assert that only libraries with matching name are returned
        assert len(result["items"]) == len(matching_libraries)
        assert result["total"] == total_matching
        
        # Check that the returned library matches the filter
        assert name_filter in result["items"][0]["name"]


def test_get_libraries_with_user_filter(db_session):
    """Test retrieval of libraries filtered by user."""
    # Create test libraries with different user IDs
    user1_id = uuid.uuid4()
    user2_id = uuid.uuid4()
    
    user1_libraries = [
        Library(
            id=uuid.uuid4(),
            name=f"User1 Library {i}",
            description=f"Description {i}",
            created_by=user1_id,
            created_at=datetime.utcnow()
        ) for i in range(2)
    ]
    
    user2_libraries = [
        Library(
            id=uuid.uuid4(),
            name=f"User2 Library {i}",
            description=f"Description {i}",
            created_by=user2_id,
            created_at=datetime.utcnow()
        ) for i in range(1)
    ]
    
    # Create a filter with created_by parameter
    filter_obj = LibraryFilter(created_by=user1_id)
    
    # Mock library.filter_libraries to return user1's libraries
    # Mock library.count_molecules for each library
    with patch.object(library, 'filter_libraries', return_value=(user1_libraries, len(user1_libraries))), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries(filter_obj)
        
        # Assert that only libraries created by the user are returned
        assert len(result["items"]) == len(user1_libraries)
        assert result["total"] == len(user1_libraries)
        
        # Check that each library has the correct user ID
        for lib in result["items"]:
            assert lib["created_by"] == user1_id


def test_get_libraries_with_molecule_filter(db_session):
    """Test retrieval of libraries filtered by contained molecule."""
    # Create test libraries and molecules
    molecule_id = uuid.uuid4()
    
    libraries_with_molecule = [
        Library(
            id=uuid.uuid4(),
            name=f"Molecule Library {i}",
            description=f"Contains the molecule",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ) for i in range(2)
    ]
    
    # Create a filter with contains_molecule_id parameter
    filter_obj = LibraryFilter(contains_molecule_id=molecule_id)
    
    # Mock molecule.get to return a valid molecule
    # Mock library.filter_libraries to return libraries containing the molecule
    with patch.object(molecule, 'get', return_value=MagicMock(id=molecule_id)), \
         patch.object(library, 'filter_libraries', return_value=(libraries_with_molecule, len(libraries_with_molecule))), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries(filter_obj)
        
        # Assert that only libraries containing the molecule are returned
        assert len(result["items"]) == len(libraries_with_molecule)
        assert result["total"] == len(libraries_with_molecule)


def test_get_libraries_with_sorting(db_session):
    """Test retrieval of libraries with sorting."""
    # Create test libraries with different names and creation dates
    libraries = [
        Library(
            id=uuid.uuid4(),
            name="A Library",
            description="Description 1",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ),
        Library(
            id=uuid.uuid4(),
            name="B Library",
            description="Description 2",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ),
        Library(
            id=uuid.uuid4(),
            name="C Library",
            description="Description 3",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        )
    ]
    
    # Create a filter with sort_by and sort_desc parameters
    filter_obj = LibraryFilter(sort_by="name", sort_desc=True)
    
    # Libraries should be returned in reverse alphabetical order
    sorted_libraries = sorted(libraries, key=lambda x: x.name, reverse=True)
    
    # Mock library.filter_libraries to return sorted libraries
    with patch.object(library, 'filter_libraries', return_value=(sorted_libraries, len(sorted_libraries))), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries(filter_obj)
        
        # Assert that libraries are returned in correct order
        assert len(result["items"]) == len(sorted_libraries)
        
        # Check order matches expected reverse alphabetical
        for i in range(len(result["items"]) - 1):
            assert result["items"][i]["name"] > result["items"][i+1]["name"]


def test_get_libraries_with_pagination(db_session):
    """Test retrieval of libraries with pagination."""
    # Create multiple test libraries (more than default limit)
    all_libraries = [
        Library(
            id=uuid.uuid4(),
            name=f"Library {i}",
            description=f"Description {i}",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow()
        ) for i in range(10)
    ]
    
    # Set pagination parameters
    skip = 2
    limit = 3
    
    # Expected subset of libraries
    expected_subset = all_libraries[skip:skip+limit]
    
    # Mock library.filter_libraries to return paginated subset but total count of all
    with patch.object(library, 'filter_libraries', return_value=(expected_subset, len(all_libraries))), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_libraries(skip=skip, limit=limit)
        
        # Assert that correct subset of libraries is returned
        assert len(result["items"]) == len(expected_subset)
        assert result["total"] == len(all_libraries)
        assert result["skip"] == skip
        assert result["limit"] == limit


def test_get_user_libraries(db_session):
    """Test retrieval of libraries for a specific user."""
    # Create test libraries with different user IDs
    user_id = uuid.uuid4()
    
    user_libraries = [
        Library(
            id=uuid.uuid4(),
            name=f"User Library {i}",
            description=f"Description {i}",
            created_by=user_id,
            created_at=datetime.utcnow()
        ) for i in range(3)
    ]
    
    total_user_libraries = len(user_libraries)
    
    # Mock library.get_by_user to return our user's libraries
    # Mock library.count_by_user to return the total count
    # Mock library.count_molecules for molecule counts
    with patch.object(library, 'get_by_user', return_value=user_libraries), \
         patch.object(library, 'count_by_user', return_value=total_user_libraries), \
         patch.object(library, 'count_molecules', return_value=5):
        
        # Call the function being tested
        result = get_user_libraries(user_id)
        
        # Assert that only libraries created by the user are returned
        assert len(result["items"]) == total_user_libraries
        assert result["total"] == total_user_libraries
        
        # Check that each library has the correct user ID
        for lib in result["items"]:
            assert lib["created_by"] == user_id


def test_get_library_with_molecules(db_session):
    """Test retrieval of a library with its molecules."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules
    test_molecules = [
        MagicMock(
            id=uuid.uuid4(),
            smiles=f"C{i}O",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow(),
            updated_at=None,
            flag_status=None
        ) for i in range(3)
    ]
    
    total_molecules = len(test_molecules)
    
    # Mock library.get_with_molecules to return our library and molecules
    # Mock library.count_molecules to return the total count
    with patch.object(library, 'get_with_molecules', return_value=(test_library, test_molecules)), \
         patch.object(library, 'count_molecules', return_value=total_molecules):
        
        # Call the function being tested
        result = get_library_with_molecules(library_id)
        
        # Assert that library details are correct
        assert result["library"]["id"] == test_library.id
        assert result["library"]["name"] == test_library.name
        assert result["library"]["description"] == test_library.description
        
        # Assert that all molecules in the library are returned
        assert len(result["molecules"]) == total_molecules
        assert result["total_molecules"] == total_molecules
        
        # Check each molecule has expected data
        for i, mol in enumerate(result["molecules"]):
            assert mol["id"] == test_molecules[i].id
            assert mol["smiles"] == test_molecules[i].smiles


def test_get_library_with_molecules_not_found(db_session):
    """Test retrieval of non-existent library with molecules."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Mock library.get_with_molecules to return None (not found)
    with patch.object(library, 'get_with_molecules', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            get_library_with_molecules(non_existent_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)


def test_add_molecules_to_library_success(db_session):
    """Test successful addition of molecules to a library."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Success response from add_molecules
    success_response = {
        "success": len(molecule_ids),
        "failures": []
    }
    
    # Mock library.get to return our test library
    # Mock molecule.get to return a valid molecule (for each molecule ID)
    # Mock library.add_molecules to return success
    with patch.object(library, 'get', return_value=test_library), \
         patch.object(molecule, 'get', return_value=MagicMock()), \
         patch.object(library, 'add_molecules', return_value=success_response):
        
        # Call the function being tested
        result = add_molecules_to_library(library_id, molecule_ids, user_id)
        
        # Assert that success count matches number of molecules
        assert result["success"] == len(molecule_ids)
        assert len(result["failures"]) == 0
        
        # Assert that molecules were associated with the library
        library.add_molecules.assert_called_once_with(db_session, library_id, molecule_ids)


def test_add_molecules_to_library_not_found(db_session):
    """Test molecule addition to non-existent library."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    
    # Mock library.get to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            add_molecules_to_library(non_existent_id, molecule_ids, user_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)


def test_add_molecules_to_library_unauthorized(db_session):
    """Test molecule addition failure with unauthorized user."""
    # Create a test library with a specific user ID
    library_id = uuid.uuid4()
    creator_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=creator_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Generate a different user ID
    different_user_id = uuid.uuid4()
    
    # Mock library.get to return our test library
    with patch.object(library, 'get', return_value=test_library):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            add_molecules_to_library(library_id, molecule_ids, different_user_id)
        
        # Assert that the exception message contains "unauthorized"
        assert "unauthorized" in str(excinfo.value)


def test_add_molecules_to_library_invalid_molecules(db_session):
    """Test addition of non-existent molecules to library."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Generate random UUIDs for non-existent molecules
    invalid_molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Mock library.get to return our test library
    # Mock molecule.get to return None (molecules don't exist)
    with patch.object(library, 'get', return_value=test_library), \
         patch.object(molecule, 'get', return_value=None):
        
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            add_molecules_to_library(library_id, invalid_molecule_ids, user_id)
        
        # Assert that the exception message contains information about the invalid molecules
        assert "not exist" in str(excinfo.value)


def test_add_molecules_to_library_already_in_library(db_session):
    """Test addition of molecules already in library."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Response indicating molecules already in library
    response = {
        "success": 0,
        "failures": [
            {"molecule_id": str(molecule_id), "error": "Molecule already in library"}
            for molecule_id in molecule_ids
        ]
    }
    
    # Mock library.get to return our test library
    # Mock molecule.get to return valid molecules
    # Mock library.add_molecules to return the "already in library" response
    with patch.object(library, 'get', return_value=test_library), \
         patch.object(molecule, 'get', return_value=MagicMock()), \
         patch.object(library, 'add_molecules', return_value=response):
        
        # Call the function being tested
        result = add_molecules_to_library(library_id, molecule_ids, user_id)
        
        # Assert that result indicates molecules were already in library
        assert result["success"] == 0
        assert len(result["failures"]) == len(molecule_ids)
        for failure in result["failures"]:
            assert "already in library" in failure["error"]


def test_remove_molecules_from_library_success(db_session):
    """Test successful removal of molecules from a library."""
    # Create a test library
    library_id = uuid.uuid4()
    user_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=user_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Success response from remove_molecules
    success_response = {
        "success": len(molecule_ids),
        "failures": []
    }
    
    # Mock library.get to return our test library
    # Mock library.remove_molecules to return success
    with patch.object(library, 'get', return_value=test_library), \
         patch.object(library, 'remove_molecules', return_value=success_response):
        
        # Call the function being tested
        result = remove_molecules_from_library(library_id, molecule_ids, user_id)
        
        # Assert that success count matches number of molecules
        assert result["success"] == len(molecule_ids)
        assert len(result["failures"]) == 0
        
        # Assert that molecules are no longer associated with the library
        library.remove_molecules.assert_called_once_with(db_session, library_id, molecule_ids)


def test_remove_molecules_from_library_not_found(db_session):
    """Test molecule removal from non-existent library."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Create test molecules
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Generate a random UUID for a user
    user_id = uuid.uuid4()
    
    # Mock library.get to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            remove_molecules_from_library(non_existent_id, molecule_ids, user_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)


def test_remove_molecules_from_library_unauthorized(db_session):
    """Test molecule removal failure with unauthorized user."""
    # Create a test library with a specific user ID
    library_id = uuid.uuid4()
    creator_id = uuid.uuid4()
    test_library = Library(
        id=library_id,
        name="Test Library",
        description="Test Description",
        created_by=creator_id,
        created_at=datetime.utcnow()
    )
    
    # Create test molecules and add them to the library
    molecule_ids = [uuid.uuid4(), uuid.uuid4(), uuid.uuid4()]
    
    # Generate a different user ID
    different_user_id = uuid.uuid4()
    
    # Mock library.get to return our test library
    with patch.object(library, 'get', return_value=test_library):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            remove_molecules_from_library(library_id, molecule_ids, different_user_id)
        
        # Assert that the exception message contains "unauthorized"
        assert "unauthorized" in str(excinfo.value)
        
        # Verify that molecules are still associated with the library
        library.remove_molecules.assert_not_called()


def test_check_molecule_in_library(db_session):
    """Test checking if a molecule is in a library."""
    # Create test IDs
    library_id = uuid.uuid4()
    molecule_in_library_id = uuid.uuid4()
    molecule_not_in_library_id = uuid.uuid4()
    
    # Mock library.is_molecule_in_library to return True for one molecule, False for the other
    with patch.object(library, 'is_molecule_in_library', side_effect=lambda db, lib_id, mol_id: mol_id == molecule_in_library_id):
        # Test with molecule that is in the library
        result1 = check_molecule_in_library(library_id, molecule_in_library_id)
        assert result1 is True
        
        # Test with molecule that is not in the library
        result2 = check_molecule_in_library(library_id, molecule_not_in_library_id)
        assert result2 is False


def test_get_library_molecules(db_session):
    """Test retrieval of molecules in a library."""
    # Create test library ID
    library_id = uuid.uuid4()
    
    # Create mock molecules in the library
    test_molecules = [
        MagicMock(
            id=uuid.uuid4(),
            smiles=f"C{i}O",
            created_by=uuid.uuid4(),
            created_at=datetime.utcnow(),
            updated_at=None,
            flag_status=None
        ) for i in range(3)
    ]
    
    total_molecules = len(test_molecules)
    
    # Mock library.get to return a valid library
    # Mock library.get_molecules to return our test molecules
    # Mock library.count_molecules to return the molecule count
    with patch.object(library, 'get', return_value=MagicMock()), \
         patch.object(library, 'get_molecules', return_value=test_molecules), \
         patch.object(library, 'count_molecules', return_value=total_molecules):
        
        # Call the function being tested
        result = get_library_molecules(library_id)
        
        # Assert that all molecules in the library are returned
        assert len(result["items"]) == total_molecules
        assert result["total"] == total_molecules
        
        # Test pagination parameters
        assert result["skip"] == 0
        assert result["limit"] == 100
        
        # Check each molecule has expected data
        for i, mol in enumerate(result["items"]):
            assert mol["id"] == test_molecules[i].id
            assert mol["smiles"] == test_molecules[i].smiles


def test_get_library_molecules_not_found(db_session):
    """Test molecule retrieval from non-existent library."""
    # Generate a random UUID for a non-existent library
    non_existent_id = uuid.uuid4()
    
    # Mock library.get to return None (not found)
    with patch.object(library, 'get', return_value=None):
        # Call the function being tested and expect an exception
        with pytest.raises(LibraryServiceException) as excinfo:
            get_library_molecules(non_existent_id)
        
        # Assert that the exception message contains "not found"
        assert "not found" in str(excinfo.value)