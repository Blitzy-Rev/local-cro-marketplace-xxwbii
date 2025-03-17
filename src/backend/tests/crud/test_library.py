"""
Unit tests for the library CRUD operations in the Molecular Data Management and CRO Integration Platform.

This test module verifies the functionality of library database operations including creation, 
retrieval, filtering, and molecule association management.
"""

import pytest
import random
from typing import List, Dict, Optional
import uuid

from app.crud.crud_library import library
from app.models.library import Library
from app.models.library_molecule import LibraryMolecule
from app.schemas.library import LibraryCreate


def test_get_by_name(db_session, create_test_library, test_pharma_user):
    """Test retrieving a library by its name and creator ID"""
    # Create a test library
    test_name = f"Test Library {random.randint(1, 1000)}"
    db_library = create_test_library(test_name, test_pharma_user.id)
    
    # Test retrieval by name and user ID
    result = library.get_by_name(db_session, test_name, test_pharma_user.id)
    assert result is not None
    assert result.name == test_name
    
    # Test with non-existent name
    non_existent = library.get_by_name(db_session, "Non-existent Library", test_pharma_user.id)
    assert non_existent is None
    
    # Test with existing name but wrong user ID
    wrong_user_id = uuid.uuid4()
    wrong_user = library.get_by_name(db_session, test_name, wrong_user_id)
    assert wrong_user is None


def test_get_with_molecules(db_session, create_test_library, create_test_molecule):
    """Test retrieving a library with its associated molecules"""
    # Create test library
    test_library = create_test_library(f"Molecule Library {random.randint(1, 1000)}", uuid.uuid4())
    
    # Create test molecules
    molecules = []
    for i in range(5):
        molecule = create_test_molecule(f"C{i}CO", test_library.created_by)
        molecules.append(molecule)
    
    # Add molecules to library
    db_session.flush()
    library.add_molecules(db_session, test_library.id, [m.id for m in molecules])
    
    # Retrieve library with molecules
    result = library.get_with_molecules(db_session, test_library.id)
    assert result is not None
    lib, mols = result
    assert lib.id == test_library.id
    assert len(mols) == 5
    
    # Test with pagination
    result_paginated = library.get_with_molecules(db_session, test_library.id, skip=2, limit=2)
    assert result_paginated is not None
    lib, mols = result_paginated
    assert len(mols) == 2
    
    # Test with non-existent library ID
    non_existent = library.get_with_molecules(db_session, uuid.uuid4())
    assert non_existent is None


def test_get_by_user(db_session, create_test_library, test_pharma_user):
    """Test retrieving libraries created by a specific user"""
    # Create multiple test libraries for the user
    num_libraries = 5
    created_libraries = []
    for i in range(num_libraries):
        lib = create_test_library(f"User Library {i}", test_pharma_user.id)
        created_libraries.append(lib)
    
    # Retrieve libraries for the user
    result = library.get_by_user(db_session, test_pharma_user.id)
    assert len(result) >= num_libraries  # May have libraries from other tests
    
    # Check that all retrieved libraries belong to the user
    for lib in result:
        assert lib.created_by == test_pharma_user.id
    
    # Test with pagination
    paginated = library.get_by_user(db_session, test_pharma_user.id, skip=1, limit=2)
    assert len(paginated) == 2
    
    # Test with a user that has no libraries
    non_user_id = uuid.uuid4()
    empty_result = library.get_by_user(db_session, non_user_id)
    assert len(empty_result) == 0


def test_count_by_user(db_session, create_test_library, test_pharma_user):
    """Test counting libraries created by a specific user"""
    # Create multiple test libraries for the user
    initial_count = library.count_by_user(db_session, test_pharma_user.id)
    num_libraries = 3
    for i in range(num_libraries):
        create_test_library(f"Count Library {i}", test_pharma_user.id)
    
    # Count libraries for the user
    count = library.count_by_user(db_session, test_pharma_user.id)
    assert count == initial_count + num_libraries
    
    # Test with a user that has no libraries
    non_user_id = uuid.uuid4()
    empty_count = library.count_by_user(db_session, non_user_id)
    assert empty_count == 0


def test_add_molecules(db_session, create_test_library, create_test_molecule):
    """Test adding molecules to a library"""
    # Create a test library
    test_library = create_test_library("Add Molecules Test", uuid.uuid4())
    
    # Create test molecules
    num_molecules = 5
    molecules = []
    for i in range(num_molecules):
        molecule = create_test_molecule(f"CC{i}O", test_library.created_by)
        molecules.append(molecule)
    
    # Add molecules to library
    molecule_ids = [m.id for m in molecules]
    result = library.add_molecules(db_session, test_library.id, molecule_ids)
    
    # Check results
    assert result["success"] == num_molecules
    assert len(result["failures"]) == 0
    
    # Verify molecules were added
    for molecule_id in molecule_ids:
        assert library.is_molecule_in_library(db_session, test_library.id, molecule_id)
    
    # Test adding same molecules again (should handle duplicates)
    duplicate_result = library.add_molecules(db_session, test_library.id, molecule_ids)
    assert duplicate_result["success"] == 0
    assert len(duplicate_result["failures"]) == num_molecules
    
    # Test with non-existent molecule IDs
    fake_id = uuid.uuid4()
    fake_result = library.add_molecules(db_session, test_library.id, [fake_id])
    assert fake_result["success"] == 0
    assert len(fake_result["failures"]) == 1
    
    # Test with non-existent library ID
    non_lib_result = library.add_molecules(db_session, uuid.uuid4(), molecule_ids)
    assert non_lib_result["success"] == 0
    assert len(non_lib_result["failures"]) == 1


def test_remove_molecules(db_session, create_test_library, create_test_molecule):
    """Test removing molecules from a library"""
    # Create a test library
    test_library = create_test_library("Remove Molecules Test", uuid.uuid4())
    
    # Create test molecules and add to library
    num_molecules = 5
    molecules = []
    for i in range(num_molecules):
        molecule = create_test_molecule(f"CC{i}C", test_library.created_by)
        molecules.append(molecule)
    
    molecule_ids = [m.id for m in molecules]
    library.add_molecules(db_session, test_library.id, molecule_ids)
    
    # Remove some molecules
    to_remove = molecule_ids[:2]
    result = library.remove_molecules(db_session, test_library.id, to_remove)
    
    # Check results
    assert result["success"] == 2
    assert len(result["failures"]) == 0
    
    # Verify molecules were removed
    for molecule_id in to_remove:
        assert not library.is_molecule_in_library(db_session, test_library.id, molecule_id)
    
    # Verify remaining molecules are still there
    for molecule_id in molecule_ids[2:]:
        assert library.is_molecule_in_library(db_session, test_library.id, molecule_id)
    
    # Test removing molecules that aren't in the library
    non_member_result = library.remove_molecules(db_session, test_library.id, to_remove)
    assert non_member_result["success"] < len(to_remove)
    
    # Test with non-existent molecule IDs
    fake_id = uuid.uuid4()
    fake_result = library.remove_molecules(db_session, test_library.id, [fake_id])
    assert fake_result["success"] == 0
    
    # Test with non-existent library ID
    non_lib_result = library.remove_molecules(db_session, uuid.uuid4(), molecule_ids)
    assert non_lib_result["success"] == 0
    assert len(non_lib_result["failures"]) == 1


def test_get_molecules(db_session, create_test_library, create_test_molecule):
    """Test retrieving molecules in a library with pagination"""
    # Create a test library
    test_library = create_test_library("Get Molecules Test", uuid.uuid4())
    
    # Create test molecules and add to library
    num_molecules = 10
    molecules = []
    for i in range(num_molecules):
        molecule = create_test_molecule(f"CC{i}N", test_library.created_by)
        molecules.append(molecule)
    
    molecule_ids = [m.id for m in molecules]
    library.add_molecules(db_session, test_library.id, molecule_ids)
    
    # Retrieve all molecules
    all_molecules = library.get_molecules(db_session, test_library.id)
    assert len(all_molecules) == num_molecules
    
    # Test with pagination
    paginated = library.get_molecules(db_session, test_library.id, skip=3, limit=4)
    assert len(paginated) == 4
    
    # Test with non-existent library ID
    non_existent = library.get_molecules(db_session, uuid.uuid4())
    assert len(non_existent) == 0


def test_count_molecules(db_session, create_test_library, create_test_molecule):
    """Test counting molecules in a library"""
    # Create a test library
    test_library = create_test_library("Count Molecules Test", uuid.uuid4())
    
    # Create test molecules and add to library
    num_molecules = 7
    molecules = []
    for i in range(num_molecules):
        molecule = create_test_molecule(f"CC{i}P", test_library.created_by)
        molecules.append(molecule)
    
    library.add_molecules(db_session, test_library.id, [m.id for m in molecules])
    
    # Count molecules
    count = library.count_molecules(db_session, test_library.id)
    assert count == num_molecules
    
    # Test with non-existent library ID
    non_existent_count = library.count_molecules(db_session, uuid.uuid4())
    assert non_existent_count == 0


def test_is_molecule_in_library(db_session, create_test_library, create_test_molecule):
    """Test checking if a molecule is in a library"""
    # Create a test library
    test_library = create_test_library("Molecule Check Test", uuid.uuid4())
    
    # Create test molecules
    molecules = []
    for i in range(5):
        molecule = create_test_molecule(f"CC{i}S", test_library.created_by)
        molecules.append(molecule)
    
    # Add only the first 3 molecules to the library
    library.add_molecules(db_session, test_library.id, [m.id for m in molecules[:3]])
    
    # Test molecules that are in the library
    for molecule in molecules[:3]:
        assert library.is_molecule_in_library(db_session, test_library.id, molecule.id)
    
    # Test molecules that are not in the library
    for molecule in molecules[3:]:
        assert not library.is_molecule_in_library(db_session, test_library.id, molecule.id)
    
    # Test with non-existent molecule ID
    assert not library.is_molecule_in_library(db_session, test_library.id, uuid.uuid4())
    
    # Test with non-existent library ID
    assert not library.is_molecule_in_library(db_session, uuid.uuid4(), molecules[0].id)


def test_filter_libraries(db_session, create_test_library, test_pharma_user):
    """Test filtering libraries based on criteria"""
    # Create test libraries with different names and creators
    user_id = test_pharma_user.id
    other_user_id = uuid.uuid4()
    
    # Create libraries for main test user
    create_test_library("Alpha Chemistry", user_id)
    create_test_library("Beta Chemistry", user_id)
    create_test_library("Alpha Physics", user_id)
    
    # Create libraries for other user
    create_test_library("Alpha Biology", other_user_id)
    create_test_library("Gamma Molecules", other_user_id)
    
    # Test filtering by name substring
    alpha_libs, alpha_count = library.filter_libraries(db_session, name_contains="Alpha")
    assert alpha_count == 3  # Alpha libraries across both users
    assert len(alpha_libs) == 3
    
    # Test filtering by user ID
    user_libs, user_count = library.filter_libraries(db_session, created_by=user_id)
    assert user_count >= 3  # Could include libraries from other tests
    for lib in user_libs:
        assert lib.created_by == user_id
    
    # Test filtering by both name and user ID
    alpha_user_libs, alpha_user_count = library.filter_libraries(
        db_session, name_contains="Alpha", created_by=user_id
    )
    assert alpha_user_count == 2  # Alpha Chemistry and Alpha Physics for the user
    assert len(alpha_user_libs) == 2
    
    # Test with sorting (ascending)
    sorted_libs, _ = library.filter_libraries(db_session, created_by=user_id, sort_by="name")
    if len(sorted_libs) >= 2:
        assert sorted_libs[0].name < sorted_libs[1].name
    
    # Test with sorting (descending)
    desc_libs, _ = library.filter_libraries(db_session, created_by=user_id, sort_by="name", sort_desc=True)
    if len(desc_libs) >= 2:
        assert desc_libs[0].name > desc_libs[1].name
    
    # Test with pagination
    paginated, total = library.filter_libraries(db_session, skip=1, limit=2)
    assert len(paginated) == 2
    assert total > 2  # We created at least 5 libraries
    
    # Test with criteria that match no libraries
    empty_libs, empty_count = library.filter_libraries(db_session, name_contains="NonExistentName")
    assert empty_count == 0
    assert len(empty_libs) == 0