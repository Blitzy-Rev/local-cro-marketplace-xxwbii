"""
Test module for library-related API endpoints in the Molecular Data Management and CRO Integration Platform.

This file contains comprehensive tests for creating, retrieving, updating, and deleting libraries,
as well as managing molecules within libraries.
"""

import pytest
import json
import uuid
from fastapi import status
from unittest.mock import MagicMock, patch

from app.services.library_service import library_service
from app.exceptions import LibraryServiceException


def test_create_library_success(client, db_session, pharma_token_headers):
    """Tests successful library creation with valid data."""
    # Create library data with valid name and description
    library_data = {
        "name": "Test Library",
        "description": "A test library for unit testing"
    }
    
    # Send POST request to /api/v1/libraries/ with library data and pharma token headers
    response = client.post(
        "/api/v1/libraries/",
        json=library_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 201 Created
    assert response.status_code == status.HTTP_201_CREATED
    
    # Assert response contains library data with correct name and description
    response_data = response.json()
    assert response_data["name"] == library_data["name"]
    assert response_data["description"] == library_data["description"]
    
    # Verify library exists in database with correct data
    created_id = response_data["id"]
    db_library = library_service.get_by_id(db_session, created_id)
    assert db_library is not None
    assert db_library.name == library_data["name"]
    assert db_library.description == library_data["description"]


@pytest.mark.parametrize("invalid_name", ['', ' ', 'a' * 101])
def test_create_library_invalid_name(client, pharma_token_headers):
    """Tests library creation failure with invalid name."""
    # Create library data with invalid name
    library_data = {
        "name": invalid_name,
        "description": "A test library with invalid name"
    }
    
    # Send POST request to /api/v1/libraries/ with library data and pharma token headers
    response = client.post(
        "/api/v1/libraries/",
        json=library_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 422 Unprocessable Entity
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    
    # Assert response contains validation error about invalid name
    response_data = response.json()
    assert "name" in str(response_data["detail"]).lower()


def test_create_library_duplicate_name(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests library creation failure with duplicate name."""
    # Create a test library with specific name
    existing_library = create_test_library("Existing Library", test_pharma_user.id)
    
    # Create library data with the same name as existing library
    library_data = {
        "name": "Existing Library",
        "description": "This should fail due to duplicate name"
    }
    
    # Send POST request to /api/v1/libraries/ with library data and pharma token headers
    response = client.post(
        "/api/v1/libraries/",
        json=library_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 400 Bad Request
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    
    # Assert response contains error about duplicate library name
    response_data = response.json()
    assert "already exists" in response_data["detail"]["message"].lower()


def test_create_library_unauthorized(client):
    """Tests library creation failure without proper authentication."""
    # Create library data with valid name and description
    library_data = {
        "name": "Unauthorized Library",
        "description": "This should fail due to missing authentication"
    }
    
    # Send POST request to /api/v1/libraries/ with library data and no authentication
    response = client.post(
        "/api/v1/libraries/",
        json=library_data
    )
    
    # Assert response status code is 401 Unauthorized
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    
    # Assert response contains error about missing token
    response_data = response.json()
    assert "not authenticated" in response_data["detail"].lower()


def test_get_library_success(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests successful retrieval of a library by ID."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Test Library", test_pharma_user.id)
    
    # Send GET request to /api/v1/libraries/{library_id} with pharma token headers
    response = client.get(
        f"/api/v1/libraries/{test_library.id}",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains library data with correct ID, name, and description
    response_data = response.json()
    assert str(test_library.id) == response_data["id"]
    assert test_library.name == response_data["name"]
    assert test_library.description == response_data["description"]
    
    # Assert response contains molecules list (empty for new library)
    assert "molecules" in response_data
    assert isinstance(response_data["molecules"], list)


def test_get_library_not_found(client, pharma_token_headers):
    """Tests library retrieval failure with non-existent ID."""
    # Generate a random UUID for non-existent library
    random_id = str(uuid.uuid4())
    
    # Send GET request to /api/v1/libraries/{random_id} with pharma token headers
    response = client.get(
        f"/api/v1/libraries/{random_id}",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND
    
    # Assert response contains error about library not found
    response_data = response.json()
    assert "not found" in response_data["detail"]["message"].lower()


def test_get_library_unauthorized_access(client, db_session, pharma_token_headers, test_admin_user, create_test_library):
    """Tests library retrieval failure when accessing another user's library."""
    # Create a test library owned by the admin user
    admin_library = create_test_library("Admin Library", test_admin_user.id)
    
    # Send GET request to /api/v1/libraries/{library_id} with pharma token headers
    response = client.get(
        f"/api/v1/libraries/{admin_library.id}",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 403 Forbidden
    assert response.status_code == status.HTTP_403_FORBIDDEN
    
    # Assert response contains error about unauthorized access
    response_data = response.json()
    assert "not authorized" in response_data["detail"]["message"].lower()


def test_update_library_success(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests successful update of a library."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Original Name", test_pharma_user.id)
    
    # Create update data with new name and description
    update_data = {
        "name": "Updated Name",
        "description": "Updated description"
    }
    
    # Send PUT request to /api/v1/libraries/{library_id} with update data and pharma token headers
    response = client.put(
        f"/api/v1/libraries/{test_library.id}",
        json=update_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains updated library data with new name and description
    response_data = response.json()
    assert response_data["name"] == update_data["name"]
    assert response_data["description"] == update_data["description"]
    
    # Verify library is updated in database with correct data
    db_library = library_service.get_by_id(db_session, test_library.id)
    assert db_library.name == update_data["name"]
    assert db_library.description == update_data["description"]


def test_update_library_not_owner(client, db_session, pharma_token_headers, test_admin_user, create_test_library):
    """Tests library update failure when user is not the owner."""
    # Create a test library owned by the admin user
    admin_library = create_test_library("Admin Library", test_admin_user.id)
    
    # Create update data with new name and description
    update_data = {
        "name": "Unauthorized Update",
        "description": "This should fail"
    }
    
    # Send PUT request to /api/v1/libraries/{library_id} with update data and pharma token headers
    response = client.put(
        f"/api/v1/libraries/{admin_library.id}",
        json=update_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 403 Forbidden
    assert response.status_code == status.HTTP_403_FORBIDDEN
    
    # Assert response contains error about not being the owner
    response_data = response.json()
    assert "not authorized" in response_data["detail"]["message"].lower()


def test_update_library_duplicate_name(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests library update failure with duplicate name."""
    # Create two test libraries owned by the pharma user
    library1 = create_test_library("Library One", test_pharma_user.id)
    library2 = create_test_library("Library Two", test_pharma_user.id)
    
    # Create update data with the name of the second library
    update_data = {
        "name": "Library Two"
    }
    
    # Send PUT request to update the first library with the name of the second library
    response = client.put(
        f"/api/v1/libraries/{library1.id}",
        json=update_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 400 Bad Request
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    
    # Assert response contains error about duplicate library name
    response_data = response.json()
    assert "already exists" in response_data["detail"]["message"].lower()


def test_delete_library_success(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests successful deletion of a library."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Library to Delete", test_pharma_user.id)
    
    # Send DELETE request to /api/v1/libraries/{library_id} with pharma token headers
    response = client.delete(
        f"/api/v1/libraries/{test_library.id}",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains success message
    response_data = response.json()
    assert response_data["success"] is True
    
    # Verify library no longer exists in database
    db_library = library_service.get_by_id(db_session, test_library.id)
    assert db_library is None


def test_delete_library_not_owner(client, db_session, pharma_token_headers, test_admin_user, create_test_library):
    """Tests library deletion failure when user is not the owner."""
    # Create a test library owned by the admin user
    admin_library = create_test_library("Admin Library", test_admin_user.id)
    
    # Send DELETE request to /api/v1/libraries/{library_id} with pharma token headers
    response = client.delete(
        f"/api/v1/libraries/{admin_library.id}",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 403 Forbidden
    assert response.status_code == status.HTTP_403_FORBIDDEN
    
    # Assert response contains error about not being the owner
    response_data = response.json()
    assert "not authorized" in response_data["detail"]["message"].lower()
    
    # Verify library still exists in database
    db_library = library_service.get_by_id(db_session, admin_library.id)
    assert db_library is not None


@pytest.mark.parametrize("filter_params", [{"name_contains": "Test"}, {"sort_by": "created_at", "sort_desc": "true"}])
def test_get_libraries_with_filters(client, db_session, pharma_token_headers, create_test_library, test_pharma_user, filter_params):
    """Tests library listing with various filters."""
    # Create multiple test libraries with different names and creation times
    create_test_library("Test Library 1", test_pharma_user.id)
    create_test_library("Test Library 2", test_pharma_user.id)
    create_test_library("Other Library", test_pharma_user.id)
    
    # Send GET request to /api/v1/libraries/ with filter parameters and pharma token headers
    response = client.get(
        "/api/v1/libraries/",
        params=filter_params,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains libraries list with pagination info
    response_data = response.json()
    assert "items" in response_data
    assert "total" in response_data
    
    # Verify filtered results match the expected libraries based on filter criteria
    if "name_contains" in filter_params and filter_params["name_contains"] == "Test":
        assert all("test" in item["name"].lower() for item in response_data["items"])


def test_get_libraries_pagination(client, db_session, pharma_token_headers, create_test_library, test_pharma_user):
    """Tests library listing with pagination."""
    # Create multiple test libraries (more than the default page size)
    for i in range(25):  # Create 25 libraries
        create_test_library(f"Test Library {i}", test_pharma_user.id)
    
    # Send GET request to /api/v1/libraries/ with skip and limit parameters and pharma token headers
    response1 = client.get(
        "/api/v1/libraries/",
        params={"skip": 0, "limit": 10},
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response1.status_code == status.HTTP_200_OK
    
    # Assert response contains libraries list with correct pagination info
    data1 = response1.json()
    assert len(data1["items"]) == 10
    assert data1["total"] >= 25
    
    # Verify returned libraries match the expected page of results
    first_page_ids = {item["id"] for item in data1["items"]}
    
    # Send another request for a different page and verify correct results
    response2 = client.get(
        "/api/v1/libraries/",
        params={"skip": 10, "limit": 10},
        headers=pharma_token_headers
    )
    
    assert response2.status_code == status.HTTP_200_OK
    data2 = response2.json()
    second_page_ids = {item["id"] for item in data2["items"]}
    
    # Ensure no overlap between pages
    assert not first_page_ids.intersection(second_page_ids)


def test_get_library_molecules(client, db_session, pharma_token_headers, create_test_library, create_test_molecule, test_pharma_user):
    """Tests retrieval of molecules in a library."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Test Molecule Library", test_pharma_user.id)
    
    # Create multiple test molecules
    molecule1 = create_test_molecule("CCO", test_pharma_user.id)
    molecule2 = create_test_molecule("CCC", test_pharma_user.id)
    
    # Add molecules to the library using library_service.add_molecules_to_library
    library_service.add_molecules_to_library(test_library.id, [molecule1.id, molecule2.id], test_pharma_user.id)
    
    # Send GET request to /api/v1/libraries/{library_id}/molecules with pharma token headers
    response = client.get(
        f"/api/v1/libraries/{test_library.id}/molecules",
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains molecules list with pagination info
    response_data = response.json()
    assert "items" in response_data
    assert "total" in response_data
    assert response_data["total"] == 2
    
    # Verify returned molecules match those added to the library
    molecule_ids = {item["id"] for item in response_data["items"]}
    assert str(molecule1.id) in molecule_ids
    assert str(molecule2.id) in molecule_ids


def test_get_library_molecules_pagination(client, db_session, pharma_token_headers, create_test_library, create_test_molecule, test_pharma_user):
    """Tests retrieval of molecules in a library with pagination."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Pagination Library", test_pharma_user.id)
    
    # Create multiple test molecules (more than the default page size)
    molecule_ids = []
    for i in range(25):  # Create 25 molecules
        mol = create_test_molecule(f"C{'C' * i}O", test_pharma_user.id)
        molecule_ids.append(mol.id)
    
    # Add all molecules to the library using library_service.add_molecules_to_library
    library_service.add_molecules_to_library(test_library.id, molecule_ids, test_pharma_user.id)
    
    # Send GET request to /api/v1/libraries/{library_id}/molecules with skip and limit parameters
    response = client.get(
        f"/api/v1/libraries/{test_library.id}/molecules",
        params={"skip": 0, "limit": 10},
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains molecules list with correct pagination info
    response_data = response.json()
    assert len(response_data["items"]) == 10
    assert response_data["total"] == 25
    
    # Verify returned molecules match the expected page of results
    first_page = response_data["items"]
    
    # Request second page to verify pagination works correctly
    response2 = client.get(
        f"/api/v1/libraries/{test_library.id}/molecules",
        params={"skip": 10, "limit": 10},
        headers=pharma_token_headers
    )
    
    assert response2.status_code == status.HTTP_200_OK
    second_page = response2.json()["items"]
    
    # Ensure pages contain different molecules
    first_page_ids = {item["id"] for item in first_page}
    second_page_ids = {item["id"] for item in second_page}
    assert not first_page_ids.intersection(second_page_ids)


def test_add_molecules_to_library(client, db_session, pharma_token_headers, create_test_library, create_test_molecule, test_pharma_user):
    """Tests adding molecules to a library."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Add Molecules Library", test_pharma_user.id)
    
    # Create multiple test molecules
    molecule1 = create_test_molecule("CCO", test_pharma_user.id)
    molecule2 = create_test_molecule("CCC", test_pharma_user.id)

    # Create molecule operation data with 'add' operation and molecule IDs
    operation_data = {
        "operation": "add",
        "molecule_ids": [str(molecule1.id), str(molecule2.id)]
    }
    
    # Send POST request to /api/v1/libraries/{library_id}/molecules with operation data
    response = client.post(
        f"/api/v1/libraries/{test_library.id}/molecules",
        json=operation_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains success count matching number of molecules
    response_data = response.json()
    assert response_data["success"] == 2
    
    # Verify molecules are added to the library in database
    library_molecules = library_service.get_library_molecules(db_session, test_library.id)
    assert library_molecules["total"] == 2


def test_add_molecules_to_library_not_owner(client, db_session, pharma_token_headers, test_admin_user, create_test_library, create_test_molecule, test_pharma_user):
    """Tests adding molecules to a library failure when user is not the owner."""
    # Create a test library owned by the admin user
    admin_library = create_test_library("Admin Library", test_admin_user.id)
    
    # Create test molecules
    molecule = create_test_molecule("CCO", test_pharma_user.id)

    # Create molecule operation data with 'add' operation and molecule IDs
    operation_data = {
        "operation": "add",
        "molecule_ids": [str(molecule.id)]
    }
    
    # Send POST request to /api/v1/libraries/{library_id}/molecules with operation data and pharma token headers
    response = client.post(
        f"/api/v1/libraries/{admin_library.id}/molecules",
        json=operation_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 403 Forbidden
    assert response.status_code == status.HTTP_403_FORBIDDEN
    
    # Assert response contains error about not being the owner
    response_data = response.json()
    assert "not authorized" in response_data["detail"]["message"].lower()
    
    # Verify molecules are not added to the library in database
    library_molecules = library_service.get_library_molecules(db_session, admin_library.id)
    assert library_molecules["total"] == 0


def test_remove_molecules_from_library(client, db_session, pharma_token_headers, create_test_library, create_test_molecule, test_pharma_user):
    """Tests removing molecules from a library."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Remove Molecules Library", test_pharma_user.id)
    
    # Create multiple test molecules
    molecule1 = create_test_molecule("CCO", test_pharma_user.id)
    molecule2 = create_test_molecule("CCC", test_pharma_user.id)
    
    # Add molecules to the library using library_service.add_molecules_to_library
    library_service.add_molecules_to_library(test_library.id, [molecule1.id, molecule2.id], test_pharma_user.id)
    
    # Create molecule operation data with 'remove' operation and molecule IDs
    operation_data = {
        "operation": "remove",
        "molecule_ids": [str(molecule1.id)]
    }
    
    # Send POST request to /api/v1/libraries/{library_id}/molecules with operation data
    response = client.post(
        f"/api/v1/libraries/{test_library.id}/molecules",
        json=operation_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK
    
    # Assert response contains success count matching number of molecules
    response_data = response.json()
    assert response_data["success"] == 1
    
    # Verify molecules are removed from the library in database
    library_molecules = library_service.get_library_molecules(db_session, test_library.id)
    assert library_molecules["total"] == 1
    remaining_molecule_id = library_molecules["items"][0]["id"]
    assert str(molecule2.id) == remaining_molecule_id


def test_invalid_molecule_operation(client, db_session, pharma_token_headers, create_test_library, create_test_molecule, test_pharma_user):
    """Tests molecule operation failure with invalid operation type."""
    # Create a test library owned by the pharma user
    test_library = create_test_library("Invalid Operation Library", test_pharma_user.id)
    
    # Create test molecules
    molecule = create_test_molecule("CCO", test_pharma_user.id)
    
    # Create molecule operation data with invalid operation type and molecule IDs
    operation_data = {
        "operation": "invalid_operation",
        "molecule_ids": [str(molecule.id)]
    }
    
    # Send POST request to /api/v1/libraries/{library_id}/molecules with operation data
    response = client.post(
        f"/api/v1/libraries/{test_library.id}/molecules",
        json=operation_data,
        headers=pharma_token_headers
    )
    
    # Assert response status code is 422 Unprocessable Entity
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    
    # Assert response contains validation error about invalid operation type
    response_data = response.json()
    assert "operation" in str(response_data["detail"]).lower()


def test_error_handling(client, pharma_token_headers):
    """Tests error handling for library service exceptions."""
    # Mock library_service.create to raise LibraryServiceException with specific error message
    error_message = "Test error message"
    with patch("app.services.library_service.create", side_effect=LibraryServiceException(error_message, {})):
        # Create library data with valid name and description
        library_data = {
            "name": "Error Test Library",
            "description": "Test library for error handling"
        }
        
        # Send POST request to /api/v1/libraries/ with library data and pharma token headers
        response = client.post(
            "/api/v1/libraries/",
            json=library_data,
            headers=pharma_token_headers
        )
        
        # Assert response status code is 400 Bad Request
        assert response.status_code == status.HTTP_400_BAD_REQUEST
        
        # Assert response contains error message from the exception
        response_data = response.json()
        assert error_message in response_data["detail"]["message"]