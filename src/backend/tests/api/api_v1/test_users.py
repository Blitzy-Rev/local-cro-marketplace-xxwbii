"""
Test module for user management API endpoints in the Molecular Data Management and CRO Integration Platform.

This file contains comprehensive tests for retrieving, updating, and managing user profiles
with role-based access control.
"""
import pytest
from fastapi import status
from unittest.mock import MagicMock, patch

from app.models.user import User
from app.constants import UserRole, UserStatus
from app.services.auth_service import change_password


@pytest.mark.parametrize('token_headers_fixture, user_fixture', [
    ('pharma_token_headers', 'test_pharma_user'),
    ('cro_token_headers', 'test_cro_user'),
    ('admin_token_headers', 'test_admin_user')
])
def test_get_current_user_profile(client, pharma_token_headers, cro_token_headers, admin_token_headers,
                                test_pharma_user, test_cro_user, test_admin_user, token_headers_fixture, user_fixture):
    """Tests retrieving current user profile information"""
    # Map fixture names to actual fixtures
    token_headers_map = {
        'pharma_token_headers': pharma_token_headers,
        'cro_token_headers': cro_token_headers,
        'admin_token_headers': admin_token_headers
    }
    
    user_map = {
        'test_pharma_user': test_pharma_user,
        'test_cro_user': test_cro_user,
        'test_admin_user': test_admin_user
    }
    
    token_headers = token_headers_map[token_headers_fixture]
    user = user_map[user_fixture]
    
    # Send request to get current user profile
    response = client.get("/api/v1/users/me", headers=token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert data["email"] == user.email
    assert data["id"] == user.id
    assert data["role"] == user.role.name.lower()


@pytest.mark.parametrize('token_headers_fixture, user_fixture', [
    ('pharma_token_headers', 'test_pharma_user'),
    ('cro_token_headers', 'test_cro_user'),
    ('admin_token_headers', 'test_admin_user')
])
def test_update_current_user_profile_success(client, pharma_token_headers, cro_token_headers, admin_token_headers,
                                          test_pharma_user, test_cro_user, test_admin_user, db_session,
                                          token_headers_fixture, user_fixture):
    """Tests successful update of current user profile information"""
    # Map fixture names to actual fixtures
    token_headers_map = {
        'pharma_token_headers': pharma_token_headers,
        'cro_token_headers': cro_token_headers,
        'admin_token_headers': admin_token_headers
    }
    
    user_map = {
        'test_pharma_user': test_pharma_user,
        'test_cro_user': test_cro_user,
        'test_admin_user': test_admin_user
    }
    
    token_headers = token_headers_map[token_headers_fixture]
    user = user_map[user_fixture]
    
    # Create update data
    new_email = f"updated_{user.email}"
    update_data = {"email": new_email}
    
    # Send request to update user profile
    response = client.put("/api/v1/users/me", json=update_data, headers=token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert data["email"] == new_email
    
    # Verify database update
    updated_user = db_session.query(User).filter(User.id == user.id).first()
    assert updated_user.email == new_email


def test_update_current_user_profile_email_exists(client, pharma_token_headers, test_pharma_user, 
                                               test_cro_user, db_session):
    """Tests update failure when email already exists"""
    # Try to update with an email that already exists
    update_data = {"email": test_cro_user.email}
    
    # Send request to update user profile
    response = client.put("/api/v1/users/me", json=update_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    data = response.json()
    assert "email already in use" in str(data).lower()
    
    # Verify database not updated
    user = db_session.query(User).filter(User.id == test_pharma_user.id).first()
    assert user.email == test_pharma_user.email


def test_change_password_success(client, pharma_token_headers, test_pharma_user, test_password, db_session):
    """Tests successful password change with valid current password"""
    # Create password change data
    new_password = "NewPassword123!"
    password_data = {
        "current_password": test_password,
        "new_password": new_password
    }
    
    # Send request to change password
    response = client.post("/api/v1/users/me/change-password", json=password_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert "password changed successfully" in str(data).lower()
    
    # Check that the password was changed by mocking auth_service
    with patch('app.services.auth_service.change_password') as mock_change_password:
        mock_change_password.return_value = test_pharma_user
        
        # Call the service directly to verify - this would normally be done in service tests
        # but we want to ensure the API endpoint actually called the service
        result = change_password(db_session, test_pharma_user.id, new_password, "AnotherNewPassword123!")
        assert result is not None


def test_change_password_incorrect_current(client, pharma_token_headers, db_session):
    """Tests password change failure with incorrect current password"""
    # Create password change data with incorrect current password
    password_data = {
        "current_password": "WrongPassword123!",
        "new_password": "NewPassword123!"
    }
    
    # Send request to change password
    response = client.post("/api/v1/users/me/change-password", json=password_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    data = response.json()
    assert "incorrect current password" in str(data).lower()


@pytest.mark.parametrize('weak_password', [
    'short',
    'nouppercase123!',
    'NOLOWERCASE123!',
    'NoSpecialChar123',
    'NoNumber!'
])
def test_change_password_weak_new(client, pharma_token_headers, test_password, weak_password, db_session):
    """Tests password change failure with weak new password"""
    # Create password change data with weak new password
    password_data = {
        "current_password": test_password,
        "new_password": weak_password
    }
    
    # Send request to change password
    response = client.post("/api/v1/users/me/change-password", json=password_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    data = response.json()
    assert "password requirements" in str(data).lower()


def test_get_user_by_id_admin(client, admin_token_headers, test_pharma_user, db_session):
    """Tests retrieving user by ID as admin"""
    # Send request to get user by ID
    response = client.get(f"/api/v1/users/{test_pharma_user.id}", headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert data["email"] == test_pharma_user.email
    assert data["id"] == test_pharma_user.id
    assert data["role"] == test_pharma_user.role.name.lower()


def test_get_user_by_id_non_admin(client, pharma_token_headers, test_cro_user, db_session):
    """Tests failure to retrieve user by ID as non-admin"""
    # Send request to get user by ID as non-admin
    response = client.get(f"/api/v1/users/{test_cro_user.id}", headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_403_FORBIDDEN
    data = response.json()
    assert "insufficient permissions" in str(data).lower()


def test_get_user_by_id_not_found(client, admin_token_headers, db_session):
    """Tests failure to retrieve non-existent user by ID"""
    # Send request to get non-existent user
    response = client.get("/api/v1/users/999999", headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_404_NOT_FOUND
    data = response.json()
    assert "user not found" in str(data).lower()


def test_update_user_admin(client, admin_token_headers, test_pharma_user, db_session):
    """Tests updating user information as admin"""
    # Create update data
    new_email = f"admin_updated_{test_pharma_user.email}"
    update_data = {"email": new_email}
    
    # Send request to update user as admin
    response = client.put(f"/api/v1/users/{test_pharma_user.id}", json=update_data, headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert data["email"] == new_email
    
    # Verify database update
    updated_user = db_session.query(User).filter(User.id == test_pharma_user.id).first()
    assert updated_user.email == new_email


def test_update_user_non_admin(client, pharma_token_headers, test_cro_user, db_session):
    """Tests failure to update user information as non-admin"""
    # Create update data
    new_email = f"updated_{test_cro_user.email}"
    update_data = {"email": new_email}
    
    # Send request to update user as non-admin
    response = client.put(f"/api/v1/users/{test_cro_user.id}", json=update_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_403_FORBIDDEN
    data = response.json()
    assert "insufficient permissions" in str(data).lower()
    
    # Verify database not updated
    user = db_session.query(User).filter(User.id == test_cro_user.id).first()
    assert user.email == test_cro_user.email


def test_update_user_email_exists(client, admin_token_headers, test_pharma_user, test_cro_user, db_session):
    """Tests update failure when email already exists"""
    # Try to update with an email that already exists
    update_data = {"email": test_cro_user.email}
    
    # Send request to update user
    response = client.put(f"/api/v1/users/{test_pharma_user.id}", json=update_data, headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    data = response.json()
    assert "email already in use" in str(data).lower()
    
    # Verify database not updated
    user = db_session.query(User).filter(User.id == test_pharma_user.id).first()
    assert user.email == test_pharma_user.email


@pytest.mark.parametrize('new_status', ['active', 'inactive', 'pending'])
def test_update_user_status_admin(client, admin_token_headers, test_pharma_user, new_status, db_session):
    """Tests updating user status as admin"""
    # Create status update data
    status_data = {"status": new_status}
    
    # Send request to update user status
    response = client.put(f"/api/v1/users/{test_pharma_user.id}/status", json=status_data, headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert data["status"] == new_status
    
    # Verify database update
    updated_user = db_session.query(User).filter(User.id == test_pharma_user.id).first()
    assert updated_user.status.name.lower() == new_status


def test_update_user_status_non_admin(client, pharma_token_headers, test_cro_user, db_session):
    """Tests failure to update user status as non-admin"""
    # Create status update data
    status_data = {"status": "inactive"}
    
    # Send request to update user status as non-admin
    response = client.put(f"/api/v1/users/{test_cro_user.id}/status", json=status_data, headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_403_FORBIDDEN
    data = response.json()
    assert "insufficient permissions" in str(data).lower()


def test_update_user_status_invalid(client, admin_token_headers, test_pharma_user, db_session):
    """Tests update failure with invalid status"""
    # Create status update data with invalid status
    status_data = {"status": "invalid_status"}
    
    # Send request to update user status
    response = client.put(f"/api/v1/users/{test_pharma_user.id}/status", json=status_data, headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    data = response.json()
    assert "invalid status" in str(data).lower()


def test_get_users_admin(client, admin_token_headers, db_session):
    """Tests retrieving list of users as admin"""
    # Send request to get users as admin
    response = client.get("/api/v1/users", headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert "items" in data
    assert isinstance(data["items"], list)
    assert len(data["items"]) > 0
    assert "total" in data
    assert "page" in data
    assert "size" in data


@pytest.mark.parametrize('filter_param, filter_value, expected_count', [
    ('role', 'pharma', 1),
    ('status', 'active', 3),
    ('email', 'pharma', 1)
])
def test_get_users_with_filters(client, admin_token_headers, filter_param, filter_value, expected_count, db_session):
    """Tests retrieving filtered list of users as admin"""
    # Send request to get users with filter
    response = client.get(f"/api/v1/users?{filter_param}={filter_value}", headers=admin_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_200_OK
    data = response.json()
    assert "items" in data
    assert isinstance(data["items"], list)
    assert len(data["items"]) == expected_count
    
    # Verify filter was applied correctly
    for item in data["items"]:
        if filter_param == 'role':
            assert item["role"] == filter_value
        elif filter_param == 'status':
            assert item["status"] == filter_value
        elif filter_param == 'email':
            assert filter_value in item["email"]


def test_get_users_pagination(client, admin_token_headers, db_session):
    """Tests user list pagination as admin"""
    # Create additional test users to ensure pagination
    for i in range(5):
        user = User(
            email=f"test_user_{i}@example.com",
            password_hash="hashed_password",
            role=UserRole.PHARMA,
            status=UserStatus.ACTIVE,
            is_active=True,
            email_verified=True,
            password_history=[]
        )
        db_session.add(user)
    db_session.commit()
    
    # Send request for first page
    response_page1 = client.get("/api/v1/users?page=1&size=2", headers=admin_token_headers)
    
    # Assert first page response
    assert response_page1.status_code == status.HTTP_200_OK
    data_page1 = response_page1.json()
    assert len(data_page1["items"]) == 2
    assert data_page1["page"] == 1
    assert data_page1["size"] == 2
    
    # Store IDs from first page
    page1_ids = [item["id"] for item in data_page1["items"]]
    
    # Send request for second page
    response_page2 = client.get("/api/v1/users?page=2&size=2", headers=admin_token_headers)
    
    # Assert second page response
    assert response_page2.status_code == status.HTTP_200_OK
    data_page2 = response_page2.json()
    assert len(data_page2["items"]) == 2
    
    # Verify different users on each page
    page2_ids = [item["id"] for item in data_page2["items"]]
    assert not any(id in page1_ids for id in page2_ids)


def test_get_users_non_admin(client, pharma_token_headers, db_session):
    """Tests failure to retrieve list of users as non-admin"""
    # Send request to get users as non-admin
    response = client.get("/api/v1/users", headers=pharma_token_headers)
    
    # Assert response
    assert response.status_code == status.HTTP_403_FORBIDDEN
    data = response.json()
    assert "insufficient permissions" in str(data).lower()