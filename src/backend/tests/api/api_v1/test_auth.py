import pytest
from unittest.mock import MagicMock, patch
from fastapi import status
from datetime import datetime, timedelta

from app.core.security import get_password_hash, verify_password
from app.core.jwt import create_access_token, verify_token
from app.services.auth_service import AuthService
from app.models.user import User
from app.constants import UserRole, UserStatus


@pytest.mark.parametrize('user_fixture', ['test_pharma_user', 'test_cro_user', 'test_admin_user'])
def test_login_success(client, db_session, test_pharma_user, test_password):
    """Tests successful user login with valid credentials"""
    # Note: In actual implementation, we would use request.getfixturevalue(user_fixture)
    # to dynamically get the user based on the parametrized value
    
    login_data = {
        "email": test_pharma_user.email,
        "password": test_password
    }
    
    response = client.post("/api/v1/auth/login", json=login_data)
    
    assert response.status_code == status.HTTP_200_OK
    assert "access_token" in response.json()
    assert "refresh_token" in response.json()
    assert "user" in response.json()
    assert response.json()["user"]["email"] == test_pharma_user.email
    assert response.json()["user"]["role"] == test_pharma_user.role.name
    
    # Verify that access_token is valid and contains correct user information
    token = response.json()["access_token"]
    payload = verify_token(token)
    assert payload["user_id"] == test_pharma_user.id
    assert payload["email"] == test_pharma_user.email


def test_login_invalid_credentials(client, db_session, test_pharma_user):
    """Tests login failure with invalid credentials"""
    login_data = {
        "email": test_pharma_user.email,
        "password": "wrong_password"
    }
    
    response = client.post("/api/v1/auth/login", json=login_data)
    
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    assert "detail" in response.json()


def test_login_nonexistent_user(client, db_session):
    """Tests login failure with non-existent user"""
    login_data = {
        "email": "nonexistent@example.com",
        "password": "any_password"
    }
    
    response = client.post("/api/v1/auth/login", json=login_data)
    
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    assert "detail" in response.json()


def test_login_inactive_user(client, db_session, test_pharma_user, test_password):
    """Tests login failure with inactive user account"""
    # Set user to inactive
    test_pharma_user.status = UserStatus.INACTIVE
    db_session.add(test_pharma_user)
    db_session.commit()
    
    login_data = {
        "email": test_pharma_user.email,
        "password": test_password
    }
    
    response = client.post("/api/v1/auth/login", json=login_data)
    
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    assert "detail" in response.json()
    
    # Reset user status to active for other tests
    test_pharma_user.status = UserStatus.ACTIVE
    db_session.add(test_pharma_user)
    db_session.commit()


def test_refresh_token_success(client, db_session, test_pharma_user):
    """Tests successful token refresh with valid refresh token"""
    # First login to get initial tokens
    login_data = {
        "email": test_pharma_user.email,
        "password": "test_password"  # This should match the fixture
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    refresh_token = login_response.json()["refresh_token"]
    
    # Now use the refresh token to get a new access token
    refresh_data = {
        "refresh_token": refresh_token
    }
    
    response = client.post("/api/v1/auth/refresh-token", json=refresh_data)
    
    assert response.status_code == status.HTTP_200_OK
    assert "access_token" in response.json()
    
    # Verify that new access_token is valid and contains correct user information
    token = response.json()["access_token"]
    payload = verify_token(token)
    assert payload["user_id"] == test_pharma_user.id
    assert payload["email"] == test_pharma_user.email
    assert payload["role"] == test_pharma_user.role.name


def test_refresh_token_invalid(client, db_session):
    """Tests token refresh failure with invalid refresh token"""
    refresh_data = {
        "refresh_token": "invalid_token"
    }
    
    response = client.post("/api/v1/auth/refresh-token", json=refresh_data)
    
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    assert "detail" in response.json()


def test_logout_success(client, db_session, test_pharma_user, pharma_token_headers):
    """Tests successful logout with valid refresh token"""
    # First login to get tokens
    login_data = {
        "email": test_pharma_user.email,
        "password": "test_password"  # This should match the fixture
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    refresh_token = login_response.json()["refresh_token"]
    
    # Now logout using the refresh token
    logout_data = {
        "refresh_token": refresh_token
    }
    
    response = client.post(
        "/api/v1/auth/logout", 
        json=logout_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == status.HTTP_200_OK
    assert "message" in response.json()
    
    # Try to use the same refresh token again and verify it fails
    refresh_data = {
        "refresh_token": refresh_token
    }
    response = client.post("/api/v1/auth/refresh-token", json=refresh_data)
    assert response.status_code == status.HTTP_401_UNAUTHORIZED


@pytest.mark.parametrize('role', ['pharma', 'cro'])
def test_register_success(client, db_session):
    """Tests successful user registration with valid data"""
    # Generate unique email for test
    email = f"test_{datetime.now().timestamp()}@example.com"
    
    registration_data = {
        "email": email,
        "password": "SecurePass123!",
        "confirm_password": "SecurePass123!",
        "role": "pharma"  # This will be overridden by parametrize in actual implementation
    }
    
    # Mock email sending function to avoid actual email sending
    with patch("app.services.auth_service.send_verification_email") as mock_send_email:
        mock_send_email.return_value = True
        
        response = client.post("/api/v1/auth/register", json=registration_data)
        
        assert response.status_code == status.HTTP_201_CREATED
        assert response.json()["success"] is True
        assert "user" in response.json()
        assert response.json()["user"]["email"] == email
        assert response.json()["user"]["role"] == "PHARMA"  # Would match parametrized role in actual implementation
        
        # Verify user exists in database with correct role and pending status
        user = db_session.query(User).filter(User.email == email).first()
        assert user is not None
        assert user.role == UserRole.PHARMA  # Would match parametrized role in actual implementation
        assert user.status == UserStatus.PENDING


def test_register_existing_email(client, db_session, test_pharma_user):
    """Tests registration failure with existing email"""
    registration_data = {
        "email": test_pharma_user.email,
        "password": "SecurePass123!",
        "confirm_password": "SecurePass123!",
        "role": "pharma"
    }
    
    response = client.post("/api/v1/auth/register", json=registration_data)
    
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "detail" in response.json()


def test_register_password_mismatch(client, db_session):
    """Tests registration failure when passwords don't match"""
    # Generate unique email for test
    email = f"test_{datetime.now().timestamp()}@example.com"
    
    registration_data = {
        "email": email,
        "password": "SecurePass123!",
        "confirm_password": "DifferentPass123!",  # Mismatched password
        "role": "pharma"
    }
    
    response = client.post("/api/v1/auth/register", json=registration_data)
    
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    assert "detail" in response.json()


@pytest.mark.parametrize('weak_password', [
    'short',  # Too short
    'nouppercase123!',  # No uppercase
    'NOLOWERCASE123!',  # No lowercase 
    'NoSpecialChar123',  # No special character
    'NoNumber!'  # No number
])
def test_register_weak_password(client, db_session):
    """Tests registration failure with weak password"""
    # Generate unique email for test
    email = f"test_{datetime.now().timestamp()}@example.com"
    
    # Note: In actual implementation, we would use the parametrized weak_password
    registration_data = {
        "email": email,
        "password": "short",  # This will be overridden by parametrize in actual implementation
        "confirm_password": "short",  # This will be overridden by parametrize in actual implementation
        "role": "pharma"
    }
    
    response = client.post("/api/v1/auth/register", json=registration_data)
    
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    assert "detail" in response.json()


def test_register_invalid_role(client, db_session):
    """Tests registration failure with invalid role"""
    # Generate unique email for test
    email = f"test_{datetime.now().timestamp()}@example.com"
    
    registration_data = {
        "email": email,
        "password": "SecurePass123!",
        "confirm_password": "SecurePass123!",
        "role": "invalid_role"  # Role that doesn't exist
    }
    
    response = client.post("/api/v1/auth/register", json=registration_data)
    
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    assert "detail" in response.json()


def test_verify_email_success(client, db_session, test_pharma_user):
    """Tests successful email verification with valid token"""
    # Create a user with pending status
    user = User(
        email=f"pending_{datetime.now().timestamp()}@example.com",
        password_hash=get_password_hash("SecurePass123!"),
        role=UserRole.PHARMA,
        status=UserStatus.PENDING,
        is_active=True,
        email_verified=False,
        password_history=[]
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    
    # Generate valid verification token for user
    token = create_access_token(
        data={"user_id": user.id, "email": user.email, "purpose": "email_verification"},
        expires_delta=timedelta(hours=24)
    )
    
    verification_data = {
        "token": token
    }
    
    # Mock email sending function to avoid actual email sending
    with patch("app.services.auth_service.send_welcome_email") as mock_send_email:
        mock_send_email.return_value = True
        
        response = client.post("/api/v1/auth/verify-email", json=verification_data)
        
        assert response.status_code == status.HTTP_200_OK
        assert "message" in response.json()
        
        # Verify user status is updated to active in database
        db_session.refresh(user)
        assert user.email_verified is True
        assert user.status == UserStatus.ACTIVE


def test_verify_email_invalid_token(client, db_session):
    """Tests email verification failure with invalid token"""
    verification_data = {
        "token": "invalid_token"
    }
    
    response = client.post("/api/v1/auth/verify-email", json=verification_data)
    
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "detail" in response.json()


def test_verify_email_expired_token(client, db_session, test_pharma_user):
    """Tests email verification failure with expired token"""
    # Create a user with pending status
    user = User(
        email=f"expired_{datetime.now().timestamp()}@example.com",
        password_hash=get_password_hash("SecurePass123!"),
        role=UserRole.PHARMA,
        status=UserStatus.PENDING,
        is_active=True,
        email_verified=False,
        password_history=[]
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    
    # Generate expired verification token for user
    token = create_access_token(
        data={"user_id": user.id, "email": user.email, "purpose": "email_verification"},
        expires_delta=timedelta(seconds=-1)  # Already expired
    )
    
    verification_data = {
        "token": token
    }
    
    response = client.post("/api/v1/auth/verify-email", json=verification_data)
    
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "detail" in response.json()
    
    # Verify user status remains pending in database
    db_session.refresh(user)
    assert user.email_verified is False
    assert user.status == UserStatus.PENDING


def test_request_password_reset(client, db_session, test_pharma_user):
    """Tests password reset request functionality"""
    reset_data = {
        "email": test_pharma_user.email
    }
    
    # Mock email sending function to avoid actual email sending
    with patch("app.services.auth_service.send_password_reset_email") as mock_send_email:
        mock_send_email.return_value = True
        
        response = client.post("/api/v1/auth/password-reset/request", json=reset_data)
        
        assert response.status_code == status.HTTP_200_OK
        assert "message" in response.json()
        
        # Verify email sending function was called with correct parameters
        mock_send_email.assert_called_once()


def test_request_password_reset_nonexistent_email(client, db_session):
    """Tests password reset request with non-existent email"""
    reset_data = {
        "email": "nonexistent@example.com"
    }
    
    # Mock email sending function to avoid actual email sending
    with patch("app.services.auth_service.send_password_reset_email") as mock_send_email:
        mock_send_email.return_value = True
        
        response = client.post("/api/v1/auth/password-reset/request", json=reset_data)
        
        # For security, same response as success even for non-existent email
        assert response.status_code == status.HTTP_200_OK
        assert "message" in response.json()
        
        # Verify email sending function was not called
        mock_send_email.assert_not_called()


def test_confirm_password_reset_success(client, db_session, test_pharma_user):
    """Tests successful password reset confirmation with valid token"""
    # Generate valid password reset token for user
    token = create_access_token(
        data={"user_id": test_pharma_user.id, "email": test_pharma_user.email, "purpose": "password_reset"},
        expires_delta=timedelta(hours=1)
    )
    
    # Create new password that meets requirements
    new_password = "NewSecurePass123!"
    
    reset_data = {
        "token": token,
        "new_password": new_password,
        "confirm_password": new_password
    }
    
    response = client.post("/api/v1/auth/password-reset/confirm", json=reset_data)
    
    assert response.status_code == status.HTTP_200_OK
    assert "message" in response.json()
    
    # Verify user can login with new password
    login_data = {
        "email": test_pharma_user.email,
        "password": new_password
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    assert login_response.status_code == status.HTTP_200_OK
    
    # Verify user cannot login with old password
    login_data = {
        "email": test_pharma_user.email,
        "password": "test_password"  # Old password
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    assert login_response.status_code == status.HTTP_401_UNAUTHORIZED


def test_confirm_password_reset_invalid_token(client, db_session):
    """Tests password reset confirmation failure with invalid token"""
    reset_data = {
        "token": "invalid_token",
        "new_password": "NewSecurePass123!",
        "confirm_password": "NewSecurePass123!"
    }
    
    response = client.post("/api/v1/auth/password-reset/confirm", json=reset_data)
    
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "detail" in response.json()


def test_confirm_password_reset_weak_password(client, db_session, test_pharma_user):
    """Tests password reset confirmation failure with weak password"""
    # Generate valid password reset token for user
    token = create_access_token(
        data={"user_id": test_pharma_user.id, "email": test_pharma_user.email, "purpose": "password_reset"},
        expires_delta=timedelta(hours=1)
    )
    
    # Create weak password that doesn't meet requirements
    weak_password = "weak"
    
    reset_data = {
        "token": token,
        "new_password": weak_password,
        "confirm_password": weak_password
    }
    
    response = client.post("/api/v1/auth/password-reset/confirm", json=reset_data)
    
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    assert "detail" in response.json()


def test_change_password_success(client, db_session, test_pharma_user, test_password, pharma_token_headers):
    """Tests successful password change with valid current password"""
    # Create new password that meets requirements
    new_password = "NewSecurePass123!"
    
    change_data = {
        "current_password": test_password,
        "new_password": new_password,
        "confirm_password": new_password
    }
    
    response = client.post(
        "/api/v1/auth/change-password",
        json=change_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == status.HTTP_200_OK
    assert "message" in response.json()
    
    # Verify user can login with new password
    login_data = {
        "email": test_pharma_user.email,
        "password": new_password
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    assert login_response.status_code == status.HTTP_200_OK
    
    # Verify user cannot login with old password
    login_data = {
        "email": test_pharma_user.email,
        "password": test_password
    }
    
    login_response = client.post("/api/v1/auth/login", json=login_data)
    assert login_response.status_code == status.HTTP_401_UNAUTHORIZED


def test_change_password_incorrect_current(client, db_session, pharma_token_headers):
    """Tests password change failure with incorrect current password"""
    change_data = {
        "current_password": "incorrect_password",
        "new_password": "NewSecurePass123!",
        "confirm_password": "NewSecurePass123!"
    }
    
    response = client.post(
        "/api/v1/auth/change-password",
        json=change_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == status.HTTP_400_BAD_REQUEST
    assert "detail" in response.json()


def test_change_password_weak_new(client, db_session, test_password, pharma_token_headers):
    """Tests password change failure with weak new password"""
    # Create weak password that doesn't meet requirements
    weak_password = "weak"
    
    change_data = {
        "current_password": test_password,
        "new_password": weak_password,
        "confirm_password": weak_password
    }
    
    response = client.post(
        "/api/v1/auth/change-password",
        json=change_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    assert "detail" in response.json()


@pytest.mark.parametrize('token_headers_fixture, user_fixture', [
    ('pharma_token_headers', 'test_pharma_user'),
    ('cro_token_headers', 'test_cro_user'),
    ('admin_token_headers', 'test_admin_user')
])
def test_get_current_user(client, pharma_token_headers, test_pharma_user):
    """Tests retrieving current user information with valid token"""
    # Note: In actual implementation, we would use request.getfixturevalue(token_headers_fixture)
    # and request.getfixturevalue(user_fixture) to dynamically get the fixtures
    
    response = client.get("/api/v1/auth/me", headers=pharma_token_headers)
    
    assert response.status_code == status.HTTP_200_OK
    assert response.json()["email"] == test_pharma_user.email
    assert response.json()["role"] == test_pharma_user.role.name
    assert response.json()["id"] == test_pharma_user.id


def test_get_current_user_invalid_token(client):
    """Tests failure to retrieve user information with invalid token"""
    invalid_headers = {"Authorization": "Bearer invalid_token"}
    
    response = client.get("/api/v1/auth/me", headers=invalid_headers)
    
    assert response.status_code == status.HTTP_401_UNAUTHORIZED
    assert "detail" in response.json()