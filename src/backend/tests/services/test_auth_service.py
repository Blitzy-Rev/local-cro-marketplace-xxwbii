import pytest
from unittest.mock import MagicMock, patch
from datetime import datetime, timedelta

from ../../app/services/auth_service import (
    authenticate_user,
    register_user,
    verify_email_token,
    create_password_reset_token,
    reset_password,
    refresh_access_token,
    change_password,
    get_current_user,
    create_tokens_for_user
)
from ../../app/crud/crud_user import user
from ../../app/core/jwt import verify_token, create_access_token, create_refresh_token
from ../../app/core/security import verify_password, get_password_hash
from ../../app/schemas/user import UserCreate
from ../../app/models/user import User
from ../../app/exceptions import AuthenticationException, ValidationException, ResourceNotFoundException
from ../../app/constants import UserRole, UserStatus


@pytest.mark.parametrize('role', ['pharma', 'cro', 'admin'])
def test_authenticate_user_success(db_session, role):
    """Test successful user authentication with valid credentials."""
    # Create a test user with the given role
    email = "test@example.com"
    password = "SecureP@ssw0rd123"
    password_hash = get_password_hash(password)
    
    mock_user = MagicMock()
    mock_user.id = 1
    mock_user.email = email
    mock_user.role = getattr(UserRole, role.upper())
    mock_user.status = UserStatus.ACTIVE
    mock_user.password_hash = password_hash
    
    # Mock the necessary functions
    with patch.object(user, 'authenticate', return_value=mock_user) as mock_authenticate, \
         patch.object(user, 'update_last_login', return_value=mock_user) as mock_update_last_login, \
         patch('../../app/core/jwt.create_token_for_user', return_value="access_token") as mock_create_access, \
         patch('../../app/core/jwt.create_refresh_token_for_user', return_value="refresh_token") as mock_create_refresh:
        
        # Call authenticate_user with valid credentials
        result_user, access_token, refresh_token = authenticate_user(db_session, email, password)
        
        # Assert that the returned user matches the created user
        assert result_user == mock_user
        # Assert that access and refresh tokens are returned
        assert access_token == "access_token"
        assert refresh_token == "refresh_token"
        # Assert that user's last login timestamp is updated
        mock_update_last_login.assert_called_once_with(db_session, mock_user)


def test_authenticate_user_invalid_email(db_session):
    """Test authentication failure with non-existent email."""
    # Call authenticate_user with non-existent email
    with patch.object(user, 'authenticate', return_value=None):
        with pytest.raises(AuthenticationException) as exc_info:
            authenticate_user(db_session, "nonexistent@example.com", "anypassword")
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid email or password" in str(exc_info.value)


def test_authenticate_user_invalid_password(db_session):
    """Test authentication failure with incorrect password."""
    # Create a test user
    email = "test@example.com"
    
    # Mock authenticate to return None (password not valid)
    with patch.object(user, 'authenticate', return_value=None):
        with pytest.raises(AuthenticationException) as exc_info:
            authenticate_user(db_session, email, "wrongpassword")
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid email or password" in str(exc_info.value)


def test_authenticate_user_inactive(db_session):
    """Test authentication failure with inactive user."""
    # Create a test user with is_active=False
    email = "inactive@example.com"
    password = "SecureP@ssw0rd123"
    
    mock_user = MagicMock()
    mock_user.email = email
    mock_user.status = UserStatus.PENDING
    
    # Call authenticate_user with valid credentials for inactive user
    with patch.object(user, 'authenticate', return_value=mock_user):
        with pytest.raises(AuthenticationException) as exc_info:
            authenticate_user(db_session, email, password)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "User account is not active" in str(exc_info.value)


@pytest.mark.parametrize('role', ['pharma', 'cro', 'admin'])
def test_register_user_success(db_session, role):
    """Test successful user registration with valid data."""
    # Create UserCreate object with valid data
    user_data = UserCreate(
        email="newuser@example.com",
        password="SecureP@ssw0rd123",
        role=role
    )
    base_url = "http://example.com"
    
    # Create a mock user object that would be returned after creation
    mock_user = MagicMock()
    mock_user.id = 1
    mock_user.email = user_data.email
    mock_user.role = getattr(UserRole, role.upper())
    mock_user.status = UserStatus.PENDING
    mock_user.email_verified = False
    
    # Mock the necessary functions
    with patch.object(user, 'get_by_email', return_value=None) as mock_get_by_email, \
         patch.object(user, 'create_with_password', return_value=mock_user) as mock_create, \
         patch('../../app/services/auth_service.generate_verification_token', return_value="token123") as mock_gen_token, \
         patch('../../app/services/auth_service.is_email_enabled', return_value=True) as mock_email_enabled, \
         patch('../../app/utils/email_utils.send_verification_email', return_value=True) as mock_send_email:
        
        # Call register_user with the user data
        result = register_user(db_session, user_data, base_url)
        
        # Assert that a user is created with correct attributes
        mock_create.assert_called_once_with(db_session, obj_in=user_data)
        # Assert that password is properly hashed (validated by create_with_password call)
        # Assert that user status is set to PENDING
        assert mock_user.status == UserStatus.PENDING
        # Assert that email_verified is False
        assert mock_user.email_verified is False
        assert result == mock_user


def test_register_user_existing_email(db_session):
    """Test registration failure with existing email."""
    # Create a test user
    existing_email = "existing@example.com"
    
    # Create UserCreate object with the same email
    user_data = UserCreate(
        email=existing_email,
        password="SecureP@ssw0rd123",
        role="pharma"
    )
    base_url = "http://example.com"
    
    # Mock existing user
    mock_existing_user = MagicMock()
    mock_existing_user.email = existing_email
    
    # Call register_user with the user data
    with patch.object(user, 'get_by_email', return_value=mock_existing_user):
        with pytest.raises(ValidationException) as exc_info:
            register_user(db_session, user_data, base_url)
        
        # Assert that ValidationException is raised with appropriate message
        assert "Email already registered" in str(exc_info.value)


@pytest.mark.parametrize('password,expected_errors', [
    ('short', ['Password must be at least 10 characters']),
    ('nouppercase123!', ['Password must contain at least one uppercase letter']),
    ('NOLOWERCASE123!', ['Password must contain at least one lowercase letter']),
    ('NoDigits!', ['Password must contain at least one digit']),
    ('NoSpecial123', ['Password must contain at least one special character'])
])
def test_register_user_weak_password(db_session, password, expected_errors):
    """Test registration failure with weak password."""
    # Create UserCreate object with weak password
    user_data = UserCreate(
        email="newuser@example.com",
        password=password,
        role="pharma"
    )
    base_url = "http://example.com"
    
    # Call register_user with the user data
    with patch.object(user, 'get_by_email', return_value=None), \
         patch('../../app/core/security.validate_password', return_value=False), \
         patch('../../app/core/security.get_validation_errors', return_value=expected_errors):
        
        with pytest.raises(ValidationException) as exc_info:
            register_user(db_session, user_data, base_url)
        
        # Assert that ValidationException is raised
        assert "Password does not meet security requirements" in str(exc_info.value)
        # Assert that error message contains expected validation errors
        assert exc_info.value.details["errors"] == expected_errors


def test_verify_email_token_success(db_session):
    """Test successful email verification with valid token."""
    # Create a test user with email_verified=False
    user_id = 1
    email = "user@example.com"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.email_verified = False
    mock_user.status = UserStatus.PENDING
    
    # Create a valid email verification token for the user
    token = "valid_token"
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "email_verification",
        "user_id": user_id
    }) as mock_verify_token, \
         patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch.object(user, 'verify_email', return_value=mock_user) as mock_verify_email, \
         patch.object(user, 'update', return_value=mock_user) as mock_update_user, \
         patch('../../app/services/auth_service.is_email_enabled', return_value=True) as mock_email_enabled, \
         patch('../../app/utils/email_utils.send_welcome_email', return_value=True) as mock_send_welcome:
        
        # Call verify_email_token with the token
        result = verify_email_token(db_session, token)
        
        # Assert that the user's email_verified is set to True
        mock_verify_email.assert_called_once_with(db_session, mock_user)
        # Assert that the user's status is updated to ACTIVE
        mock_update_user.assert_called_once()
        assert mock_update_user.call_args[1]['obj_in'] == {"status": UserStatus.ACTIVE}
        assert result == mock_user


def test_verify_email_token_invalid(db_session):
    """Test email verification failure with invalid token."""
    # Call verify_email_token with an invalid token
    token = "invalid_token"
    
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "not_email_verification",
        "user_id": 1
    }):
        with pytest.raises(AuthenticationException) as exc_info:
            verify_email_token(db_session, token)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid verification token" in str(exc_info.value)


def test_verify_email_token_expired(db_session):
    """Test email verification failure with expired token."""
    # Create a test user
    token = "expired_token"
    
    # Create an expired email verification token
    with patch('../../app/core/jwt.verify_token', side_effect=AuthenticationException("Token expired")):
        with pytest.raises(AuthenticationException) as exc_info:
            verify_email_token(db_session, token)
        
        # Assert that AuthenticationException is raised with 'Token expired' message
        assert "Token expired" in str(exc_info.value)


def test_verify_email_token_user_not_found(db_session):
    """Test email verification failure with non-existent user."""
    # Create a token with non-existent user_id
    token = "valid_token_nonexistent_user"
    non_existent_user_id = 999
    
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "email_verification",
        "user_id": non_existent_user_id
    }), \
         patch.object(user, 'get', return_value=None):
        
        with pytest.raises(ResourceNotFoundException) as exc_info:
            verify_email_token(db_session, token)
        
        # Assert that ResourceNotFoundException is raised with appropriate message
        assert "User not found" in str(exc_info.value)


def test_create_password_reset_token_success(db_session):
    """Test successful creation of password reset token."""
    # Create a test user
    email = "user@example.com"
    base_url = "http://example.com"
    
    mock_user = MagicMock()
    mock_user.id = 1
    mock_user.email = email
    
    # Mock email_utils.send_password_reset_email
    with patch.object(user, 'get_by_email', return_value=mock_user) as mock_get_by_email, \
         patch('../../app/services/auth_service.generate_password_reset_token', return_value="token123") as mock_gen_token, \
         patch('../../app/services/auth_service.is_email_enabled', return_value=True) as mock_email_enabled, \
         patch('../../app/utils/email_utils.send_password_reset_email', return_value=True) as mock_send_email:
        
        # Call create_password_reset_token with user's email
        result = create_password_reset_token(db_session, email, base_url)
        
        # Assert that the function returns True
        assert result is True
        # Assert that send_password_reset_email was called with correct arguments
        mock_send_email.assert_called_once_with(
            email=mock_user.email,
            username=mock_user.email.split('@')[0],
            token="token123",
            reset_url=f"{base_url}/reset-password"
        )


def test_create_password_reset_token_user_not_found(db_session):
    """Test password reset token creation failure with non-existent user."""
    # Call create_password_reset_token with non-existent email
    email = "nonexistent@example.com"
    base_url = "http://example.com"
    
    with patch.object(user, 'get_by_email', return_value=None):
        with pytest.raises(ResourceNotFoundException) as exc_info:
            create_password_reset_token(db_session, email, base_url)
        
        # Assert that ResourceNotFoundException is raised with appropriate message
        assert "User not found" in str(exc_info.value)


def test_reset_password_success(db_session):
    """Test successful password reset with valid token."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    old_password_hash = get_password_hash("OldP@ssw0rd")
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.password_hash = old_password_hash
    
    # Create a valid password reset token
    token = "valid_reset_token"
    new_password = "NewSecureP@ssw0rd123"
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "password_reset",
        "user_id": user_id
    }) as mock_verify_token, \
         patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch('../../app/core/security.validate_password', return_value=True) as mock_validate, \
         patch.object(user, 'update_with_password', return_value=mock_user) as mock_update_password:
        
        # Call reset_password with token and new password
        result = reset_password(db_session, token, new_password)
        
        # Assert that user's password is updated
        mock_update_password.assert_called_once_with(db_session, db_obj=mock_user, obj_in={"password": new_password})
        # Assert that new password works for authentication (validated by other functions)
        assert result == mock_user


def test_reset_password_invalid_token(db_session):
    """Test password reset failure with invalid token."""
    # Call reset_password with invalid token and new password
    token = "invalid_token"
    new_password = "NewSecureP@ssw0rd123"
    
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "not_password_reset",
        "user_id": 1
    }):
        with pytest.raises(AuthenticationException) as exc_info:
            reset_password(db_session, token, new_password)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid password reset token" in str(exc_info.value)


def test_reset_password_weak_password(db_session):
    """Test password reset failure with weak password."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    
    # Create a valid password reset token
    token = "valid_reset_token"
    weak_password = "weak"
    validation_errors = ["Password must be at least 10 characters"]
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.verify_token', return_value={
        "purpose": "password_reset",
        "user_id": user_id
    }), \
         patch.object(user, 'get', return_value=mock_user), \
         patch('../../app/core/security.validate_password', return_value=False), \
         patch('../../app/core/security.get_validation_errors', return_value=validation_errors):
        
        with pytest.raises(ValidationException) as exc_info:
            reset_password(db_session, token, weak_password)
        
        # Assert that ValidationException is raised with appropriate message
        assert "Password does not meet security requirements" in str(exc_info.value)
        assert exc_info.value.details["errors"] == validation_errors


def test_refresh_access_token_success(db_session):
    """Test successful access token refresh with valid refresh token."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    role = "pharma"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.is_active = True
    mock_user.status = UserStatus.ACTIVE
    
    # Create a valid refresh token for the user
    refresh_token = "valid_refresh_token"
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.verify_token', return_value={
        "token_type": "refresh",
        "user_id": user_id,
        "email": email,
        "role": role
    }) as mock_verify_token, \
         patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch('../../app/core/jwt.create_token_for_user', return_value="new_access_token") as mock_create_token:
        
        # Call refresh_access_token with the refresh token
        result = refresh_access_token(db_session, refresh_token)
        
        # Assert that a new access token is returned
        assert result == "new_access_token"
        # Verify that the new token contains correct user information
        mock_create_token.assert_called_once_with(
            user_id=user_id,
            email=email,
            role=role
        )


def test_refresh_access_token_invalid_token(db_session):
    """Test token refresh failure with invalid token."""
    # Call refresh_access_token with invalid token
    invalid_token = "invalid_token"
    
    with patch('../../app/core/jwt.verify_token', side_effect=AuthenticationException("Invalid token")):
        with pytest.raises(AuthenticationException) as exc_info:
            refresh_access_token(db_session, invalid_token)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid token" in str(exc_info.value)


def test_refresh_access_token_not_refresh_token(db_session):
    """Test token refresh failure with access token instead of refresh token."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    role = "pharma"
    
    # Create an access token (not refresh token)
    access_token = "access_token_not_refresh"
    
    with patch('../../app/core/jwt.verify_token', return_value={
        "token_type": "access",
        "user_id": user_id,
        "email": email,
        "role": role
    }):
        with pytest.raises(AuthenticationException) as exc_info:
            refresh_access_token(db_session, access_token)
        
        # Assert that AuthenticationException is raised with 'Not a refresh token' message
        assert "Invalid refresh token" in str(exc_info.value)


def test_change_password_success(db_session):
    """Test successful password change with valid credentials."""
    # Create a test user with known password
    user_id = 1
    email = "user@example.com"
    current_password = "CurrentP@ssw0rd"
    current_password_hash = get_password_hash(current_password)
    new_password = "NewSecureP@ssw0rd123"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.password_hash = current_password_hash
    
    # Mock the necessary functions
    with patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch('../../app/core/security.verify_password', return_value=True) as mock_verify, \
         patch('../../app/core/security.validate_password', return_value=True) as mock_validate, \
         patch.object(user, 'update_with_password', return_value=mock_user) as mock_update_password:
        
        # Call change_password with user_id, current password, and new password
        result = change_password(db_session, user_id, current_password, new_password)
        
        # Assert that user's password is updated
        mock_update_password.assert_called_once_with(db_session, db_obj=mock_user, obj_in={"password": new_password})
        # Assert that new password works for authentication (validated by other functions)
        assert result == mock_user
        # Assert that old password no longer works (handled by update_with_password)


def test_change_password_invalid_current(db_session):
    """Test password change failure with incorrect current password."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    current_password = "CurrentP@ssw0rd"
    wrong_current_password = "WrongP@ssw0rd"
    current_password_hash = get_password_hash(current_password)
    new_password = "NewSecureP@ssw0rd123"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.password_hash = current_password_hash
    
    # Call change_password with user_id, wrong current password, and new password
    with patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch('../../app/core/security.verify_password', return_value=False) as mock_verify:
        
        with pytest.raises(AuthenticationException) as exc_info:
            change_password(db_session, user_id, wrong_current_password, new_password)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Current password is incorrect" in str(exc_info.value)


def test_change_password_weak_new(db_session):
    """Test password change failure with weak new password."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    current_password = "CurrentP@ssw0rd"
    current_password_hash = get_password_hash(current_password)
    weak_new_password = "weak"
    validation_errors = ["Password must be at least 10 characters"]
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.password_hash = current_password_hash
    
    # Call change_password with user_id, correct current password, and weak new password
    with patch.object(user, 'get', return_value=mock_user) as mock_get_user, \
         patch('../../app/core/security.verify_password', return_value=True) as mock_verify, \
         patch('../../app/core/security.validate_password', return_value=False) as mock_validate, \
         patch('../../app/core/security.get_validation_errors', return_value=validation_errors) as mock_get_errors:
        
        with pytest.raises(ValidationException) as exc_info:
            change_password(db_session, user_id, current_password, weak_new_password)
        
        # Assert that ValidationException is raised with appropriate message
        assert "Password does not meet security requirements" in str(exc_info.value)
        assert exc_info.value.details["errors"] == validation_errors


def test_get_current_user_success(db_session):
    """Test successful retrieval of current user from valid token."""
    # Create a test user
    user_id = 1
    email = "user@example.com"
    
    mock_user = MagicMock()
    mock_user.id = user_id
    mock_user.email = email
    mock_user.is_active = True
    mock_user.status = UserStatus.ACTIVE
    
    # Create a valid access token for the user
    token = "valid_access_token"
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.verify_token', return_value={
        "token_type": "access",
        "user_id": user_id
    }) as mock_verify_token, \
         patch.object(user, 'get', return_value=mock_user) as mock_get_user:
        
        # Call get_current_user with the token
        result = get_current_user(db_session, token)
        
        # Assert that the correct user is returned
        assert result == mock_user
        mock_get_user.assert_called_once_with(db_session, id=user_id)


def test_get_current_user_invalid_token(db_session):
    """Test current user retrieval failure with invalid token."""
    # Call get_current_user with invalid token
    invalid_token = "invalid_token"
    
    with patch('../../app/core/jwt.verify_token', side_effect=AuthenticationException("Invalid token")):
        with pytest.raises(AuthenticationException) as exc_info:
            get_current_user(db_session, invalid_token)
        
        # Assert that AuthenticationException is raised with appropriate message
        assert "Invalid token" in str(exc_info.value)


def test_get_current_user_not_access_token(db_session):
    """Test current user retrieval failure with refresh token instead of access token."""
    # Create a test user
    user_id = 1
    
    # Create a refresh token (not access token)
    refresh_token = "refresh_token_not_access"
    
    with patch('../../app/core/jwt.verify_token', return_value={
        "token_type": "refresh",
        "user_id": user_id
    }):
        with pytest.raises(AuthenticationException) as exc_info:
            get_current_user(db_session, refresh_token)
        
        # Assert that AuthenticationException is raised with 'Not an access token' message
        assert "Invalid access token" in str(exc_info.value)


def test_create_tokens_for_user(db_session):
    """Test creation of access and refresh tokens for a user."""
    # Create a test user
    mock_user = MagicMock()
    mock_user.id = 1
    mock_user.email = "user@example.com"
    mock_user.role = UserRole.PHARMA
    
    # Mock the necessary functions
    with patch('../../app/core/jwt.create_token_for_user', return_value="access_token") as mock_create_access, \
         patch('../../app/core/jwt.create_refresh_token_for_user', return_value="refresh_token") as mock_create_refresh:
        
        # Call create_tokens_for_user with the user
        result = create_tokens_for_user(mock_user)
        
        # Assert that both access_token and refresh_token are returned
        assert result == {
            "access_token": "access_token",
            "refresh_token": "refresh_token"
        }
        
        # Verify that tokens contain correct user information
        mock_create_access.assert_called_once_with(
            user_id=mock_user.id,
            email=mock_user.email,
            role=mock_user.role.name
        )
        mock_create_refresh.assert_called_once_with(
            user_id=mock_user.id,
            email=mock_user.email,
            role=mock_user.role.name
        )