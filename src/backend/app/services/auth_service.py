"""
Service module that implements authentication and user management functionality for the 
Molecular Data Management and CRO Integration Platform. This module provides methods for 
user registration, authentication, token management, password reset, and email verification 
without external dependencies, supporting the platform's requirement for fully local deployment.
"""

from datetime import datetime, timedelta
from typing import Optional, Dict, Any, Tuple, List
import uuid
import secrets

from sqlalchemy.orm import Session

from ..crud.crud_user import user
from ..core.security import (
    get_password_hash,
    verify_password,
    validate_password,
    get_validation_errors
)
from ..core.jwt import (
    create_access_token,
    create_refresh_token,
    verify_token,
    create_token_for_user,
    create_refresh_token_for_user
)
from ..db.session import get_db
from ..models.user import User
from ..schemas.user import UserCreate, UserUpdate, UserResponse
from ..schemas.token import TokenResponse, TokenData
from ..exceptions import AuthenticationException, ValidationException, ResourceNotFoundException
from ..utils.email_utils import (
    send_verification_email,
    send_password_reset_email,
    send_welcome_email,
    is_email_enabled
)
from ..constants import (
    UserRole,
    UserStatus,
    EMAIL_VERIFICATION_TOKEN_EXPIRY_HOURS,
    PASSWORD_RESET_TOKEN_EXPIRY_HOURS
)


def authenticate_user(db: Session, email: str, password: str) -> Tuple[User, str, str]:
    """
    Authenticates a user with email and password, returning user and tokens if successful.
    
    Args:
        db: Database session
        email: User's email address
        password: User's password
        
    Returns:
        Tuple containing user object, access token, and refresh token
        
    Raises:
        AuthenticationException: If authentication fails
    """
    # Authenticate user with credentials
    authenticated_user = user.authenticate(db, email=email, password=password)
    if not authenticated_user:
        raise AuthenticationException("Invalid email or password")
    
    # Ensure user is active
    if authenticated_user.status != UserStatus.ACTIVE:
        raise AuthenticationException("User account is not active")
    
    # Update last login timestamp
    authenticated_user = user.update_last_login(db, authenticated_user)
    
    # Create tokens
    access_token = create_token_for_user(
        user_id=authenticated_user.id,
        email=authenticated_user.email,
        role=authenticated_user.role.name
    )
    
    refresh_token = create_refresh_token_for_user(
        user_id=authenticated_user.id,
        email=authenticated_user.email,
        role=authenticated_user.role.name
    )
    
    return authenticated_user, access_token, refresh_token


def register_user(db: Session, user_data: UserCreate, base_url: str) -> User:
    """
    Registers a new user with validation and sends verification email.
    
    Args:
        db: Database session
        user_data: User creation data
        base_url: Base URL for verification link
        
    Returns:
        Created user object
        
    Raises:
        ValidationException: If email already exists or password is invalid
    """
    # Check if user with this email already exists
    existing_user = user.get_by_email(db, email=user_data.email)
    if existing_user:
        raise ValidationException("Email already registered")
    
    # Validate password
    if not validate_password(user_data.password):
        errors = get_validation_errors(user_data.password)
        raise ValidationException(
            "Password does not meet security requirements",
            {"errors": errors}
        )
    
    # Create user with hashed password
    new_user = user.create_with_password(db, obj_in=user_data)
    
    # Generate verification token
    token = generate_verification_token(new_user)
    
    # Send verification email if email is enabled
    if is_email_enabled():
        send_verification_email(
            email=new_user.email,
            username=new_user.email.split('@')[0],  # Use part before @ as username
            token=token,
            verification_url=f"{base_url}/verify-email"
        )
    
    return new_user


def verify_email_token(db: Session, token: str) -> User:
    """
    Verifies an email verification token and activates the user account.
    
    Args:
        db: Database session
        token: Email verification token
        
    Returns:
        Verified user object
        
    Raises:
        AuthenticationException: If token is invalid
        ResourceNotFoundException: If user not found
    """
    # Verify token
    payload = verify_token(token)
    
    # Check token purpose
    if payload.get("purpose") != "email_verification":
        raise AuthenticationException("Invalid verification token")
    
    # Get user ID from token
    user_id = payload.get("user_id")
    if not user_id:
        raise AuthenticationException("Invalid token payload")
    
    # Get user by ID
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise ResourceNotFoundException("User not found")
    
    # If already verified, just return the user
    if db_user.email_verified:
        return db_user
    
    # Mark email as verified
    verified_user = user.verify_email(db, db_user)
    
    # Update user status to active
    if verified_user.status == UserStatus.PENDING:
        update_data = {"status": UserStatus.ACTIVE}
        verified_user = user.update(db, db_obj=verified_user, obj_in=update_data)
    
    # Send welcome email if email is enabled
    if is_email_enabled():
        send_welcome_email(
            email=verified_user.email,
            username=verified_user.email.split('@')[0]  # Use part before @ as username
        )
    
    return verified_user


def create_password_reset_token(db: Session, email: str, base_url: str) -> bool:
    """
    Creates a password reset token for a user and sends reset email.
    
    Args:
        db: Database session
        email: User's email address
        base_url: Base URL for reset link
        
    Returns:
        True if token was created and email sent, False otherwise
        
    Raises:
        ResourceNotFoundException: If user not found
    """
    # Get user by email
    db_user = user.get_by_email(db, email=email)
    if not db_user:
        raise ResourceNotFoundException("User not found")
    
    # Generate password reset token
    token = generate_password_reset_token(db_user)
    
    # Send password reset email if email is enabled
    if is_email_enabled():
        success = send_password_reset_email(
            email=db_user.email,
            username=db_user.email.split('@')[0],  # Use part before @ as username
            token=token,
            reset_url=f"{base_url}/reset-password"
        )
        return success
    
    return False


def reset_password(db: Session, token: str, new_password: str) -> User:
    """
    Resets a user's password using a valid reset token.
    
    Args:
        db: Database session
        token: Password reset token
        new_password: New password
        
    Returns:
        User object with updated password
        
    Raises:
        AuthenticationException: If token is invalid
        ValidationException: If password is invalid
        ResourceNotFoundException: If user not found
    """
    # Verify token
    payload = verify_token(token)
    
    # Check token purpose
    if payload.get("purpose") != "password_reset":
        raise AuthenticationException("Invalid password reset token")
    
    # Get user ID from token
    user_id = payload.get("user_id")
    if not user_id:
        raise AuthenticationException("Invalid token payload")
    
    # Get user by ID
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise ResourceNotFoundException("User not found")
    
    # Validate new password
    if not validate_password(new_password):
        errors = get_validation_errors(new_password)
        raise ValidationException(
            "Password does not meet security requirements",
            {"errors": errors}
        )
    
    # Update password
    update_data = {"password": new_password}
    updated_user = user.update_with_password(db, db_obj=db_user, obj_in=update_data)
    
    return updated_user


def refresh_access_token(db: Session, refresh_token: str) -> str:
    """
    Creates a new access token using a valid refresh token.
    
    Args:
        db: Database session
        refresh_token: Refresh token
        
    Returns:
        New access token
        
    Raises:
        AuthenticationException: If token is invalid or user is not active
    """
    # Verify refresh token
    payload = verify_token(refresh_token)
    
    # Check token type
    if payload.get("token_type") != "refresh":
        raise AuthenticationException("Invalid refresh token")
    
    # Get user details from token
    user_id = payload.get("user_id")
    email = payload.get("email")
    role = payload.get("role")
    
    if not user_id or not email or not role:
        raise AuthenticationException("Invalid token payload")
    
    # Get user by ID to ensure they still exist and are active
    db_user = user.get(db, id=user_id)
    if not db_user or not db_user.is_active or db_user.status != UserStatus.ACTIVE:
        raise AuthenticationException("User not active or not found")
    
    # Create new access token
    access_token = create_token_for_user(
        user_id=user_id,
        email=email,
        role=role
    )
    
    return access_token


def change_password(db: Session, user_id: int, current_password: str, new_password: str) -> User:
    """
    Changes a user's password after verifying current password.
    
    Args:
        db: Database session
        user_id: User ID
        current_password: Current password
        new_password: New password
        
    Returns:
        User object with updated password
        
    Raises:
        ResourceNotFoundException: If user not found
        AuthenticationException: If current password is invalid
        ValidationException: If new password is invalid
    """
    # Get user by ID
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise ResourceNotFoundException("User not found")
    
    # Verify current password
    if not verify_password(current_password, db_user.password_hash):
        raise AuthenticationException("Current password is incorrect")
    
    # Validate new password
    if not validate_password(new_password):
        errors = get_validation_errors(new_password)
        raise ValidationException(
            "Password does not meet security requirements",
            {"errors": errors}
        )
    
    # Update password
    update_data = {"password": new_password}
    updated_user = user.update_with_password(db, db_obj=db_user, obj_in=update_data)
    
    return updated_user


def generate_verification_token(user: User) -> str:
    """
    Generates an email verification token for a user.
    
    Args:
        user: User object
        
    Returns:
        Email verification token
    """
    # Create token payload
    payload = {
        "user_id": user.id,
        "email": user.email,
        "purpose": "email_verification"
    }
    
    # Set expiration time
    expires_delta = timedelta(hours=EMAIL_VERIFICATION_TOKEN_EXPIRY_HOURS)
    
    # Create token
    token = create_access_token(payload, expires_delta)
    
    return token


def generate_password_reset_token(user: User) -> str:
    """
    Generates a password reset token for a user.
    
    Args:
        user: User object
        
    Returns:
        Password reset token
    """
    # Create token payload
    payload = {
        "user_id": user.id,
        "email": user.email,
        "purpose": "password_reset"
    }
    
    # Set expiration time
    expires_delta = timedelta(hours=PASSWORD_RESET_TOKEN_EXPIRY_HOURS)
    
    # Create token
    token = create_access_token(payload, expires_delta)
    
    return token


def get_current_user(db: Session, token: str) -> User:
    """
    Gets the current user from a valid access token.
    
    Args:
        db: Database session
        token: Access token
        
    Returns:
        Current user object
        
    Raises:
        AuthenticationException: If token is invalid or user is not active
        ResourceNotFoundException: If user not found
    """
    # Verify token
    payload = verify_token(token)
    
    # Check token type
    if payload.get("token_type") != "access":
        raise AuthenticationException("Invalid access token")
    
    # Get user ID from token
    user_id = payload.get("user_id")
    if not user_id:
        raise AuthenticationException("Invalid token payload")
    
    # Get user by ID
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise ResourceNotFoundException("User not found")
    
    # Check if user is active
    if not db_user.is_active or db_user.status != UserStatus.ACTIVE:
        raise AuthenticationException("User account is not active")
    
    return db_user


def create_tokens_for_user(user: User) -> Dict[str, str]:
    """
    Creates access and refresh tokens for a user.
    
    Args:
        user: User object
        
    Returns:
        Dictionary containing access_token and refresh_token
    """
    # Create access token
    access_token = create_token_for_user(
        user_id=user.id,
        email=user.email,
        role=user.role.name
    )
    
    # Create refresh token
    refresh_token = create_refresh_token_for_user(
        user_id=user.id,
        email=user.email,
        role=user.role.name
    )
    
    # Return both tokens
    return {
        "access_token": access_token,
        "refresh_token": refresh_token
    }