"""
Core authentication module that implements user authentication, token validation, and role-based 
access control for the Molecular Data Management and CRO Integration Platform. This module provides
a comprehensive authentication framework without external dependencies, supporting the platform's 
requirement for fully local deployment.
"""

from typing import Dict, List, Optional, Any, Union
from datetime import datetime, timedelta

from fastapi import Depends, HTTPException, status
from fastapi.security import SecurityScopes, OAuth2PasswordBearer
from sqlalchemy.orm import Session

from .config import get_settings
from .security import verify_password
from .jwt import (
    create_access_token,
    create_refresh_token,
    verify_token,
    get_token_data
)
from ..schemas.token import TokenData
from ..crud.crud_user import user
from ..exceptions import AuthenticationException, AuthorizationException
from ..constants import ROLE_PERMISSIONS, UserRole, API_V1_PREFIX

# OAuth2 password bearer scheme for token extraction
oauth2_scheme = OAuth2PasswordBearer(tokenUrl=f"{API_V1_PREFIX}/auth/login", auto_error=False)


def authenticate_user(db: Session, email: str, password: str) -> Dict[str, Any]:
    """
    Authenticates a user with email and password.
    
    Args:
        db: Database session
        email: User's email
        password: User's password
        
    Returns:
        User data dictionary if authentication succeeds
        
    Raises:
        AuthenticationException: If authentication fails
    """
    # Authenticate user with provided credentials
    authenticated_user = user.authenticate(db, email, password)
    
    # If authentication fails, raise exception
    if not authenticated_user:
        raise AuthenticationException("Invalid email or password")
    
    # Check if user is active
    if not user.is_active(authenticated_user):
        raise AuthenticationException("Inactive user")
    
    # Update last login timestamp
    user.update_last_login(db, authenticated_user)
    
    # Return user data as dictionary
    return {
        "user_id": authenticated_user.id,
        "email": authenticated_user.email,
        "role": authenticated_user.role.name  # Convert Enum to string
    }


def create_user_tokens(user_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Creates access and refresh tokens for a user.
    
    Args:
        user_data: User data dictionary with user_id, email, and role
        
    Returns:
        Dictionary containing access and refresh tokens
    """
    # Extract user data
    user_id = user_data["user_id"]
    email = user_data["email"]
    role = user_data["role"]
    
    # Create token data
    token_data = {
        "user_id": user_id,
        "email": email,
        "role": role
    }
    
    # Create access and refresh tokens
    access_token = create_access_token(token_data)
    refresh_token = create_refresh_token(token_data)
    
    # Return tokens
    return {
        "access_token": access_token,
        "refresh_token": refresh_token
    }


def get_current_user(
    security_scopes: SecurityScopes,
    token: str = Depends(oauth2_scheme)
) -> Dict[str, Any]:
    """
    Extracts and validates the current user from a JWT token.
    
    Args:
        security_scopes: Security scopes for OAuth2
        token: JWT token string
        
    Returns:
        Current user information from token
        
    Raises:
        AuthenticationException: If token is invalid or user is not found
    """
    if token is None:
        raise AuthenticationException("Not authenticated")
    
    try:
        # Validate token and extract data
        token_data = get_token_data(token)
        
        # Check if token contains required scopes
        if security_scopes.scopes:
            # Since our token doesn't contain scopes directly, we would need to
            # map user roles to scopes or check permissions here
            raise AuthenticationException(
                "Not enough permissions",
                {"required": security_scopes.scopes, "provided": []}
            )
        
        # Return user data from token
        return {
            "user_id": token_data.user_id,
            "email": token_data.email,
            "role": token_data.role
        }
    except Exception as e:
        raise AuthenticationException(f"Invalid authentication credentials: {str(e)}")


def get_current_active_user(current_user: Dict[str, Any], db: Session) -> Dict[str, Any]:
    """
    Ensures the current user is active.
    
    Args:
        current_user: Current user information from token
        db: Database session
        
    Returns:
        Current active user information
        
    Raises:
        AuthenticationException: If user is not active
    """
    # Get user from database
    db_user = user.get(db, id=current_user["user_id"])
    
    # Check if user is active
    if not db_user or not user.is_active(db_user):
        raise AuthenticationException("Inactive user")
    
    return current_user


def check_permissions(current_user: Dict[str, Any], required_permissions: List[str]) -> bool:
    """
    Checks if a user has the required permissions.
    
    Args:
        current_user: Current user information
        required_permissions: List of required permission strings
        
    Returns:
        True if user has all required permissions, False otherwise
    """
    # Extract role from user data
    role = current_user.get("role")
    
    # Get permissions for this role
    user_permissions = ROLE_PERMISSIONS.get(UserRole[role], [])
    
    # Check if user has all required permissions
    return all(permission in user_permissions for permission in required_permissions)


def require_permissions(required_permissions: List[str]):
    """
    Decorator that requires specific permissions for a function.
    
    Args:
        required_permissions: List of required permission strings
        
    Returns:
        Decorated function that checks permissions
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Get current_user from function arguments
            current_user = None
            for arg in args:
                if isinstance(arg, dict) and "user_id" in arg and "role" in arg:
                    current_user = arg
                    break
            
            if not current_user:
                for key, value in kwargs.items():
                    if key == "current_user" or (isinstance(value, dict) and "user_id" in value and "role" in value):
                        current_user = value
                        break
            
            if not current_user:
                raise AuthorizationException("User information not found in request")
            
            # Check permissions
            if not check_permissions(current_user, required_permissions):
                raise AuthorizationException(
                    "Not enough permissions",
                    {"required": required_permissions, "role": current_user.get("role")}
                )
            
            # If permissions check passes, call the function
            return func(*args, **kwargs)
        
        return wrapper
    
    return decorator


def refresh_access_token(refresh_token: str) -> Dict[str, Any]:
    """
    Creates a new access token from a valid refresh token.
    
    Args:
        refresh_token: Refresh token string
        
    Returns:
        Dictionary containing new access token
        
    Raises:
        AuthenticationException: If refresh token is invalid
    """
    # Verify refresh token
    token_data = verify_token(refresh_token)
    
    # Check if it's a refresh token
    if token_data.get("token_type") != "refresh":
        raise AuthenticationException("Invalid refresh token")
    
    # Extract user data
    user_id = token_data.get("user_id")
    email = token_data.get("email")
    role = token_data.get("role")
    
    # Create token data
    token_data = {
        "user_id": user_id,
        "email": email,
        "role": role
    }
    
    # Create new access token
    access_token = create_access_token(token_data)
    
    # Return new access token
    return {
        "access_token": access_token
    }


class RoleChecker:
    """
    Utility class for checking user roles.
    """
    
    @staticmethod
    def is_admin(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has admin role.
        
        Args:
            user: User information dictionary
            
        Returns:
            True if user is admin, False otherwise
        """
        return user.get("role") == UserRole.ADMIN.name
    
    @staticmethod
    def is_pharma(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has pharma role.
        
        Args:
            user: User information dictionary
            
        Returns:
            True if user is pharma, False otherwise
        """
        return user.get("role") == UserRole.PHARMA.name
    
    @staticmethod
    def is_cro(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has CRO role.
        
        Args:
            user: User information dictionary
            
        Returns:
            True if user is CRO, False otherwise
        """
        return user.get("role") == UserRole.CRO.name
    
    @staticmethod
    def has_role(user: Dict[str, Any], role: UserRole) -> bool:
        """
        Checks if a user has a specific role.
        
        Args:
            user: User information dictionary
            role: Role to check
            
        Returns:
            True if user has the specified role, False otherwise
        """
        return user.get("role") == role.name