"""
Central security module for the Molecular Data Management and CRO Integration Platform.

This module re-exports essential security functions from core modules while maintaining 
backward compatibility, ensuring secure authentication, authorization, and data protection 
without external dependencies.
"""

import warnings  # standard library
from typing import Dict, List, Optional, Any, Callable  # standard library

from fastapi import HTTPException, Depends  # version 0.95+
from fastapi.security import Security, HTTPBearer, HTTPAuthorizationCredentials  # version 0.95+

# Import core security functions
from .core.security import verify_password, get_password_hash
from .core.jwt import (
    create_access_token,
    create_refresh_token,
    decode_token,
    verify_token,
    oauth2_scheme,
    generate_password_reset_token,
    verify_password_reset_token,
    generate_email_verification_token,
    verify_email_verification_token,
)
from .core.auth import (
    get_current_user,
    get_role_permissions,
    check_permissions,
    require_permissions,
    ROLE_PERMISSIONS,
    RoleChecker,
)
from .utils.security_utils import (
    validate_password_strength,
    get_password_validation_errors,
    sanitize_input,
    mask_sensitive_data,
)
from .constants import UserRole
from .exceptions import AuthenticationException, AuthorizationException

# Track whether deprecation warning has been shown
_deprecation_warning_shown = False

def show_deprecation_warning():
    """
    Shows a deprecation warning for direct imports from this module.
    """
    global _deprecation_warning_shown
    if not _deprecation_warning_shown:
        warnings.warn(
            "Direct import from app.security is deprecated. "
            "Please import from app.core.security, app.core.jwt, app.core.auth, "
            "or app.utils.security_utils instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        _deprecation_warning_shown = True

def get_user_permissions(user: Dict[str, Any]) -> List[str]:
    """
    Gets the permissions for a user based on their role.
    
    Args:
        user: User data dictionary containing role information
        
    Returns:
        List of permission strings for the user's role
    """
    show_deprecation_warning()
    role = user.get("role")
    return get_role_permissions(role)

def has_permission(user: Dict[str, Any], permission: str) -> bool:
    """
    Checks if a user has a specific permission.
    
    Args:
        user: User data dictionary
        permission: Permission string to check
        
    Returns:
        True if user has the permission, False otherwise
    """
    show_deprecation_warning()
    permissions = get_user_permissions(user)
    return permission in permissions

def create_user_tokens(user_data: Dict[str, Any]) -> Dict[str, str]:
    """
    Creates access and refresh tokens for a user.
    
    Args:
        user_data: User data dictionary with user_id, email, and role
        
    Returns:
        Dictionary with access_token and refresh_token
    """
    show_deprecation_warning()
    user_id = user_data.get("user_id")
    email = user_data.get("email")
    role = user_data.get("role")
    
    token_data = {
        "user_id": user_id,
        "email": email,
        "role": role
    }
    
    access_token = create_access_token(token_data)
    refresh_token = create_refresh_token(token_data)
    
    return {
        "access_token": access_token,
        "refresh_token": refresh_token
    }

class SecurityUtils:
    """
    Utility class providing security helper methods.
    """
    
    @staticmethod
    def is_pharma_user(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has pharma role.
        
        Args:
            user: User data dictionary
            
        Returns:
            True if user has pharma role, False otherwise
        """
        if "role" not in user:
            return False
        return user["role"] == UserRole.PHARMA.name
    
    @staticmethod
    def is_cro_user(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has CRO role.
        
        Args:
            user: User data dictionary
            
        Returns:
            True if user has CRO role, False otherwise
        """
        if "role" not in user:
            return False
        return user["role"] == UserRole.CRO.name
    
    @staticmethod
    def is_admin_user(user: Dict[str, Any]) -> bool:
        """
        Checks if a user has admin role.
        
        Args:
            user: User data dictionary
            
        Returns:
            True if user has admin role, False otherwise
        """
        if "role" not in user:
            return False
        return user["role"] == UserRole.ADMIN.name