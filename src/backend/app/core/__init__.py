"""
Core initialization module for the Molecular Data Management and CRO Integration Platform.

This module exports essential components for authentication, configuration, security,
and JWT handling, providing a clean interface for importing core functionality
throughout the application. It implements security mechanisms without external
dependencies, supporting the platform's requirement for fully local deployment.
"""

# Configuration components
from .config import (
    get_settings,
    Settings,
    get_database_url,
    get_redis_settings,
    get_minio_settings,
)

# Security components
from .security import (
    verify_password,
    get_password_hash,
    validate_password,
    PasswordPolicy,
)

# JWT handling components
from .jwt import (
    create_access_token,
    create_refresh_token,
    verify_token,
    get_token_data,
)

# Authentication and authorization components
from .auth import (
    authenticate_user,
    get_current_user,
    get_current_active_user,
    check_permissions,
    require_permissions,
    RoleChecker,
)

__all__ = [
    # Configuration
    "get_settings",
    "Settings",
    "get_database_url",
    "get_redis_settings",
    "get_minio_settings",
    
    # Security
    "verify_password",
    "get_password_hash",
    "validate_password",
    "PasswordPolicy",
    
    # JWT
    "create_access_token",
    "create_refresh_token",
    "verify_token",
    "get_token_data",
    
    # Authentication and Authorization
    "authenticate_user",
    "get_current_user",
    "get_current_active_user",
    "check_permissions",
    "require_permissions",
    "RoleChecker",
]