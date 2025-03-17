"""
Configuration module for the Molecular Data Management and CRO Integration Platform.

This module provides a simplified interface for accessing application configuration
settings throughout the system. It acts as a compatibility layer that redirects to 
the core configuration module while maintaining backward compatibility with existing imports.

This configuration approach supports fully local deployment without external dependencies
as required by the platform specifications.
"""

import warnings  # standard library

# Import functions and classes from core.config
from .core.config import (
    get_settings,
    initialize_settings,
    get_database_url,
    get_redis_settings,
    get_minio_settings,
    get_bucket_names
)

# Import Settings class for type annotations
from .core.settings import Settings

# Global variable to track if deprecation warning has been shown
_deprecation_warning_shown = False


def show_deprecation_warning() -> None:
    """
    Shows a deprecation warning for direct imports from this module.
    
    This function tracks whether the warning has already been shown to avoid
    repetitive warnings during application runtime.
    """
    global _deprecation_warning_shown
    
    if not _deprecation_warning_shown:
        warnings.warn(
            "Direct imports from app.config are deprecated. "
            "Please use app.core.config instead.",
            DeprecationWarning,
            stacklevel=3
        )
        _deprecation_warning_shown = True


def get_app_settings() -> Settings:
    """
    Returns the application settings instance.
    
    This is a wrapper around get_settings() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        Settings: Application settings instance
    """
    show_deprecation_warning()
    return get_settings()


def init_app_settings() -> Settings:
    """
    Initializes the application settings.
    
    This is a wrapper around initialize_settings() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        Settings: Initialized settings instance
    """
    show_deprecation_warning()
    return initialize_settings()


def get_db_url() -> str:
    """
    Returns the database connection URL.
    
    This is a wrapper around get_database_url() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        str: Database connection URL
    """
    show_deprecation_warning()
    return get_database_url()


def get_redis_config() -> dict:
    """
    Returns Redis connection configuration.
    
    This is a wrapper around get_redis_settings() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        dict: Redis connection parameters
    """
    show_deprecation_warning()
    return get_redis_settings()


def get_minio_config() -> dict:
    """
    Returns MinIO connection configuration.
    
    This is a wrapper around get_minio_settings() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        dict: MinIO connection parameters
    """
    show_deprecation_warning()
    return get_minio_settings()


def get_storage_buckets() -> dict:
    """
    Returns storage bucket names for different file types.
    
    This is a wrapper around get_bucket_names() from core.config that provides
    backward compatibility with existing code.
    
    Returns:
        dict: Dictionary mapping file types to bucket names
    """
    show_deprecation_warning()
    return get_bucket_names()

# Re-export original functions for backward compatibility
__all__ = [
    # Wrapper functions
    "get_app_settings",
    "init_app_settings",
    "get_db_url",
    "get_redis_config",
    "get_minio_config",
    "get_storage_buckets",
    
    # Original functions
    "get_settings",
    "initialize_settings",
    "get_database_url",
    "get_redis_settings",
    "get_minio_settings",
    "get_bucket_names"
]