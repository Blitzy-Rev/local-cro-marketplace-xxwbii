"""
Configuration module for the Molecular Data Management and CRO Integration Platform.

This module provides application settings, environment variable handling, and configuration
management. It implements a singleton pattern to ensure consistent configuration across
the application while supporting fully local deployment without external dependencies.
"""

import os  # standard library
import pathlib  # standard library
from pathlib import Path  # standard library
from typing import Dict, Any, Optional, List  # standard library
import functools  # standard library
from functools import lru_cache  # standard library

import pydantic  # v2.0+
from pydantic import BaseSettings, Field  # v2.0+
from dotenv import load_dotenv  # v1.0+

from ..constants import (
    PROJECT_NAME,
    API_V1_PREFIX,
    TOKEN_EXPIRY_MINUTES,
    REFRESH_TOKEN_EXPIRY_DAYS,
    ALGORITHM,
    BUCKET_NAMES,
    LOG_LEVELS,
    DEFAULT_LOG_LEVEL,
    DEFAULT_LOG_FORMAT,
)

# Base directory and environment file path
BASE_DIR = Path(__file__).resolve().parent.parent.parent
ENV_FILE = os.path.join(BASE_DIR, '.env')

# Singleton instance
_settings_instance = None


def get_env_variable(name: str, default: str) -> str:
    """
    Gets an environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: Default value if environment variable is not set
        
    Returns:
        Value of the environment variable or default
    """
    # Load environment variables from .env file if it exists
    if os.path.exists(ENV_FILE):
        load_dotenv(ENV_FILE)
    
    return os.getenv(name, default)


def get_bool_env_variable(name: str, default: bool) -> bool:
    """
    Gets a boolean environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: Default value if environment variable is not set or conversion fails
        
    Returns:
        Boolean value of the environment variable or default
    """
    value = get_env_variable(name, str(default)).lower()
    
    if value in ('true', 'yes', '1', 'y', 't'):
        return True
    elif value in ('false', 'no', '0', 'n', 'f'):
        return False
    else:
        return default


def get_int_env_variable(name: str, default: int) -> int:
    """
    Gets an integer environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: Default value if environment variable is not set or conversion fails
        
    Returns:
        Integer value of the environment variable or default
    """
    value = get_env_variable(name, str(default))
    
    try:
        return int(value)
    except ValueError:
        return default


class Settings(BaseSettings):
    """
    Application settings loaded from environment variables with defaults.
    
    This class defines all configuration settings for the application, with defaults that
    can be overridden by environment variables. It follows Pydantic's BaseSettings pattern
    for environment variable parsing and validation.
    """
    
    # Application settings
    PROJECT_NAME: str = PROJECT_NAME
    API_V1_PREFIX: str = API_V1_PREFIX
    
    # Database settings
    DATABASE_URL: str = Field(
        default="postgresql://postgres:postgres@localhost:5432/molecular_platform",
        description="Database connection URL"
    )
    DATABASE_POOL_SIZE: int = Field(
        default=5,
        description="Database connection pool size"
    )
    DATABASE_MAX_OVERFLOW: int = Field(
        default=10,
        description="Maximum database connection overflow"
    )
    
    # Redis settings
    REDIS_HOST: str = Field(
        default="localhost",
        description="Redis host"
    )
    REDIS_PORT: int = Field(
        default=6379,
        description="Redis port"
    )
    REDIS_PASSWORD: str = Field(
        default="",
        description="Redis password"
    )
    REDIS_DB: int = Field(
        default=0,
        description="Redis database index"
    )
    
    # MinIO settings
    MINIO_HOST: str = Field(
        default="localhost",
        description="MinIO host"
    )
    MINIO_PORT: int = Field(
        default=9000,
        description="MinIO port"
    )
    MINIO_ACCESS_KEY: str = Field(
        default="minioadmin",
        description="MinIO access key"
    )
    MINIO_SECRET_KEY: str = Field(
        default="minioadmin",
        description="MinIO secret key"
    )
    MINIO_SECURE: bool = Field(
        default=False,
        description="Use secure connection for MinIO"
    )
    MINIO_REGION: str = Field(
        default="us-east-1",
        description="MinIO region"
    )
    
    # JWT Settings
    JWT_PRIVATE_KEY_PATH: str = Field(
        default=os.path.join(BASE_DIR, "keys", "jwt-private.pem"),
        description="Path to JWT private key file"
    )
    JWT_PUBLIC_KEY_PATH: str = Field(
        default=os.path.join(BASE_DIR, "keys", "jwt-public.pem"),
        description="Path to JWT public key file"
    )
    JWT_ALGORITHM: str = Field(
        default=ALGORITHM,
        description="JWT signing algorithm"
    )
    ACCESS_TOKEN_EXPIRE_MINUTES: int = Field(
        default=TOKEN_EXPIRY_MINUTES,
        description="Access token expiration time in minutes"
    )
    REFRESH_TOKEN_EXPIRE_DAYS: int = Field(
        default=REFRESH_TOKEN_EXPIRY_DAYS,
        description="Refresh token expiration time in days"
    )
    
    # CORS settings
    CORS_ORIGINS: str = Field(
        default="http://localhost:3000,http://localhost:8000",
        description="Comma-separated list of allowed origins"
    )
    CORS_METHODS: str = Field(
        default="GET,POST,PUT,DELETE,OPTIONS",
        description="Comma-separated list of allowed methods"
    )
    CORS_HEADERS: str = Field(
        default="*",
        description="Comma-separated list of allowed headers"
    )
    CORS_CREDENTIALS: bool = Field(
        default=True,
        description="Allow credentials"
    )
    
    # Logging settings
    LOG_LEVEL: str = Field(
        default=DEFAULT_LOG_LEVEL,
        description="Log level"
    )
    LOG_FORMAT: str = Field(
        default=DEFAULT_LOG_FORMAT,
        description="Log format"
    )
    
    # Debug mode
    DEBUG: bool = Field(
        default=False,
        description="Debug mode"
    )
    
    # Server settings
    SERVER_HOST: str = Field(
        default="0.0.0.0",
        description="Server host"
    )
    SERVER_PORT: int = Field(
        default=8000,
        description="Server port"
    )
    SERVER_WORKERS: int = Field(
        default=1,
        description="Number of worker processes"
    )
    
    # Bucket names
    BUCKET_NAMES: Dict[str, str] = Field(
        default=BUCKET_NAMES,
        description="Storage bucket names"
    )
    
    def __init__(self, **kwargs):
        """
        Initializes the Settings class with values from environment variables.
        """
        # Load environment variables from .env file if it exists
        if os.path.exists(ENV_FILE):
            load_dotenv(ENV_FILE)
            
        super().__init__(**kwargs)
    
    @property
    def get_jwt_private_key(self) -> str:
        """
        Loads the JWT private key from file.
        
        Returns:
            Private key as string
            
        Raises:
            FileNotFoundError: If the private key file does not exist
        """
        if not os.path.exists(self.JWT_PRIVATE_KEY_PATH):
            raise FileNotFoundError(f"JWT private key file not found: {self.JWT_PRIVATE_KEY_PATH}")
        
        with open(self.JWT_PRIVATE_KEY_PATH, "r") as f:
            private_key = f.read()
            
        return private_key
    
    @property
    def get_jwt_public_key(self) -> str:
        """
        Loads the JWT public key from file.
        
        Returns:
            Public key as string
            
        Raises:
            FileNotFoundError: If the public key file does not exist
        """
        if not os.path.exists(self.JWT_PUBLIC_KEY_PATH):
            raise FileNotFoundError(f"JWT public key file not found: {self.JWT_PUBLIC_KEY_PATH}")
        
        with open(self.JWT_PUBLIC_KEY_PATH, "r") as f:
            public_key = f.read()
            
        return public_key
    
    @property
    def get_cors_origins(self) -> List[str]:
        """
        Parses CORS origins from string to list.
        
        Returns:
            List of allowed origins
        """
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]
    
    @property
    def get_cors_methods(self) -> List[str]:
        """
        Parses CORS methods from string to list.
        
        Returns:
            List of allowed methods
        """
        return [method.strip() for method in self.CORS_METHODS.split(",")]
    
    @property
    def get_cors_headers(self) -> List[str]:
        """
        Parses CORS headers from string to list.
        
        Returns:
            List of allowed headers
        """
        return [header.strip() for header in self.CORS_HEADERS.split(",")]


def initialize_settings() -> Settings:
    """
    Initializes the application settings singleton.
    
    Returns:
        Settings instance
    """
    global _settings_instance
    
    if _settings_instance is None:
        _settings_instance = Settings()
        
    return _settings_instance


@lru_cache()
def get_settings() -> Settings:
    """
    Gets the application settings singleton instance.
    
    Using lru_cache ensures that the settings are only loaded once and reused.
    
    Returns:
        Settings instance
    """
    return initialize_settings()


def get_database_url() -> str:
    """
    Gets the database connection URL from settings.
    
    Returns:
        Database connection URL
    """
    settings = get_settings()
    return settings.DATABASE_URL


def get_redis_settings() -> Dict[str, Any]:
    """
    Gets Redis connection settings from application settings.
    
    Returns:
        Redis connection settings as a dictionary
    """
    settings = get_settings()
    return {
        "host": settings.REDIS_HOST,
        "port": settings.REDIS_PORT,
        "password": settings.REDIS_PASSWORD,
        "db": settings.REDIS_DB
    }


def get_minio_settings() -> Dict[str, Any]:
    """
    Gets MinIO connection settings from application settings.
    
    Returns:
        MinIO connection settings as a dictionary
    """
    settings = get_settings()
    return {
        "host": settings.MINIO_HOST,
        "port": settings.MINIO_PORT,
        "access_key": settings.MINIO_ACCESS_KEY,
        "secret_key": settings.MINIO_SECRET_KEY,
        "secure": settings.MINIO_SECURE,
        "region": settings.MINIO_REGION
    }


def get_bucket_names() -> Dict[str, str]:
    """
    Gets storage bucket names from application settings.
    
    Returns:
        Storage bucket names as a dictionary
    """
    settings = get_settings()
    return settings.BUCKET_NAMES