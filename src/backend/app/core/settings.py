"""
Settings module for the Molecular Data Management and CRO Integration Platform.

This module defines application settings and configuration parameters loaded
from environment variables with sensible defaults. It supports the requirement
for fully local deployment without external dependencies.
"""

import os  # standard library
from pathlib import Path  # standard library
from typing import Dict, List, Optional, Any  # standard library

from pydantic import BaseSettings  # version 2.0+
from dotenv import load_dotenv  # version 1.0+

from ..constants import (
    PROJECT_NAME,
    API_V1_PREFIX,
    TOKEN_EXPIRY_MINUTES,
    REFRESH_TOKEN_EXPIRY_DAYS,
    ALGORITHM,
    BUCKET_NAMES,
)

# Base directory and environment file path
BASE_DIR = Path(__file__).resolve().parent.parent.parent
ENV_FILE = os.path.join(BASE_DIR, '.env')

# Load environment variables from .env file if it exists
if os.path.exists(ENV_FILE):
    load_dotenv(ENV_FILE)


def get_env_variable(name: str, default: str = "") -> str:
    """Gets an environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: The default value if the variable is not found
        
    Returns:
        The value of the environment variable or the default value
    """
    return os.environ.get(name, default)


def get_bool_env_variable(name: str, default: bool = False) -> bool:
    """Gets a boolean environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: The default boolean value if the variable is not found or invalid
        
    Returns:
        The boolean value of the environment variable or the default value
    """
    value = get_env_variable(name, "").lower()
    
    if value in ("true", "yes", "1", "t", "y"):
        return True
    elif value in ("false", "no", "0", "f", "n"):
        return False
    
    return default


def get_int_env_variable(name: str, default: int = 0) -> int:
    """Gets an integer environment variable with a default value.
    
    Args:
        name: The name of the environment variable
        default: The default integer value if the variable is not found or invalid
        
    Returns:
        The integer value of the environment variable or the default value
    """
    value = get_env_variable(name, "")
    
    try:
        return int(value)
    except (ValueError, TypeError):
        return default


class Settings(BaseSettings):
    """Application settings loaded from environment variables with defaults."""
    
    # Application settings
    PROJECT_NAME: str = PROJECT_NAME
    API_V1_PREFIX: str = API_V1_PREFIX
    
    # Database settings
    DATABASE_URL: str = get_env_variable(
        "DATABASE_URL", 
        "postgresql://postgres:postgres@postgres:5432/molecular_data"
    )
    DATABASE_POOL_SIZE: int = get_int_env_variable("DATABASE_POOL_SIZE", 20)
    DATABASE_MAX_OVERFLOW: int = get_int_env_variable("DATABASE_MAX_OVERFLOW", 10)
    
    # Redis settings
    REDIS_HOST: str = get_env_variable("REDIS_HOST", "redis")
    REDIS_PORT: int = get_int_env_variable("REDIS_PORT", 6379)
    REDIS_PASSWORD: str = get_env_variable("REDIS_PASSWORD", "")
    REDIS_DB: int = get_int_env_variable("REDIS_DB", 0)
    
    # MinIO settings (S3-compatible storage)
    MINIO_HOST: str = get_env_variable("MINIO_HOST", "minio")
    MINIO_PORT: int = get_int_env_variable("MINIO_PORT", 9000)
    MINIO_ACCESS_KEY: str = get_env_variable("MINIO_ACCESS_KEY", "minioadmin")
    MINIO_SECRET_KEY: str = get_env_variable("MINIO_SECRET_KEY", "minioadmin")
    MINIO_SECURE: bool = get_bool_env_variable("MINIO_SECURE", False)
    MINIO_REGION: str = get_env_variable("MINIO_REGION", "us-east-1")
    
    # JWT settings
    JWT_PRIVATE_KEY_PATH: str = get_env_variable(
        "JWT_PRIVATE_KEY_PATH", 
        os.path.join(BASE_DIR, "keys", "jwt-private.pem")
    )
    JWT_PUBLIC_KEY_PATH: str = get_env_variable(
        "JWT_PUBLIC_KEY_PATH", 
        os.path.join(BASE_DIR, "keys", "jwt-public.pem")
    )
    JWT_ALGORITHM: str = ALGORITHM
    ACCESS_TOKEN_EXPIRE_MINUTES: int = get_int_env_variable(
        "ACCESS_TOKEN_EXPIRE_MINUTES", 
        TOKEN_EXPIRY_MINUTES
    )
    REFRESH_TOKEN_EXPIRE_DAYS: int = get_int_env_variable(
        "REFRESH_TOKEN_EXPIRE_DAYS", 
        REFRESH_TOKEN_EXPIRY_DAYS
    )
    
    # CORS settings
    CORS_ORIGINS: str = get_env_variable("CORS_ORIGINS", "*")
    CORS_METHODS: str = get_env_variable("CORS_METHODS", "GET,POST,PUT,DELETE,OPTIONS")
    CORS_HEADERS: str = get_env_variable("CORS_HEADERS", "*")
    CORS_CREDENTIALS: bool = get_bool_env_variable("CORS_CREDENTIALS", True)
    
    # Logging settings
    LOG_LEVEL: str = get_env_variable("LOG_LEVEL", "INFO")
    LOG_FORMAT: str = get_env_variable(
        "LOG_FORMAT", 
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    
    # Debug mode
    DEBUG: bool = get_bool_env_variable("DEBUG", False)
    
    # Server settings
    SERVER_HOST: str = get_env_variable("SERVER_HOST", "0.0.0.0")
    SERVER_PORT: int = get_int_env_variable("SERVER_PORT", 8000)
    SERVER_WORKERS: int = get_int_env_variable("SERVER_WORKERS", 1)
    
    # Storage buckets from constants
    BUCKET_NAMES: Dict[str, str] = BUCKET_NAMES
    
    class Config:
        """Pydantic configuration class."""
        env_file = ENV_FILE if os.path.exists(ENV_FILE) else None
        case_sensitive = True
    
    @property
    def get_jwt_private_key(self) -> str:
        """Loads the JWT private key from file.
        
        Returns:
            Private key as string
            
        Raises:
            FileNotFoundError: If the private key file doesn't exist
        """
        if not os.path.exists(self.JWT_PRIVATE_KEY_PATH):
            raise FileNotFoundError(
                f"JWT private key not found at {self.JWT_PRIVATE_KEY_PATH}"
            )
        
        with open(self.JWT_PRIVATE_KEY_PATH, "r") as f:
            return f.read()
    
    @property
    def get_jwt_public_key(self) -> str:
        """Loads the JWT public key from file.
        
        Returns:
            Public key as string
            
        Raises:
            FileNotFoundError: If the public key file doesn't exist
        """
        if not os.path.exists(self.JWT_PUBLIC_KEY_PATH):
            raise FileNotFoundError(
                f"JWT public key not found at {self.JWT_PUBLIC_KEY_PATH}"
            )
        
        with open(self.JWT_PUBLIC_KEY_PATH, "r") as f:
            return f.read()
    
    @property
    def get_cors_origins(self) -> List[str]:
        """Parses CORS origins from string to list.
        
        Returns:
            List of allowed origins
        """
        if self.CORS_ORIGINS == "*":
            return ["*"]
        return [origin.strip() for origin in self.CORS_ORIGINS.split(",")]
    
    @property
    def get_cors_methods(self) -> List[str]:
        """Parses CORS methods from string to list.
        
        Returns:
            List of allowed methods
        """
        return [method.strip() for method in self.CORS_METHODS.split(",")]
    
    @property
    def get_cors_headers(self) -> List[str]:
        """Parses CORS headers from string to list.
        
        Returns:
            List of allowed headers
        """
        if self.CORS_HEADERS == "*":
            return ["*"]
        return [header.strip() for header in self.CORS_HEADERS.split(",")]


# Export the Settings class and utility functions
__all__ = ["Settings", "get_env_variable", "get_bool_env_variable", "get_int_env_variable"]