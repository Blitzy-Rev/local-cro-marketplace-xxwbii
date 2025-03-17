"""
Pydantic schema models for user data validation, serialization, and API responses.

This module defines the data models used for user management operations, with
appropriate validation rules, type annotations, and documentation.
"""

from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel, Field, EmailStr, constr

from ..constants import UserRole, UserStatus


class UserBase(BaseModel):
    """Base Pydantic model for user data with common fields."""
    email: EmailStr = Field(..., description="User's email address", example="user@example.com")


class UserCreate(UserBase):
    """Pydantic model for user creation with password validation."""
    # Password must be at least 10 characters and include lowercase, uppercase, 
    # digit, and special character
    password: constr(min_length=10, regex=r'^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[@$!%*?&])[A-Za-z\d@$!%*?&]{10,}$') = Field(
        ..., 
        description="User's password, must be at least 10 characters and include lowercase, uppercase, digit, and special character",
        example="SecureP@ssw0rd"
    )
    role: str = Field(..., description="User role (pharma, cro, admin)", example="pharma")

    def validate_role(self, values: dict) -> dict:
        """Validate that role is one of the allowed values."""
        role = values.get('role')
        allowed_roles = [role.name.lower() for role in UserRole]
        if role.lower() not in allowed_roles:
            raise ValueError(f"Role must be one of {allowed_roles}")
        return values


class UserUpdate(BaseModel):
    """Pydantic model for updating user information."""
    email: Optional[EmailStr] = Field(None, description="User's email address", example="user@example.com")
    password: Optional[str] = Field(None, description="User's password", example="NewSecureP@ssw0rd")
    role: Optional[str] = Field(None, description="User role (pharma, cro, admin)", example="pharma")
    status: Optional[str] = Field(None, description="User status (active, inactive, pending, locked)", example="active")

    def validate_role(self, values: dict) -> dict:
        """Validate that role is one of the allowed values if provided."""
        role = values.get('role')
        if role is not None:
            allowed_roles = [role.name.lower() for role in UserRole]
            if role.lower() not in allowed_roles:
                raise ValueError(f"Role must be one of {allowed_roles}")
        return values

    def validate_status(self, values: dict) -> dict:
        """Validate that status is one of the allowed values if provided."""
        status = values.get('status')
        if status is not None:
            allowed_statuses = [status.name.lower() for status in UserStatus]
            if status.lower() not in allowed_statuses:
                raise ValueError(f"Status must be one of {allowed_statuses}")
        return values


class UserInDB(UserBase):
    """Pydantic model for user data as stored in the database, including password hash."""
    id: int = Field(..., description="Unique user ID", example=1)
    password_hash: str = Field(..., description="Hashed password", example="$2b$12$EixZaYVK1fsbw1ZfbX3OXePaWxn96p36WQoeG6Lruj3vjPGga31lW")
    role: str = Field(..., description="User role (pharma, cro, admin)", example="pharma")
    status: str = Field(..., description="User status (active, inactive, pending, locked)", example="active")
    is_active: bool = Field(..., description="Whether the user is active", example=True)
    email_verified: bool = Field(..., description="Whether the email address has been verified", example=True)
    created_at: datetime = Field(..., description="Timestamp of user creation", example="2023-01-01T00:00:00Z")
    updated_at: Optional[datetime] = Field(None, description="Timestamp of last update", example="2023-01-02T00:00:00Z")
    last_login: Optional[datetime] = Field(None, description="Timestamp of last login", example="2023-01-03T00:00:00Z")
    password_history: Optional[List[str]] = Field(None, description="List of previously used password hashes", 
                                                example=["$2b$12$1ab2c3d4e5f6g7h8i9j0k.", "$2b$12$9i8h7g6f5e4d3c2b1a0z."])


class UserResponse(UserBase):
    """Pydantic model for user data in API responses, excluding sensitive information."""
    id: int = Field(..., description="Unique user ID", example=1)
    role: str = Field(..., description="User role (pharma, cro, admin)", example="pharma")


class UserDetailResponse(UserBase):
    """Pydantic model for detailed user data in API responses."""
    id: int = Field(..., description="Unique user ID", example=1)
    role: str = Field(..., description="User role (pharma, cro, admin)", example="pharma")
    status: str = Field(..., description="User status (active, inactive, pending, locked)", example="active")
    is_active: bool = Field(..., description="Whether the user is active", example=True)
    email_verified: bool = Field(..., description="Whether the email address has been verified", example=True)
    created_at: datetime = Field(..., description="Timestamp of user creation", example="2023-01-01T00:00:00Z")
    updated_at: Optional[datetime] = Field(None, description="Timestamp of last update", example="2023-01-02T00:00:00Z")
    last_login: Optional[datetime] = Field(None, description="Timestamp of last login", example="2023-01-03T00:00:00Z")