"""
Pydantic schema models for user-related API requests and responses in the API v1 namespace.

These schemas are used for data validation, serialization, and documentation in the
user management endpoints.
"""

from datetime import datetime
from typing import Optional, List

from pydantic import BaseModel, Field, EmailStr

from ...schemas.user import UserResponse, UserUpdate
from .auth import MessageResponse


class UserProfileResponse(BaseModel):
    """Schema for user profile information in API responses."""
    id: int = Field(..., description="Unique identifier for the user")
    email: str = Field(..., description="User's email address")
    role: str = Field(..., description="User's role (pharma, cro, admin)")
    status: str = Field(..., description="User's status (active, inactive, pending)")
    email_verified: bool = Field(..., description="Whether the user's email is verified")
    created_at: datetime = Field(..., description="When the user account was created")
    last_login: Optional[datetime] = Field(None, description="When the user last logged in")


class UserUpdateRequest(BaseModel):
    """Schema for updating user profile information."""
    email: Optional[EmailStr] = Field(None, description="New email address")
    password: Optional[str] = Field(None, description="New password")


class UserFilterParams(BaseModel):
    """Schema for filtering users in list endpoints."""
    email: Optional[str] = Field(None, description="Filter by email address")
    role: Optional[str] = Field(None, description="Filter by role")
    status: Optional[str] = Field(None, description="Filter by status")
    email_verified: Optional[bool] = Field(None, description="Filter by email verification status")
    created_after: Optional[datetime] = Field(None, description="Filter users created after this date")
    created_before: Optional[datetime] = Field(None, description="Filter users created before this date")

    def validate_role(self, values: dict) -> dict:
        """Validate that role is one of the allowed values if provided."""
        role = values.get('role')
        if role is not None:
            allowed_roles = ['pharma', 'cro', 'admin']
            if role not in allowed_roles:
                raise ValueError(f"Role must be one of {allowed_roles}")
        return values

    def validate_status(self, values: dict) -> dict:
        """Validate that status is one of the allowed values if provided."""
        status = values.get('status')
        if status is not None:
            allowed_statuses = ['active', 'inactive', 'pending']
            if status not in allowed_statuses:
                raise ValueError(f"Status must be one of {allowed_statuses}")
        return values


class UserListResponse(BaseModel):
    """Schema for paginated list of users in API responses."""
    items: List[UserProfileResponse] = Field(..., description="List of user profiles")
    total: int = Field(..., description="Total number of users matching the filters")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")


class UserStatusUpdateRequest(BaseModel):
    """Schema for updating user status."""
    status: str = Field(..., description="New user status")

    def validate_status(self, values: dict) -> dict:
        """Validate that status is one of the allowed values."""
        status = values.get('status')
        allowed_statuses = ['active', 'inactive', 'pending']
        if status not in allowed_statuses:
            raise ValueError(f"Status must be one of {allowed_statuses}")
        return values


class ChangePasswordRequest(BaseModel):
    """Schema for changing user password."""
    current_password: str = Field(..., description="Current password")
    new_password: str = Field(..., description="New password")
    confirm_password: str = Field(..., description="Confirm new password")

    def validate_passwords_match(self, values: dict) -> dict:
        """Validate that new_password and confirm_password match."""
        new_password = values.get('new_password')
        confirm_password = values.get('confirm_password')
        if new_password != confirm_password:
            raise ValueError("Passwords do not match")
        return values