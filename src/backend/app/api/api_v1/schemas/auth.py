"""
Pydantic schema models for authentication-related requests and responses.

This module defines the data models used for authentication operations, including
login, registration, token refresh, password management, and more. These models
include appropriate validation rules, type annotations, and documentation.
"""

from pydantic import BaseModel, Field, EmailStr, constr
from ...schemas.token import Token, TokenResponse
from ...schemas.user import UserResponse

class LoginRequest(BaseModel):
    """Schema for login request data."""
    email: EmailStr = Field(..., description="User's email address")
    password: str = Field(..., description="User's password")

class LoginResponse(BaseModel):
    """Schema for login response data."""
    tokens: TokenResponse = Field(..., description="JWT tokens for authentication")
    user: UserResponse = Field(..., description="User information")

class RefreshTokenRequest(BaseModel):
    """Schema for refresh token request data."""
    refresh_token: str = Field(..., description="JWT refresh token")

class RegistrationRequest(BaseModel):
    """Schema for user registration request data."""
    email: EmailStr = Field(..., description="User's email address")
    password: constr(min_length=10, regex='^(?=.*[a-z])(?=.*[A-Z])(?=.*\\d)(?=.*[@$!%*?&])[A-Za-z\\d@$!%*?&]{10,}$') = Field(
        ..., 
        description="User's password, must be at least 10 characters and include lowercase, uppercase, digit, and special character"
    )
    confirm_password: str = Field(..., description="Confirmation of user's password")
    role: str = Field(..., description="User role (pharma, cro, admin)")

    def validate_passwords_match(self, values):
        """Validate that password and confirm_password match."""
        password = values.get('password')
        confirm_password = values.get('confirm_password')
        if password != confirm_password:
            raise ValueError('Passwords do not match')
        return values

    def validate_role(self, values):
        """Validate that role is one of the allowed values."""
        role = values.get('role')
        allowed_roles = ['pharma', 'cro', 'admin']
        if role not in allowed_roles:
            raise ValueError(f"Role must be one of: {', '.join(allowed_roles)}")
        return values

class RegistrationResponse(BaseModel):
    """Schema for user registration response data."""
    success: bool = Field(..., description="Whether the registration was successful")
    user: UserResponse = Field(..., description="User information")
    message: str = Field(..., description="Registration result message")

class EmailVerificationRequest(BaseModel):
    """Schema for email verification request data."""
    token: str = Field(..., description="Email verification token")

class PasswordResetRequest(BaseModel):
    """Schema for password reset request data."""
    email: EmailStr = Field(..., description="User's email address")

class PasswordResetConfirmRequest(BaseModel):
    """Schema for password reset confirmation request data."""
    token: str = Field(..., description="Password reset token")
    new_password: constr(min_length=10, regex='^(?=.*[a-z])(?=.*[A-Z])(?=.*\\d)(?=.*[@$!%*?&])[A-Za-z\\d@$!%*?&]{10,}$') = Field(
        ..., 
        description="User's new password, must be at least 10 characters and include lowercase, uppercase, digit, and special character"
    )
    confirm_password: str = Field(..., description="Confirmation of user's new password")

    def validate_passwords_match(self, values):
        """Validate that new_password and confirm_password match."""
        new_password = values.get('new_password')
        confirm_password = values.get('confirm_password')
        if new_password != confirm_password:
            raise ValueError('Passwords do not match')
        return values

class ChangePasswordRequest(BaseModel):
    """Schema for change password request data."""
    current_password: str = Field(..., description="User's current password")
    new_password: constr(min_length=10, regex='^(?=.*[a-z])(?=.*[A-Z])(?=.*\\d)(?=.*[@$!%*?&])[A-Za-z\\d@$!%*?&]{10,}$') = Field(
        ..., 
        description="User's new password, must be at least 10 characters and include lowercase, uppercase, digit, and special character"
    )
    confirm_password: str = Field(..., description="Confirmation of user's new password")

    def validate_passwords_match(self, values):
        """Validate that new_password and confirm_password match."""
        new_password = values.get('new_password')
        confirm_password = values.get('confirm_password')
        if new_password != confirm_password:
            raise ValueError('Passwords do not match')
        return values

class MessageResponse(BaseModel):
    """Schema for generic message response data."""
    success: bool = Field(..., description="Whether the operation was successful")
    message: str = Field(..., description="Operation result message")