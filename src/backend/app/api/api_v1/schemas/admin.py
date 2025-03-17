"""
Pydantic schema models for admin-related API requests and responses.

This module defines schema models for administrative operations including
user management, system statistics and monitoring, activity logs, and system alerts.
"""

from datetime import datetime
from typing import Optional, List, Dict, Any

from pydantic import BaseModel, Field, EmailStr, constr

from ....schemas.user import UserBase, UserCreate, UserUpdate, UserResponse
from .auth import MessageResponse


class UserCreateRequest(BaseModel):
    """Schema for admin user creation request."""
    email: EmailStr = Field(..., description="User's email address", example="user@example.com")
    password: str = Field(..., description="User's password", example="SecureP@ssw0rd")
    role: str = Field(..., description="User role (pharma, cro, admin)", example="pharma")
    status: Optional[str] = Field(None, description="User status (active, inactive, pending)", example="active")
    email_verified: Optional[bool] = Field(None, description="Whether the email is verified", example=True)
    
    def validate_role(self, values: dict) -> dict:
        """Validate that role is one of the allowed values."""
        role = values.get('role')
        allowed_roles = ['pharma', 'cro', 'admin']
        if role.lower() not in allowed_roles:
            raise ValueError(f"Role must be one of {allowed_roles}")
        return values
    
    def validate_status(self, values: dict) -> dict:
        """Validate that status is one of the allowed values if provided."""
        status = values.get('status')
        if status is not None:
            allowed_statuses = ['active', 'inactive', 'pending']
            if status.lower() not in allowed_statuses:
                raise ValueError(f"Status must be one of {allowed_statuses}")
        return values


class UserAdminUpdateRequest(BaseModel):
    """Schema for admin user update request with additional fields."""
    email: Optional[EmailStr] = Field(None, description="User's email address", example="user@example.com")
    role: Optional[str] = Field(None, description="User role (pharma, cro, admin)", example="pharma")
    status: Optional[str] = Field(None, description="User status (active, inactive, pending)", example="active")
    email_verified: Optional[bool] = Field(None, description="Whether the email is verified", example=True)
    is_active: Optional[bool] = Field(None, description="Whether the user is active", example=True)
    
    def validate_role(self, values: dict) -> dict:
        """Validate that role is one of the allowed values if provided."""
        role = values.get('role')
        if role is not None:
            allowed_roles = ['pharma', 'cro', 'admin']
            if role.lower() not in allowed_roles:
                raise ValueError(f"Role must be one of {allowed_roles}")
        return values
    
    def validate_status(self, values: dict) -> dict:
        """Validate that status is one of the allowed values if provided."""
        status = values.get('status')
        if status is not None:
            allowed_statuses = ['active', 'inactive', 'pending']
            if status.lower() not in allowed_statuses:
                raise ValueError(f"Status must be one of {allowed_statuses}")
        return values


class UserPasswordResetRequest(BaseModel):
    """Schema for admin password reset request."""
    new_password: str = Field(..., description="New password for the user", example="NewSecureP@ssw0rd")
    
    def validate_password(self, values: dict) -> dict:
        """Validate that password meets complexity requirements."""
        password = values.get('new_password')
        # Check minimum length
        if len(password) < 10:
            raise ValueError("Password must be at least 10 characters long")
        
        # Check for lowercase, uppercase, digit, and special character
        has_lower = any(c.islower() for c in password)
        has_upper = any(c.isupper() for c in password)
        has_digit = any(c.isdigit() for c in password)
        has_special = any(c in "@$!%*?&" for c in password)
        
        if not (has_lower and has_upper and has_digit and has_special):
            raise ValueError("Password must include lowercase, uppercase, digit, and special character")
        
        return values


class SystemStatsResponse(BaseModel):
    """Schema for system statistics response."""
    user_counts: Dict[str, int] = Field(..., description="User counts by role", example={"pharma": 35, "cro": 8, "admin": 2})
    active_users: int = Field(..., description="Number of active users", example=12)
    total_molecules: int = Field(..., description="Total number of molecules in the system", example=5432)
    total_libraries: int = Field(..., description="Total number of libraries", example=48)
    total_experiments: int = Field(..., description="Total number of experiments", example=120)
    total_submissions: int = Field(..., description="Total number of submissions", example=95)
    resource_usage: Dict[str, float] = Field(..., description="System resource usage", example={"cpu": 32.5, "memory": 45.8, "disk": 68.2})
    timestamp: datetime = Field(..., description="Timestamp of the statistics", example="2023-06-05T14:30:00Z")


class SystemResourcesResponse(BaseModel):
    """Schema for detailed system resources response."""
    cpu: Dict[str, float] = Field(..., description="CPU usage statistics", example={"usage_percent": 32.5, "load_avg_1m": 1.2, "load_avg_5m": 1.1, "load_avg_15m": 0.9})
    memory: Dict[str, float] = Field(..., description="Memory usage statistics", example={"total_gb": 16.0, "used_gb": 7.3, "used_percent": 45.8})
    disk: Dict[str, float] = Field(..., description="Disk usage statistics", example={"total_gb": 500.0, "used_gb": 341.0, "used_percent": 68.2})
    database: Dict[str, Any] = Field(..., description="Database statistics", example={"connections": 12, "size_mb": 2340, "active_queries": 3})
    storage: Dict[str, Any] = Field(..., description="Storage statistics", example={"total_gb": 1000.0, "used_gb": 423.5, "used_percent": 42.35, "file_count": 12450})
    timestamp: datetime = Field(..., description="Timestamp of the resource statistics", example="2023-06-05T14:30:00Z")


class ActivityLogResponse(BaseModel):
    """Schema for activity log entry."""
    id: int = Field(..., description="Activity log entry ID", example=12345)
    user_id: Optional[int] = Field(None, description="ID of the user who performed the action", example=42)
    username: Optional[str] = Field(None, description="Username of the user who performed the action", example="jsmith")
    action: str = Field(..., description="Action performed", example="create_molecule")
    resource_type: str = Field(..., description="Type of resource affected", example="molecule")
    resource_id: Optional[int] = Field(None, description="ID of the resource affected", example=789)
    details: Optional[Dict[str, Any]] = Field(None, description="Additional details about the action", example={"smiles": "CCO", "name": "Ethanol"})
    ip_address: str = Field(..., description="IP address from which the action was performed", example="192.168.1.100")
    timestamp: datetime = Field(..., description="Timestamp of the action", example="2023-06-05T14:32:15Z")


class ActivityLogFilterParams(BaseModel):
    """Schema for activity log filtering parameters."""
    user_id: Optional[int] = Field(None, description="Filter by user ID", example=42)
    action: Optional[str] = Field(None, description="Filter by action", example="create_molecule")
    resource_type: Optional[str] = Field(None, description="Filter by resource type", example="molecule")
    resource_id: Optional[int] = Field(None, description="Filter by resource ID", example=789)
    start_date: Optional[datetime] = Field(None, description="Filter by start date", example="2023-06-01T00:00:00Z")
    end_date: Optional[datetime] = Field(None, description="Filter by end date", example="2023-06-05T23:59:59Z")


class ActivityLogListResponse(BaseModel):
    """Schema for paginated list of activity logs."""
    items: List[ActivityLogResponse] = Field(..., description="List of activity log entries")
    total: int = Field(..., description="Total number of activity log entries matching the filter", example=1245)
    page: int = Field(..., description="Current page number", example=1)
    size: int = Field(..., description="Number of items per page", example=20)


class SystemAlertResponse(BaseModel):
    """Schema for system alert."""
    id: int = Field(..., description="Alert ID", example=42)
    title: str = Field(..., description="Alert title", example="Storage usage high")
    message: str = Field(..., description="Alert message", example="Storage usage is approaching 75% capacity")
    severity: str = Field(..., description="Alert severity", example="warning")
    alert_type: str = Field(..., description="Alert type", example="storage")
    resolved: bool = Field(..., description="Whether the alert has been resolved", example=False)
    created_at: datetime = Field(..., description="Timestamp when the alert was created", example="2023-06-05T14:30:00Z")
    resolved_at: Optional[datetime] = Field(None, description="Timestamp when the alert was resolved", example="2023-06-05T15:45:00Z")
    resolved_by: Optional[int] = Field(None, description="ID of the user who resolved the alert", example=1)
    details: Optional[Dict[str, Any]] = Field(None, description="Additional alert details", example={"current_usage": 75.2, "threshold": 75.0})


class SystemAlertListResponse(BaseModel):
    """Schema for list of system alerts."""
    items: List[SystemAlertResponse] = Field(..., description="List of system alerts")
    total: int = Field(..., description="Total number of alerts", example=5)
    counts_by_severity: Dict[str, int] = Field(..., description="Count of alerts by severity", example={"critical": 1, "warning": 3, "info": 1})


class SystemAlertUpdateRequest(BaseModel):
    """Schema for updating system alert."""
    resolved: bool = Field(..., description="Whether the alert should be marked as resolved", example=True)
    resolution_notes: Optional[str] = Field(None, description="Notes about how the alert was resolved", example="Cleaned up unused files to reduce storage usage")