"""
Notification schema module for the Molecular Data Management and CRO Integration Platform.

This module defines Pydantic schema models for notification data validation, 
serialization, and API responses. These schemas support the notification system
for informing users about system events and updates.
"""

from datetime import datetime
from typing import Optional, List, Dict, Any

from pydantic import BaseModel, Field  # version 2.0+

from ..constants import NotificationType

# Generate a list of valid notification types from the enum for validation
NOTIFICATION_TYPES = [str(item.name).lower() for item in NotificationType]


def validate_notification_type(notification_type: str) -> str:
    """
    Validate that notification type is one of the allowed values.
    
    Args:
        notification_type: The notification type to validate
        
    Returns:
        The validated notification type
        
    Raises:
        ValueError: If the notification type is not valid
    """
    notification_type = notification_type.lower()
    if notification_type not in NOTIFICATION_TYPES:
        raise ValueError(
            f"Invalid notification type. Must be one of: {', '.join(NOTIFICATION_TYPES)}"
        )
    return notification_type


class NotificationBase(BaseModel):
    """Base Pydantic model for notification data with common fields."""
    user_id: int = Field(..., description="ID of the user who will receive the notification")
    type: str = Field(..., description="Type of notification")
    message: str = Field(..., description="Notification message content")
    data: Optional[Dict[str, Any]] = Field(None, description="Additional data related to the notification")
    
    def validate_type(self, values: dict) -> dict:
        """
        Validate that notification type is one of the allowed values.
        
        Args:
            values: The values dictionary passed to the model
            
        Returns:
            The validated values dictionary
        """
        if values.get("type"):
            values["type"] = validate_notification_type(values["type"])
        return values


class NotificationCreate(NotificationBase):
    """Pydantic model for notification creation."""
    pass


class NotificationRead(NotificationBase):
    """Pydantic model for notification data in API responses."""
    id: int = Field(..., description="Unique identifier for the notification")
    read_status: bool = Field(..., description="Whether the notification has been read")
    created_at: datetime = Field(..., description="Timestamp when the notification was created")
    read_at: Optional[datetime] = Field(None, description="Timestamp when the notification was read")


class NotificationUpdate(BaseModel):
    """Pydantic model for updating notification status."""
    read_status: bool = Field(..., description="Whether the notification has been read")


class NotificationFilter(BaseModel):
    """Pydantic model for filtering notifications in API requests."""
    read_status: Optional[bool] = Field(None, description="Filter by read status (true/false)")
    type: Optional[str] = Field(None, description="Filter by notification type")
    from_date: Optional[datetime] = Field(None, description="Filter notifications created after this date")
    to_date: Optional[datetime] = Field(None, description="Filter notifications created before this date")
    
    def validate_type(self, values: dict) -> dict:
        """
        Validate that notification type is one of the allowed values if provided.
        
        Args:
            values: The values dictionary passed to the model
            
        Returns:
            The validated values dictionary
        """
        if values.get("type") is not None:
            values["type"] = validate_notification_type(values["type"])
        return values


class NotificationList(BaseModel):
    """Pydantic model for paginated list of notifications in API responses."""
    items: List[NotificationRead] = Field(..., description="List of notification items")
    total: int = Field(..., description="Total number of items matching the query")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")
    pages: int = Field(..., description="Total number of pages")


class NotificationBulkUpdate(BaseModel):
    """Pydantic model for bulk updating notification status."""
    notification_ids: List[int] = Field(..., description="List of notification IDs to update")
    read_status: bool = Field(..., description="Whether the notifications have been read")


class NotificationCount(BaseModel):
    """Pydantic model for notification count in API responses."""
    total: int = Field(..., description="Total number of notifications")
    unread: int = Field(..., description="Number of unread notifications")