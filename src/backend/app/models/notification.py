"""
Notification model representing system notifications to users.

This module defines the SQLAlchemy ORM model for the notification entity,
which represents notifications sent to users for various system events
such as experiment status changes, CRO submissions, quotes, and results.
"""

from datetime import datetime  # standard library
from sqlalchemy import Column, String, Boolean, DateTime, Enum, Integer, JSON, ForeignKey  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+

from ..db.base_class import Base
from ..constants import NotificationType


class Notification(Base):
    """
    Notification model representing a system notification sent to a user.
    
    Notifications are created for various system events such as experiment status changes,
    CRO submissions, quotes, and result uploads. They help keep users informed about
    relevant activities in the system.
    
    Attributes:
        id (int): Unique identifier for the notification (inherited from Base)
        user_id (int): Foreign key to the user receiving the notification
        type (NotificationType): Type of notification (e.g., EXPERIMENT_STATUS_CHANGE)
        message (str): Human-readable notification message
        data (JSON): Additional structured data related to the notification
        read_status (bool): Whether the notification has been read by the user
        created_at (datetime): When the notification was created
        read_at (datetime): When the notification was marked as read (if applicable)
        user (User): Relationship to the user receiving the notification
    """
    
    user_id = Column(Integer, ForeignKey("users.id"), index=True, nullable=False)
    type = Column(Enum(NotificationType), nullable=False)
    message = Column(String(255), nullable=False)
    data = Column(JSON, nullable=True)
    read_status = Column(Boolean, default=False, nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    read_at = Column(DateTime, nullable=True)
    
    # Relationship to User model
    user = relationship("User", back_populates="notifications")