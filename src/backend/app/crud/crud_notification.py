from typing import Optional, List, Dict, Any, Union
from datetime import datetime
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select, update, and_, or_, func

from .base import CRUDBase
from ..models.notification import Notification
from ..schemas.notification import NotificationCreate, NotificationUpdate, NotificationFilter


class CRUDNotification(CRUDBase[Notification, NotificationCreate, NotificationUpdate]):
    """
    CRUD operations for Notification model with specialized methods for notification management.
    
    This class extends the base CRUD operations with notification-specific functionality
    such as filtering by user, marking notifications as read, and handling notification
    status updates.
    """
    
    def __init__(self):
        """Initialize the CRUD notification object with the Notification model."""
        super().__init__(Notification)
    
    def create_notification(self, db: Session, obj_in: NotificationCreate) -> Notification:
        """
        Create a new notification.
        
        Args:
            db: Database session
            obj_in: Notification data to create
            
        Returns:
            Created notification instance
        """
        # Convert obj_in to dict and set default values if needed
        obj_in_data = obj_in.dict()
        
        # Set read_status to False if not provided
        if "read_status" not in obj_in_data:
            obj_in_data["read_status"] = False
            
        # Set created_at to current datetime
        obj_in_data["created_at"] = datetime.utcnow()
        
        # Call parent create method
        return super().create(db, obj_in_data)
    
    def get_by_user(
        self, 
        db: Session, 
        user_id: int, 
        skip: int = 0, 
        limit: int = 100,
        filters: Optional[NotificationFilter] = None,
        count: bool = False
    ) -> Union[List[Notification], int]:
        """
        Get notifications for a specific user with pagination and optional filtering.
        
        Args:
            db: Database session
            user_id: ID of the user whose notifications to retrieve
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            filters: Optional filter criteria
            count: If True, return count instead of results
            
        Returns:
            List of notifications or count if count=True
        """
        # Create base query filtering by user_id
        query = select(Notification).where(Notification.user_id == user_id)
        
        # Apply additional filters if provided
        if filters:
            filter_conditions = []
            
            if filters.read_status is not None:
                filter_conditions.append(Notification.read_status == filters.read_status)
                
            if filters.type:
                filter_conditions.append(Notification.type == filters.type)
                
            if filters.from_date:
                filter_conditions.append(Notification.created_at >= filters.from_date)
                
            if filters.to_date:
                filter_conditions.append(Notification.created_at <= filters.to_date)
                
            if filter_conditions:
                query = query.where(and_(*filter_conditions))
        
        # If count requested, return count only
        if count:
            count_query = select(func.count()).select_from(Notification).where(query.whereclause)
            return db.execute(count_query).scalar_one()
        
        # Otherwise, apply ordering and pagination
        query = query.order_by(Notification.created_at.desc())
        query = query.offset(skip).limit(limit)
        
        # Execute query and return results
        result = db.execute(query).scalars().all()
        return list(result)
    
    def get_unread_by_user(
        self, 
        db: Session, 
        user_id: int, 
        skip: int = 0, 
        limit: int = 100,
        filters: Optional[NotificationFilter] = None,
        count: bool = False
    ) -> Union[List[Notification], int]:
        """
        Get unread notifications for a specific user with pagination.
        
        Args:
            db: Database session
            user_id: ID of the user whose unread notifications to retrieve
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            filters: Optional filter criteria
            count: If True, return count instead of results
            
        Returns:
            List of unread notifications or count if count=True
        """
        # Create filter for unread notifications
        unread_filter = NotificationFilter(read_status=False)
        
        # Merge with provided filters if any
        if filters:
            merged_filters = filters.dict(exclude_unset=True)
            merged_filters["read_status"] = False
            filters = NotificationFilter(**merged_filters)
        else:
            filters = unread_filter
        
        # Call get_by_user with unread filter
        return self.get_by_user(db, user_id, skip, limit, filters, count)
    
    def count_unread_by_user(self, db: Session, user_id: int) -> int:
        """
        Count unread notifications for a specific user.
        
        Args:
            db: Database session
            user_id: ID of the user whose unread notifications to count
            
        Returns:
            Count of unread notifications
        """
        return self.get_unread_by_user(db, user_id, count=True)
    
    def mark_as_read(self, db: Session, notification: Notification) -> Notification:
        """
        Mark a notification as read.
        
        Args:
            db: Database session
            notification: Notification instance to mark as read
            
        Returns:
            Updated notification instance
        """
        notification.read_status = True
        notification.read_at = datetime.utcnow()
        
        db.add(notification)
        db.commit()
        db.refresh(notification)
        
        return notification
    
    def mark_all_as_read(self, db: Session, user_id: int) -> int:
        """
        Mark all unread notifications for a user as read.
        
        Args:
            db: Database session
            user_id: ID of the user whose notifications to mark as read
            
        Returns:
            Number of notifications marked as read
        """
        # Get all unread notifications for the user
        notifications = self.get_unread_by_user(db, user_id)
        
        # Mark each notification as read
        count = 0
        for notification in notifications:
            self.mark_as_read(db, notification)
            count += 1
            
        return count
    
    def get_by_type(
        self, 
        db: Session, 
        notification_type: str, 
        skip: int = 0, 
        limit: int = 100
    ) -> List[Notification]:
        """
        Get notifications of a specific type with pagination.
        
        Args:
            db: Database session
            notification_type: Type of notifications to retrieve
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of notifications of the specified type
        """
        query = select(Notification).where(Notification.type == notification_type)
        query = query.order_by(Notification.created_at.desc())
        query = query.offset(skip).limit(limit)
        
        result = db.execute(query).scalars().all()
        return list(result)
    
    def delete_old_notifications(self, db: Session, older_than: datetime) -> int:
        """
        Delete notifications older than a specified date.
        
        Args:
            db: Database session
            older_than: Delete notifications created before this date
            
        Returns:
            Number of notifications deleted
        """
        # Create delete query for old notifications
        query = delete(Notification).where(Notification.created_at < older_than)
        
        # Execute delete operation
        result = db.execute(query)
        db.commit()
        
        # Return number of deleted rows
        return result.rowcount


# Create singleton instance for application-wide use
notification = CRUDNotification()