"""
Notification service for the Molecular Data Management and CRO Integration Platform.

This module provides high-level functions for creating, managing, and delivering
notifications to users about system events such as experiment status changes,
CRO submissions, quotes, and result uploads. It acts as a facade over the
notification CRUD operations and background tasks.
"""

import logging  # standard library
from typing import Dict, List, Any, Optional, Union  # standard library
from datetime import datetime  # standard library

import redis  # redis 4.5+
from fastapi import HTTPException  # fastapi 0.95+

from ..crud.crud_notification import notification
from ..constants import NotificationType
from ..schemas.notification import (
    NotificationCreate, 
    NotificationUpdate, 
    NotificationFilter, 
    NotificationList, 
    NotificationCount
)
from ..worker.tasks.notification_tasks import (
    deliver_notification,
    batch_deliver_notifications,
    send_experiment_status_notification,
    send_submission_notification,
    send_quote_notification,
    send_results_notification,
    send_system_notification,
    mark_notification_read,
    cleanup_old_notifications as cleanup_old_notifications_task
)
from ..db.session import get_db

# Set up logger for this module
logger = logging.getLogger(__name__)

def create_notification(user_id: int, notification_type: str, message: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Create a new notification in the database.
    
    Args:
        user_id: ID of the user to receive the notification
        notification_type: Type of notification (from NotificationType enum)
        message: Notification message content
        data: Optional additional structured data related to the notification
        
    Returns:
        Created notification data
    """
    logger.info(f"Creating notification of type {notification_type} for user {user_id}")
    
    with get_db() as db:
        # Create notification object using Pydantic schema
        notification_data = NotificationCreate(
            user_id=user_id,
            type=notification_type,
            message=message,
            data=data
        )
        
        # Create notification in database
        db_notification = notification.create_notification(db, notification_data)
        
        # Convert DB model to dictionary
        notification_dict = {
            "id": db_notification.id,
            "user_id": db_notification.user_id,
            "type": db_notification.type,
            "message": db_notification.message,
            "data": db_notification.data,
            "read_status": db_notification.read_status,
            "created_at": db_notification.created_at.isoformat() if db_notification.created_at else None
        }
        
        logger.info(f"Notification created with ID {db_notification.id}")
        
        return notification_dict

def get_notifications(user_id: int, skip: int = 0, limit: int = 20, filters: Optional[NotificationFilter] = None) -> NotificationList:
    """
    Get notifications for a user with pagination and optional filtering.
    
    Args:
        user_id: ID of the user whose notifications to retrieve
        skip: Number of records to skip (for pagination)
        limit: Maximum number of records to return
        filters: Optional filter criteria
        
    Returns:
        Paginated list of notifications
    """
    logger.debug(f"Getting notifications for user {user_id} (skip={skip}, limit={limit})")
    
    with get_db() as db:
        # Get notifications
        notifications = notification.get_by_user(db, user_id, skip, limit, filters)
        
        # Get total count for pagination
        total = notification.get_by_user(db, user_id, count=True, filters=filters)
        
        # Calculate total pages
        pages = (total + limit - 1) // limit if limit > 0 else 0
        
        # Create response object
        result = NotificationList(
            items=notifications,
            total=total,
            page=(skip // limit) + 1 if limit > 0 else 1,
            size=limit,
            pages=pages
        )
        
        return result

def get_unread_notifications(user_id: int, skip: int = 0, limit: int = 20, filters: Optional[NotificationFilter] = None) -> NotificationList:
    """
    Get unread notifications for a user with pagination and optional filtering.
    
    Args:
        user_id: ID of the user whose unread notifications to retrieve
        skip: Number of records to skip (for pagination)
        limit: Maximum number of records to return
        filters: Optional filter criteria
        
    Returns:
        Paginated list of unread notifications
    """
    logger.debug(f"Getting unread notifications for user {user_id} (skip={skip}, limit={limit})")
    
    with get_db() as db:
        # Get unread notifications
        notifications = notification.get_unread_by_user(db, user_id, skip, limit, filters)
        
        # Get total count for pagination
        total = notification.get_unread_by_user(db, user_id, count=True, filters=filters)
        
        # Calculate total pages
        pages = (total + limit - 1) // limit if limit > 0 else 0
        
        # Create response object
        result = NotificationList(
            items=notifications,
            total=total,
            page=(skip // limit) + 1 if limit > 0 else 1,
            size=limit,
            pages=pages
        )
        
        return result

def count_notifications(user_id: int) -> NotificationCount:
    """
    Count total and unread notifications for a user.
    
    Args:
        user_id: ID of the user whose notifications to count
        
    Returns:
        Count of total and unread notifications
    """
    logger.debug(f"Counting notifications for user {user_id}")
    
    with get_db() as db:
        # Get total count
        total = notification.get_by_user(db, user_id, count=True)
        
        # Get unread count
        unread = notification.count_unread_by_user(db, user_id)
        
        # Create response object
        result = NotificationCount(
            total=total,
            unread=unread
        )
        
        return result

def mark_notification_as_read(notification_id: int, user_id: int) -> Dict[str, Any]:
    """
    Mark a notification as read.
    
    Args:
        notification_id: ID of the notification to mark as read
        user_id: ID of the user who owns the notification
        
    Returns:
        Updated notification data
    """
    logger.info(f"Marking notification {notification_id} as read for user {user_id}")
    
    # Use the task to mark the notification as read asynchronously
    mark_notification_read.delay(notification_id, user_id)
    
    return {
        "notification_id": notification_id,
        "status": "success"
    }

def mark_all_notifications_as_read(user_id: int) -> Dict[str, Any]:
    """
    Mark all unread notifications for a user as read.
    
    Args:
        user_id: ID of the user whose notifications to mark as read
        
    Returns:
        Result with count of notifications marked as read
    """
    logger.info(f"Marking all notifications as read for user {user_id}")
    
    with get_db() as db:
        # Mark all notifications as read
        count = notification.mark_all_as_read(db, user_id)
        
        return {
            "user_id": user_id,
            "count": count,
            "status": "success"
        }

def send_experiment_status_change(user_id: int, experiment_id: int, status: str, experiment_name: str) -> Dict[str, Any]:
    """
    Send a notification about experiment status change.
    
    Args:
        user_id: ID of the user to notify
        experiment_id: ID of the experiment
        status: New status of the experiment
        experiment_name: Name of the experiment
        
    Returns:
        Task information
    """
    logger.info(f"Sending experiment status notification for experiment {experiment_id} to user {user_id}")
    
    # Use the task to send the notification asynchronously
    task = send_experiment_status_notification.delay(user_id, experiment_id, status, experiment_name)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }

def send_submission_created(cro_user_id: int, submission_id: int, experiment_name: str, pharma_company: str) -> Dict[str, Any]:
    """
    Send a notification about new submission to CRO.
    
    Args:
        cro_user_id: ID of the CRO user to notify
        submission_id: ID of the submission
        experiment_name: Name of the experiment
        pharma_company: Name of the pharmaceutical company
        
    Returns:
        Task information
    """
    logger.info(f"Sending submission notification for submission {submission_id} to CRO user {cro_user_id}")
    
    # Use the task to send the notification asynchronously
    task = send_submission_notification.delay(cro_user_id, submission_id, experiment_name, pharma_company)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }

def send_quote_provided(pharma_user_id: int, submission_id: int, experiment_name: str, cro_name: str, price: float) -> Dict[str, Any]:
    """
    Send a notification about quote provided by CRO.
    
    Args:
        pharma_user_id: ID of the pharma user to notify
        submission_id: ID of the submission
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        price: Quoted price
        
    Returns:
        Task information
    """
    logger.info(f"Sending quote notification for submission {submission_id} to pharma user {pharma_user_id}")
    
    # Use the task to send the notification asynchronously
    task = send_quote_notification.delay(pharma_user_id, submission_id, experiment_name, cro_name, price)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }

def send_results_uploaded(pharma_user_id: int, submission_id: int, result_id: int, experiment_name: str, cro_name: str) -> Dict[str, Any]:
    """
    Send a notification about results uploaded by CRO.
    
    Args:
        pharma_user_id: ID of the pharma user to notify
        submission_id: ID of the submission
        result_id: ID of the result
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        
    Returns:
        Task information
    """
    logger.info(f"Sending results notification for submission {submission_id} to pharma user {pharma_user_id}")
    
    # Use the task to send the notification asynchronously
    task = send_results_notification.delay(pharma_user_id, submission_id, result_id, experiment_name, cro_name)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }

def send_system_alert(message: str, user_ids: List[int], data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Send a system alert notification to multiple users.
    
    Args:
        message: System alert message
        user_ids: List of user IDs to notify
        data: Optional additional structured data
        
    Returns:
        Task information
    """
    logger.info(f"Sending system notification to {len(user_ids)} users")
    
    # Use the task to send the notification asynchronously
    task = send_system_notification.delay(message, user_ids, data)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }

def cleanup_old_notifications(days: int = 90) -> Dict[str, Any]:
    """
    Delete notifications older than the specified days.
    
    Args:
        days: Number of days to keep notifications (default: 90)
        
    Returns:
        Task information
    """
    logger.info(f"Cleaning up notifications older than {days} days")
    
    # Use the task to clean up old notifications asynchronously
    task = cleanup_old_notifications_task.delay(days)
    
    return {
        "task_id": task.id,
        "status": "pending"
    }