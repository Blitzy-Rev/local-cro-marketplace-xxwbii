# standard library
import logging
from typing import List, Dict, Optional, Any, Tuple, Union

# redis 4.5+
import redis

# fastapi 0.95+
from fastapi import HTTPException

from ..crud import crud_notification as notification
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
from ..exceptions import ValidationException, ResourceNotFoundException

# Set up logger for this module
logger = logging.getLogger(__name__)

# Redis client - will be initialized on first use
_redis_client = None

def get_redis_client():
    """
    Get a Redis client instance.
    
    This function initializes the Redis client lazily on first use.
    In a real implementation, it would use configuration from settings.
    """
    global _redis_client
    if _redis_client is None:
        # In a real implementation, this would use settings from config
        # For now, using default values that match those in celery_config.py
        _redis_client = redis.Redis(
            host='localhost',
            port=6379,
            db=0
        )
    return _redis_client

def publish_notification(user_id: int, notification_data: Dict[str, Any]) -> bool:
    """
    Publishes a notification to a Redis channel for real-time delivery.
    
    Args:
        user_id: User ID to publish notification for
        notification_data: Notification data to publish
        
    Returns:
        True if published successfully, False otherwise
    """
    try:
        redis_client = get_redis_client()
        channel = f"user:{user_id}:notifications"
        redis_client.publish(channel, json.dumps(notification_data))
        logger.debug(f"Published notification to channel {channel}")
        return True
    except Exception as e:
        logger.error(f"Error publishing notification to Redis: {str(e)}")
        # Don't let Redis failures stop the notification process
        return False

@app.task(name='notification_tasks.deliver_notification', queue='notifications')
def deliver_notification(user_id: int, notification_type: str, message: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Delivers a notification to a user through the appropriate channels.
    
    Args:
        user_id: ID of the user to deliver the notification to
        notification_type: Type of notification (from NotificationType enum)
        message: Notification message content
        data: Optional additional structured data related to the notification
        
    Returns:
        Created notification data
    """
    logger.info(f"Delivering notification of type {notification_type} to user {user_id}")
    
    # Create notification in database
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
        
        # Publish to Redis for real-time delivery
        publish_notification(user_id, notification_dict)
        
        logger.info(f"Notification delivered successfully to user {user_id}")
        
        # Return notification data
        return notification_dict

@app.task(name='notification_tasks.batch_deliver_notifications', queue='notifications')
def batch_deliver_notifications(user_ids: List[int], notification_type: str, message: str, data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Delivers the same notification to multiple users.
    
    Args:
        user_ids: List of user IDs to deliver the notification to
        notification_type: Type of notification (from NotificationType enum)
        message: Notification message content
        data: Optional additional structured data related to the notification
        
    Returns:
        Result with count of notifications delivered
    """
    logger.info(f"Batch delivering notification of type {notification_type} to {len(user_ids)} users")
    
    results = []
    
    # Deliver notification to each user asynchronously
    for user_id in user_ids:
        # Use .delay() to create an asynchronous task for each user
        task = deliver_notification.delay(user_id, notification_type, message, data)
        results.append(task.id)
    
    logger.info(f"Batch notification delivery initiated for {len(user_ids)} users")
    
    return {
        "count": len(user_ids),
        "task_ids": results
    }

@app.task(name='notification_tasks.mark_notification_read', queue='notifications')
def mark_notification_read(notification_id: int, user_id: int) -> Dict[str, Any]:
    """
    Marks a notification as read.
    
    Args:
        notification_id: ID of the notification to mark as read
        user_id: ID of the user who owns the notification
        
    Returns:
        Updated notification data
    """
    logger.info(f"Marking notification {notification_id} as read for user {user_id}")
    
    with get_db() as db:
        # Get notification by ID
        db_notification = notification.get(db, notification_id)
        
        # Check if notification exists and belongs to the user
        if db_notification and db_notification.user_id == user_id:
            # Mark notification as read
            updated_notification = notification.mark_as_read(db, db_notification)
            
            # Convert DB model to dictionary
            result = {
                "id": updated_notification.id,
                "user_id": updated_notification.user_id,
                "type": updated_notification.type,
                "message": updated_notification.message,
                "data": updated_notification.data,
                "read_status": updated_notification.read_status,
                "created_at": updated_notification.created_at.isoformat() if updated_notification.created_at else None,
                "read_at": updated_notification.read_at.isoformat() if updated_notification.read_at else None
            }
            
            logger.info(f"Notification {notification_id} marked as read")
            return result
        else:
            logger.error(f"Notification {notification_id} not found or does not belong to user {user_id}")
            return {"error": "Notification not found or access denied"}

@app.task(name='notification_tasks.cleanup_old_notifications', queue='notifications')
def cleanup_old_notifications(days: int = 90) -> Dict[str, Any]:
    """
    Deletes notifications older than the specified days.
    
    Args:
        days: Number of days to keep notifications (default: 90)
        
    Returns:
        Result with count of notifications deleted
    """
    logger.info(f"Cleaning up notifications older than {days} days")
    
    # Calculate cutoff date
    cutoff_date = datetime.utcnow() - timedelta(days=days)
    
    with get_db() as db:
        # Delete old notifications
        deleted_count = notification.delete_old_notifications(db, cutoff_date)
        
        logger.info(f"Deleted {deleted_count} old notifications")
        
        return {
            "count": deleted_count,
            "cutoff_date": cutoff_date.isoformat()
        }

@app.task(name='notification_tasks.send_experiment_status_notification', queue='notifications')
def send_experiment_status_notification(user_id: int, experiment_id: int, status: str, experiment_name: str) -> Dict[str, Any]:
    """
    Sends a notification about experiment status change.
    
    Args:
        user_id: ID of the user to notify
        experiment_id: ID of the experiment
        status: New status of the experiment
        experiment_name: Name of the experiment
        
    Returns:
        Created notification data
    """
    logger.info(f"Sending experiment status notification for experiment {experiment_id} to user {user_id}")
    
    # Format message
    message = f"Experiment '{experiment_name}' status changed to {status}"
    
    # Create data dictionary
    data = {
        "experiment_id": experiment_id,
        "status": status
    }
    
    # Deliver notification
    return deliver_notification(user_id, NotificationType.EXPERIMENT_STATUS_CHANGE.name, message, data)

@app.task(name='notification_tasks.send_submission_notification', queue='notifications')
def send_submission_notification(cro_user_id: int, submission_id: int, experiment_name: str, pharma_company: str) -> Dict[str, Any]:
    """
    Sends a notification about new submission to CRO.
    
    Args:
        cro_user_id: ID of the CRO user to notify
        submission_id: ID of the submission
        experiment_name: Name of the experiment
        pharma_company: Name of the pharmaceutical company
        
    Returns:
        Created notification data
    """
    logger.info(f"Sending submission notification for submission {submission_id} to CRO user {cro_user_id}")
    
    # Format message
    message = f"New submission from {pharma_company} for experiment '{experiment_name}'"
    
    # Create data dictionary
    data = {
        "submission_id": submission_id,
        "experiment_name": experiment_name
    }
    
    # Deliver notification
    return deliver_notification(cro_user_id, NotificationType.SUBMISSION_CREATED.name, message, data)

@app.task(name='notification_tasks.send_quote_notification', queue='notifications')
def send_quote_notification(pharma_user_id: int, submission_id: int, experiment_name: str, cro_name: str, price: float) -> Dict[str, Any]:
    """
    Sends a notification about quote provided by CRO.
    
    Args:
        pharma_user_id: ID of the pharma user to notify
        submission_id: ID of the submission
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        price: Quoted price
        
    Returns:
        Created notification data
    """
    logger.info(f"Sending quote notification for submission {submission_id} to pharma user {pharma_user_id}")
    
    # Format message
    message = f"Quote received from {cro_name} for experiment '{experiment_name}': ${price:.2f}"
    
    # Create data dictionary
    data = {
        "submission_id": submission_id,
        "experiment_name": experiment_name,
        "price": price
    }
    
    # Deliver notification
    return deliver_notification(pharma_user_id, NotificationType.QUOTE_PROVIDED.name, message, data)

@app.task(name='notification_tasks.send_results_notification', queue='notifications')
def send_results_notification(pharma_user_id: int, submission_id: int, result_id: int, experiment_name: str, cro_name: str) -> Dict[str, Any]:
    """
    Sends a notification about results uploaded by CRO.
    
    Args:
        pharma_user_id: ID of the pharma user to notify
        submission_id: ID of the submission
        result_id: ID of the result
        experiment_name: Name of the experiment
        cro_name: Name of the CRO
        
    Returns:
        Created notification data
    """
    logger.info(f"Sending results notification for submission {submission_id} to pharma user {pharma_user_id}")
    
    # Format message
    message = f"Results uploaded by {cro_name} for experiment '{experiment_name}'"
    
    # Create data dictionary
    data = {
        "submission_id": submission_id,
        "result_id": result_id,
        "experiment_name": experiment_name
    }
    
    # Deliver notification
    return deliver_notification(pharma_user_id, NotificationType.RESULTS_UPLOADED.name, message, data)

@app.task(name='notification_tasks.send_system_notification', queue='notifications')
def send_system_notification(message: str, user_ids: List[int], data: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Sends a system alert notification to multiple users.
    
    Args:
        message: System alert message
        user_ids: List of user IDs to notify
        data: Optional additional structured data
        
    Returns:
        Result with count of notifications sent
    """
    logger.info(f"Sending system notification to {len(user_ids)} users")
    
    # Deliver notification to multiple users
    return batch_deliver_notifications(user_ids, NotificationType.SYSTEM_ALERT.name, message, data)