"""
Unit tests for the notification task functions in the Celery worker module.

This test file verifies the correct behavior of asynchronous notification delivery,
batch notifications, notification status updates, and specialized notification types
for different system events.
"""

import pytest
from unittest.mock import MagicMock, patch, Mock
from datetime import datetime, timedelta

from ...app.worker.tasks.notification_tasks import (
    deliver_notification,
    batch_deliver_notifications,
    mark_notification_read,
    cleanup_old_notifications,
    send_experiment_status_notification,
    send_submission_notification,
    send_quote_notification,
    send_results_notification,
    send_system_notification
)
from ...app.crud.crud_notification import notification
from ...app.constants import NotificationType
from ...app.schemas.notification import NotificationCreate


@pytest.mark.worker
def test_deliver_notification(db_session, test_pharma_user):
    """Tests the deliver_notification function creates and stores a notification correctly"""
    # Set up test data
    user_id = test_pharma_user.id
    notification_type = NotificationType.EXPERIMENT_STATUS_CHANGE.name
    message = "Test notification message"
    data = {"test_key": "test_value"}

    # Call the function
    result = deliver_notification(user_id, notification_type, message, data)

    # Verify notification was created
    assert result["user_id"] == user_id
    assert result["type"] == notification_type
    assert result["message"] == message
    assert result["data"] == data
    assert result["read_status"] is False
    assert "created_at" in result
    
    # Verify it's in the database
    db_notification = notification.get(db_session, result["id"])
    assert db_notification is not None
    assert db_notification.user_id == user_id
    assert db_notification.type.name == notification_type
    assert db_notification.message == message
    assert db_notification.data == data
    assert db_notification.read_status is False
    assert db_notification.created_at is not None


@pytest.mark.worker
def test_batch_deliver_notifications(db_session, test_pharma_user, test_cro_user):
    """Tests the batch_deliver_notifications function delivers notifications to multiple users"""
    # Set up test data
    user_ids = [test_pharma_user.id, test_cro_user.id]
    notification_type = NotificationType.SYSTEM_ALERT.name
    message = "Test batch notification"
    data = {"batch": True}

    # Mock deliver_notification.delay to track calls
    with patch('...app.worker.tasks.notification_tasks.deliver_notification.delay') as mock_deliver:
        # Set up the mock
        mock_deliver.return_value = Mock()
        mock_deliver.return_value.id = "mock-task-id"
        
        # Call the function
        result = batch_deliver_notifications(user_ids, notification_type, message, data)
        
        # Verify deliver_notification.delay was called for each user
        assert mock_deliver.call_count == len(user_ids)
        
        # Verify correct parameters were passed
        calls = mock_deliver.call_args_list
        for i, user_id in enumerate(user_ids):
            args, kwargs = calls[i]
            assert args[0] == user_id
            assert args[1] == notification_type
            assert args[2] == message
            assert args[3] == data
            
        # Verify result contains correct count
        assert result["count"] == len(user_ids)
        assert len(result["task_ids"]) == len(user_ids)


@pytest.mark.worker
def test_mark_notification_read(db_session, test_pharma_user):
    """Tests the mark_notification_read function updates notification status correctly"""
    # Create a test notification
    notification_data = NotificationCreate(
        user_id=test_pharma_user.id,
        type=NotificationType.EXPERIMENT_STATUS_CHANGE.name,
        message="Test notification to mark as read",
        data=None
    )
    db_notification = notification.create_notification(db_session, notification_data)
    
    # Verify it starts with read_status=False
    assert db_notification.read_status is False
    assert db_notification.read_at is None
    
    # Call the function
    result = mark_notification_read(db_notification.id, test_pharma_user.id)
    
    # Verify notification is now marked as read
    updated_notification = notification.get(db_session, db_notification.id)
    assert updated_notification.read_status is True
    assert updated_notification.read_at is not None
    
    # Verify result contains updated status
    assert result["read_status"] is True
    assert "read_at" in result and result["read_at"] is not None


@pytest.mark.worker
def test_mark_notification_read_wrong_user(db_session, test_pharma_user, test_cro_user):
    """Tests the mark_notification_read function handles wrong user correctly"""
    # Create a test notification for pharma user
    notification_data = NotificationCreate(
        user_id=test_pharma_user.id,
        type=NotificationType.EXPERIMENT_STATUS_CHANGE.name,
        message="Test notification for wrong user",
        data=None
    )
    db_notification = notification.create_notification(db_session, notification_data)
    
    # Call function with CRO user ID (wrong user)
    result = mark_notification_read(db_notification.id, test_cro_user.id)
    
    # Verify notification is still unread
    updated_notification = notification.get(db_session, db_notification.id)
    assert updated_notification.read_status is False
    assert updated_notification.read_at is None
    
    # Verify result contains error
    assert "error" in result


@pytest.mark.worker
def test_cleanup_old_notifications(db_session, test_pharma_user):
    """Tests the cleanup_old_notifications function deletes old notifications"""
    # Create old notifications (91 days old)
    old_date = datetime.utcnow() - timedelta(days=91)
    
    # Manually set created_at to old date for test
    old_notification = notification.create_notification(
        db_session, 
        NotificationCreate(
            user_id=test_pharma_user.id,
            type=NotificationType.SYSTEM_ALERT.name,
            message="Old notification to clean up",
            data=None
        )
    )
    old_notification.created_at = old_date
    db_session.add(old_notification)
    db_session.commit()
    
    # Create recent notification
    recent_notification = notification.create_notification(
        db_session, 
        NotificationCreate(
            user_id=test_pharma_user.id,
            type=NotificationType.SYSTEM_ALERT.name,
            message="Recent notification to keep",
            data=None
        )
    )
    
    # Call cleanup function with 90 days
    result = cleanup_old_notifications(90)
    
    # Verify old notification is gone
    assert notification.get(db_session, old_notification.id) is None
    
    # Verify recent notification is still there
    assert notification.get(db_session, recent_notification.id) is not None
    
    # Verify result contains correct count
    assert result["count"] == 1


@pytest.mark.worker
def test_send_experiment_status_notification(db_session, test_pharma_user):
    """Tests the send_experiment_status_notification function creates correct notification"""
    # Set up test data
    user_id = test_pharma_user.id
    experiment_id = 123
    status = "IN_PROGRESS"
    experiment_name = "Test Experiment"
    
    # Mock deliver_notification to verify calls
    with patch('...app.worker.tasks.notification_tasks.deliver_notification') as mock_deliver:
        # Set up mock return value
        mock_notification = {
            "id": 1,
            "user_id": user_id,
            "type": NotificationType.EXPERIMENT_STATUS_CHANGE.name,
            "message": f"Experiment '{experiment_name}' status changed to {status}",
            "data": {"experiment_id": experiment_id, "status": status},
            "read_status": False,
            "created_at": datetime.utcnow().isoformat()
        }
        mock_deliver.return_value = mock_notification
        
        # Call the function
        result = send_experiment_status_notification(user_id, experiment_id, status, experiment_name)
        
        # Verify deliver_notification was called with correct parameters
        mock_deliver.assert_called_once_with(
            user_id, 
            NotificationType.EXPERIMENT_STATUS_CHANGE.name,
            f"Experiment '{experiment_name}' status changed to {status}",
            {"experiment_id": experiment_id, "status": status}
        )
        
        # Verify result matches mock return value
        assert result == mock_notification


@pytest.mark.worker
def test_send_submission_notification(db_session, test_cro_user):
    """Tests the send_submission_notification function creates correct notification"""
    # Set up test data
    cro_user_id = test_cro_user.id
    submission_id = 456
    experiment_name = "Test Experiment"
    pharma_company = "Test Pharma"
    
    # Mock deliver_notification to verify calls
    with patch('...app.worker.tasks.notification_tasks.deliver_notification') as mock_deliver:
        # Set up mock return value
        mock_notification = {
            "id": 2,
            "user_id": cro_user_id,
            "type": NotificationType.SUBMISSION_CREATED.name,
            "message": f"New submission from {pharma_company} for experiment '{experiment_name}'",
            "data": {"submission_id": submission_id, "experiment_name": experiment_name},
            "read_status": False,
            "created_at": datetime.utcnow().isoformat()
        }
        mock_deliver.return_value = mock_notification
        
        # Call the function
        result = send_submission_notification(cro_user_id, submission_id, experiment_name, pharma_company)
        
        # Verify deliver_notification was called with correct parameters
        mock_deliver.assert_called_once_with(
            cro_user_id, 
            NotificationType.SUBMISSION_CREATED.name,
            f"New submission from {pharma_company} for experiment '{experiment_name}'",
            {"submission_id": submission_id, "experiment_name": experiment_name}
        )
        
        # Verify result matches mock return value
        assert result == mock_notification


@pytest.mark.worker
def test_send_quote_notification(db_session, test_pharma_user):
    """Tests the send_quote_notification function creates correct notification"""
    # Set up test data
    pharma_user_id = test_pharma_user.id
    submission_id = 789
    experiment_name = "Test Experiment"
    cro_name = "Test CRO"
    price = 1250.00
    
    # Mock deliver_notification to verify calls
    with patch('...app.worker.tasks.notification_tasks.deliver_notification') as mock_deliver:
        # Set up mock return value
        mock_notification = {
            "id": 3,
            "user_id": pharma_user_id,
            "type": NotificationType.QUOTE_PROVIDED.name,
            "message": f"Quote received from {cro_name} for experiment '{experiment_name}': ${price:.2f}",
            "data": {"submission_id": submission_id, "experiment_name": experiment_name, "price": price},
            "read_status": False,
            "created_at": datetime.utcnow().isoformat()
        }
        mock_deliver.return_value = mock_notification
        
        # Call the function
        result = send_quote_notification(pharma_user_id, submission_id, experiment_name, cro_name, price)
        
        # Verify deliver_notification was called with correct parameters
        mock_deliver.assert_called_once_with(
            pharma_user_id, 
            NotificationType.QUOTE_PROVIDED.name,
            f"Quote received from {cro_name} for experiment '{experiment_name}': ${price:.2f}",
            {"submission_id": submission_id, "experiment_name": experiment_name, "price": price}
        )
        
        # Verify result matches mock return value
        assert result == mock_notification


@pytest.mark.worker
def test_send_results_notification(db_session, test_pharma_user):
    """Tests the send_results_notification function creates correct notification"""
    # Set up test data
    pharma_user_id = test_pharma_user.id
    submission_id = 101
    result_id = 202
    experiment_name = "Test Experiment"
    cro_name = "Test CRO"
    
    # Mock deliver_notification to verify calls
    with patch('...app.worker.tasks.notification_tasks.deliver_notification') as mock_deliver:
        # Set up mock return value
        mock_notification = {
            "id": 4,
            "user_id": pharma_user_id,
            "type": NotificationType.RESULTS_UPLOADED.name,
            "message": f"Results uploaded by {cro_name} for experiment '{experiment_name}'",
            "data": {"submission_id": submission_id, "result_id": result_id, "experiment_name": experiment_name},
            "read_status": False,
            "created_at": datetime.utcnow().isoformat()
        }
        mock_deliver.return_value = mock_notification
        
        # Call the function
        result = send_results_notification(pharma_user_id, submission_id, result_id, experiment_name, cro_name)
        
        # Verify deliver_notification was called with correct parameters
        mock_deliver.assert_called_once_with(
            pharma_user_id, 
            NotificationType.RESULTS_UPLOADED.name,
            f"Results uploaded by {cro_name} for experiment '{experiment_name}'",
            {"submission_id": submission_id, "result_id": result_id, "experiment_name": experiment_name}
        )
        
        # Verify result matches mock return value
        assert result == mock_notification


@pytest.mark.worker
def test_send_system_notification(db_session, test_pharma_user, test_cro_user):
    """Tests the send_system_notification function sends notifications to multiple users"""
    # Set up test data
    user_ids = [test_pharma_user.id, test_cro_user.id]
    message = "System maintenance scheduled"
    data = {"maintenance_time": "2023-07-01T12:00:00Z"}
    
    # Mock batch_deliver_notifications to verify calls
    with patch('...app.worker.tasks.notification_tasks.batch_deliver_notifications') as mock_batch:
        # Set up mock return value
        mock_result = {
            "count": len(user_ids),
            "task_ids": ["task-1", "task-2"]
        }
        mock_batch.return_value = mock_result
        
        # Call the function
        result = send_system_notification(message, user_ids, data)
        
        # Verify batch_deliver_notifications was called with correct parameters
        mock_batch.assert_called_once_with(
            user_ids,
            NotificationType.SYSTEM_ALERT.name,
            message,
            data
        )
        
        # Verify result matches mock return value
        assert result == mock_result