import pytest
from unittest.mock import MagicMock, patch
from datetime import datetime

# Import notification service functions for testing
from ../../app/services/notification_service import (
    create_notification, get_notifications, get_unread_notifications,
    count_notifications, mark_notification_as_read, mark_all_notifications_as_read,
    send_experiment_status_change, send_submission_created, send_quote_provided,
    send_results_uploaded, send_system_alert, cleanup_old_notifications
)

# Import notification CRUD operations for mocking
from ../../app/crud/crud_notification import notification

# Import notification type enumeration for testing
from ../../app/constants import NotificationType

# Import notification schemas for testing
from ../../app/schemas/notification import NotificationCreate, NotificationFilter, NotificationList, NotificationCount

# Import notification tasks for mocking
from ../../app/worker/tasks/notification_tasks import (
    deliver_notification, batch_deliver_notifications, send_experiment_status_notification,
    send_submission_notification, send_quote_notification, send_results_notification,
    send_system_notification, mark_notification_read, cleanup_old_notifications as cleanup_old_notifications_task
)


def test_create_notification(db_session):
    """Test creating a notification."""
    # Set up test data with user_id, notification_type, message, and data
    user_id = 1
    notification_type = NotificationType.SYSTEM_ALERT.name
    message = "Test notification"
    data = {"key": "value"}
    
    # Mock notification.create_notification to return a mock notification
    mock_notification = MagicMock()
    mock_notification.id = 1
    mock_notification.user_id = user_id
    mock_notification.type = notification_type
    mock_notification.message = message
    mock_notification.data = data
    mock_notification.read_status = False
    mock_notification.created_at = datetime.utcnow()
    
    with patch.object(notification, 'create_notification', return_value=mock_notification) as mock_create:
        # Call create_notification with test data
        result = create_notification(user_id, notification_type, message, data)
        
        # Assert notification.create_notification was called with correct parameters
        mock_create.assert_called_once()
        create_call_args = mock_create.call_args[0]
        assert len(create_call_args) == 2  # db and notification_data
        notification_data = create_call_args[1]
        assert isinstance(notification_data, NotificationCreate)
        assert notification_data.user_id == user_id
        assert notification_data.type == notification_type
        assert notification_data.message == message
        assert notification_data.data == data
        
        # Assert the function returns the expected notification data
        assert result["id"] == mock_notification.id
        assert result["user_id"] == user_id
        assert result["type"] == notification_type
        assert result["message"] == message
        assert result["data"] == data
        assert result["read_status"] == False
        assert "created_at" in result


def test_get_notifications(db_session):
    """Test retrieving notifications for a user."""
    # Set up test data with user_id, skip, limit, and filters
    user_id = 1
    skip = 0
    limit = 20
    filters = NotificationFilter(read_status=False)
    
    # Mock notification.get_by_user to return a list of mock notifications
    mock_notifications = [MagicMock() for _ in range(3)]
    total_count = len(mock_notifications)
    
    with patch.object(notification, 'get_by_user', side_effect=[mock_notifications, total_count]) as mock_get:
        # Call get_notifications with test data
        result = get_notifications(user_id, skip, limit, filters)
        
        # Assert notification.get_by_user was called with correct parameters
        assert mock_get.call_count == 2
        
        # First call - get notifications
        first_call_args = mock_get.call_args_list[0][0]
        assert first_call_args[1] == user_id
        assert first_call_args[2] == skip
        assert first_call_args[3] == limit
        assert first_call_args[4] == filters
        
        # Second call - get count
        second_call_kwargs = mock_get.call_args_list[1][1]
        assert 'count' in second_call_kwargs
        assert second_call_kwargs['count'] == True
        
        # Assert the function returns a NotificationList with correct items, total, page, and size
        assert isinstance(result, NotificationList)
        assert result.items == mock_notifications
        assert result.total == total_count
        assert result.page == 1
        assert result.size == limit
        assert result.pages == 1  # Since total is 3 and limit is 20, we have 1 page


def test_get_unread_notifications(db_session):
    """Test retrieving unread notifications for a user."""
    # Set up test data with user_id, skip, limit, and filters
    user_id = 1
    skip = 0
    limit = 20
    filters = NotificationFilter(read_status=False)
    
    # Mock notification.get_unread_by_user to return a list of mock notifications
    mock_notifications = [MagicMock() for _ in range(3)]
    total_count = len(mock_notifications)
    
    with patch.object(notification, 'get_unread_by_user', side_effect=[mock_notifications, total_count]) as mock_get_unread:
        # Call get_unread_notifications with test data
        result = get_unread_notifications(user_id, skip, limit, filters)
        
        # Assert notification.get_unread_by_user was called with correct parameters
        assert mock_get_unread.call_count == 2
        
        # First call - get notifications
        first_call_args = mock_get_unread.call_args_list[0][0]
        assert first_call_args[1] == user_id
        assert first_call_args[2] == skip
        assert first_call_args[3] == limit
        assert first_call_args[4] == filters
        
        # Second call - get count
        second_call_kwargs = mock_get_unread.call_args_list[1][1]
        assert 'count' in second_call_kwargs
        assert second_call_kwargs['count'] == True
        
        # Assert the function returns a NotificationList with correct items, total, page, and size
        assert isinstance(result, NotificationList)
        assert result.items == mock_notifications
        assert result.total == total_count
        assert result.page == 1
        assert result.size == limit
        assert result.pages == 1  # Since total is 3 and limit is 20, we have 1 page


def test_count_notifications(db_session):
    """Test counting total and unread notifications for a user."""
    # Set up test data with user_id
    user_id = 1
    
    # Mock notification.get_by_user with count=True to return total count
    # Mock notification.count_unread_by_user to return unread count
    total_count = 10
    unread_count = 5
    
    with patch.object(notification, 'get_by_user', return_value=total_count) as mock_get:
        with patch.object(notification, 'count_unread_by_user', return_value=unread_count) as mock_count_unread:
            # Call count_notifications with test data
            result = count_notifications(user_id)
            
            # Assert notification.get_by_user was called with correct parameters
            mock_get.assert_called_once()
            get_call_args = mock_get.call_args[0]
            assert get_call_args[1] == user_id
            assert mock_get.call_args[1]['count'] == True
            
            # Assert notification.count_unread_by_user was called with correct parameters
            mock_count_unread.assert_called_once()
            count_call_args = mock_count_unread.call_args[0]
            assert count_call_args[1] == user_id
            
            # Assert the function returns a NotificationCount with correct total and unread counts
            assert isinstance(result, NotificationCount)
            assert result.total == total_count
            assert result.unread == unread_count


def test_mark_notification_as_read(db_session):
    """Test marking a notification as read."""
    # Set up test data with notification_id and user_id
    notification_id = 1
    user_id = 1
    
    # Mock mark_notification_read.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.mark_notification_read.delay', return_value=mock_task) as mock_mark_read:
        # Call mark_notification_as_read with test data
        result = mark_notification_as_read(notification_id, user_id)
        
        # Assert mark_notification_read.delay was called with correct parameters
        mock_mark_read.assert_called_once_with(notification_id, user_id)
        
        # Assert the function returns a success response with notification_id
        assert result["notification_id"] == notification_id
        assert result["status"] == "success"


def test_mark_all_notifications_as_read(db_session):
    """Test marking all notifications as read for a user."""
    # Set up test data with user_id
    user_id = 1
    
    # Mock notification.mark_all_as_read to return a count of marked notifications
    count = 5
    
    with patch.object(notification, 'mark_all_as_read', return_value=count) as mock_mark_all:
        # Call mark_all_notifications_as_read with test data
        result = mark_all_notifications_as_read(user_id)
        
        # Assert notification.mark_all_as_read was called with correct parameters
        mock_mark_all.assert_called_once()
        mark_all_args = mock_mark_all.call_args[0]
        assert mark_all_args[1] == user_id
        
        # Assert the function returns a success response with count of notifications marked as read
        assert result["user_id"] == user_id
        assert result["count"] == count
        assert result["status"] == "success"


def test_send_experiment_status_change(db_session):
    """Test sending a notification about experiment status change."""
    # Set up test data with user_id, experiment_id, status, and experiment_name
    user_id = 1
    experiment_id = 1
    status = "COMPLETED"
    experiment_name = "Test Experiment"
    
    # Mock send_experiment_status_notification.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.send_experiment_status_notification.delay', return_value=mock_task) as mock_send:
        # Call send_experiment_status_change with test data
        result = send_experiment_status_change(user_id, experiment_id, status, experiment_name)
        
        # Assert send_experiment_status_notification.delay was called with correct parameters
        mock_send.assert_called_once_with(user_id, experiment_id, status, experiment_name)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"


def test_send_submission_created(db_session):
    """Test sending a notification about new submission to CRO."""
    # Set up test data with cro_user_id, submission_id, experiment_name, and pharma_company
    cro_user_id = 1
    submission_id = 1
    experiment_name = "Test Experiment"
    pharma_company = "Test Pharma"
    
    # Mock send_submission_notification.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.send_submission_notification.delay', return_value=mock_task) as mock_send:
        # Call send_submission_created with test data
        result = send_submission_created(cro_user_id, submission_id, experiment_name, pharma_company)
        
        # Assert send_submission_notification.delay was called with correct parameters
        mock_send.assert_called_once_with(cro_user_id, submission_id, experiment_name, pharma_company)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"


def test_send_quote_provided(db_session):
    """Test sending a notification about quote provided by CRO."""
    # Set up test data with pharma_user_id, submission_id, experiment_name, cro_name, and price
    pharma_user_id = 1
    submission_id = 1
    experiment_name = "Test Experiment"
    cro_name = "Test CRO"
    price = 1000.0
    
    # Mock send_quote_notification.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.send_quote_notification.delay', return_value=mock_task) as mock_send:
        # Call send_quote_provided with test data
        result = send_quote_provided(pharma_user_id, submission_id, experiment_name, cro_name, price)
        
        # Assert send_quote_notification.delay was called with correct parameters
        mock_send.assert_called_once_with(pharma_user_id, submission_id, experiment_name, cro_name, price)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"


def test_send_results_uploaded(db_session):
    """Test sending a notification about results uploaded by CRO."""
    # Set up test data with pharma_user_id, submission_id, result_id, experiment_name, and cro_name
    pharma_user_id = 1
    submission_id = 1
    result_id = 1
    experiment_name = "Test Experiment"
    cro_name = "Test CRO"
    
    # Mock send_results_notification.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.send_results_notification.delay', return_value=mock_task) as mock_send:
        # Call send_results_uploaded with test data
        result = send_results_uploaded(pharma_user_id, submission_id, result_id, experiment_name, cro_name)
        
        # Assert send_results_notification.delay was called with correct parameters
        mock_send.assert_called_once_with(pharma_user_id, submission_id, result_id, experiment_name, cro_name)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"


def test_send_system_alert(db_session):
    """Test sending a system alert notification to multiple users."""
    # Set up test data with message, user_ids, and data
    message = "System maintenance scheduled"
    user_ids = [1, 2, 3]
    data = {"start_time": "2023-06-01T10:00:00", "end_time": "2023-06-01T12:00:00"}
    
    # Mock send_system_notification.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.send_system_notification.delay', return_value=mock_task) as mock_send:
        # Call send_system_alert with test data
        result = send_system_alert(message, user_ids, data)
        
        # Assert send_system_notification.delay was called with correct parameters
        mock_send.assert_called_once_with(message, user_ids, data)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"


def test_cleanup_old_notifications(db_session):
    """Test cleaning up old notifications."""
    # Set up test data with days parameter
    days = 30
    
    # Mock cleanup_old_notifications.delay to return a mock task
    mock_task = MagicMock()
    mock_task.id = "task-123"
    
    with patch('../../app/worker/tasks/notification_tasks.cleanup_old_notifications.delay', return_value=mock_task) as mock_cleanup:
        # Call cleanup_old_notifications with test data
        result = cleanup_old_notifications(days)
        
        # Assert cleanup_old_notifications.delay was called with correct parameters
        mock_cleanup.assert_called_once_with(days)
        
        # Assert the function returns task information with task_id
        assert result["task_id"] == mock_task.id
        assert result["status"] == "pending"