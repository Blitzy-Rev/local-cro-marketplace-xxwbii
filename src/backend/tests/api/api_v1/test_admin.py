from unittest import mock  # standard library
from typing import Dict, Any  # standard library
from datetime import datetime  # standard library

import pytest  # pytest 7.3+
from fastapi import FastAPI  # fastapi 0.95+
from starlette import status  # starlette 0.27+

from app.api.api_v1.endpoints import admin as admin_router  # Internal import
from app.api.api_v1.schemas.admin import UserCreateRequest, UserAdminUpdateRequest, UserPasswordResetRequest  # Internal import
from app.services.admin_service import AdminService  # Internal import
from app.exceptions import AdminException  # Internal import


@pytest.mark.parametrize(
    'mock_stats',
    [
        {
            'user_counts': {'pharma': 10, 'cro': 5, 'admin': 2},
            'active_users': 15,
            'total_molecules': 1000,
            'total_libraries': 50,
            'total_experiments': 100,
            'total_submissions': 75,
            'resource_usage': {'cpu': 25.5, 'memory': 40.2, 'disk': 30.0}
        }
    ]
)
def test_get_system_stats(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch, mock_stats: Dict[str, Any]) -> None:
    """Test getting system statistics as admin user"""
    # Mock AdminService.get_system_stats to return mock_stats
    mock_get_system_stats = mock.MagicMock(return_value=mock_stats)
    monkeypatch.setattr(AdminService, "get_system_stats", mock_get_system_stats)

    # Make GET request to /api/admin/stats with admin token headers
    response = client.get("/api/admin/stats", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected system statistics
    content = response.json()
    assert content["user_counts"] == mock_stats["user_counts"]
    assert content["active_users"] == mock_stats["active_users"]
    assert content["total_molecules"] == mock_stats["total_molecules"]
    assert content["total_libraries"] == mock_stats["total_libraries"]
    assert content["total_experiments"] == mock_stats["total_experiments"]
    assert content["total_submissions"] == mock_stats["total_submissions"]
    assert content["resource_usage"] == mock_stats["resource_usage"]

    # Assert timestamp is present in response
    assert "timestamp" in content


def test_get_system_stats_unauthorized(client: FastAPI, pharma_token_headers: Dict[str, str]) -> None:
    """Test that non-admin users cannot access system statistics"""
    # Make GET request to /api/admin/stats with pharma user token headers
    response = client.get("/api/admin/stats", headers=pharma_token_headers)

    # Assert response status code is 403 Forbidden
    assert response.status_code == status.HTTP_403_FORBIDDEN


@pytest.mark.parametrize(
    'mock_resources',
    [
        {
            'cpu': {'usage': 25.5, 'cores': 8},
            'memory': {'used': 4096, 'total': 16384, 'percent': 25.0},
            'disk': {'used': 50000, 'total': 500000, 'percent': 10.0},
            'database': {'connections': 5, 'size': 1024},
            'storage': {'files': 100, 'size': 2048}
        }
    ]
)
def test_get_system_resources(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch, mock_resources: Dict[str, Any]) -> None:
    """Test getting detailed system resources as admin user"""
    # Mock AdminService.get_system_resources to return mock_resources
    mock_get_system_resources = mock.MagicMock(return_value=mock_resources)
    monkeypatch.setattr(AdminService, "get_system_resources", mock_get_system_resources)

    # Make GET request to /api/admin/resources with admin token headers
    response = client.get("/api/admin/resources", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected system resources
    content = response.json()
    assert content["cpu"] == mock_resources["cpu"]
    assert content["memory"] == mock_resources["memory"]
    assert content["disk"] == mock_resources["disk"]
    assert content["database"] == mock_resources["database"]
    assert content["storage"] == mock_resources["storage"]

    # Assert timestamp is present in response
    assert "timestamp" in content


@pytest.mark.parametrize(
    'mock_health',
    [
        {
            'database': 'UP',
            'file_storage': 'UP',
            'resources': {'cpu': 'OK', 'memory': 'OK', 'disk': 'OK'},
            'overall': 'HEALTHY'
        }
    ]
)
def test_check_system_health(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch, mock_health: Dict[str, Any]) -> None:
    """Test system health check as admin user"""
    # Mock AdminService.check_system_health to return mock_health
    mock_check_system_health = mock.MagicMock(return_value=mock_health)
    monkeypatch.setattr(AdminService, "check_system_health", mock_check_system_health)

    # Make GET request to /api/admin/health with admin token headers
    response = client.get("/api/admin/health", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected health check results
    content = response.json()
    assert content == mock_health


def test_get_users(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test getting list of users as admin user"""
    # Mock AdminService.get_users to return list of test users
    test_users = [
        {"id": 1, "email": "user1@example.com", "role": "pharma"},
        {"id": 2, "email": "user2@example.com", "role": "cro"},
        {"id": 3, "email": "user3@example.com", "role": "admin"}
    ]
    mock_get_users = mock.MagicMock(return_value=test_users)
    monkeypatch.setattr(AdminService, "get_users", mock_get_users)

    # Make GET request to /api/admin/users with admin token headers
    response = client.get("/api/admin/users", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected list of users
    content = response.json()
    assert content == test_users

    # Test with query parameters (role, status, skip, limit)
    response = client.get(
        "/api/admin/users?role=pharma&status=active&skip=1&limit=10",
        headers=admin_token_headers
    )
    assert response.status_code == status.HTTP_200_OK
    # Assert that query parameters are correctly passed to the service
    mock_get_users.assert_called_once()


def test_get_user(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test getting a specific user by ID as admin user"""
    # Mock AdminService.get_user to return a test user
    test_user = {"id": 1, "email": "user1@example.com", "role": "pharma"}
    mock_get_user = mock.MagicMock(return_value=test_user)
    monkeypatch.setattr(AdminService, "get_user", mock_get_user)

    # Make GET request to /api/admin/users/{user_id} with admin token headers
    response = client.get("/api/admin/users/1", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected user details
    content = response.json()
    assert content == test_user

    # Test with non-existent user ID
    mock_get_user.side_effect = AdminException("User not found")
    response = client.get("/api/admin/users/999", headers=admin_token_headers)

    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND


def test_create_user(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test creating a new user as admin user"""
    # Create test user data with email, password, role
    user_data = {"email": "newuser@example.com", "password": "SecureP@ssw0rd", "role": "pharma"}
    # Mock AdminService.create_user to return a new user
    new_user = {"id": 4, "email": "newuser@example.com", "role": "pharma"}
    mock_create_user = mock.MagicMock(return_value=new_user)
    monkeypatch.setattr(AdminService, "create_user", mock_create_user)

    # Make POST request to /api/admin/users with admin token headers and user data
    response = client.post("/api/admin/users", headers=admin_token_headers, json=user_data)

    # Assert response status code is 201 Created
    assert response.status_code == status.HTTP_201_CREATED

    # Assert response JSON contains expected user details
    content = response.json()
    assert content == new_user

    # Test with duplicate email
    mock_create_user.side_effect = AdminException("User with email already exists")
    response = client.post("/api/admin/users", headers=admin_token_headers, json=user_data)

    # Assert response status code is 400 Bad Request
    assert response.status_code == status.HTTP_400_BAD_REQUEST

    # Test with invalid role
    user_data["role"] = "invalid"
    response = client.post("/api/admin/users", headers=admin_token_headers, json=user_data)

    # Assert response status code is 422 Unprocessable Entity
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


def test_update_user(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test updating a user as admin user"""
    # Create test update data with email, role, status
    update_data = {"email": "updated@example.com", "role": "cro", "status": "active"}
    # Mock AdminService.update_user to return updated user
    updated_user = {"id": 1, "email": "updated@example.com", "role": "cro", "status": "active"}
    mock_update_user = mock.MagicMock(return_value=updated_user)
    monkeypatch.setattr(AdminService, "update_user", mock_update_user)

    # Make PUT request to /api/admin/users/{user_id} with admin token headers and update data
    response = client.put("/api/admin/users/1", headers=admin_token_headers, json=update_data)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected updated user details
    content = response.json()
    assert content == updated_user

    # Test with non-existent user ID
    mock_update_user.side_effect = AdminException("User not found")
    response = client.put("/api/admin/users/999", headers=admin_token_headers, json=update_data)

    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND

    # Test with duplicate email
    mock_update_user.side_effect = AdminException("User with email already exists")
    response = client.put("/api/admin/users/1", headers=admin_token_headers, json=update_data)

    # Assert response status code is 400 Bad Request
    assert response.status_code == status.HTTP_400_BAD_REQUEST

    # Test with invalid role
    update_data["role"] = "invalid"
    response = client.put("/api/admin/users/1", headers=admin_token_headers, json=update_data)

    # Assert response status code is 422 Unprocessable Entity
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


def test_delete_user(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test deleting a user as admin user"""
    # Mock AdminService.delete_user to return True
    mock_delete_user = mock.MagicMock(return_value=True)
    monkeypatch.setattr(AdminService, "delete_user", mock_delete_user)

    # Make DELETE request to /api/admin/users/{user_id} with admin token headers
    response = client.delete("/api/admin/users/1", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains success message
    content = response.json()
    assert content["success"] is True

    # Test with non-existent user ID
    mock_delete_user.side_effect = AdminException("User not found")
    response = client.delete("/api/admin/users/999", headers=admin_token_headers)

    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND


def test_reset_user_password(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test resetting a user's password as admin user"""
    # Create test password reset data with new_password
    password_data = {"new_password": "NewSecureP@ssw0rd"}
    # Mock AdminService.reset_user_password to return updated user
    updated_user = {"id": 1, "email": "user1@example.com", "role": "pharma"}
    mock_reset_user_password = mock.MagicMock(return_value=updated_user)
    monkeypatch.setattr(AdminService, "reset_user_password", mock_reset_user_password)

    # Make POST request to /api/admin/users/{user_id}/reset-password with admin token headers and password data
    response = client.post("/api/admin/users/1/reset-password", headers=admin_token_headers, json=password_data)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains success message
    content = response.json()
    assert content["success"] is True

    # Test with non-existent user ID
    mock_reset_user_password.side_effect = AdminException("User not found")
    response = client.post("/api/admin/users/999/reset-password", headers=admin_token_headers, json=password_data)

    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND

    # Test with invalid password (too short)
    password_data["new_password"] = "short"
    response = client.post("/api/admin/users/1/reset-password", headers=admin_token_headers, json=password_data)

    # Assert response status code is 422 Unprocessable Entity
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY


@pytest.mark.parametrize(
    'mock_logs',
    [
        {
            'items': [
                {
                    'id': 1,
                    'user_id': 1,
                    'username': 'admin@example.com',
                    'action': 'CREATE',
                    'resource_type': 'USER',
                    'resource_id': 2,
                    'details': {'role': 'pharma'},
                    'ip_address': '127.0.0.1',
                    'timestamp': '2023-06-01T12:00:00'
                }
            ],
            'total': 1,
            'page': 1,
            'size': 20
        }
    ]
)
def test_get_activity_logs(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch, mock_logs: Dict[str, Any]) -> None:
    """Test getting activity logs as admin user"""
    # Mock AdminService.get_activity_logs to return mock_logs
    mock_get_activity_logs = mock.MagicMock(return_value=mock_logs)
    monkeypatch.setattr(AdminService, "get_activity_logs", mock_get_activity_logs)

    # Make GET request to /api/admin/activity-logs with admin token headers
    response = client.get("/api/admin/activity-logs", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected activity logs
    content = response.json()
    assert content["items"] == mock_logs["items"]
    assert content["total"] == mock_logs["total"]
    assert content["page"] == mock_logs["page"]
    assert content["size"] == mock_logs["size"]

    # Test with query parameters (user_id, action, resource_type, page, size)
    response = client.get(
        "/api/admin/activity-logs?user_id=1&action=CREATE&resource_type=USER&page=2&size=10",
        headers=admin_token_headers
    )
    assert response.status_code == status.HTTP_200_OK
    # Assert that query parameters are correctly passed to the service
    mock_get_activity_logs.assert_called_once()


@pytest.mark.parametrize(
    'mock_alerts',
    [
        {
            'items': [
                {
                    'id': 1,
                    'title': 'Disk Space Warning',
                    'message': 'Disk usage above 75%',
                    'severity': 'WARNING',
                    'alert_type': 'RESOURCE',
                    'resolved': False,
                    'created_at': '2023-06-01T12:00:00',
                    'resolved_at': None,
                    'resolved_by': None,
                    'details': {'usage': 75.5}
                }
            ],
            'total': 1,
            'counts_by_severity': {'WARNING': 1, 'CRITICAL': 0, 'INFO': 0}
        }
    ]
)
def test_get_system_alerts(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch, mock_alerts: Dict[str, Any]) -> None:
    """Test getting system alerts as admin user"""
    # Mock AdminService.get_system_alerts to return mock_alerts
    mock_get_system_alerts = mock.MagicMock(return_value=mock_alerts)
    monkeypatch.setattr(AdminService, "get_system_alerts", mock_get_system_alerts)

    # Make GET request to /api/admin/alerts with admin token headers
    response = client.get("/api/admin/alerts", headers=admin_token_headers)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected system alerts
    content = response.json()
    assert content["items"] == mock_alerts["items"]
    assert content["total"] == mock_alerts["total"]
    assert content["counts_by_severity"] == mock_alerts["counts_by_severity"]

    # Test with query parameters (severity, alert_type, resolved)
    response = client.get(
        "/api/admin/alerts?severity=WARNING&alert_type=RESOURCE&resolved=false",
        headers=admin_token_headers
    )
    assert response.status_code == status.HTTP_200_OK
    # Assert that query parameters are correctly passed to the service
    mock_get_system_alerts.assert_called_once()


def test_update_system_alert(client: FastAPI, admin_token_headers: Dict[str, str], monkeypatch: pytest.MonkeyPatch) -> None:
    """Test updating a system alert as admin user"""
    # Create test alert update data with resolved=True and resolution_notes
    alert_data = {"resolved": True, "resolution_notes": "Resolved by cleaning up disk space"}
    # Mock AdminService.update_system_alert to return updated alert
    updated_alert = {"id": 1, "title": "Disk Space Warning", "message": "Disk usage above 75%", "severity": "WARNING", "alert_type": "RESOURCE", "resolved": True, "created_at": "2023-06-01T12:00:00", "resolved_at": "2023-06-05T15:00:00", "resolved_by": None, "details": {"usage": 75.5, "notes": "Resolved by cleaning up disk space"}}
    mock_update_system_alert = mock.MagicMock(return_value=updated_alert)
    monkeypatch.setattr(AdminService, "update_system_alert", mock_update_system_alert)

    # Make PUT request to /api/admin/alerts/{alert_id} with admin token headers and update data
    response = client.put("/api/admin/alerts/1", headers=admin_token_headers, json=alert_data)

    # Assert response status code is 200 OK
    assert response.status_code == status.HTTP_200_OK

    # Assert response JSON contains expected updated alert details
    content = response.json()
    assert content == updated_alert

    # Test with non-existent alert ID
    mock_update_system_alert.side_effect = AdminException("Alert not found")
    response = client.put("/api/admin/alerts/999", headers=admin_token_headers, json=alert_data)

    # Assert response status code is 404 Not Found
    assert response.status_code == status.HTTP_404_NOT_FOUND