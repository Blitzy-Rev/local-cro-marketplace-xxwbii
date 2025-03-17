"""
Test module for the health check endpoints of the Molecular Data Management and CRO Integration Platform API.

This file contains tests to verify the functionality of various health check endpoints including liveness,
readiness, deep health checks, database connectivity, and system metrics.
"""

import pytest
from unittest.mock import MagicMock
from fastapi import status

from ..conftest import client, db_session, mock_redis, mock_minio


def test_health_live_endpoint(client):
    """Tests the /health/live endpoint to ensure it returns a 200 OK status and UP status"""
    response = client.get("/api/v1/health/live")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "status" in response_data
    assert response_data["status"] == "UP"


def test_health_ready_endpoint(client):
    """Tests the /health/ready endpoint to ensure it returns a 200 OK status and component health information"""
    response = client.get("/api/v1/health/ready")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "status" in response_data
    
    assert "details" in response_data
    details = response_data["details"]
    
    assert "database" in details
    assert "status" in details["database"]
    assert "latency_ms" in details["database"]
    
    assert "redis" in details
    assert "status" in details["redis"]
    assert "latency_ms" in details["redis"]


def test_health_deep_endpoint(client):
    """Tests the /health/deep endpoint to ensure it returns a 200 OK status and comprehensive component health information"""
    response = client.get("/api/v1/health/deep")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "status" in response_data
    
    assert "components" in response_data
    components = response_data["components"]
    
    assert "database" in components
    assert "redis" in components
    assert "minio" in components
    assert "celery" in components
    
    for component_name, component_data in components.items():
        assert "status" in component_data
        assert component_data["status"] in ["UP", "DOWN"]


def test_health_db_endpoint(client):
    """Tests the /health/db endpoint to ensure it returns database connectivity status"""
    response = client.get("/api/v1/health/db")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "status" in response_data
    assert "latency_ms" in response_data
    assert "details" in response_data


@pytest.mark.parametrize("psutil_available", [True, False])
def test_health_system_endpoint(client, mocker, psutil_available):
    """Tests the /health/system endpoint to ensure it returns system metrics"""
    # Mock psutil availability based on test parameter
    if psutil_available:
        # Create mock with expected return values
        mock_psutil = MagicMock()
        mock_psutil.virtual_memory.return_value.percent = 45.5
        mock_psutil.cpu_percent.return_value = 32.1
        mock_psutil.disk_usage.return_value.percent = 67.8
        
        # Patch the module in the health check module
        mocker.patch("app.api.api_v1.endpoints.health.psutil", mock_psutil)
    else:
        # Mock that psutil is not available by raising ImportError when imported
        mocker.patch("app.api.api_v1.endpoints.health.psutil", None)
    
    response = client.get("/api/v1/health/system")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "status" in response_data
    assert "metrics" in response_data
    
    if psutil_available:
        metrics = response_data["metrics"]
        assert "memory" in metrics
        assert "cpu" in metrics
        assert "disk" in metrics
    else:
        assert "message" in response_data["metrics"]
        assert "psutil not available" in response_data["metrics"]["message"].lower()


def test_database_health_check_failure(client, mocker):
    """Tests the database health check when the database connection fails"""
    # Mock the db connection to fail
    mocker.patch("app.api.api_v1.endpoints.health.check_database_connection", 
                 side_effect=Exception("Database connection failed"))
    
    response = client.get("/api/v1/health/db")
    
    # Even with failure, status code should be 200 OK, with status DOWN in the body
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert response_data["status"] == "DOWN"
    assert "error" in response_data
    assert "database connection failed" in response_data["error"].lower()


def test_redis_health_check_failure(client, mocker):
    """Tests the Redis health check when the Redis connection fails"""
    # Mock the Redis health check to fail
    mocker.patch("app.api.api_v1.endpoints.health.check_redis_connection", 
                 side_effect=Exception("Redis connection failed"))
    
    response = client.get("/api/v1/health/ready")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "details" in response_data
    assert "redis" in response_data["details"]
    assert response_data["details"]["redis"]["status"] == "DOWN"
    assert "error" in response_data["details"]["redis"]
    assert "redis connection failed" in response_data["details"]["redis"]["error"].lower()


def test_minio_health_check_failure(client, mocker):
    """Tests the MinIO health check when the MinIO connection fails"""
    # Mock the MinIO health check to fail
    mocker.patch("app.api.api_v1.endpoints.health.check_minio_connection", 
                 side_effect=Exception("MinIO connection failed"))
    
    response = client.get("/api/v1/health/deep")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "components" in response_data
    assert "minio" in response_data["components"]
    assert response_data["components"]["minio"]["status"] == "DOWN"
    assert "error" in response_data["components"]["minio"]
    assert "minio connection failed" in response_data["components"]["minio"]["error"].lower()


def test_celery_health_check_failure(client, mocker):
    """Tests the Celery health check when the Celery connection fails"""
    # Mock the Celery health check to fail
    mocker.patch("app.api.api_v1.endpoints.health.check_celery_connection", 
                 side_effect=Exception("Celery connection failed"))
    
    response = client.get("/api/v1/health/deep")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert "components" in response_data
    assert "celery" in response_data["components"]
    assert response_data["components"]["celery"]["status"] == "DOWN"
    assert "error" in response_data["components"]["celery"]
    assert "celery connection failed" in response_data["components"]["celery"]["error"].lower()


@pytest.mark.parametrize(
    "component_statuses,expected_overall", 
    [
        ({"database": "UP", "redis": "UP"}, "UP"),
        ({"database": "DOWN", "redis": "UP"}, "DOWN"),
        ({"database": "UP", "redis": "DOWN"}, "DOWN")
    ]
)
def test_overall_health_status(client, mocker, component_statuses, expected_overall):
    """Tests that the overall health status is correctly determined based on component statuses"""
    # Mock the component health checks to return specified statuses
    for component, status_value in component_statuses.items():
        if component == "database":
            if status_value == "UP":
                mocker.patch("app.api.api_v1.endpoints.health.check_database_connection", 
                             return_value={"status": "UP", "latency_ms": 1.5})
            else:
                mocker.patch("app.api.api_v1.endpoints.health.check_database_connection", 
                             side_effect=Exception("Database connection failed"))
        elif component == "redis":
            if status_value == "UP":
                mocker.patch("app.api.api_v1.endpoints.health.check_redis_connection", 
                             return_value={"status": "UP", "latency_ms": 0.8})
            else:
                mocker.patch("app.api.api_v1.endpoints.health.check_redis_connection", 
                             side_effect=Exception("Redis connection failed"))
    
    response = client.get("/api/v1/health/ready")
    
    assert response.status_code == status.HTTP_200_OK
    
    response_data = response.json()
    assert response_data["status"] == expected_overall