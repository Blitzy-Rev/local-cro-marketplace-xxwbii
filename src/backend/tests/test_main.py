import pytest  # -- 7.3+ Testing framework for Python
from fastapi.testclient import TestClient  # -- 0.95+ Test client for FastAPI applications
from fastapi import status  # -- 0.95+ HTTP status codes for response validation

from app import app  # src/backend/app/__init__.py Import the main FastAPI application instance for testing

# Define a test client to be used in the tests
client = TestClient(app)

def test_read_main_health():
    """Tests the main health endpoint to verify the application is running"""
    # Send GET request to /health/live endpoint
    response = client.get("/health/live")
    # Verify response status code is 200
    assert response.status_code == 200
    # Verify response JSON contains status: UP
    assert response.json() == {"status": "UP"}

def test_read_readiness():
    """Tests the readiness health endpoint to verify dependencies are available"""
    # Send GET request to /health/ready endpoint
    response = client.get("/health/ready")
    # Verify response status code is 200
    assert response.status_code == 200
    # Verify response JSON contains status field
    assert "status" in response.json()
    # Verify response JSON contains details with database and redis status
    assert "details" in response.json()
    assert "database" in response.json()["details"]
    assert "redis" in response.json()["details"]

def test_api_docs_available():
    """Tests that the API documentation endpoints are available"""
    # Send GET request to /api/docs endpoint
    response = client.get("/api/docs")
    # Verify response status code is 200
    assert response.status_code == 200

    # Send GET request to /api/redoc endpoint
    response = client.get("/api/redoc")
    # Verify response status code is 200
    assert response.status_code == 200

    # Send GET request to /api/openapi.json endpoint
    response = client.get("/api/openapi.json")
    # Verify response status code is 200
    assert response.status_code == 200

def test_cors_headers():
    """Tests that CORS headers are properly set in responses"""
    # Send OPTIONS request to /health/live endpoint with Origin header
    response = client.options("/health/live", headers={"Origin": "http://example.com"})
    # Verify response status code is 200
    assert response.status_code == 200
    # Verify response contains Access-Control-Allow-Origin header
    assert "access-control-allow-origin" in response.headers
    # Verify response contains Access-Control-Allow-Methods header
    assert "access-control-allow-methods" in response.headers
    # Verify response contains Access-Control-Allow-Headers header
    assert "access-control-allow-headers" in response.headers

def test_validation_error_handler():
    """Tests the custom validation error handler for invalid requests"""
    # Send POST request to an endpoint with invalid data
    response = client.post("/api/auth/login", json={"email": "invalid-email", "password": ""})
    # Verify response status code is 422
    assert response.status_code == status.HTTP_422_UNPROCESSABLE_ENTITY
    # Verify response JSON contains detail field with validation errors
    assert "detail" in response.json()

def test_http_error_handler():
    """Tests the custom HTTP error handler for HTTP exceptions"""
    # Send GET request to a non-existent endpoint
    response = client.get("/api/non-existent-endpoint")
    # Verify response status code is 404
    assert response.status_code == status.HTTP_404_NOT_FOUND
    # Verify response JSON contains detail field with error message
    assert "detail" in response.json()

def test_api_prefix_configuration():
    """Tests that the API prefix is correctly configured"""
    # Send GET request to /api/v1/health/live endpoint
    response = client.get("/api/v1/health/live")
    # Verify response status code is 200
    assert response.status_code == 200
    # Verify response JSON contains status: UP
    assert response.json() == {"status": "UP"}