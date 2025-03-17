"""
Test module for the CSV upload, mapping, and processing API endpoints in the Molecular Data Management and CRO Integration Platform. Validates the functionality of CSV file handling, column mapping, processing status tracking, and error handling.
"""

import pytest
import os
import io
import json
from typing import Dict, List, Any

from ....app.api.api_v1.schemas.csv import (
    CSVUploadResponse,
    CSVMappingResponse,
    CSVStatusResponse,
    AvailablePropertiesResponse,
    CSV_PROCESSING_STATUS,
)
from ....app.exceptions import ValidationException

# Test data
TEST_CSV_CONTENT = """SMILES,MW,LogP,Activity
CCO,46.07,-0.14,78.5
CCCCO,74.12,0.88,45.2
c1ccccc1,78.11,1.90,92.3"""

TEST_CSV_HEADERS = ["SMILES", "MW", "LogP", "Activity"]

TEST_CSV_MAPPING = {
    "SMILES": "smiles",
    "MW": "molecular_weight",
    "LogP": "logp",
    "Activity": "activity"
}

API_PREFIX = "/api/v1/csv"


def create_test_csv_file():
    """Creates a test CSV file in memory for testing."""
    file = io.BytesIO(TEST_CSV_CONTENT.encode())
    file.seek(0)
    return file


@pytest.mark.parametrize('file_id', ['test_file_123'])
def test_upload_csv_success(client, pharma_token_headers, mocker, file_id):
    """Tests successful CSV file upload."""
    test_file = create_test_csv_file()
    
    # Mock file storage service to save the CSV file
    mocker.patch(
        "app.services.file_storage_service.save_csv_file",
        return_value=file_id
    )
    
    # Mock CSV service to validate the file
    mocker.patch(
        "app.services.csv_service.validate_csv_file",
        return_value=(TEST_CSV_HEADERS, 3)
    )
    
    # Mock CSV service to get headers
    mocker.patch(
        "app.services.csv_service.get_csv_headers",
        return_value=TEST_CSV_HEADERS
    )
    
    response = client.post(
        f"{API_PREFIX}/upload",
        files={"file": ("test.csv", test_file, "text/csv")},
        headers=pharma_token_headers
    )
    
    assert response.status_code == 200
    data = response.json()
    assert data["file_id"] == file_id
    assert data["headers"] == TEST_CSV_HEADERS
    assert data["row_count"] == 3
    assert "message" in data


def test_upload_csv_invalid_file(client, pharma_token_headers):
    """Tests CSV upload with invalid file."""
    empty_file = io.BytesIO(b"")
    
    response = client.post(
        f"{API_PREFIX}/upload",
        files={"file": ("test.csv", empty_file, "text/csv")},
        headers=pharma_token_headers
    )
    
    assert response.status_code == 400
    assert "message" in response.json()


def test_upload_csv_unauthorized(client):
    """Tests CSV upload without authentication."""
    test_file = create_test_csv_file()
    
    response = client.post(
        f"{API_PREFIX}/upload",
        files={"file": ("test.csv", test_file, "text/csv")}
    )
    
    assert response.status_code == 401


@pytest.mark.parametrize('job_id', ['job_123'])
def test_map_csv_columns_success(client, pharma_token_headers, mocker, test_pharma_user, job_id):
    """Tests successful CSV column mapping."""
    # Mock CSV service to get headers
    mocker.patch(
        "app.services.csv_service.get_csv_headers",
        return_value=TEST_CSV_HEADERS
    )
    
    # Mock validation of column mapping
    mocker.patch(
        "app.services.csv_service.validate_column_mapping",
        return_value=True
    )
    
    # Mock creation of processing job
    mocker.patch(
        "app.services.csv_service.create_processing_job",
        return_value=job_id
    )
    
    request_data = {
        "file_id": "test_file_123",
        "mapping": TEST_CSV_MAPPING
    }
    
    response = client.post(
        f"{API_PREFIX}/map",
        json=request_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == 202
    data = response.json()
    assert data["job_id"] == job_id
    assert data["status"] == CSV_PROCESSING_STATUS.PENDING.name
    assert "message" in data


def test_map_csv_columns_invalid_mapping(client, pharma_token_headers, mocker):
    """Tests CSV mapping with invalid mapping configuration."""
    # Mock CSV service to get headers
    mocker.patch(
        "app.services.csv_service.get_csv_headers",
        return_value=TEST_CSV_HEADERS
    )
    
    # Mock validation to raise exception
    mocker.patch(
        "app.services.csv_service.validate_column_mapping",
        side_effect=ValidationException("Invalid mapping: SMILES column required")
    )
    
    request_data = {
        "file_id": "test_file_123",
        "mapping": {
            "MW": "molecular_weight",
            "LogP": "logp",
            "Activity": "activity"
        }
    }
    
    response = client.post(
        f"{API_PREFIX}/map",
        json=request_data,
        headers=pharma_token_headers
    )
    
    assert response.status_code == 422
    assert "message" in response.json()


@pytest.mark.parametrize('status', [
    CSV_PROCESSING_STATUS.PENDING,
    CSV_PROCESSING_STATUS.PROCESSING,
    CSV_PROCESSING_STATUS.COMPLETED
])
def test_get_processing_status(client, pharma_token_headers, mocker, test_pharma_user, status):
    """Tests retrieving CSV processing status."""
    job_id = "job_123"
    
    # Create mock job status based on test parameter
    mock_status = {
        "job_id": job_id,
        "user_id": test_pharma_user.id,
        "status": status.name,
        "progress": 50 if status == CSV_PROCESSING_STATUS.PROCESSING else None,
        "summary": {
            "total_rows": 3,
            "processed_rows": 3,
            "successful_rows": 3,
            "failed_rows": 0,
            "errors": []
        } if status == CSV_PROCESSING_STATUS.COMPLETED else None
    }
    
    # Mock job status retrieval
    mocker.patch(
        "app.services.csv_service.get_job_status",
        return_value=mock_status
    )
    
    response = client.get(
        f"{API_PREFIX}/status/{job_id}",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 200
    data = response.json()
    assert data["job_id"] == job_id
    assert data["status"] == status.name
    
    if status == CSV_PROCESSING_STATUS.PROCESSING:
        assert data["progress"] == 50
    elif status == CSV_PROCESSING_STATUS.COMPLETED:
        assert "summary" in data
        assert data["summary"]["total_rows"] == 3
        assert data["summary"]["successful_rows"] == 3


def test_get_processing_status_not_found(client, pharma_token_headers, mocker):
    """Tests retrieving status for non-existent job."""
    job_id = "non_existent_job"
    
    # Mock job status to raise KeyError
    mocker.patch(
        "app.services.csv_service.get_job_status",
        side_effect=KeyError("Job not found")
    )
    
    response = client.get(
        f"{API_PREFIX}/status/{job_id}",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 404
    assert "message" in response.json()


def test_cancel_processing(client, pharma_token_headers, mocker, test_pharma_user):
    """Tests cancelling a CSV processing job."""
    job_id = "job_123"
    
    # Mock job cancellation
    mocker.patch(
        "app.services.csv_service.cancel_job",
        return_value=True
    )
    
    response = client.post(
        f"{API_PREFIX}/cancel/{job_id}",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 200
    assert "msg" in response.json()


def test_cancel_processing_not_found(client, pharma_token_headers, mocker):
    """Tests cancelling a non-existent processing job."""
    job_id = "non_existent_job"
    
    # Mock job cancellation to raise KeyError
    mocker.patch(
        "app.services.csv_service.cancel_job",
        side_effect=KeyError("Job not found")
    )
    
    response = client.post(
        f"{API_PREFIX}/cancel/{job_id}",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 404
    assert "message" in response.json()


def test_get_available_properties(client, pharma_token_headers, mocker):
    """Tests retrieving available system properties for mapping."""
    mock_properties = [
        {
            "name": "smiles",
            "description": "SMILES representation of molecule",
            "data_type": "string",
            "required": True
        },
        {
            "name": "molecular_weight",
            "description": "Molecular weight in g/mol",
            "data_type": "number",
            "required": False
        },
        {
            "name": "logp",
            "description": "Partition coefficient",
            "data_type": "number",
            "required": False
        }
    ]
    
    # Mock available properties retrieval
    mocker.patch(
        "app.services.csv_service.get_available_properties",
        return_value=mock_properties
    )
    
    response = client.get(
        f"{API_PREFIX}/properties",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 200
    data = response.json()
    assert "properties" in data
    assert len(data["properties"]) == 3
    assert data["properties"][0]["name"] == "smiles"
    assert data["properties"][1]["data_type"] == "number"


def test_generate_import_report(client, pharma_token_headers, mocker):
    """Tests generating a CSV import report."""
    job_id = "job_123"
    
    response = client.post(
        f"{API_PREFIX}/report/{job_id}",
        headers=pharma_token_headers
    )
    
    assert response.status_code == 202
    assert "msg" in response.json()


def test_end_to_end_csv_workflow(client, pharma_token_headers, mocker, test_pharma_user):
    """Tests the complete CSV workflow from upload to processing."""
    file_id = "test_file_123"
    job_id = "job_123"
    
    test_file = create_test_csv_file()
    
    # Mock file storage service
    mocker.patch(
        "app.services.file_storage_service.save_csv_file",
        return_value=file_id
    )
    
    # Mock CSV validation and header retrieval
    mocker.patch(
        "app.services.csv_service.validate_csv_file",
        return_value=(TEST_CSV_HEADERS, 3)
    )
    
    mocker.patch(
        "app.services.csv_service.get_csv_headers",
        return_value=TEST_CSV_HEADERS
    )
    
    # Mock column mapping validation
    mocker.patch(
        "app.services.csv_service.validate_column_mapping",
        return_value=True
    )
    
    # Mock job creation
    mocker.patch(
        "app.services.csv_service.create_processing_job",
        return_value=job_id
    )
    
    # Mock job status progression
    job_status_sequence = [
        {
            "job_id": job_id,
            "user_id": test_pharma_user.id,
            "status": CSV_PROCESSING_STATUS.PROCESSING.name,
            "progress": 50
        },
        {
            "job_id": job_id,
            "user_id": test_pharma_user.id,
            "status": CSV_PROCESSING_STATUS.COMPLETED.name,
            "progress": 100,
            "summary": {
                "total_rows": 3,
                "processed_rows": 3,
                "successful_rows": 3,
                "failed_rows": 0,
                "errors": []
            }
        }
    ]
    
    mock_get_job_status = mocker.patch("app.services.csv_service.get_job_status")
    mock_get_job_status.side_effect = job_status_sequence
    
    # Step 1: Upload CSV file
    upload_response = client.post(
        f"{API_PREFIX}/upload",
        files={"file": ("test.csv", test_file, "text/csv")},
        headers=pharma_token_headers
    )
    
    assert upload_response.status_code == 200
    upload_data = upload_response.json()
    assert upload_data["file_id"] == file_id
    
    # Step 2: Map CSV columns
    mapping_request = {
        "file_id": file_id,
        "mapping": TEST_CSV_MAPPING
    }
    
    mapping_response = client.post(
        f"{API_PREFIX}/map",
        json=mapping_request,
        headers=pharma_token_headers
    )
    
    assert mapping_response.status_code == 202
    mapping_data = mapping_response.json()
    assert mapping_data["job_id"] == job_id
    
    # Step 3: Check processing status
    status_response = client.get(
        f"{API_PREFIX}/status/{job_id}",
        headers=pharma_token_headers
    )
    
    assert status_response.status_code == 200
    status_data = status_response.json()
    assert status_data["status"] == CSV_PROCESSING_STATUS.COMPLETED.name
    assert "summary" in status_data
    assert status_data["summary"]["successful_rows"] == 3