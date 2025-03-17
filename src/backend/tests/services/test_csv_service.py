# src/backend/tests/services/test_csv_service.py
import os
import io
import tempfile
from typing import List, Dict, Any
import pandas as pd
import pytest
from unittest.mock import MagicMock

from src.backend.app.services.csv_service import CSVService, CSVProcessingException, ValidationException
from src.backend.app.api.api_v1.schemas import csv as csv_schemas
from src.backend.app.constants import DEFAULT_MOLECULE_PROPERTIES

# Mock the FileStorageService and MoleculeService for testing purposes
@pytest.fixture
def mock_file_storage_service():
    """
    Fixture to mock the FileStorageService.
    """
    return MagicMock()

@pytest.fixture
def mock_molecule_service():
    """
    Fixture to mock the MoleculeService.
    """
    return MagicMock()

@pytest.fixture
def csv_service(mock_file_storage_service, mock_molecule_service):
    """
    Fixture to create a CSVService instance with mocked dependencies.
    """
    return CSVService(mock_file_storage_service, mock_molecule_service)

def create_test_csv(filename: str, data: List[Dict[str, Any]]) -> str:
    """
    Creates a test CSV file with sample molecule data.

    Args:
        filename: The name of the CSV file to create.
        data: A list of dictionaries representing the data to write to the CSV file.

    Returns:
        The path to the created CSV file.
    """
    df = pd.DataFrame(data)
    temp_dir = tempfile.mkdtemp()
    file_path = os.path.join(temp_dir, filename)
    df.to_csv(file_path, index=False)
    return file_path

def test_validate_csv_file(csv_service, mock_file_storage_service):
    """
    Test CSV file validation functionality.
    """
    # Create a test CSV file
    test_data = [{'smiles': 'CCO', 'MW': 46.07, 'LogP': -0.14}]
    file_path = create_test_csv('test.csv', test_data)

    # Mock the file storage service to return the test CSV file
    with open(file_path, 'rb') as f:
        mock_file_storage_service.get_csv_file.return_value = (io.BytesIO(f.read()), 'test.csv', 'text/csv')

    # Call the validate_csv_file method
    validation_results = csv_service.validate_csv_file('test_file_id')

    # Assert that the validation results are correct
    assert validation_results['is_valid'] is True
    assert validation_results['headers'] == ['smiles', 'MW', 'LogP']
    assert validation_results['row_count'] == 1

    # Clean up the test CSV file
    os.remove(file_path)

def test_get_csv_headers(csv_service, mock_file_storage_service):
    """
    Test extraction of headers from CSV files.
    """
    # Create a test CSV file
    test_data = [{'smiles': 'CCO', 'MW': 46.07, 'LogP': -0.14}]
    file_path = create_test_csv('test.csv', test_data)

    # Mock the file storage service to return the test CSV file
    with open(file_path, 'rb') as f:
        mock_file_storage_service.get_csv_file.return_value = (io.BytesIO(f.read()), 'test.csv', 'text/csv')

    # Call the get_csv_headers method
    headers = csv_service.get_csv_headers('test_file_id')

    # Assert that the headers are correct
    assert headers == ['smiles', 'MW', 'LogP']

    # Clean up the test CSV file
    os.remove(file_path)

def test_validate_column_mapping(csv_service):
    """
    Test validation of column mapping configuration.
    """
    # Define test headers and mapping
    headers = ['smiles', 'MW', 'LogP']
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}

    # Call the validate_column_mapping method
    is_valid = csv_service.validate_column_mapping(headers, mapping)

    # Assert that the mapping is valid
    assert is_valid is True

    # Test with invalid mapping (missing SMILES)
    mapping = {'MW': 'MW', 'LogP': 'LogP'}
    with pytest.raises(ValidationException) as exc_info:
        csv_service.validate_column_mapping(headers, mapping)
    assert "SMILES column must be mapped" in str(exc_info.value)

    # Test with invalid mapping (column not in headers)
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'Invalid': 'Invalid'}
    with pytest.raises(ValidationException) as exc_info:
        csv_service.validate_column_mapping(headers, mapping)
    assert "Mapped column 'Invalid' does not exist in CSV headers" in str(exc_info.value)

def test_get_available_properties(csv_service, mock_molecule_service):
    """
    Test retrieval of available system properties.
    """
    # Mock the molecule service to return a list of available properties
    mock_molecule_service.get_available_properties.return_value = ['MW', 'LogP', 'HBA', 'HBD']

    # Call the get_available_properties method
    properties = csv_service.get_available_properties()

    # Assert that the properties are correct
    assert properties == [{'name': 'MW', 'description': 'MW', 'data_type': 'string'},
                          {'name': 'LogP', 'description': 'LogP', 'data_type': 'string'},
                          {'name': 'HBA', 'description': 'HBA', 'data_type': 'string'},
                          {'name': 'HBD', 'description': 'HBD', 'data_type': 'string'}]

def test_create_processing_job(csv_service):
    """
    Test creation of CSV processing jobs.
    """
    # Define test parameters
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1

    # Call the create_processing_job method
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Assert that the job ID is a UUID
    assert isinstance(job_id, str)
    # Assert that the job status is PENDING
    assert csv_service.JOB_STATUS_CACHE[job_id]['status'] == csv_schemas.CSV_PROCESSING_STATUS.PENDING.name
    assert csv_service.JOB_STATUS_CACHE[job_id]['file_id'] == file_id
    assert csv_service.JOB_STATUS_CACHE[job_id]['mapping'] == mapping
    assert csv_service.JOB_STATUS_CACHE[job_id]['user_id'] == user_id

def test_update_job_status(csv_service):
    """
    Test updating of job status.
    """
    # Create a test processing job
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the update_job_status method
    status = csv_schemas.CSV_PROCESSING_STATUS.PROCESSING.name
    is_updated = csv_service.update_job_status(job_id, status, progress=50)

    # Assert that the job status is updated
    assert is_updated is True
    assert csv_service.JOB_STATUS_CACHE[job_id]['status'] == status
    assert csv_service.JOB_STATUS_CACHE[job_id]['progress'] == 50

def test_get_job_status(csv_service):
    """
    Test retrieval of job status.
    """
    # Create a test processing job
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the get_job_status method
    job_status = csv_service.get_job_status(job_id)

    # Assert that the job status is correct
    assert job_status['status'] == csv_schemas.CSV_PROCESSING_STATUS.PENDING.name
    assert job_status['file_id'] == file_id
    assert job_status['mapping'] == mapping
    assert job_status['user_id'] == user_id

    # Test with invalid job ID
    with pytest.raises(Exception) as exc_info:
        csv_service.get_job_status('invalid_job_id')
    assert "Job with ID invalid_job_id not found" in str(exc_info.value)

def test_cancel_job(csv_service):
    """
    Test cancellation of processing jobs.
    """
    # Create a test processing job
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the cancel_job method
    is_cancelled = csv_service.cancel_job(job_id, user_id)

    # Assert that the job is cancelled
    assert is_cancelled is True
    assert csv_service.JOB_STATUS_CACHE[job_id]['status'] == csv_schemas.CSV_PROCESSING_STATUS.CANCELLED.name

    # Test with invalid job ID
    is_cancelled = csv_service.cancel_job('invalid_job_id', user_id)
    assert is_cancelled is False

def test_process_csv_file(csv_service, mock_file_storage_service, mock_molecule_service):
    """
    Test end-to-end CSV processing.
    """
    # Create a test CSV file
    test_data = [{'smiles': 'CCO', 'MW': 46.07, 'LogP': -0.14},
                 {'smiles': 'c1ccccc1', 'MW': 78.11, 'LogP': 1.90}]
    file_path = create_test_csv('test.csv', test_data)

    # Mock the file storage service to return the test CSV file
    with open(file_path, 'rb') as f:
        mock_file_storage_service.get_csv_file.return_value = (io.BytesIO(f.read()), 'test.csv', 'text/csv')

    # Mock the molecule service to return a successful bulk creation result
    mock_molecule_service.bulk_create.return_value = {'created_molecules': [], 'total_created': 2, 'total_duplicates': 0, 'total_processed': 2, 'total_failed': 0}

    # Define test parameters
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the process_csv_file method
    processing_results = csv_service.process_csv_file(job_id, file_id, mapping, user_id)

    # Assert that the processing results are correct
    assert processing_results['status'] == csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name
    assert processing_results['summary']['total_rows'] == 2
    assert processing_results['summary']['processed_rows'] == 2
    assert processing_results['summary']['successful_rows'] == 2
    assert processing_results['summary']['failed_rows'] == 0

    # Clean up the test CSV file
    os.remove(file_path)

def test_cleanup_job_resources(csv_service, mock_file_storage_service):
    """
    Test cleanup of job resources.
    """
    # Create a test processing job
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Mock the file storage service to return a successful deletion
    mock_file_storage_service.delete_file.return_value = True

    # Call the cleanup_job_resources method
    is_cleaned = csv_service.cleanup_job_resources(job_id)

    # Assert that the resources are cleaned up
    assert is_cleaned is True
    mock_file_storage_service.delete_file.assert_called_once_with(file_id, 'csv-uploads')
    assert job_id not in csv_service.JOB_STATUS_CACHE

def test_get_job_statistics(csv_service):
    """
    Test retrieval of job statistics.
    """
    # Create test processing jobs
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id1 = csv_service.create_processing_job(file_id, mapping, user_id)
    job_id2 = csv_service.create_processing_job(file_id, mapping, user_id)
    csv_service.update_job_status(job_id1, csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name)
    csv_service.update_job_status(job_id2, csv_schemas.CSV_PROCESSING_STATUS.FAILED.name)

    # Call the get_job_statistics method
    statistics = csv_service.get_job_statistics()

    # Assert that the statistics are correct
    assert statistics['total_jobs'] == 2
    assert statistics['total_pending'] == 0
    assert statistics['total_processing'] == 0
    assert statistics['total_completed'] == 1
    assert statistics['total_failed'] == 1
    assert statistics['total_cancelled'] == 0

def test_error_handling(csv_service, mock_file_storage_service):
    """
    Test error handling in CSV processing.
    """
    # Mock the file storage service to raise an exception
    mock_file_storage_service.get_csv_file.side_effect = Exception('File not found')

    # Define test parameters
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles', 'MW': 'MW', 'LogP': 'LogP'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the process_csv_file method and assert that it raises an exception
    with pytest.raises(Exception) as exc_info:
        csv_service.process_csv_file(job_id, file_id, mapping, user_id)
    assert "File not found" in str(exc_info.value)

def test_batch_processor(csv_service, mock_file_storage_service, mock_molecule_service):
    """
    Test batch processing functionality.
    """
    # Create a test CSV file
    test_data = [{'smiles': 'CCO', 'MW': 46.07, 'LogP': -0.14},
                 {'smiles': 'c1ccccc1', 'MW': 78.11, 'LogP': 1.90}]
    file_path = create_test_csv('test.csv', test_data)

    # Mock the file storage service to return the test CSV file
    with open(file_path, 'rb') as f:
        mock_file_storage_service.get_csv_file.return_value = (io.BytesIO(f.read()), 'test.csv', 'text/csv')

    # Mock the molecule service to return a successful bulk creation result
    mock_molecule_service.bulk_create.return_value = {'created_molecules': [], 'total_created': 2, 'total_duplicates': 0, 'total_processed': 2, 'total_failed': 0}

    # Define test parameters
    file_id = 'test_file_id'
    mapping = {'smiles': 'smiles'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the process_csv_file method
    processing_results = csv_service.process_csv_file(job_id, file_id, mapping, user_id)

    # Assert that the processing results are correct
    assert processing_results['status'] == csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name
    assert processing_results['summary']['total_rows'] == 2
    assert processing_results['summary']['processed_rows'] == 2
    assert processing_results['summary']['successful_rows'] == 2
    assert processing_results['summary']['failed_rows'] == 0

    # Clean up the test CSV file
    os.remove(file_path)

def test_performance_large_csv(csv_service, mock_file_storage_service, mock_molecule_service):
    """
    Test performance with large CSV files.
    """
    # Create a large test CSV file with 10,000 rows
    num_rows = 10000
    test_data = [{'smiles': 'CCO', 'MW': 46.07, 'LogP': -0.14} for _ in range(num_rows)]
    file_path = create_test_csv('large_test.csv', test_data)

    # Mock the file storage service to return the test CSV file
    with open(file_path, 'rb') as f:
        mock_file_storage_service.get_csv_file.return_value = (io.BytesIO(f.read()), 'large_test.csv', 'text/csv')

    # Mock the molecule service to return a successful bulk creation result
    mock_molecule_service.bulk_create.return_value = {'created_molecules': [], 'total_created': num_rows, 'total_duplicates': 0, 'total_processed': num_rows, 'total_failed': 0}

    # Define test parameters
    file_id = 'large_test_file_id'
    mapping = {'smiles': 'smiles'}
    user_id = 1
    job_id = csv_service.create_processing_job(file_id, mapping, user_id)

    # Call the process_csv_file method and measure the execution time
    import time
    start_time = time.time()
    processing_results = csv_service.process_csv_file(job_id, file_id, mapping, user_id)
    end_time = time.time()

    # Assert that the processing results are correct
    assert processing_results['status'] == csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name
    assert processing_results['summary']['total_rows'] == num_rows
    assert processing_results['summary']['processed_rows'] == num_rows
    assert processing_results['summary']['successful_rows'] == num_rows
    assert processing_results['summary']['failed_rows'] == 0

    # Assert that the processing time is less than 30 seconds
    processing_time = end_time - start_time
    assert processing_time < 30

    # Clean up the test CSV file
    os.remove(file_path)