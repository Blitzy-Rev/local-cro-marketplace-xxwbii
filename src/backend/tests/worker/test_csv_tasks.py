# src/backend/tests/worker/test_csv_tasks.py
import pytest  # pytest 7.0+
from unittest.mock import patch, MagicMock, Mock  # standard library
from datetime import datetime, timedelta  # standard library
import io  # standard library
from io import BytesIO  # standard library
from typing import Dict, List, Any  # standard library

import celery  # celery 5.2+

from app.worker.tasks.csv_tasks import process_csv, cleanup_csv_files, cancel_csv_processing, retry_csv_processing, get_csv_processing_stats  # Import Celery tasks
from app.services.csv_service import CSVService  # Import CSVService
from app.services.file_storage_service import FileStorageService  # Import FileStorageService
from app.services.molecule_service import MoleculeService  # Import MoleculeService
from app.api.api_v1.schemas.csv import CSV_PROCESSING_STATUS  # Import CSV_PROCESSING_STATUS enum


@pytest.mark.parametrize('job_id,file_id,mapping,user_id', [('job-123', 'file-123', {'col1': 'SMILES', 'col2': 'LogP'}, 1)])
def test_process_csv_success(job_id: str, file_id: str, mapping: Dict[str, str], user_id: int) -> None:
    """Tests successful CSV processing task execution"""
    # Mock CSVService.process_csv_file to return success result
    with patch.object(CSVService, 'process_csv_file', return_value={'status': 'COMPLETED', 'message': 'Success'}) as mock_process_csv:
        # Call process_csv task with test parameters
        result = process_csv(job_id=job_id, file_id=file_id, mapping=mapping, user_id=user_id)

        # Verify CSVService.process_csv_file was called with correct parameters
        mock_process_csv.assert_called_once_with(job_id, file_id, mapping, user_id)

        # Verify task returns expected success result
        assert result == {'status': 'COMPLETED', 'message': 'Success'}


@pytest.mark.parametrize('job_id,file_id,mapping,user_id', [('job-123', 'file-123', {'col1': 'SMILES', 'col2': 'LogP'}, 1)])
def test_process_csv_failure(job_id: str, file_id: str, mapping: Dict[str, str], user_id: int) -> None:
    """Tests CSV processing task error handling"""
    # Mock CSVService.process_csv_file to raise an exception
    with patch.object(CSVService, 'process_csv_file', side_effect=Exception("Test Exception")) as mock_process_csv,\
         patch.object(CSVService, 'update_job_status') as mock_update_job_status,\
         patch('app.worker.tasks.csv_tasks.process_csv.retry') as mock_task_retry:
        mock_task_retry.side_effect = Exception("Retry Exception")  # Simulate retry failure

        # Call process_csv task with test parameters
        with pytest.raises(Exception, match="Retry Exception"):
            process_csv(job_id=job_id, file_id=file_id, mapping=mapping, user_id=user_id)

        # Verify CSVService.update_job_status was called with FAILED status
        mock_update_job_status.assert_called_once_with(job_id, CSV_PROCESSING_STATUS.FAILED.name, errors={'error': 'Test Exception'})

        # Verify task.retry was called with appropriate parameters
        mock_task_retry.assert_called_once()

def test_cleanup_csv_files_with_ids() -> None:
    """Tests CSV file cleanup with specific file IDs"""
    # Create list of test file IDs
    file_ids = ['file-1', 'file-2', 'file-3']

    # Mock FileStorageService.delete_file to return success
    with patch.object(FileStorageService, 'delete_file', return_value=True) as mock_delete_file:
        # Call cleanup_csv_files task with file IDs
        result = cleanup_csv_files(file_ids=file_ids)

        # Verify FileStorageService.delete_file was called for each file ID
        assert mock_delete_file.call_count == len(file_ids)
        for file_id in file_ids:
            mock_delete_file.assert_any_call(file_id, 'csv-uploads')

        # Verify task returns expected cleanup statistics
        assert result == {'deleted_count': len(file_ids), 'failed_count': 0}


def test_cleanup_csv_files_by_age() -> None:
    """Tests CSV file cleanup based on file age"""
    # Mock FileStorageService.list_files to return test files with timestamps
    now = datetime.utcnow()
    old_file_time = now - timedelta(days=8)
    recent_file_time = now - timedelta(days=2)
    test_files = [
        {'name': 'old-file.csv', 'last_modified': old_file_time},
        {'name': 'recent-file.csv', 'last_modified': recent_file_time}
    ]
    with patch.object(FileStorageService, 'list_files', return_value=test_files) as mock_list_files,\
         patch.object(FileStorageService, 'delete_file', return_value=True) as mock_delete_file:
        # Call cleanup_csv_files task with empty file IDs list
        result = cleanup_csv_files(file_ids=[])

        # Verify FileStorageService.list_files was called with correct bucket
        mock_list_files.assert_called_once_with('csv-uploads')

        # Verify FileStorageService.delete_file was called for files older than retention period
        mock_delete_file.assert_called_once_with('old-file.csv', 'csv-uploads')

        # Verify task returns expected cleanup statistics
        assert result == {'deleted_count': 1, 'failed_count': 0}


@pytest.mark.parametrize('job_id,user_id', [('job-123', 1)])
def test_cancel_csv_processing(job_id: str, user_id: int) -> None:
    """Tests cancellation of CSV processing job"""
    # Mock CSVService.update_job_status to return success
    with patch.object(CSVService, 'cancel_job', return_value=True) as mock_cancel_job:
        # Call cancel_csv_processing task with test parameters
        result = cancel_csv_processing(job_id=job_id, user_id=user_id)

        # Verify CSVService.update_job_status was called with CANCELLED status
        mock_cancel_job.assert_called_once_with(job_id, user_id)

        # Verify task returns expected cancellation result
        assert result == {'job_id': job_id, 'status': 'cancelled', 'success': True}


@pytest.mark.parametrize('job_id,user_id', [('job-123', 1)])
def test_retry_csv_processing_success(job_id: str, user_id: int) -> None:
    """Tests successful retry of failed CSV processing job"""
    # Mock CSVService.get_job_status to return FAILED status
    with patch.object(CSVService, 'get_job_status', return_value={'status': CSV_PROCESSING_STATUS.FAILED.name, 'file_id': 'file-123', 'mapping': {'col1': 'SMILES'}}),\
         patch.object(CSVService, 'create_processing_job', return_value='new-job-456') as mock_create_processing_job,\
         patch('app.worker.tasks.csv_tasks.process_csv.delay') as mock_process_csv_delay:
        # Call retry_csv_processing task with test parameters
        result = retry_csv_processing(job_id=job_id, user_id=user_id)

        # Verify CSVService.get_job_status was called with correct job_id
        # Verify process_csv.delay was called with new job parameters
        mock_process_csv_delay.assert_called_once_with('new-job-456', 'file-123', {'col1': 'SMILES'}, user_id)

        # Verify task returns expected retry result with new job_id
        assert result == {'old_job_id': job_id, 'new_job_id': 'new-job-456', 'status': 'retried', 'success': True}


@pytest.mark.parametrize('job_id,user_id', [('job-123', 1)])
def test_retry_csv_processing_not_failed(job_id: str, user_id: int) -> None:
    """Tests retry attempt on a job that is not in FAILED state"""
    # Mock CSVService.get_job_status to return non-FAILED status
    with patch.object(CSVService, 'get_job_status', return_value={'status': CSV_PROCESSING_STATUS.COMPLETED.name}),\
         patch('app.worker.tasks.csv_tasks.process_csv.delay') as mock_process_csv_delay:
        # Call retry_csv_processing task with test parameters
        result = retry_csv_processing(job_id=job_id, user_id=user_id)

        # Verify CSVService.get_job_status was called with correct job_id
        # Verify process_csv.delay was not called
        mock_process_csv_delay.assert_not_called()

        # Verify task returns error result indicating job is not in FAILED state
        assert result == {'job_id': job_id, 'status': 'retry_failed', 'success': False, 'error': 'Job status is not FAILED'}


@pytest.mark.parametrize('user_id', [1, None])
def test_get_csv_processing_stats(user_id: int) -> None:
    """Tests retrieval of CSV processing statistics"""
    # Mock CSVService.get_job_statistics to return test statistics
    with patch.object(CSVService, 'get_job_statistics', return_value={'total_jobs': 10, 'total_pending': 2}) as mock_get_job_statistics:
        # Call get_csv_processing_stats task with test user_id
        result = get_csv_processing_stats(user_id=user_id)

        # Verify CSVService.get_job_statistics was called with correct user_id
        mock_get_job_statistics.assert_called_once_with(user_id)

        # Verify task returns expected statistics
        assert result == {'total_jobs': 10, 'total_pending': 2}