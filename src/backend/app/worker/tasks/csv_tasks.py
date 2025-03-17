"""
Celery task module for asynchronous CSV file processing in the Molecular Data Management and CRO Integration Platform.
Implements background tasks for processing large CSV files containing molecular data, handling validation, property calculation, and database storage in an efficient manner.
"""

import os  # standard library
import logging  # standard library
from typing import Dict, List, Any, Optional  # standard library

from celery import Celery  # celery 5.2+

from ..celery_app import app  # Import Celery application instance for task registration
from ..logging_config import logger  # Application logger for task logging
from ...services.csv_service import CSVService  # Service for CSV file processing operations
from ...services.file_storage_service import FileStorageService  # Service for file storage operations
from ...services.molecule_service import MoleculeService  # Service for molecule operations
from ...api.api_v1.schemas.csv import CSV_PROCESSING_STATUS  # Enumeration of CSV processing status values
from ...constants import BUCKET_NAMES  # Dictionary of bucket names for different file types
from ...constants import CSV_CHUNK_SIZE  # Number of rows to process in each chunk for CSV processing


# Initialize services
csv_service = CSVService(file_storage_service=FileStorageService(), molecule_service=MoleculeService())
file_storage_service = FileStorageService()
molecule_service = MoleculeService()

# Constants
CSV_FILE_RETENTION_DAYS = 7  # Number of days to retain CSV files


@app.task(bind=True, name='csv_tasks.process_csv', max_retries=3, acks_late=True)
def process_csv(self, job_id: str, file_id: str, mapping: Dict[str, str], user_id: int) -> Dict[str, Any]:
    """
    Celery task for processing a CSV file asynchronously.

    Args:
        job_id (str): Unique identifier for the processing job.
        file_id (str): Identifier of the CSV file in storage.
        mapping (Dict[str, str]): Mapping of CSV columns to system properties.
        user_id (int): Identifier of the user initiating the process.

    Returns:
        Dict[str, Any]: Processing results with statistics.
    """
    logger.info(f"Starting CSV processing task for job_id: {job_id}, file_id: {file_id}")
    try:
        # Process the CSV file using csv_service
        processing_results = csv_service.process_csv_file(job_id, file_id, mapping, user_id)
        logger.info(f"CSV processing completed successfully for job_id: {job_id}")
        return processing_results
    except Exception as e:
        logger.error(f"CSV processing failed for job_id: {job_id} with error: {str(e)}")
        # Update job status to FAILED with error information
        csv_service.update_job_status(job_id, CSV_PROCESSING_STATUS.FAILED.name, errors={"error": str(e)})
        # Retry the task if appropriate (up to max_retries)
        if self.request.retries < self.max_retries:
            logger.info(f"Retrying task for job_id: {job_id}, attempt: {self.request.retries + 1}")
            raise self.retry(exc=e, countdown=60)  # Retry after 60 seconds
        else:
            logger.error(f"Max retries reached for job_id: {job_id}, task failed permanently")
            raise  # Re-raise exception if retry limit reached


@app.task(name='csv_tasks.cleanup_csv_files')
def cleanup_csv_files(file_ids: List[str] = None) -> Dict[str, Any]:
    """
    Celery task for cleaning up old CSV files based on retention policy.

    Args:
        file_ids (List[str], optional): List of specific file IDs to clean up. If None, cleans up all files older than retention period.

    Returns:
        Dict[str, Any]: Cleanup results with statistics (deleted count, failed count).
    """
    logger.info("Starting CSV file cleanup task")
    deleted_count = 0
    failed_count = 0

    try:
        # Calculate cutoff date based on CSV_FILE_RETENTION_DAYS
        cutoff_date = (os.path.utcnow() - os.path.timedelta(days=CSV_FILE_RETENTION_DAYS)).isoformat()

        # List files in CSV uploads bucket
        files = file_storage_service.list_files(BUCKET_NAMES["CSV_UPLOADS"])

        # Filter files by last modified date older than cutoff date
        files_to_delete = [f for f in files if f["last_modified"].isoformat() < cutoff_date]

        # If specific file_ids are provided, filter the list further
        if file_ids:
            files_to_delete = [f for f in files_to_delete if f["name"] in file_ids]

        # Delete files
        for file in files_to_delete:
            try:
                file_storage_service.delete_file(file["name"], BUCKET_NAMES["CSV_UPLOADS"])
                logger.info(f"Deleted CSV file: {file['name']}")
                deleted_count += 1
            except Exception as e:
                logger.error(f"Failed to delete CSV file: {file['name']} with error: {str(e)}")
                failed_count += 1

        logger.info(f"CSV file cleanup completed: {deleted_count} deleted, {failed_count} failed")
        return {"deleted_count": deleted_count, "failed_count": failed_count}
    except Exception as e:
        logger.error(f"CSV file cleanup task failed with error: {str(e)}")
        return {"deleted_count": deleted_count, "failed_count": failed_count, "error": str(e)}


@app.task(name='csv_tasks.cancel_csv_processing')
def cancel_csv_processing(job_id: str, user_id: int) -> Dict[str, Any]:
    """
    Celery task for cancelling an ongoing CSV processing job.

    Args:
        job_id (str): Unique identifier for the processing job.
        user_id (int): Identifier of the user requesting the cancellation.

    Returns:
        Dict[str, Any]: Cancellation result.
    """
    logger.info(f"Starting cancellation request for job_id: {job_id}")
    try:
        # Update job status to CANCELLED using csv_service
        cancellation_result = csv_service.cancel_job(job_id, user_id)
        logger.info(f"Cancellation result for job_id: {job_id} - {cancellation_result}")
        return {"job_id": job_id, "status": "cancelled", "success": cancellation_result}
    except Exception as e:
        logger.error(f"Cancellation request failed for job_id: {job_id} with error: {str(e)}")
        return {"job_id": job_id, "status": "cancellation_failed", "success": False, "error": str(e)}


@app.task(name='csv_tasks.retry_csv_processing')
def retry_csv_processing(job_id: str, user_id: int) -> Dict[str, Any]:
    """
    Celery task for retrying a failed CSV processing job.

    Args:
        job_id (str): Unique identifier for the failed processing job.
        user_id (int): Identifier of the user requesting the retry.

    Returns:
        Dict[str, Any]: Retry result with new job_id.
    """
    logger.info(f"Starting retry request for job_id: {job_id}")
    try:
        # Get job status information from csv_service
        job_status = csv_service.get_job_status(job_id, user_id)

        # If job status is not FAILED, return error
        if job_status["status"] != CSV_PROCESSING_STATUS.FAILED.name:
            logger.error(f"Retry request failed for job_id: {job_id} - Job status is not FAILED")
            return {"job_id": job_id, "status": "retry_failed", "success": False, "error": "Job status is not FAILED"}

        # Create a new processing job with the same parameters
        file_id = job_status["file_id"]
        mapping = job_status["mapping"]
        new_job_id = csv_service.create_processing_job(file_id, mapping, user_id)

        # Trigger process_csv task with new job_id
        process_csv.delay(new_job_id, file_id, mapping, user_id)

        logger.info(f"Retry request succeeded for job_id: {job_id} - New job_id: {new_job_id}")
        return {"old_job_id": job_id, "new_job_id": new_job_id, "status": "retried", "success": True}
    except Exception as e:
        logger.error(f"Retry request failed for job_id: {job_id} with error: {str(e)}")
        return {"job_id": job_id, "status": "retry_failed", "success": False, "error": str(e)}


@app.task(name='csv_tasks.get_csv_processing_stats')
def get_csv_processing_stats(user_id: Optional[int] = None) -> Dict[str, Any]:
    """
    Celery task for retrieving CSV processing statistics.

    Args:
        user_id (Optional[int], optional): Identifier of the user to filter by. Defaults to None.

    Returns:
        Dict[str, Any]: Processing statistics.
    """
    logger.info(f"Starting statistics request for user_id: {user_id}")
    try:
        # Get job statistics from csv_service
        statistics = csv_service.get_job_statistics(user_id)
        logger.info(f"Statistics request succeeded for user_id: {user_id} - Statistics: {statistics}")
        return statistics
    except Exception as e:
        logger.error(f"Statistics request failed for user_id: {user_id} with error: {str(e)}")
        return {"error": str(e)}