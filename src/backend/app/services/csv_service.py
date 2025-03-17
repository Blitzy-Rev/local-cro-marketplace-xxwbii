# src/backend/app/services/csv_service.py
"""
Service layer for CSV file processing in the Molecular Data Management and CRO Integration Platform.
This service provides high-level business logic for uploading, validating, mapping, and processing CSV files containing molecular data. It coordinates between API endpoints, file storage, and background processing tasks.
"""

import os  # standard library
import io  # standard library
import uuid  # standard library
import tempfile  # standard library
import datetime  # standard library
import typing  # standard library
from typing import Dict, List, Optional, Any, Tuple, BinaryIO, Union  # standard library

import pandas as pd  # pandas 2.0+
from sqlalchemy import exc  # sqlalchemy 2.0+
from sqlalchemy.orm import Session  # sqlalchemy.orm 2.0+

from ..logging_config import logger  # Application logger for error and info logging
from ..exceptions import CSVProcessingException, ValidationException  # Exception for CSV processing errors
from .file_storage_service import FileStorageService  # Service for file storage operations
from .molecule_service import MoleculeService  # Service for molecule operations
from ..db.session import get_db  # Database session context manager
from ..utils import csv_utils  # Utility functions to read CSV files
from ..api.api_v1.schemas import csv as csv_schemas  # Enumeration of CSV processing status values
from ..constants import DEFAULT_MOLECULE_PROPERTIES, CSV_CHUNK_SIZE  # List of default molecular properties


class CSVService:
    """
    Service class for CSV file processing operations
    """

    def __init__(self, file_storage_service: FileStorageService, molecule_service: Optional[MoleculeService] = None):
        """
        Initializes the CSVService with file storage service

        Args:
            file_storage_service: The FileStorageService instance to use for file operations
            molecule_service: The MoleculeService instance to use for molecule operations
        """
        self.file_storage_service = file_storage_service
        self.molecule_service = molecule_service
        self.logger = logger

    def validate_csv_file(self, file_id: str) -> Dict[str, Any]:
        """
        Validates a CSV file structure and content

        Args:
            file_id: The ID of the CSV file to validate

        Returns:
            Validation results with headers and row count
        """
        try:
            # Get CSV file from storage
            file_data, filename, _ = self.file_storage_service.get_csv_file(file_id)

            # Save file to temporary location
            with tempfile.NamedTemporaryFile(delete=False, suffix=filename) as temp_file:
                temp_file.write(file_data.read())
                temp_file_path = temp_file.name

            # Validate CSV structure
            validation_results = csv_utils.validate_csv_structure(temp_file_path)

            # Clean up temporary file
            os.remove(temp_file_path)

            return validation_results
        except Exception as e:
            raise CSVProcessingException(f"Error validating CSV file: {str(e)}")

    def get_csv_headers(self, file_id: str) -> List[str]:
        """
        Gets headers from a CSV file

        Args:
            file_id: The ID of the CSV file

        Returns:
            List of CSV header column names
        """
        try:
            # Get CSV file from storage
            file_data, filename, _ = self.file_storage_service.get_csv_file(file_id)

            # Save file to temporary location
            with tempfile.NamedTemporaryFile(delete=False, suffix=filename) as temp_file:
                temp_file.write(file_data.read())
                temp_file_path = temp_file.name

            # Extract headers
            headers = csv_utils.get_csv_headers(temp_file_path)

            # Clean up temporary file
            os.remove(temp_file_path)

            return headers
        except Exception as e:
            raise CSVProcessingException(f"Error getting CSV headers: {str(e)}")

    def validate_column_mapping(self, headers: List[str], mapping: Dict[str, str]) -> bool:
        """
        Validates mapping of CSV columns to system properties

        Args:
            headers: List of CSV header column names
            mapping: Dictionary mapping CSV column names to system properties

        Returns:
            True if mapping is valid
        """
        try:
            is_valid, error_messages = csv_utils.validate_mapping(headers, mapping)
            if not is_valid:
                raise ValidationException(f"Invalid column mapping: {', '.join(error_messages)}")
            return True
        except Exception as e:
            raise ValidationException(f"Error validating column mapping: {str(e)}")

    def get_available_properties(self) -> List[Dict[str, Any]]:
        """
        Gets list of available system properties for mapping

        Returns:
            List of property definitions with name, description, and type
        """
        try:
            # Get available properties from molecule service
            properties = self.molecule_service.get_available_properties()

            # Format properties as dictionaries
            formatted_properties = [{"name": prop, "description": prop, "data_type": "string"} for prop in properties]

            return formatted_properties
        except Exception as e:
            raise Exception(f"Error getting available properties: {str(e)}")

    def create_processing_job(self, file_id: str, mapping: Dict[str, str], user_id: int) -> str:
        """
        Creates a new CSV processing job

        Args:
            file_id: The ID of the CSV file to process
            mapping: The mapping of CSV columns to system properties
            user_id: The ID of the user who initiated the job

        Returns:
            Job ID for the created processing job
        """
        try:
            job_id = str(uuid.uuid4())
            # Store job status in a global dictionary (replace with a proper cache)
            JOB_STATUS_CACHE[job_id] = {
                "status": csv_schemas.CSV_PROCESSING_STATUS.PENDING.name,
                "file_id": file_id,
                "mapping": mapping,
                "user_id": user_id,
                "created_at": datetime.datetime.now().isoformat(),
                "last_updated": datetime.datetime.now().isoformat()
            }
            return job_id
        except Exception as e:
            raise Exception(f"Error creating processing job: {str(e)}")

    def update_job_status(self, job_id: str, status: str, progress: Optional[int] = None, summary: Optional[Dict[str, Any]] = None, errors: Optional[Dict[str, Any]] = None) -> bool:
        """
        Updates the status of a CSV processing job

        Args:
            job_id: The ID of the job to update
            status: The new status of the job
            progress: The processing progress as a percentage (optional)
            summary: A summary of the processing results (optional)
            errors: A dictionary of errors encountered during processing (optional)

        Returns:
            True if update was successful
        """
        try:
            if job_id not in JOB_STATUS_CACHE:
                return False

            JOB_STATUS_CACHE[job_id]["status"] = status
            JOB_STATUS_CACHE[job_id]["last_updated"] = datetime.datetime.now().isoformat()

            if progress is not None:
                JOB_STATUS_CACHE[job_id]["progress"] = progress
            if summary is not None:
                JOB_STATUS_CACHE[job_id]["summary"] = summary
            if errors is not None:
                JOB_STATUS_CACHE[job_id]["errors"] = errors

            return True
        except Exception as e:
            raise Exception(f"Error updating job status: {str(e)}")

    def get_job_status(self, job_id: str, user_id: Optional[int] = None) -> Dict[str, Any]:
        """
        Gets the current status of a CSV processing job

        Args:
            job_id: The ID of the job to retrieve
            user_id: The ID of the user requesting the status (optional)

        Returns:
            Job status information with status, progress, and summary
        """
        try:
            if job_id not in JOB_STATUS_CACHE:
                raise Exception(f"Job with ID {job_id} not found")

            job = JOB_STATUS_CACHE[job_id]

            # Verify job belongs to user if user_id is provided
            if user_id is not None and job["user_id"] != user_id:
                raise Exception("Job does not belong to user")

            return job
        except Exception as e:
            raise Exception(f"Error getting job status: {str(e)}")

    def cancel_job(self, job_id: str, user_id: int) -> bool:
        """
        Cancels an ongoing CSV processing job

        Args:
            job_id: The ID of the job to cancel
            user_id: The ID of the user requesting the cancellation

        Returns:
            True if cancellation was successful
        """
        try:
            if job_id not in JOB_STATUS_CACHE:
                return False

            job = JOB_STATUS_CACHE[job_id]

            # Verify job belongs to user
            if job["user_id"] != user_id:
                return False

            # Check if job is in a cancellable state
            if job["status"] not in [csv_schemas.CSV_PROCESSING_STATUS.PENDING.name, csv_schemas.CSV_PROCESSING_STATUS.PROCESSING.name]:
                return False

            # Update job status to CANCELLED
            job["status"] = csv_schemas.CSV_PROCESSING_STATUS.CANCELLED.name
            job["last_updated"] = datetime.datetime.now().isoformat()

            return True
        except Exception as e:
            raise Exception(f"Error cancelling job: {str(e)}")

    def process_csv_file(self, job_id: str, file_id: str, mapping: Dict[str, str], user_id: int) -> Dict[str, Any]:
        """
        Processes a CSV file and extracts molecules

        Args:
            job_id: The ID of the processing job
            file_id: The ID of the CSV file to process
            mapping: The mapping of CSV columns to system properties
            user_id: The ID of the user who initiated the job

        Returns:
            Processing results with statistics and status
        """
        try:
            # Update job status to PROCESSING
            self.update_job_status(job_id, csv_schemas.CSV_PROCESSING_STATUS.PROCESSING.name)

            # Get CSV file from storage
            file_data, filename, _ = self.file_storage_service.get_csv_file(file_id)

            # Save file to temporary location
            with tempfile.NamedTemporaryFile(delete=False, suffix=filename) as temp_file:
                temp_file.write(file_data.read())
                temp_file_path = temp_file.name

            # Process CSV in batches
            def batch_processor(df: pd.DataFrame, batch_index: int, total_row_count: int) -> Dict[str, Any]:
                """Processes a batch of CSV data and extracts molecules"""
                # Map CSV columns to system properties
                mapped_df = csv_utils.map_csv_columns(df, mapping)

                # Extract molecules from DataFrame
                extraction_results = csv_utils.extract_molecules_from_dataframe(mapped_df)

                # Calculate progress
                progress = int((batch_index * CSV_CHUNK_SIZE / total_row_count) * 100)

                # Return batch processing results
                return {
                    'batch_index': batch_index,
                    'processed_rows': len(df),
                    'successful_rows': extraction_results['valid_count'],
                    'failed_rows': extraction_results['error_count'],
                    'molecules': extraction_results['valid_molecules'],
                    'errors': extraction_results['errors'],
                    'progress': progress
                }

            processing_results = csv_utils.process_csv_in_batches(temp_file_path, batch_processor)

            # Extract molecules from processing results
            molecules = []
            for batch_result in processing_results['batch_results']:
                if 'molecules' in batch_result:
                    molecules.extend(batch_result['molecules'])

            # Create molecules in database
            bulk_creation_results = self.molecule_service.bulk_create(molecules, user_id)

            # Generate import summary
            import_summary = csv_utils.generate_csv_import_summary(processing_results)

            # Update job status to COMPLETED with summary
            self.update_job_status(job_id, csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name, summary=import_summary)

            # Clean up temporary file
            os.remove(temp_file_path)

            return {
                "status": csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name,
                "summary": import_summary,
                "bulk_creation_results": bulk_creation_results
            }
        except Exception as e:
            # Update job status to FAILED
            self.update_job_status(job_id, csv_schemas.CSV_PROCESSING_STATUS.FAILED.name, errors={"error": str(e)})
            raise

    def cleanup_job_resources(self, job_id: str) -> bool:
        """
        Cleans up resources associated with a completed job

        Args:
            job_id: The ID of the job to cleanup

        Returns:
            True if cleanup was successful
        """
        try:
            if job_id not in JOB_STATUS_CACHE:
                return False

            # Get file_id from job status
            file_id = JOB_STATUS_CACHE[job_id]["file_id"]

            # Delete CSV file from storage
            self.file_storage_service.delete_file(file_id, BUCKET_NAMES["CSV_UPLOADS"])

            # Remove job from JOB_STATUS_CACHE if older than retention period
            del JOB_STATUS_CACHE[job_id]

            return True
        except Exception as e:
            raise Exception(f"Error cleaning up job resources: {str(e)}")

    def get_job_statistics(self, user_id: Optional[int] = None, status: Optional[str] = None) -> Dict[str, Any]:
        """
        Gets statistics for completed CSV processing jobs

        Args:
            user_id: The ID of the user to filter by (optional)
            status: The status to filter by (optional)

        Returns:
            Job statistics with counts by status and user
        """
        try:
            total_pending = 0
            total_processing = 0
            total_completed = 0
            total_failed = 0
            total_cancelled = 0
            total_jobs = 0

            for job_id, job in JOB_STATUS_CACHE.items():
                # Filter by user_id if provided
                if user_id is not None and job["user_id"] != user_id:
                    continue

                # Filter by status if provided
                if status is not None and job["status"] != status:
                    continue

                total_jobs += 1

                if job["status"] == csv_schemas.CSV_PROCESSING_STATUS.PENDING.name:
                    total_pending += 1
                elif job["status"] == csv_schemas.CSV_PROCESSING_STATUS.PROCESSING.name:
                    total_processing += 1
                elif job["status"] == csv_schemas.CSV_PROCESSING_STATUS.COMPLETED.name:
                    total_completed += 1
                elif job["status"] == csv_schemas.CSV_PROCESSING_STATUS.FAILED.name:
                    total_failed += 1
                elif job["status"] == csv_schemas.CSV_PROCESSING_STATUS.CANCELLED.name:
                    total_cancelled += 1

            return {
                "total_jobs": total_jobs,
                "total_pending": total_pending,
                "total_processing": total_processing,
                "total_completed": total_completed,
                "total_failed": total_failed,
                "total_cancelled": total_cancelled
            }
        except Exception as e:
            raise Exception(f"Error getting job statistics: {str(e)}")


# Global variable to store job status (replace with a proper cache)
JOB_STATUS_CACHE: Dict[str, Dict[str, Any]] = {}

def generate_job_id() -> str:
    """
    Generates a unique job ID for CSV processing

    Returns:
        Unique job ID
    """
    return str(uuid.uuid4())

def batch_processor(df: pd.DataFrame, batch_index: int, total_batches: int, mapping: Dict[str, str]) -> Dict[str, Any]:
    """
    Processes a batch of CSV data and extracts molecules

    Args:
        df: DataFrame containing CSV data
        batch_index: Index of the current batch
        total_batches: Total number of batches
        mapping: Mapping of CSV columns to system properties

    Returns:
        Batch processing results with molecules and statistics
    """
    # Map CSV columns to system properties
    mapped_df = csv_utils.map_csv_columns(df, mapping)

    # Extract molecules from DataFrame
    extraction_results = csv_utils.extract_molecules_from_dataframe(mapped_df)

    # Calculate progress
    progress = int((batch_index / total_batches) * 100)

    # Return batch processing results
    return {
        'batch_index': batch_index,
        'processed_rows': len(df),
        'successful_rows': extraction_results['valid_count'],
        'failed_rows': extraction_results['error_count'],
        'molecules': extraction_results['valid_molecules'],
        'errors': extraction_results['errors'],
        'progress': progress
    }