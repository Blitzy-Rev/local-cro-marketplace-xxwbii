"""
Implements Celery tasks for asynchronous processing of experimental results in the Molecular Data Management and CRO Integration Platform. This module handles background operations related to result file processing, data extraction, validation, and notification delivery for experimental results uploaded by CROs.
"""

import logging
from typing import Dict, List, Any, Optional, Union, Tuple, UUID
import io
from io import BytesIO

import pandas  # pandas 2.0+

from celery import Celery  # celery 5.2+

from ..celery_app import app
from ...services.result_service import result_service
from ...services.file_storage_service import file_storage_service
from ...services.notification_service import send_results_uploaded
from ...schemas.result import ResultDataCreate, ResultUpdate
from ...constants import ResultStatus
from ...exceptions import FileStorageException, ValidationException, ResourceNotFoundException

# Initialize logger
logger = logging.getLogger(__name__)


@app.task(name='result_tasks.process_result_file', queue='results')
def process_result_file(result_id: UUID, file_id: UUID, cro_user_id: int) -> Dict[str, Any]:
    """
    Processes a result file and extracts structured data.

    Args:
        result_id: ID of the result
        file_id: ID of the file
        cro_user_id: ID of the CRO user

    Returns:
        Dict[str, Any]: Processing result with extracted data count
    """
    logger.info(f"Attempting to process result file for result ID: {result_id}, file ID: {file_id}")
    try:
        # Retrieve result file
        file_data, filename, content_type = file_storage_service.get_result_file(file_id)

        # Determine file type and process accordingly
        if filename.endswith('.csv'):
            extracted_data = process_csv_result(file_data, result_id, cro_user_id)
        elif filename.endswith(('.xls', '.xlsx')):
            extracted_data = process_excel_result(file_data, result_id, cro_user_id)
        else:
            raise ValidationException(f"Unsupported file type: {filename}", details={"filename": filename})

        # Add data to result
        for data_point in extracted_data:
            result_data_create = ResultDataCreate(**data_point)
            result_service.add_result_data(result_data_create)

        # Update result status to UPLOADED if currently PENDING
        db_result = result_service.get_result_by_id(result_id)
        if db_result and db_result.status == ResultStatus.PENDING.name:
            result_update = ResultUpdate(status=ResultStatus.UPLOADED.name)
            result_service.update_result(result_id, result_update)

        data_point_count = len(extracted_data)
        logger.info(f"Successfully processed result file for result ID: {result_id}, extracted {data_point_count} data points")
        return {"status": "success", "data_point_count": data_point_count}

    except FileStorageException as e:
        logger.error(f"File storage error processing result file: {str(e)}")
        return {"status": "error", "message": str(e)}
    except ValidationException as e:
        logger.error(f"Validation error processing result file: {str(e)}")
        return {"status": "error", "message": str(e)}
    except Exception as e:
        logger.exception(f"Unexpected error processing result file: {str(e)}")
        return {"status": "error", "message": "An unexpected error occurred"}


def process_csv_result(file_data: io.BytesIO, result_id: UUID, cro_user_id: int) -> List[Dict[str, Any]]:
    """
    Processes a CSV result file and extracts data.

    Args:
        file_data: CSV file data as BytesIO object
        result_id: ID of the result
        cro_user_id: ID of the CRO user

    Returns:
        List[Dict[str, Any]]: List of extracted data points
    """
    logger.info(f"Processing CSV result file for result ID: {result_id}")
    try:
        # Read CSV data using pandas
        df = pandas.read_csv(file_data)

        # Validate CSV structure
        if 'molecule_id' not in df.columns:
            raise ValidationException("CSV must contain a 'molecule_id' column", details={"result_id": result_id})
        if len(df.columns) < 2:
            raise ValidationException("CSV must contain at least one data column", details={"result_id": result_id})

        extracted_data = []
        for index, row in df.iterrows():
            molecule_id = row['molecule_id']
            for column in df.columns:
                if column != 'molecule_id':
                    data_name = column
                    data_value = row[column]

                    # Create data point dictionary
                    data_point = {
                        "result_id": result_id,
                        "molecule_id": molecule_id,
                        "data_name": data_name,
                        "data_value": data_value,
                        "data_unit": None  # Unit is not available in CSV
                    }
                    extracted_data.append(data_point)

        logger.info(f"Successfully extracted {len(extracted_data)} data points from CSV for result ID: {result_id}")
        return extracted_data

    except ValidationException as e:
        logger.error(f"Validation error processing CSV result file: {str(e)}")
        raise
    except Exception as e:
        logger.exception(f"Unexpected error processing CSV result file: {str(e)}")
        raise


def process_excel_result(file_data: io.BytesIO, result_id: UUID, cro_user_id: int) -> List[Dict[str, Any]]:
    """
    Processes an Excel result file and extracts data.

    Args:
        file_data: Excel file data as BytesIO object
        result_id: ID of the result
        cro_user_id: ID of the CRO user

    Returns:
        List[Dict[str, Any]]: List of extracted data points
    """
    logger.info(f"Processing Excel result file for result ID: {result_id}")
    try:
        # Read Excel data using pandas
        df = pandas.read_excel(file_data)

        # Validate Excel structure
        if 'molecule_id' not in df.columns:
            raise ValidationException("Excel must contain a 'molecule_id' column", details={"result_id": result_id})
        if len(df.columns) < 2:
            raise ValidationException("Excel must contain at least one data column", details={"result_id": result_id})

        extracted_data = []
        for index, row in df.iterrows():
            molecule_id = row['molecule_id']
            for column in df.columns:
                if column != 'molecule_id':
                    data_name = column
                    data_value = row[column]

                    # Create data point dictionary
                    data_point = {
                        "result_id": result_id,
                        "molecule_id": molecule_id,
                        "data_name": data_name,
                        "data_value": data_value,
                        "data_unit": None  # Unit is not available in Excel
                    }
                    extracted_data.append(data_point)

        logger.info(f"Successfully extracted {len(extracted_data)} data points from Excel for result ID: {result_id}")
        return extracted_data

    except ValidationException as e:
        logger.error(f"Validation error processing Excel result file: {str(e)}")
        raise
    except Exception as e:
        logger.exception(f"Unexpected error processing Excel result file: {str(e)}")
        raise


@app.task(name='result_tasks.update_result_status', queue='results')
def update_result_status(result_id: UUID, status: str, notes: Optional[str] = None) -> Dict[str, Any]:
    """
    Updates the status of a result.

    Args:
        result_id: ID of the result
        status: New status
        notes: Optional notes

    Returns:
        Dict[str, Any]: Updated result data
    """
    logger.info(f"Attempting to update result status for result ID: {result_id} to status: {status}")
    try:
        # Create ResultUpdate object
        result_update = ResultUpdate(status=status, notes=notes)

        # Update result status
        updated_result = result_service.update_result(result_id, result_update)

        logger.info(f"Successfully updated result status for result ID: {result_id} to status: {status}")
        return {"status": "success", "result": updated_result}

    except ResourceNotFoundException as e:
        logger.error(f"Resource not found error updating result status: {str(e)}")
        return {"status": "error", "message": str(e)}
    except ValidationException as e:
        logger.error(f"Validation error updating result status: {str(e)}")
        return {"status": "error", "message": str(e)}
    except Exception as e:
        logger.exception(f"Unexpected error updating result status: {str(e)}")
        return {"status": "error", "message": "An unexpected error occurred"}


@app.task(name='result_tasks.batch_process_result_files', queue='results')
def batch_process_result_files(result_id: UUID, file_ids: List[UUID], cro_user_id: int) -> Dict[str, Any]:
    """
    Processes multiple result files in batch.

    Args:
        result_id: ID of the result
        file_ids: List of file IDs
        cro_user_id: ID of the CRO user

    Returns:
        Dict[str, Any]: Batch processing results with task IDs
    """
    logger.info(f"Attempting to batch process result files for result ID: {result_id}, file IDs: {file_ids}")
    try:
        results = []
        for file_id in file_ids:
            # Use .delay() to create an asynchronous task for each file
            task = process_result_file.delay(result_id, file_id, cro_user_id)
            results.append(task.id)

        logger.info(f"Successfully created batch processing tasks for {len(file_ids)} files")
        return {"status": "success", "file_count": len(file_ids), "task_ids": results}

    except Exception as e:
        logger.exception(f"Unexpected error creating batch processing tasks: {str(e)}")
        return {"status": "error", "message": "An unexpected error occurred"}


@app.task(name='result_tasks.generate_result_report', queue='results')
def generate_result_report(result_id: UUID, report_format: str) -> Dict[str, Any]:
    """
    Generates a summary report for a result.

    Args:
        result_id: ID of the result
        report_format: Report format (e.g., PDF, Excel)

    Returns:
        Dict[str, Any]: Report generation result with file path
    """
    logger.info(f"Attempting to generate result report for result ID: {result_id}, format: {report_format}")
    try:
        # Get result data
        result_data = result_service.get_result_by_id(result_id, include_data=True)

        # Generate report file
        report_file = generate_report(result_data, report_format)

        # Upload report file
        object_name = file_storage_service.upload_result_file(report_file, f"result_{result_id}.{report_format}")

        logger.info(f"Successfully generated result report for result ID: {result_id}, file path: {object_name}")
        return {"status": "success", "file_path": object_name, "format": report_format}

    except ResourceNotFoundException as e:
        logger.error(f"Resource not found error generating result report: {str(e)}")
        return {"status": "error", "message": str(e)}
    except Exception as e:
        logger.exception(f"Unexpected error generating result report: {str(e)}")
        return {"status": "error", "message": "An unexpected error occurred"}


@app.task(name='result_tasks.cleanup_temporary_result_files', queue='results')
def cleanup_temporary_result_files(file_paths: List[str]) -> Dict[str, Any]:
    """
    Cleans up temporary result files after processing.

    Args:
        file_paths: List of file paths to delete

    Returns:
        Dict[str, Any]: Cleanup result with count of files deleted
    """
    logger.info(f"Attempting to cleanup temporary result files: {file_paths}")
    deleted_count = 0
    try:
        for file_path in file_paths:
            try:
                file_storage_service.delete_file(file_path)
                deleted_count += 1
            except Exception as e:
                logger.warning(f"Failed to delete temporary file {file_path}: {str(e)}")

        logger.info(f"Successfully cleaned up {deleted_count} temporary result files")
        return {"status": "success", "deleted_count": deleted_count}

    except Exception as e:
        logger.exception(f"Unexpected error cleaning up temporary result files: {str(e)}")
        return {"status": "error", "message": "An unexpected error occurred"}


def generate_report(result_data: Dict[str, Any], report_format: str) -> io.BytesIO:
    """
    Generates a report file for a result.

    Args:
        result_data: Result data
        report_format: Report format (e.g., PDF, Excel)

    Returns:
        io.BytesIO: Report file as BytesIO object
    """
    # Placeholder for report generation logic
    # In a real implementation, this would generate a report file
    # based on the result data and report format
    # For now, return a dummy file
    report_data = "Dummy report data"
    report_file = io.BytesIO(report_data.encode('utf-8'))
    return report_file