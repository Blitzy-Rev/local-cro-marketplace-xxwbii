# src/backend/app/api/api_v1/endpoints/csv.py
"""
FastAPI endpoint implementation for CSV file upload, mapping, and processing in the Molecular Data Management and CRO Integration Platform. This module provides API routes for handling CSV files containing molecular data, including file upload, header mapping, processing status tracking, and system property information.
"""

import os  # standard library
import typing  # standard library
from typing import Dict, List, Any  # standard library

from fastapi import APIRouter  # fastapi 0.95+
from fastapi import Depends, HTTPException, status, UploadFile, File, Form, BackgroundTasks  # fastapi 0.95+
from sqlalchemy.orm import Session  # sqlalchemy.orm 2.0+

from ..deps import get_current_pharma_user, get_db_session  # Dependency to ensure user has pharma role
from ..schemas.csv import CSVUploadResponse, CSVMappingRequest, CSVMappingResponse, CSVStatusResponse, AvailablePropertiesResponse, SystemProperty, CSV_PROCESSING_STATUS  # Response schema for CSV upload endpoint
from ..schemas.msg import Msg  # Generic message response schema
from ...services.csv_service import CSVService  # Service for CSV file processing operations
from ...services.file_storage_service import FileStorageService  # Service for file storage operations
from ...worker.tasks.csv_tasks import process_csv, cancel_csv_processing, generate_csv_import_report  # Background task for processing CSV files
from ...exceptions import CSVProcessingException, ValidationException  # Custom exception for CSV processing errors
from ...logging_config import logger  # Application logger for error and info logging

# Initialize FastAPI router
router = APIRouter(prefix='/csv', tags=['csv'])

# Initialize services
file_storage_service = FileStorageService()
csv_service = CSVService(file_storage_service)

@router.post('/upload', response_model=CSVUploadResponse, status_code=status.HTTP_200_OK)
def upload_csv(
    file: UploadFile = File(...),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> CSVUploadResponse:
    """
    Upload a CSV file containing molecular data
    """
    try:
        # Log CSV upload attempt with filename and user ID
        logger.info(f"CSV upload requested for file: {file.filename} by user: {current_user['user_id']}")

        # Validate file is not empty
        if file.file.tell() == 0:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Uploaded file is empty"
            )

        # Save CSV file using file_storage_service.save_csv_file()
        file_id = file_storage_service.save_csv_file(file.file, file.filename)

        # Validate CSV file structure using csv_service.validate_csv_file()
        validation_results = csv_service.validate_csv_file(file_id)

        # Get CSV headers using csv_service.get_csv_headers()
        headers = csv_service.get_csv_headers(file_id)

        # Return CSVUploadResponse with file_id, headers, and row_count
        return CSVUploadResponse(
            file_id=file_id,
            headers=headers,
            row_count=validation_results['row_count'],
            message="File uploaded successfully"
        )
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"CSV upload failed for file: {file.filename} with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CSV upload failed: {str(e)}"
        )

@router.post('/map', response_model=CSVMappingResponse, status_code=status.HTTP_202_ACCEPTED)
def map_csv_columns(
    mapping_request: CSVMappingRequest,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db_session),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> CSVMappingResponse:
    """
    Map CSV columns to system properties and start processing
    """
    try:
        # Log CSV mapping request with file_id and user ID
        logger.info(f"CSV mapping requested for file_id: {mapping_request.file_id} by user: {current_user['user_id']}")

        # Get CSV headers using csv_service.get_csv_headers()
        headers = csv_service.get_csv_headers(mapping_request.file_id)

        # Validate column mapping using csv_service.validate_column_mapping()
        csv_service.validate_column_mapping(headers, mapping_request.mapping)

        # Create processing job using csv_service.create_processing_job()
        job_id = csv_service.create_processing_job(mapping_request.file_id, mapping_request.mapping, current_user['user_id'])

        # Add process_csv task to background_tasks with job_id, file_id, mapping, and user_id
        background_tasks.add_task(process_csv, job_id, mapping_request.file_id, mapping_request.mapping, current_user['user_id'])

        # Return CSVMappingResponse with job_id and status
        return CSVMappingResponse(
            job_id=job_id,
            status=CSV_PROCESSING_STATUS.PENDING.name,
            message="CSV mapping submitted for processing"
        )
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"CSV mapping failed for file_id: {mapping_request.file_id} with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CSV mapping failed: {str(e)}"
        )

@router.get('/status/{job_id}', response_model=CSVStatusResponse, status_code=status.HTTP_200_OK)
def get_processing_status(
    job_id: str,
    db: Session = Depends(get_db_session),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> CSVStatusResponse:
    """
    Get the status of a CSV processing job
    """
    try:
        # Log status request with job_id and user ID
        logger.info(f"CSV status requested for job_id: {job_id} by user: {current_user['user_id']}")

        # Get job status using csv_service.get_job_status()
        job_status = csv_service.get_job_status(job_id, current_user['user_id'])

        # Return CSVStatusResponse with job status information
        return CSVStatusResponse(
            job_id=job_id,
            status=job_status['status'],
            progress=job_status.get('progress'),
            summary=job_status.get('summary')
        )
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"CSV status retrieval failed for job_id: {job_id} with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CSV status retrieval failed: {str(e)}"
        )

@router.post('/cancel/{job_id}', response_model=Msg, status_code=status.HTTP_200_OK)
def cancel_processing(
    job_id: str,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db_session),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> Msg:
    """
    Cancel an ongoing CSV processing job
    """
    try:
        # Log cancellation request with job_id and user ID
        logger.info(f"CSV cancellation requested for job_id: {job_id} by user: {current_user['user_id']}")

        # Add cancel_csv_processing task to background_tasks with job_id and user_id
        background_tasks.add_task(cancel_csv_processing, job_id, current_user['user_id'])

        # Return success message
        return Msg(msg="CSV processing cancellation requested")
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"CSV cancellation failed for job_id: {job_id} with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CSV cancellation failed: {str(e)}"
        )

@router.get('/properties', response_model=AvailablePropertiesResponse, status_code=status.HTTP_200_OK)
def get_available_properties(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> AvailablePropertiesResponse:
    """
    Get list of available system properties for CSV mapping
    """
    try:
        # Log properties request with user ID
        logger.info(f"Available properties requested by user: {current_user['user_id']}")

        # Get available properties using csv_service.get_available_properties()
        available_properties = csv_service.get_available_properties()

        # Convert properties to SystemProperty objects
        properties = [SystemProperty(**prop) for prop in available_properties]

        # Return AvailablePropertiesResponse with properties list
        return AvailablePropertiesResponse(properties=properties)
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"Available properties retrieval failed with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Available properties retrieval failed: {str(e)}"
        )

@router.post('/report/{job_id}', response_model=Msg, status_code=status.HTTP_202_ACCEPTED)
def generate_import_report(
    job_id: str,
    background_tasks: BackgroundTasks,
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> Msg:
    """
    Generate a detailed report for a CSV import job
    """
    try:
        # Log report generation request with job_id and user ID
        logger.info(f"Report generation requested for job_id: {job_id} by user: {current_user['user_id']}")

        # Add generate_csv_import_report task to background_tasks with job_id and user_id
        background_tasks.add_task(generate_csv_import_report, job_id, current_user['user_id'])

        # Return success message with report generation status
        return Msg(msg="CSV import report generation requested")
    except HTTPException as http_exc:
        # Re-raise HTTPExceptions
        raise http_exc
    except Exception as e:
        # Handle exceptions and return appropriate HTTP error responses
        logger.error(f"Report generation failed for job_id: {job_id} with error: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CSV import report generation failed: {str(e)}"
        )