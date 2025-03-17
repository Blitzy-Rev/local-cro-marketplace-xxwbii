import logging
from typing import List, Dict, Any, Optional
from uuid import UUID

import io  # Import the io module
from fastapi import APIRouter, Depends, HTTPException, status, UploadFile, File, Form, Query, Path, Body
from fastapi.responses import StreamingResponse  # Import StreamingResponse

from ..deps import get_db_session, get_current_cro_user, get_current_pharma_user, get_current_admin_user
from ..schemas.results import (
    ResultDetailRead, ResultUpload, ResultStatusUpdate, ResultFileUpload, ResultDataPoint,
    ResultAnalysis, ResultExport, ResultBatchApproval
)
from ...schemas.result import ResultRead, ResultDetailedRead, ResultList, ResultFilter, ResultApproval
from ...schemas.msg import Msg
from ...services import result_service
from ...exceptions import ValidationException, ResourceNotFoundException

# Initialize logger
logger = logging.getLogger(__name__)

# Define API router for results
router = APIRouter(prefix='/results', tags=['results'])


@router.get('/', response_model=ResultList)
def get_results(
    filters: ResultFilter = Depends(),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Get a list of results with optional filtering
    """
    logger.info(f"Getting results with filters: {filters}")
    
    # Add user_id to filters if current user is a pharma user
    # filters.user_id = current_user.get("user_id")
    
    # Call result_service.get_results with the filters
    results = result_service.get_results(filters)
    
    # Return the results as a ResultList object
    return results


@router.get('/submission/{submission_id}', response_model=ResultList)
def get_results_by_submission(
    submission_id: UUID = Path(..., description="The ID of the submission"),
    skip: int = Query(0, ge=0, description="Skip N results"),
    limit: int = Query(100, ge=1, le=100, description="Limit to N results"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Get results for a specific submission
    """
    logger.info(f"Getting results for submission: {submission_id}")
    
    # Call result_service.get_results_by_submission with submission_id, skip, and limit
    results = result_service.get_results_by_submission(submission_id, skip, limit)
    
    # Return the results as a ResultList object
    return results


@router.get('/{result_id}', response_model=ResultDetailRead)
def get_result(
    result_id: UUID = Path(..., description="The ID of the result"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Get a specific result by ID
    """
    logger.info(f"Getting result by ID: {result_id}")
    
    # Call result_service.get_result_by_id with result_id
    result = result_service.get_result_by_id(result_id)
    
    # If result not found, raise HTTPException with 404 status
    if not result:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result not found"
        )
    
    # Return the result as a ResultDetailRead object
    return result


@router.post('/', response_model=ResultDetailRead, status_code=status.HTTP_201_CREATED)
def create_result(
    result_data: ResultUpload = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_cro_user)
):
    """
    Create a new result for a submission
    """
    logger.info(f"Creating new result for submission: {result_data.submission_id}")
    
    # Extract CRO user ID from current_user
    cro_user_id = current_user.get("user_id")
    
    # Try to call result_service.create_result with result_data and cro_user_id
    try:
        result = result_service.create_result(result_data, cro_user_id)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the created result as a ResultDetailRead object
    return result


@router.put('/{result_id}', response_model=ResultDetailRead)
def update_result(
    result_id: UUID = Path(..., description="The ID of the result to update"),
    result_data: ResultStatusUpdate = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_cro_user)
):
    """
    Update an existing result
    """
    logger.info(f"Updating result with ID: {result_id}")
    
    # Try to call result_service.update_result with result_id and result_data
    try:
        result = result_service.update_result(result_id, result_data)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the updated result as a ResultDetailRead object
    return result


@router.post('/{result_id}/approve', response_model=ResultDetailRead)
def approve_reject_result(
    result_id: UUID = Path(..., description="The ID of the result to approve/reject"),
    approval_data: ResultApproval = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Approve or reject a result
    """
    logger.info(f"Approving/Rejecting result with ID: {result_id}")
    
    # Extract user ID from current_user
    user_id = current_user.get("user_id")
    
    # Try to call result_service.approve_reject_result with result_id, approval_data, and user_id
    try:
        result = result_service.approve_reject_result(result_id, approval_data, user_id)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the updated result as a ResultDetailRead object
    return result

@router.post('/batch-approve', response_model=Msg)
def batch_approve_reject_results(
    batch_data: ResultBatchApproval = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Approve or reject multiple results in a batch
    """
    logger.info(f"Batch approving/rejecting results")
    
    # Extract user ID from current_user
    user_id = current_user.get("user_id")
    
    success_count = 0
    failed_count = 0
    
    # For each result_id in batch_data.result_ids:
    for result_id in batch_data.result_ids:
        # Create ResultApproval object from batch_data
        approval_data = ResultApproval(approved=batch_data.approved, notes=batch_data.notes)
        
        # Try to call result_service.approve_reject_result
        try:
            result_service.approve_reject_result(result_id, approval_data, user_id)
            success_count += 1
        except Exception as e:
            logger.error(f"Failed to approve/reject result {result_id}: {e}")
            failed_count += 1
    
    # Return message with counts of successful and failed operations
    return Msg(msg=f"Successfully processed {success_count} results. Failed to process {failed_count} results.")


@router.post('/{result_id}/files', status_code=status.HTTP_201_CREATED)
def upload_result_file(
    result_id: UUID = Path(..., description="The ID of the result"),
    file: UploadFile = File(..., description="The file to upload"),
    current_user: Dict[str, Any] = Depends(get_current_cro_user)
):
    """
    Upload a file for a result
    """
    logger.info(f"Uploading result file for result with ID: {result_id}")
    
    # Extract CRO user ID from current_user
    cro_user_id = current_user.get("user_id")
    
    # Try to call result_service.upload_result_file with result_id, file, and cro_user_id
    try:
        file_metadata = result_service.upload_result_file(result_id, file, cro_user_id)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the file metadata
    return file_metadata


@router.get('/files/{file_id}')
def get_result_file(
    file_id: UUID = Path(..., description="The ID of the file"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Download a file associated with a result
    """
    logger.info(f"Downloading result file with ID: {file_id}")
    
    # Try to call result_service.get_result_file with file_id
    try:
        file_data, content_type, filename = result_service.get_result_file(file_id)
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Create and return StreamingResponse with appropriate headers
    headers = {
        'Content-Disposition': f'attachment; filename="{filename}"',
        'Content-Type': content_type
    }
    return StreamingResponse(file_data, media_type=content_type, headers=headers)


@router.delete('/files/{file_id}', response_model=Msg)
def delete_result_file(
    file_id: UUID = Path(..., description="The ID of the file"),
    current_user: Dict[str, Any] = Depends(get_current_cro_user)
):
    """
    Delete a file associated with a result
    """
    logger.info(f"Deleting result file with ID: {file_id}")
    
    # Extract CRO user ID from current_user
    cro_user_id = current_user.get("user_id")
    
    # Try to call result_service.delete_result_file with file_id and cro_user_id
    try:
        result_service.delete_result_file(file_id, cro_user_id)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return success message
    return Msg(msg="File deleted successfully")


@router.post('/{result_id}/data', status_code=status.HTTP_201_CREATED)
def add_result_data(
    result_id: UUID = Path(..., description="The ID of the result"),
    data_point: ResultDataPoint = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_cro_user)
):
    """
    Add structured data to a result for a specific molecule
    """
    logger.info(f"Adding result data for result with ID: {result_id}")
    
    # Extract CRO user ID from current_user
    cro_user_id = current_user.get("user_id")
    
    # Try to call result_service.add_result_data with result_id, data_point parameters, and cro_user_id
    try:
        data_record = result_service.add_result_data(result_id, data_point, cro_user_id)
    except ValidationException as e:
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=e.message
        )
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the created data record
    return data_record


@router.get('/{result_id}/data')
def get_result_data(
    result_id: UUID = Path(..., description="The ID of the result"),
    molecule_id: Optional[UUID] = Query(None, description="Filter by molecule ID"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Get structured data for a result, optionally filtered by molecule
    """
    logger.info(f"Getting result data for result with ID: {result_id}")
    
    # Try to call result_service.get_result_data with result_id and molecule_id
    try:
        data_records = result_service.get_result_data(result_id, molecule_id)
    except ResourceNotFoundException as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=e.message
        )
    
    # Return the list of data records
    return data_records

@router.post('/export')
def export_result(
    export_data: ResultExport = Body(...),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Export result data in various formats
    """
    logger.info(f"Exporting result data for result with ID: {export_data.result_id}")

    # Get result data from result_service.get_result_by_id
    result = result_service.get_result_by_id(export_data.result_id)
    if not result:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result not found"
        )

    # Format data according to requested export_format (csv, xlsx, json)
    # Create appropriate content type based on format
    # Generate filename with result ID and format extension
    # Return StreamingResponse with appropriate headers
    return StreamingResponse(content="Exported data", media_type="text/csv", headers={"Content-Disposition": "attachment;filename=result.csv"})

@router.get('/{result_id}/analyze', response_model=ResultAnalysis)
def analyze_result(
    result_id: UUID = Path(..., description="The ID of the result"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
):
    """
    Generate analysis of result data
    """
    logger.info(f"Analyzing result data for result with ID: {result_id}")

    # Get result data from result_service.get_result_by_id
    result = result_service.get_result_by_id(result_id)
    if not result:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="Result not found"
        )

    # Get structured data from result_service.get_result_data
    # Calculate summary statistics for each data type
    # Organize data by molecule
    # Generate visualization metadata if requested
    # Return ResultAnalysis object with analysis results
    return ResultAnalysis(result_id=result_id, summary_statistics={}, molecule_data={})