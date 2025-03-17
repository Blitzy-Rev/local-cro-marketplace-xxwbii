"""Implements FastAPI endpoints for CRO submission management in the Molecular Data Management and CRO Integration Platform. 
This module provides API routes for creating, retrieving, updating, and managing submissions to Contract Research Organizations (CROs), 
including quote handling and status transitions."""

from typing import List, Dict, Any, Optional  # Type annotations for function parameters and return values
from uuid import UUID  # Type for submission and experiment IDs

from fastapi import APIRouter, Depends, Path, Query, Body, HTTPException, status  # FastAPI utilities for dependency injection, HTTP responses, and parameter validation

from .deps import get_db_session, get_current_pharma_user, get_current_cro_user, get_current_admin_user  # Database session dependency for API endpoints
from ..schemas.submissions import (  # Schemas for submission data
    SubmissionCreateRequest, SubmissionUpdateRequest, SubmissionResponse,
    SubmissionDetailResponse, SubmissionListResponse, SubmissionFilterParams,
    QuoteProvideRequest, QuoteResponseRequest, SubmissionStatusUpdateRequest
)
from ...services import submission_service  # Service layer for submission business logic
from ...exceptions import ValidationException, ResourceNotFoundException  # Exception classes

# Create an APIRouter instance for submission-related endpoints
submissions_router = APIRouter(prefix="/submissions", tags=["submissions"])

@submissions_router.get("/", response_model=SubmissionListResponse)
def get_submissions(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    filters: SubmissionFilterParams = Depends()
) -> SubmissionListResponse:
    """Get a paginated list of submissions with optional filtering"""
    user_id = current_user.get("user_id")  # Extract user_id from current_user
    if current_user.get("role") != "admin":
        filters.created_by = user_id  # Limit results to user's submissions
    try:
        submission_list_response = submission_service.get_submissions(filters=filters, skip=filters.skip, limit=filters.limit)  # Get submissions using submission_service
        return submission_list_response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.get("/{submission_id}", response_model=SubmissionDetailResponse)
def get_submission(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    submission_id: UUID = Path(..., description="The ID of the submission to retrieve")
) -> SubmissionDetailResponse:
    """Get detailed information about a specific submission"""
    try:
        db_submission = submission_service.get_submission_by_id(submission_id, include_results=True)  # Get submission using submission_service
        if not db_submission:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {submission_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the creator of the experiment, the assigned CRO, or an admin
        user_id = current_user.get("user_id")
        user_role = current_user.get("role")
        experiment_creator_id = db_submission["experiment"]["created_by"]
        cro_id = db_submission["cro"]["id"]

        if user_id != experiment_creator_id and user_id != cro_id and user_role != "admin":
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to view this submission"
            )  # Raise HTTPException if not authorized

        return db_submission  # Return the submission detail response
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.post("/", response_model=SubmissionResponse, status_code=status.HTTP_201_CREATED)
def create_submission(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    submission_data: SubmissionCreateRequest = Body(...)
) -> SubmissionResponse:
    """Create a new submission for an experiment to a CRO"""
    user_id = current_user.get("user_id")  # Extract user_id from current_user
    try:
        created_submission = submission_service.create_submission(submission_data=submission_data, user_id=user_id)  # Create submission using submission_service
        return created_submission  # Return the created submission response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.put("/{submission_id}", response_model=SubmissionResponse)
def update_submission(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    submission_id: UUID = Path(..., description="The ID of the submission to update"),
    submission_data: SubmissionUpdateRequest = Body(...)
) -> SubmissionResponse:
    """Update an existing submission"""
    try:
        user_id = current_user.get("user_id")
        updated_submission = submission_service.update_submission(submission_id=submission_id, submission_data=submission_data, user_id=user_id)  # Update submission using submission_service
        return updated_submission  # Return the updated submission response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.put("/{submission_id}/status", response_model=SubmissionResponse)
def update_submission_status(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    submission_id: UUID = Path(..., description="The ID of the submission to update"),
    status_data: SubmissionStatusUpdateRequest = Body(...)
) -> SubmissionResponse:
    """Update the status of a submission"""
    try:
        user_id = current_user.get("user_id")
        db_submission = submission_service.get_submission_by_id(submission_id)  # Get submission using submission_service
        if not db_submission:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {submission_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the creator of the experiment, the assigned CRO, or an admin
        user_role = current_user.get("role")
        experiment_creator_id = db_submission["experiment"]["created_by"]
        cro_id = db_submission["cro"]["id"]

        if user_id != experiment_creator_id and user_id != cro_id and user_role != "admin":
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to update this submission"
            )  # Raise HTTPException if not authorized

        updated_submission = submission_service.update_submission_status(submission_id=submission_id, status=status_data.status, user_id=user_id)  # Update submission status using submission_service
        return updated_submission  # Return the updated submission response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.post("/{submission_id}/quote", response_model=SubmissionResponse)
def provide_quote(
    current_user: Dict[str, Any] = Depends(get_current_cro_user),
    submission_id: UUID = Path(..., description="The ID of the submission"),
    quote_data: QuoteProvideRequest = Body(...)
) -> SubmissionResponse:
    """Allows a CRO to provide a quote for a submission"""
    user_id = current_user.get("user_id")  # Extract user_id from current_user
    try:
        db_submission = submission_service.get_submission_by_id(submission_id)  # Get submission using submission_service
        if not db_submission:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {submission_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the assigned CRO for the submission
        cro_id = db_submission["cro"]["id"]
        if user_id != cro_id:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to provide a quote for this submission"
            )  # Raise HTTPException if not authorized

        updated_submission = submission_service.provide_quote(submission_id=submission_id, quote_data=quote_data, cro_user_id=user_id)  # Provide quote using submission_service
        return updated_submission  # Return the updated submission response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.post("/{submission_id}/quote/response", response_model=SubmissionResponse)
def respond_to_quote(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    submission_id: UUID = Path(..., description="The ID of the submission"),
    response_data: QuoteResponseRequest = Body(...)
) -> SubmissionResponse:
    """Allows a pharma user to approve or reject a quote"""
    user_id = current_user.get("user_id")  # Extract user_id from current_user
    try:
        db_submission = submission_service.get_submission_by_id(submission_id)  # Get submission using submission_service
        if not db_submission:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {submission_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the creator of the experiment or an admin
        user_role = current_user.get("role")
        experiment_creator_id = db_submission["experiment"]["created_by"]

        if user_id != experiment_creator_id and user_role != "admin":
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to respond to the quote for this submission"
            )  # Raise HTTPException if not authorized

        updated_submission = submission_service.respond_to_quote(submission_id=submission_id, response_data=response_data, user_id=user_id)  # Respond to quote using submission_service
        return updated_submission  # Return the updated submission response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.get("/experiment/{experiment_id}", response_model=SubmissionListResponse)
def get_submissions_by_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description="The ID of the experiment"),
    skip: int = Query(0, ge=0, description="Skip N records"),
    limit: int = Query(100, ge=1, le=1000, description="Limit to N records")
) -> SubmissionListResponse:
    """Get submissions for a specific experiment"""
    try:
        user_id = current_user.get("user_id")
        db_experiment = submission_service.get_submission_by_id(experiment_id)  # Get submission using submission_service
        if not db_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {experiment_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the creator of the experiment or an admin
        user_role = current_user.get("role")
        experiment_creator_id = db_experiment["experiment"]["created_by"]

        if user_id != experiment_creator_id and user_role != "admin":
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to view submissions for this experiment"
            )  # Raise HTTPException if not authorized

        submission_list_response = submission_service.get_submissions_by_experiment(experiment_id=experiment_id, user_id=user_id, skip=skip, limit=limit)  # Get submissions using submission_service
        return submission_list_response  # Return the submission list response
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

@submissions_router.get("/cro", response_model=SubmissionListResponse)
def get_submissions_by_cro(
    current_user: Dict[str, Any] = Depends(get_current_cro_user),
    skip: int = Query(0, ge=0, description="Skip N records"),
    limit: int = Query(100, ge=1, le=1000, description="Limit to N records")
) -> SubmissionListResponse:
    """Get submissions assigned to a specific CRO"""
    user_id = current_user.get("user_id")  # Extract user_id from current_user
    try:
        submission_list_response = submission_service.get_submissions_by_cro(cro_user_id=user_id, skip=skip, limit=limit)  # Get submissions using submission_service
        return submission_list_response  # Return the submission list response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException

@submissions_router.get("/status/{status}", response_model=SubmissionListResponse)
def get_submissions_by_status(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    status: str = Path(..., description="The submission status"),
    skip: int = Query(0, ge=0, description="Skip N records"),
    limit: int = Query(100, ge=1, le=1000, description="Limit to N records")
) -> SubmissionListResponse:
    """Get submissions with a specific status"""
    try:
        user_id = current_user.get("user_id")
        user_role = current_user.get("role")
        submission_list_response = submission_service.get_submissions_by_status(status=status, skip=skip, limit=limit)  # Get submissions using submission_service
        if user_role != "admin":
            # If user is not admin, filter results to only show submissions for experiments created by the user
            submission_list_response.items = [s for s in submission_list_response.items if s.experiment.created_by == user_id]
            submission_list_response.total = len(submission_list_response.items)
        return submission_list_response  # Return the submission list response
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))  # Convert ValidationException to HTTPException

@submissions_router.get("/{submission_id}/results", response_model=Dict[str, Any])
def get_submission_results(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    submission_id: UUID = Path(..., description="The ID of the submission"),
    skip: int = Query(0, ge=0, description="Skip N records"),
    limit: int = Query(100, ge=1, le=1000, description="Limit to N records")
) -> Dict[str, Any]:
    """Get results for a specific submission"""
    try:
        db_submission = submission_service.get_submission_by_id(submission_id)  # Get submission using submission_service
        if not db_submission:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"Submission with id {submission_id} not found")  # Raise HTTPException if submission not found

        # Check if current user is the creator of the experiment, the assigned CRO, or an admin
        user_id = current_user.get("user_id")
        user_role = current_user.get("role")
        experiment_creator_id = db_submission["experiment"]["created_by"]
        cro_id = db_submission["cro"]["id"]

        if user_id != experiment_creator_id and user_id != cro_id and user_role != "admin":
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="User not authorized to view results for this submission"
            )  # Raise HTTPException if not authorized

        results_list_response = submission_service.get_submission_results(submission_id=submission_id, skip=skip, limit=limit)  # Get results using submission_service
        return results_list_response  # Return the results list response
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))  # Convert ResourceNotFoundException to HTTPException

def handle_validation_exception(exc: ValidationException):
    """Utility function to handle ValidationException and convert to HTTPException"""
    raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=exc.message)

def handle_resource_not_found_exception(exc: ResourceNotFoundException):
    """Utility function to handle ResourceNotFoundException and convert to HTTPException"""
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=exc.message)