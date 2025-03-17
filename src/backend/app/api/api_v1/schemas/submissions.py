"""
Defines Pydantic schema models for submission-related API requests and responses in the API v1 namespace.
These schemas are used for data validation, serialization, and documentation in the CRO submission and integration endpoints.
"""

from typing import Optional, List, Dict, Any  # standard library
from uuid import UUID  # standard library
from datetime import datetime  # standard library
from typing import TYPE_CHECKING  # standard library

from pydantic import BaseModel, Field, validator, root_validator  # pydantic version 2.0+

# Import base submission schemas for API schemas
from ....schemas.submission import (
    SubmissionBase, SubmissionCreate, SubmissionUpdate, SubmissionInDBBase, 
    SubmissionRead, SubmissionDetailedRead, SubmissionList, SubmissionFilter
)
from ....models.submission import SubmissionStatus  # Import submission status enumeration
from .experiments import ExperimentResponse  # Import experiment schema for submission relationships
from .users import UserProfileResponse  # Import user schema for CRO information

# Create a list of valid submission status values for validation
SUBMISSION_STATUS_VALUES = [status.name for status in SubmissionStatus]


@validator('status', pre=True)
def validate_status(cls, v):
    """
    Validates that a submission status is one of the allowed values.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        The validated status if valid
        
    Raises:
        ValueError: If the status is not a valid SubmissionStatus value
    """
    if v is None:
        return 'PENDING'  # Default status
    
    try:
        # Convert to uppercase for case-insensitive comparison
        v_upper = v.upper() if isinstance(v, str) else v
        # Check if it's a valid enum value
        if v_upper not in SUBMISSION_STATUS_VALUES:
            raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")
        return v_upper
    except (AttributeError, ValueError):
        raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")


@root_validator(pre=True)
def validate_details(cls, values):
    """
    Validates that submission details are of the correct type and format.
    
    Args:
        cls: The class being validated
        values: The values dictionary
        
    Returns:
        The validated values dictionary
        
    Raises:
        ValueError: If details are invalid
    """
    details = values.get('details')
    if details is None:
        return values
    
    if not isinstance(details, dict):
        raise ValueError("Details must be a dictionary")
    
    # Validate detail keys and values
    for detail_name, detail_value in details.items():
        if not detail_name or not isinstance(detail_name, str):
            raise ValueError(f"Detail name must be a non-empty string: {detail_name}")
        
        # Check if detail name is a valid identifier
        if not detail_name.isidentifier():
            raise ValueError(f"Detail name must be a valid identifier: {detail_name}")
    
    return values


class SubmissionCreateRequest(BaseModel):
    """Schema for creating a new submission through the API."""
    experiment_id: UUID = Field(..., description="ID of the experiment being submitted")
    cro_id: UUID = Field(..., description="ID of the CRO user assigned to this submission")
    status: Optional[str] = Field(None, description="Initial status of the submission")
    details: Optional[Dict[str, str]] = Field(None, description="Additional details about the submission")


class SubmissionUpdateRequest(BaseModel):
    """Schema for updating an existing submission through the API."""
    status: Optional[str] = Field(None, description="Updated status of the submission")
    details: Optional[Dict[str, str]] = Field(None, description="Updated submission details")


class SubmissionResponse(BaseModel):
    """Schema for submission data in API responses."""
    id: UUID = Field(..., description="Unique identifier for the submission")
    experiment_id: UUID = Field(..., description="ID of the experiment being submitted")
    cro_id: UUID = Field(..., description="ID of the CRO user assigned to this submission")
    status: str = Field(..., description="Current status of the submission")
    submitted_at: datetime = Field(..., description="Timestamp when the submission was created")
    updated_at: Optional[datetime] = Field(None, description="Timestamp when the submission was last updated")
    experiment: ExperimentResponse = Field(..., description="Details of the experiment")
    cro: UserProfileResponse = Field(..., description="Details of the assigned CRO")
    details: Dict[str, str] = Field(default_factory=dict, description="Additional details about the submission")


class SubmissionDetailResponse(BaseModel):
    """Schema for detailed submission data in API responses."""
    id: UUID = Field(..., description="Unique identifier for the submission")
    experiment_id: UUID = Field(..., description="ID of the experiment being submitted")
    cro_id: UUID = Field(..., description="ID of the CRO user assigned to this submission")
    status: str = Field(..., description="Current status of the submission")
    submitted_at: datetime = Field(..., description="Timestamp when the submission was created")
    updated_at: Optional[datetime] = Field(None, description="Timestamp when the submission was last updated")
    experiment: ExperimentResponse = Field(..., description="Details of the experiment")
    cro: UserProfileResponse = Field(..., description="Details of the assigned CRO")
    details: Dict[str, str] = Field(default_factory=dict, description="Additional details about the submission")
    results: List["ResultRead"] = Field(default_factory=list, description="List of experimental results if available")


class SubmissionListResponse(BaseModel):
    """Schema for paginated list of submissions in API responses."""
    items: List[SubmissionResponse] = Field(..., description="List of submissions")
    total: int = Field(..., description="Total number of submissions matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")
    filters: Optional[Dict[str, Any]] = Field(None, description="Applied filters")


class SubmissionFilterParams(BaseModel):
    """Schema for filtering submissions in API requests."""
    experiment_id: Optional[UUID] = Field(None, description="Filter by experiment ID")
    cro_id: Optional[UUID] = Field(None, description="Filter by CRO user ID")
    status: Optional[str] = Field(None, description="Filter by submission status")
    submitted_after: Optional[datetime] = Field(None, description="Filter by submission date (after)")
    submitted_before: Optional[datetime] = Field(None, description="Filter by submission date (before)")
    skip: Optional[int] = Field(0, description="Number of submissions to skip (pagination)")
    limit: Optional[int] = Field(100, description="Maximum number of submissions to return")
    sort_by: Optional[str] = Field("submitted_at", description="Field to sort by")
    sort_desc: Optional[bool] = Field(True, description="Sort in descending order if true")
    
    @validator('status', pre=True)
    def validate_status(cls, v):
        """
        Validates that the status filter is one of the allowed values.
        
        Args:
            cls: The class being validated
            v: The value to validate
            
        Returns:
            The validated status if valid
            
        Raises:
            ValueError: If the status is not a valid SubmissionStatus value
        """
        if v is None:
            return None
        
        try:
            # Convert to uppercase for case-insensitive comparison
            v_upper = v.upper() if isinstance(v, str) else v
            # Check if it's a valid enum value
            if v_upper not in SUBMISSION_STATUS_VALUES:
                raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")
            return v_upper
        except (AttributeError, ValueError):
            raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")


class QuoteProvideRequest(BaseModel):
    """Schema for providing a quote for a submission."""
    price: float = Field(..., description="Quoted price in USD")
    currency: Optional[str] = Field("USD", description="Currency of the quoted price")
    turnaround_time: int = Field(..., description="Estimated turnaround time in days")
    notes: Optional[str] = Field(None, description="Additional notes or conditions")


class QuoteResponseRequest(BaseModel):
    """Schema for responding to a quote."""
    approved: bool = Field(..., description="Whether the quote is approved")
    notes: Optional[str] = Field(None, description="Additional notes or explanation")


class SubmissionStatusUpdateRequest(BaseModel):
    """Schema for updating the status of a submission."""
    status: str = Field(..., description="New submission status")
    notes: Optional[str] = Field(None, description="Additional notes about the status change")
    
    @validator('status', pre=True)
    def validate_status(cls, v):
        """
        Validates that the status is one of the allowed values.
        
        Args:
            cls: The class being validated
            v: The value to validate
            
        Returns:
            The validated status if valid
            
        Raises:
            ValueError: If the status is not a valid SubmissionStatus value
        """
        if v not in SUBMISSION_STATUS_VALUES:
            raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")
        return v