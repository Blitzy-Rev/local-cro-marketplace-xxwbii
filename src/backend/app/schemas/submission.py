"""
Pydantic schema models for submission data validation, serialization, and API responses
in the Molecular Data Management and CRO Integration Platform.

These schemas represent the data structures for experiment submissions to Contract Research
Organizations (CROs), including submission details, quotes, and status tracking.
"""

from typing import Optional, List, Dict, Any  # standard library
from uuid import UUID  # standard library
from datetime import datetime  # standard library

from pydantic import BaseModel, Field, validator  # pydantic version 2.0+

from ..constants import SubmissionStatus  # Internal import for status enumeration
from .experiment import ExperimentResponse  # Internal import for experiment schema
from .user import UserResponse  # Internal import for user schema

# List of valid submission status values for validation
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
        return SubmissionStatus.PENDING.name
    
    try:
        # Convert to uppercase for case-insensitive comparison
        v_upper = v.upper() if isinstance(v, str) else v
        # Check if it's a valid enum value
        if v_upper not in SUBMISSION_STATUS_VALUES:
            raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")
        return v_upper
    except (AttributeError, ValueError):
        raise ValueError(f"Invalid status: {v}. Must be one of {SUBMISSION_STATUS_VALUES}")


@validator('details', pre=True)
def validate_details(cls, v):
    """
    Validates that submission details are properly formatted.
    
    Args:
        cls: The class being validated
        v: The value to validate (list of detail dictionaries)
        
    Returns:
        The validated details list if valid
        
    Raises:
        ValueError: If details are not properly formatted
    """
    if v is None:
        return []
    
    if not isinstance(v, list):
        raise ValueError("Details must be a list")
    
    for detail in v:
        if not isinstance(detail, dict):
            raise ValueError("Each detail must be a dictionary")
        
        if 'detail_name' not in detail or not detail['detail_name']:
            raise ValueError("Each detail must have a non-empty 'detail_name'")
        
        if 'detail_value' not in detail:
            raise ValueError("Each detail must have a 'detail_value'")
    
    return v


class SubmissionDetailBase(BaseModel):
    """
    Base Pydantic model for submission detail data with common fields.
    
    This model defines the basic structure for submission details, including
    the detail name and value fields.
    """
    detail_name: str = Field(..., description="Name of the submission detail")
    detail_value: str = Field(..., description="Value of the submission detail")


class SubmissionDetailCreate(SubmissionDetailBase):
    """
    Pydantic model for creating a new submission detail.
    
    Extends the base model with a reference to the submission it belongs to.
    """
    submission_id: UUID = Field(..., description="ID of the submission this detail belongs to")


class SubmissionDetailRead(SubmissionDetailBase):
    """
    Pydantic model for submission detail data in API responses.
    
    Extends the base model with additional fields needed for API responses,
    such as the detail ID and submission ID.
    """
    id: UUID = Field(..., description="Unique identifier for the detail")
    submission_id: UUID = Field(..., description="ID of the submission this detail belongs to")


class SubmissionBase(BaseModel):
    """
    Base Pydantic model for submission data with common fields.
    
    This model defines the basic structure for submissions, including the
    experiment ID, CRO ID, status, and details fields.
    """
    experiment_id: UUID = Field(..., description="ID of the experiment being submitted")
    cro_id: Optional[UUID] = Field(None, description="ID of the CRO user assigned to this submission")
    status: Optional[str] = Field(SubmissionStatus.PENDING.name, 
                               description="Current status of the submission")
    details: Optional[List[SubmissionDetailCreate]] = Field(default_factory=list, 
                                                        description="List of submission details")
    notes: Optional[str] = Field(None, description="Additional notes or instructions")


class SubmissionCreate(SubmissionBase):
    """
    Pydantic model for creating a new submission.
    
    Inherits all fields from SubmissionBase.
    """
    pass


class SubmissionUpdate(BaseModel):
    """
    Pydantic model for updating an existing submission.
    
    All fields are optional to allow partial updates.
    """
    cro_id: Optional[UUID] = Field(None, description="Updated CRO user ID")
    status: Optional[str] = Field(None, description="Updated submission status")
    details: Optional[List[SubmissionDetailCreate]] = Field(None, 
                                                        description="Updated submission details")
    notes: Optional[str] = Field(None, description="Updated notes")


class SubmissionInDBBase(SubmissionBase):
    """
    Base Pydantic model for submission data as stored in the database.
    
    Extends the base model with additional fields that are present in the database,
    such as ID, timestamps, and pricing information.
    """
    id: UUID = Field(..., description="Unique identifier for the submission")
    submitted_at: datetime = Field(..., description="Timestamp when the submission was created")
    updated_at: Optional[datetime] = Field(None, 
                                        description="Timestamp when the submission was last updated")
    price: Optional[float] = Field(None, description="Price quoted by the CRO in USD")
    turnaround_days: Optional[int] = Field(None, description="Estimated turnaround time in days")


class SubmissionRead(SubmissionInDBBase):
    """
    Pydantic model for submission data in API responses.
    
    Extends the database model with additional fields needed for API responses,
    such as the experiment and CRO details.
    """
    experiment: ExperimentResponse = Field(..., description="Details of the experiment")
    cro: Optional[UserResponse] = Field(None, description="Details of the assigned CRO")
    details: List[SubmissionDetailRead] = Field(default_factory=list, 
                                             description="List of submission details")


class SubmissionDetailedRead(SubmissionInDBBase):
    """
    Pydantic model for detailed submission data in API responses.
    
    Extends the SubmissionRead model with additional fields for more detailed views,
    such as experiment results.
    """
    experiment: ExperimentResponse = Field(..., description="Details of the experiment")
    cro: Optional[UserResponse] = Field(None, description="Details of the assigned CRO")
    details: List[SubmissionDetailRead] = Field(default_factory=list, 
                                             description="List of submission details")
    results: Optional[List[Dict]] = Field(default_factory=list, 
                                       description="List of experimental results if available")


class SubmissionList(BaseModel):
    """
    Pydantic model for paginated list of submissions in API responses.
    
    Used for paginated responses containing multiple submissions.
    Includes the total count, page number, and page size for pagination.
    """
    items: List[SubmissionRead] = Field(..., description="List of submissions")
    total: int = Field(..., description="Total number of submissions matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")
    filters: Optional[Dict[str, Any]] = Field(None, description="Applied filters")


class SubmissionFilter(BaseModel):
    """
    Pydantic model for filtering submissions by various criteria.
    
    This model defines the available filter parameters that can be used
    when querying submissions. All fields are optional to allow flexible filtering.
    """
    experiment_id: Optional[UUID] = Field(None, description="Filter by experiment ID")
    cro_id: Optional[UUID] = Field(None, description="Filter by CRO user ID")
    status: Optional[str] = Field(None, description="Filter by submission status")
    submitted_after: Optional[datetime] = Field(None, 
                                             description="Filter by submission date (after)")
    submitted_before: Optional[datetime] = Field(None, 
                                              description="Filter by submission date (before)")
    skip: Optional[int] = Field(0, description="Number of submissions to skip (pagination)")
    limit: Optional[int] = Field(100, description="Maximum number of submissions to return")
    sort_by: Optional[str] = Field("submitted_at", 
                                 description="Field to sort by")
    sort_desc: Optional[bool] = Field(True, 
                                   description="Sort in descending order if true")
    
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


class QuoteProvide(BaseModel):
    """
    Pydantic model for CRO to provide a quote for a submission.
    
    This model is used when a CRO responds to a submission request with pricing
    and estimated turnaround time.
    """
    price: float = Field(..., description="Quoted price in USD")
    turnaround_days: int = Field(..., description="Estimated turnaround time in days")
    notes: Optional[str] = Field(None, description="Additional notes or conditions")


class QuoteResponse(BaseModel):
    """
    Pydantic model for pharma user to respond to a quote.
    
    This model is used when a pharma user either approves or rejects a quote
    provided by a CRO.
    """
    approved: bool = Field(..., description="Whether the quote is approved")
    notes: Optional[str] = Field(None, description="Additional notes or explanation")


class SubmissionStatusUpdate(BaseModel):
    """
    Pydantic model for updating the status of a submission.
    
    This model is used to change the status of a submission (e.g., from
    IN_PROGRESS to COMPLETED).
    """
    status: str = Field(..., description="New submission status")
    notes: Optional[str] = Field(None, description="Additional notes about the status change")