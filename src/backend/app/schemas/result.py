"""
Pydantic schema models for experimental results in the Molecular Data Management and CRO Integration Platform.

This module defines the data models used for results uploaded by CRO users, including result files,
structured data points, and status tracking. These schemas are used for data validation,
serialization, and API responses throughout the application.
"""

# External imports
from pydantic import BaseModel, Field, validator  # pydantic version 2.0+
from typing import Optional, List, Dict, Any  # standard library
from uuid import UUID  # standard library
from datetime import datetime  # standard library

# Internal imports
from ..models.result import Result
from ..models.result_file import ResultFile
from ..models.result_data import ResultData
from .submission import SubmissionRead
from .molecule import MoleculeResponse
from ..constants import ResultStatus
from ..utils.validation_utils import validate_result_file

# List of valid result status values for validation
RESULT_STATUS_VALUES = [status.name for status in ResultStatus]


@validator('status', pre=True)
def validate_status(cls, v):
    """
    Validates that a result status is one of the allowed values.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        The validated status if valid
        
    Raises:
        ValueError: If the status is not a valid ResultStatus value
    """
    if v is None:
        return ResultStatus.PENDING.name
    
    try:
        # Convert to uppercase for case-insensitive comparison
        v_upper = v.upper() if isinstance(v, str) else v
        # Check if it's a valid enum value
        if v_upper not in RESULT_STATUS_VALUES:
            raise ValueError(f"Invalid status: {v}. Must be one of {RESULT_STATUS_VALUES}")
        return v_upper
    except (AttributeError, ValueError):
        raise ValueError(f"Invalid status: {v}. Must be one of {RESULT_STATUS_VALUES}")


@validator('data', pre=True)
def validate_data(cls, v):
    """
    Validates that result data points are properly formatted.
    
    Args:
        cls: The class being validated
        v: The value to validate (list of data point dictionaries)
        
    Returns:
        The validated data list if valid
        
    Raises:
        ValueError: If data points are not properly formatted
    """
    if v is None:
        return []
    
    if not isinstance(v, list):
        raise ValueError("Data points must be a list")
    
    for data_point in v:
        if not isinstance(data_point, dict):
            raise ValueError("Each data point must be a dictionary")
        
        if 'molecule_id' not in data_point or not data_point['molecule_id']:
            raise ValueError("Each data point must have a non-empty 'molecule_id'")
        
        if 'data_name' not in data_point or not data_point['data_name']:
            raise ValueError("Each data point must have a non-empty 'data_name'")
        
        if 'data_value' not in data_point:
            raise ValueError("Each data point must have a 'data_value'")
    
    return v


class ResultFileBase(BaseModel):
    """
    Base Pydantic model for result file data with common fields.
    
    This model defines the basic structure for result files, including file name,
    path, size, and type.
    """
    file_name: str = Field(..., description="Original name of the uploaded file")
    file_path: str = Field(..., description="Path to the file in the object storage system")
    file_size: int = Field(..., description="Size of the file in bytes")
    file_type: str = Field(..., description="MIME type or file extension")


class ResultFileCreate(ResultFileBase):
    """
    Pydantic model for creating a new result file.
    
    Extends the base model with a reference to the result it belongs to.
    """
    result_id: UUID = Field(..., description="ID of the result this file belongs to")


class ResultFileRead(ResultFileBase):
    """
    Pydantic model for result file data in API responses.
    
    Extends the base model with additional fields needed for API responses,
    such as the file ID, result ID, and upload timestamp.
    """
    id: UUID = Field(..., description="Unique identifier for the file")
    result_id: UUID = Field(..., description="ID of the result this file belongs to")
    uploaded_at: datetime = Field(..., description="Timestamp when the file was uploaded")


class ResultDataBase(BaseModel):
    """
    Base Pydantic model for result data point with common fields.
    
    This model defines the basic structure for structured result data points,
    linking specific values to both a result and a molecule.
    """
    molecule_id: UUID = Field(..., description="ID of the molecule this data is associated with")
    data_name: str = Field(..., description="Name or identifier of the data point (e.g., 'IC50', 'Solubility')")
    data_value: float = Field(..., description="Numerical value of the measurement")
    data_unit: Optional[str] = Field(None, description="Unit of measurement (e.g., 'nM', 'mg/mL')")


class ResultDataCreate(ResultDataBase):
    """
    Pydantic model for creating a new result data point.
    
    Extends the base model with a reference to the result it belongs to.
    """
    result_id: UUID = Field(..., description="ID of the result this data belongs to")


class ResultDataRead(ResultDataBase):
    """
    Pydantic model for result data point in API responses.
    
    Extends the base model with additional fields needed for API responses,
    such as the data point ID, result ID, and molecule details.
    """
    id: UUID = Field(..., description="Unique identifier for the data point")
    result_id: UUID = Field(..., description="ID of the result this data belongs to")
    molecule: MoleculeResponse = Field(..., description="Details of the associated molecule")


class ResultBase(BaseModel):
    """
    Base Pydantic model for result data with common fields.
    
    This model defines the basic structure for experimental results, including
    the submission ID, status, and any notes or comments.
    """
    submission_id: UUID = Field(..., description="ID of the submission this result belongs to")
    status: Optional[str] = Field(ResultStatus.PENDING.name, 
                               description="Current status of the result")
    notes: Optional[str] = Field(None, description="Additional notes or comments about the result")

    # Use the validator to check status values
    _validate_status = validator('status', pre=True)(validate_status)


class ResultCreate(ResultBase):
    """
    Pydantic model for creating a new result.
    
    Extends the base model with fields for associated files and structured data points.
    """
    files: Optional[List[ResultFileCreate]] = Field(default_factory=list, 
                                                 description="List of result files")
    data: Optional[List[ResultDataCreate]] = Field(default_factory=list,
                                               description="List of structured data points")

    # Use the validator to check data formatting
    _validate_data = validator('data', pre=True)(validate_data)


class ResultUpdate(BaseModel):
    """
    Pydantic model for updating an existing result.
    
    All fields are optional to allow partial updates. Only provided fields
    will be updated in the database.
    """
    status: Optional[str] = Field(None, description="Updated status of the result")
    notes: Optional[str] = Field(None, description="Updated notes or comments")

    # Use the validator to check status values
    _validate_status = validator('status', pre=True)(validate_status)


class ResultInDBBase(ResultBase):
    """
    Base Pydantic model for result data as stored in the database.
    
    Extends the base model with additional fields that are present in the database,
    such as ID and timestamps.
    """
    id: UUID = Field(..., description="Unique identifier for the result")
    uploaded_at: datetime = Field(..., description="Timestamp when the result was uploaded")
    approved_at: Optional[datetime] = Field(None, description="Timestamp when the result was approved (if applicable)")


class ResultRead(ResultInDBBase):
    """
    Pydantic model for result data in API responses.
    
    Extends the database model with additional fields needed for API responses,
    such as the submission details and associated files.
    """
    submission: SubmissionRead = Field(..., description="Details of the associated submission")
    files: List[ResultFileRead] = Field(default_factory=list, description="List of result files")
    file_count: int = Field(..., description="Number of files associated with this result")


class ResultDetailedRead(ResultInDBBase):
    """
    Pydantic model for detailed result data in API responses.
    
    Extends the basic result model with comprehensive details including
    all structured data points for in-depth analysis.
    """
    submission: SubmissionRead = Field(..., description="Details of the associated submission")
    files: List[ResultFileRead] = Field(default_factory=list, description="List of result files")
    data_points: List[ResultDataRead] = Field(default_factory=list, description="List of structured data points")


class ResultList(BaseModel):
    """
    Pydantic model for paginated list of results in API responses.
    
    Used for paginated responses containing multiple results.
    Includes the total count, page number, and page size for pagination.
    """
    items: List[ResultRead] = Field(..., description="List of results")
    total: int = Field(..., description="Total number of results matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")
    filters: Optional[Dict[str, Any]] = Field(None, description="Applied filters")


class ResultFilter(BaseModel):
    """
    Pydantic model for filtering results by various criteria.
    
    This model defines the available filter parameters that can be used
    when querying results. All fields are optional to allow flexible filtering.
    """
    submission_id: Optional[UUID] = Field(None, description="Filter by submission ID")
    status: Optional[str] = Field(None, description="Filter by result status")
    uploaded_after: Optional[datetime] = Field(None, description="Filter by upload date (after)")
    uploaded_before: Optional[datetime] = Field(None, description="Filter by upload date (before)")
    skip: Optional[int] = Field(0, description="Number of results to skip (pagination)")
    limit: Optional[int] = Field(100, description="Maximum number of results to return")
    sort_by: Optional[str] = Field("uploaded_at", description="Field to sort by")
    sort_desc: Optional[bool] = Field(True, description="Sort in descending order if true")
    
    # Use the validator to check status values
    _validate_status = validator('status', pre=True)(validate_status)


class ResultApproval(BaseModel):
    """
    Pydantic model for approving or rejecting a result.
    
    This model is used when a pharma user responds to uploaded results,
    either approving them or requesting changes.
    """
    approved: bool = Field(..., description="Whether the result is approved")
    notes: Optional[str] = Field(None, description="Additional notes explaining the decision")