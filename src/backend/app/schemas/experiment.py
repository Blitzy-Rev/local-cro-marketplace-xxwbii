"""
Pydantic schema models for experiment data validation, serialization, and API responses
in the Molecular Data Management and CRO Integration Platform.

This module defines the data models used for experiment creation, management, and submission
to CROs, with appropriate validation rules, type annotations, and documentation.
"""

from typing import Optional, List, Dict, Any, Union  # standard library
from datetime import datetime  # standard library
from uuid import UUID  # standard library

from pydantic import BaseModel, Field, validator  # pydantic version 2.0+

from ..constants import ExperimentStatus  # Internal import for status enumeration
from .experiment_type import ExperimentTypeRead  # Internal import for experiment type reference
from .user import UserResponse  # Internal import for user reference
from .molecule import MoleculeResponse  # Internal import for molecule reference


@validator('status', pre=True)
def validate_status(cls, v):
    """
    Validates that an experiment status value is one of the allowed values.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        The validated status if valid
        
    Raises:
        ValueError: If the status is not a valid ExperimentStatus value
    """
    if v is None:
        return None
    
    try:
        # Convert to uppercase for case-insensitive comparison
        v_upper = v.upper() if isinstance(v, str) else v
        # Check if it's a valid enum value
        valid_statuses = [status.name for status in ExperimentStatus]
        if v_upper not in valid_statuses:
            raise ValueError(f"Invalid status: {v}. Must be one of {valid_statuses}")
        return v_upper
    except (AttributeError, ValueError) as e:
        raise ValueError(f"Invalid status: {v}. Must be one of {[status.name for status in ExperimentStatus]}")


@validator('parameters', pre=True)
def validate_parameters(cls, v):
    """
    Validates that experiment parameters are properly formatted.
    
    Args:
        cls: The class being validated
        v: The value to validate (list of parameter dictionaries)
        
    Returns:
        The validated parameters list if valid
        
    Raises:
        ValueError: If parameters are not properly formatted
    """
    if v is None:
        return []
    
    if not isinstance(v, list):
        raise ValueError("Parameters must be a list")
    
    for param in v:
        if not isinstance(param, dict):
            raise ValueError("Each parameter must be a dictionary")
        
        if 'parameter_name' not in param or not param['parameter_name']:
            raise ValueError("Each parameter must have a non-empty 'parameter_name'")
        
        if 'parameter_value' not in param:
            raise ValueError("Each parameter must have a 'parameter_value'")
    
    return v


class ExperimentParameterBase(BaseModel):
    """
    Base Pydantic model for experiment parameter data with common fields.
    
    This model defines the basic structure for experiment parameters, including
    the parameter name and value fields.
    """
    parameter_name: str = Field(..., description="Name of the parameter")
    parameter_value: str = Field(..., description="Value of the parameter")


class ExperimentParameterCreate(ExperimentParameterBase):
    """
    Pydantic model for creating a new experiment parameter.
    
    Inherits all fields from ExperimentParameterBase.
    """
    pass


class ExperimentParameterResponse(ExperimentParameterBase):
    """
    Pydantic model for experiment parameter data in API responses.
    
    Extends the base model with additional fields needed for API responses.
    """
    id: UUID = Field(..., description="Unique identifier for the parameter")
    experiment_id: UUID = Field(..., description="ID of the experiment this parameter belongs to")


class ExperimentBase(BaseModel):
    """
    Base Pydantic model for experiment data with common fields.
    
    This model defines the basic structure for experiments, including the name,
    type, description, status, and parameters fields.
    """
    name: str = Field(..., description="Name of the experiment")
    type_id: UUID = Field(..., description="ID of the experiment type")
    description: Optional[str] = Field(None, description="Detailed description of the experiment")
    status: Optional[str] = Field(ExperimentStatus.DRAFT.name, description="Current status of the experiment")
    parameters: Optional[List[ExperimentParameterCreate]] = Field(default_factory=list, 
                                                               description="List of experiment parameters")


class ExperimentCreate(ExperimentBase):
    """
    Pydantic model for creating a new experiment.
    
    Extends the base model with fields needed for experiment creation, such as
    the creator ID and associated molecule IDs.
    """
    created_by: Optional[UUID] = Field(None, description="ID of the user creating the experiment")
    molecule_ids: Optional[List[UUID]] = Field(default_factory=list, 
                                           description="IDs of molecules to associate with the experiment")


class ExperimentUpdate(BaseModel):
    """
    Pydantic model for updating an existing experiment.
    
    All fields are optional to allow partial updates. Only provided fields
    will be updated in the database.
    """
    name: Optional[str] = Field(None, description="Updated name of the experiment")
    type_id: Optional[UUID] = Field(None, description="Updated experiment type ID")
    description: Optional[str] = Field(None, description="Updated description of the experiment")
    status: Optional[str] = Field(None, description="Updated status of the experiment")
    parameters: Optional[List[ExperimentParameterCreate]] = Field(None, 
                                                               description="Updated list of experiment parameters")
    molecule_ids: Optional[List[UUID]] = Field(None, 
                                           description="Updated list of molecule IDs to associate with the experiment")


class ExperimentInDBBase(ExperimentBase):
    """
    Base Pydantic model for experiment data as stored in the database.
    
    Extends the base model with additional fields that are present in the database,
    such as ID, creator ID, and timestamps.
    """
    id: UUID = Field(..., description="Unique identifier for the experiment")
    created_by: UUID = Field(..., description="ID of the user who created this experiment")
    created_at: datetime = Field(..., description="Timestamp when the experiment was created")
    updated_at: Optional[datetime] = Field(None, description="Timestamp when the experiment was last updated")

    class Config:
        orm_mode = True


class ExperimentResponse(ExperimentInDBBase):
    """
    Pydantic model for experiment data in API responses.
    
    Extends the database model with additional fields needed for API responses,
    such as the experiment type, parameters, creator details, and molecule count.
    """
    experiment_type: ExperimentTypeRead = Field(..., description="Details of the experiment type")
    parameters: List[ExperimentParameterResponse] = Field(default_factory=list, 
                                                      description="List of experiment parameters")
    creator: Optional[UserResponse] = Field(None, description="User who created this experiment")
    molecule_count: int = Field(..., description="Number of molecules associated with this experiment")

    class Config:
        orm_mode = True


class ExperimentDetailResponse(ExperimentInDBBase):
    """
    Pydantic model for detailed experiment data in API responses.
    
    Extends the database model with more comprehensive details, including
    the full list of associated molecules and submission history.
    """
    experiment_type: ExperimentTypeRead = Field(..., description="Details of the experiment type")
    parameters: List[ExperimentParameterResponse] = Field(default_factory=list, 
                                                      description="List of experiment parameters")
    creator: Optional[UserResponse] = Field(None, description="User who created this experiment")
    molecules: List[MoleculeResponse] = Field(default_factory=list, 
                                          description="List of molecules associated with this experiment")
    submissions: Optional[List[Dict]] = Field(default_factory=list, 
                                           description="List of submissions for this experiment")

    class Config:
        orm_mode = True


class ExperimentListResponse(BaseModel):
    """
    Pydantic model for paginated list of experiments in API responses.
    
    Used for paginated responses containing multiple experiments.
    Includes the total count, page number, and page size for pagination.
    """
    items: List[ExperimentResponse] = Field(..., description="List of experiments")
    total: int = Field(..., description="Total number of experiments matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")


class ExperimentFilter(BaseModel):
    """
    Pydantic model for filtering experiments by various criteria.
    
    This model defines the available filter parameters that can be used
    when querying experiments. All fields are optional to allow flexible filtering.
    """
    name: Optional[str] = Field(None, description="Filter by experiment name (partial match)")
    type_id: Optional[UUID] = Field(None, description="Filter by experiment type ID")
    status: Optional[str] = Field(None, description="Filter by experiment status")
    created_by: Optional[UUID] = Field(None, description="Filter by creator user ID")
    created_after: Optional[datetime] = Field(None, description="Filter by creation date (after)")
    created_before: Optional[datetime] = Field(None, description="Filter by creation date (before)")
    molecule_id: Optional[UUID] = Field(None, description="Filter by associated molecule ID")
    sort_by: Optional[str] = Field(None, description="Field to sort by")
    sort_desc: Optional[bool] = Field(False, description="Sort in descending order if true")