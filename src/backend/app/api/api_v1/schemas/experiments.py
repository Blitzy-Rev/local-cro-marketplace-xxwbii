"""
Defines Pydantic schema models for experiment-related API requests and responses in the API v1 namespace.
These schemas are used for data validation, serialization, and documentation in the experiment management endpoints.
"""

from typing import Optional, List, Dict, Any, Union
from uuid import UUID
from datetime import datetime

from pydantic import BaseModel, Field, validator, root_validator

# Import base experiment schemas for API schemas
from ....schemas.experiment import ExperimentBase, ExperimentCreate, ExperimentUpdate, ExperimentRead, ExperimentDetailRead, ExperimentList, ExperimentFilter, ExperimentMoleculeOperation, ExperimentStatusUpdate
from ....models.experiment import ExperimentStatus
from .users import UserProfileResponse
from .molecules import MoleculeRead

# Create a list of valid experiment status values
EXPERIMENT_STATUS_VALUES = [status.name for status in ExperimentStatus]


@validator('name', pre=True)
def validate_experiment_name(cls, v):
    """
    Validates that an experiment name meets the required criteria.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        The validated experiment name if valid
        
    Raises:
        ValueError: If the name is empty or too long
    """
    if v is None or len(v.strip()) == 0:
        raise ValueError("Experiment name cannot be empty")
    
    if len(v) > 100:
        raise ValueError("Experiment name cannot exceed 100 characters")
    
    return v


@validator('status', pre=True)
def validate_status(cls, v):
    """
    Validates that an experiment status is one of the allowed values.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        The validated status if valid
        
    Raises:
        ValueError: If the status is not a valid ExperimentStatus value
    """
    if v is None:
        return 'DRAFT'  # Default status
    
    if v not in EXPERIMENT_STATUS_VALUES:
        raise ValueError(f"Invalid status: {v}. Must be one of {EXPERIMENT_STATUS_VALUES}")
    
    return v


@root_validator(pre=True)
def validate_parameters(cls, values):
    """
    Validates that experiment parameters are of the correct type and format.
    
    Args:
        cls: The class being validated
        values: The values dictionary
        
    Returns:
        The validated values dictionary
        
    Raises:
        ValueError: If parameters are invalid
    """
    parameters = values.get('parameters')
    if parameters is None:
        return values
    
    if not isinstance(parameters, dict):
        raise ValueError("Parameters must be a dictionary")
    
    # Validate parameter names and values
    for param_name, param_value in parameters.items():
        if not param_name or not isinstance(param_name, str):
            raise ValueError(f"Parameter name must be a non-empty string: {param_name}")
        
        # Check if parameter name is a valid identifier
        if not param_name.isidentifier():
            raise ValueError(f"Parameter name must be a valid identifier: {param_name}")
    
    return values


class ExperimentCreateRequest(BaseModel):
    """Schema for creating a new experiment through the API."""
    name: str
    type_id: UUID
    status: Optional[str] = None
    parameters: Optional[Dict[str, str]] = None


class ExperimentUpdateRequest(BaseModel):
    """Schema for updating an existing experiment through the API."""
    name: Optional[str] = None
    status: Optional[str] = None
    parameters: Optional[Dict[str, str]] = None


class ExperimentResponse(BaseModel):
    """Schema for experiment data in API responses."""
    id: UUID
    name: str
    type_id: UUID
    status: str
    parameters: Dict[str, str]
    created_by: UUID
    created_at: datetime
    updated_at: Optional[datetime] = None
    experiment_type: Dict[str, Any]
    creator: UserProfileResponse
    molecule_count: int


class ExperimentDetailResponse(BaseModel):
    """Schema for detailed experiment data in API responses."""
    id: UUID
    name: str
    type_id: UUID
    status: str
    parameters: Dict[str, str]
    created_by: UUID
    created_at: datetime
    updated_at: Optional[datetime] = None
    experiment_type: Dict[str, Any]
    creator: UserProfileResponse
    molecules: List[MoleculeRead]
    submissions: List[Dict[str, Any]]


class ExperimentListResponse(BaseModel):
    """Schema for paginated list of experiments in API responses."""
    items: List[ExperimentResponse]
    total: int
    page: int
    size: int
    filters: Optional[Dict[str, Any]] = None


class ExperimentFilterParams(BaseModel):
    """Schema for filtering experiments in API requests."""
    name_contains: Optional[str] = None
    type_id: Optional[UUID] = None
    status: Optional[str] = None
    created_by: Optional[UUID] = None
    molecule_id: Optional[UUID] = None
    skip: Optional[int] = None
    limit: Optional[int] = None
    sort_by: Optional[str] = None
    sort_desc: Optional[bool] = None
    
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
            ValueError: If the status is not a valid ExperimentStatus value
        """
        if v is None:
            return None
        
        if v not in EXPERIMENT_STATUS_VALUES:
            raise ValueError(f"Invalid status filter: {v}. Must be one of {EXPERIMENT_STATUS_VALUES}")
        
        return v


class ExperimentMoleculeOperationRequest(BaseModel):
    """Schema for adding or removing molecules from an experiment."""
    molecule_ids: List[UUID]
    operation: str  # 'add' or 'remove'
    
    @validator('operation', pre=True)
    def validate_operation(cls, v):
        """
        Validates that the operation is one of the allowed values.
        
        Args:
            cls: The class being validated
            v: The value to validate
            
        Returns:
            The validated operation if valid
            
        Raises:
            ValueError: If the operation is not 'add' or 'remove'
        """
        if v not in ['add', 'remove']:
            raise ValueError(f"Invalid operation: {v}. Must be 'add' or 'remove'")
        
        return v


class ExperimentStatusUpdateRequest(BaseModel):
    """Schema for updating the status of an experiment."""
    status: str
    notes: Optional[str] = None
    
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
            ValueError: If the status is not a valid ExperimentStatus value
        """
        if v not in EXPERIMENT_STATUS_VALUES:
            raise ValueError(f"Invalid status: {v}. Must be one of {EXPERIMENT_STATUS_VALUES}")
        
        return v


class ExperimentTypeResponse(BaseModel):
    """Schema for experiment type data in API responses."""
    id: UUID
    name: str
    description: Optional[str] = None
    category: Optional[str] = None