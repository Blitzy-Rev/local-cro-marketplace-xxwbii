from typing import Optional, List
from uuid import UUID  # standard library
from pydantic import BaseModel, Field, validator  # pydantic version 2.0+

from ..models.experiment_type import ExperimentType


@validator('name', pre=True)
def validate_experiment_type_name(cls, v):
    """
    Validates that an experiment type name meets the required criteria.
    
    This validator ensures that the experiment type name is not empty
    and does not exceed the maximum allowed length.
    
    Args:
        cls: The class being validated
        v: The value to validate
        
    Returns:
        str: The validated experiment type name if valid
        
    Raises:
        ValueError: If the name is empty or too long
    """
    if v is None or len(v.strip()) == 0:
        raise ValueError("Experiment type name cannot be empty")
    
    if len(v) > 100:
        raise ValueError("Experiment type name cannot exceed 100 characters")
    
    return v


class ExperimentTypeBase(BaseModel):
    """
    Base Pydantic model for experiment type data with common fields.
    
    Provides the foundation for validating and serializing experiment type data
    throughout the application. Contains fields common to all experiment type
    representations.
    """
    name: str
    description: Optional[str] = None
    category: Optional[str] = None


class ExperimentTypeCreate(ExperimentTypeBase):
    """
    Pydantic model for creating a new experiment type.
    
    Used for validating data when creating new experiment types in the system.
    Inherits all fields from ExperimentTypeBase.
    """
    pass


class ExperimentTypeUpdate(BaseModel):
    """
    Pydantic model for updating an existing experiment type.
    
    All fields are optional to allow partial updates. Only provided fields
    will be updated in the database.
    """
    name: Optional[str] = None
    description: Optional[str] = None
    category: Optional[str] = None


class ExperimentTypeInDBBase(ExperimentTypeBase):
    """
    Base Pydantic model for experiment type data as stored in the database.
    
    Extends the base model to include the database ID field and configures
    the model to work with ORM objects.
    """
    id: UUID

    class Config:
        orm_mode = True


class ExperimentTypeRead(ExperimentTypeInDBBase):
    """
    Pydantic model for experiment type data in API responses.
    
    This model represents the complete experiment type data that will be
    returned in API responses. Inherits all fields from ExperimentTypeInDBBase.
    """
    pass


class ExperimentTypeList(BaseModel):
    """
    Pydantic model for a list of experiment types in API responses.
    
    Used for paginated responses containing multiple experiment types.
    Includes the total count for pagination purposes.
    """
    items: List[ExperimentTypeRead]
    total: int