"""
Pydantic schema models for molecular property validation, serialization, and API responses.

This module defines the data models for molecular properties in the Molecular Data Management
and CRO Integration Platform, supporting CSV import, property filtering, and data validation.
"""

from pydantic import BaseModel, Field, validator
from typing import Optional, Union, Dict, Any
from uuid import UUID
import re

from ..utils.validation_utils import validate_molecular_property

# Regex pattern for property name validation: starts with letter, followed by letters, digits, or underscores
PROPERTY_NAME_REGEX = r'^[a-zA-Z][a-zA-Z0-9_]*$'


class MoleculePropertyBase(BaseModel):
    """Base Pydantic model for molecular property data with common fields"""
    property_name: str = Field(..., description="Name of the molecular property")
    property_value: Optional[Union[float, int, str, bool]] = Field(None, description="Value of the property")
    property_unit: Optional[str] = Field(None, description="Unit of measurement for the property")
    is_calculated: Optional[bool] = Field(False, description="Whether the property was calculated or imported")

    @validator('property_name', pre=True)
    def validate_property_name(cls, v):
        """Validates that a property name follows the required format (alphanumeric with underscores, starting with a letter)"""
        if not v or not isinstance(v, str):
            raise ValueError("Property name must be a non-empty string")
        
        if not re.match(PROPERTY_NAME_REGEX, v):
            raise ValueError(f"Property name '{v}' does not match required format (must start with a letter and contain only letters, numbers, and underscores)")
        
        return v

    @validator('property_value', pre=True)
    def validate_property_value(cls, v):
        """Validates that a property value is of an appropriate type (numeric for most properties)"""
        if v is None:
            return None
        
        # Try to convert string values to appropriate numeric types
        if isinstance(v, str):
            try:
                # Try to convert to float
                float_val = float(v)
                # If it's a whole number, convert to int
                if float_val.is_integer():
                    return int(float_val)
                return float_val
            except ValueError:
                # If not a number, keep as string
                pass
        
        # Validate that the value is of an appropriate type
        if not isinstance(v, (float, int, str, bool)):
            raise ValueError(f"Property value must be a number, string, or boolean, got {type(v).__name__}")
        
        return v


class MoleculePropertyCreate(MoleculePropertyBase):
    """Pydantic model for creating a new molecular property"""
    molecule_id: Optional[UUID] = Field(None, description="ID of the molecule this property belongs to")


class MoleculePropertyUpdate(BaseModel):
    """Pydantic model for updating an existing molecular property"""
    property_value: Optional[Union[float, int, str, bool]] = Field(None, description="Updated value of the property")
    property_unit: Optional[str] = Field(None, description="Updated unit of measurement")


class MoleculePropertyInDBBase(MoleculePropertyBase):
    """Base Pydantic model for molecular property data as stored in the database"""
    id: UUID = Field(..., description="Unique identifier for the property record")
    molecule_id: UUID = Field(..., description="ID of the molecule this property belongs to")


class MoleculePropertyResponse(MoleculePropertyInDBBase):
    """Pydantic model for molecular property data in API responses"""
    pass


class PropertyFilter(BaseModel):
    """Pydantic model for filtering molecules by property values"""
    name: str = Field(..., description="Name of the property to filter by")
    value: Optional[Union[float, int, str, bool]] = Field(None, description="Exact value to match")
    min_value: Optional[Union[float, int, str, bool]] = Field(None, description="Minimum value (inclusive)")
    max_value: Optional[Union[float, int, str, bool]] = Field(None, description="Maximum value (inclusive)")
    operator: Optional[str] = Field(None, description="Comparison operator (eq, gt, lt, gte, lte, contains)")