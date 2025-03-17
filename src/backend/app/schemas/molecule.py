"""
Pydantic schema models for molecule data validation, serialization, and API responses.

This module defines the data models used for molecule-related operations in the 
Molecular Data Management and CRO Integration Platform, with appropriate validation 
rules, type annotations, and documentation.
"""

from pydantic import BaseModel, Field, validator
from typing import Optional, List, Dict, Any
from uuid import UUID
from datetime import datetime

from ..models.molecule import Molecule
from .molecule_property import MoleculePropertyBase, MoleculePropertyCreate, MoleculePropertyResponse
from .user import UserResponse
from ..constants import MOLECULE_FLAG_STATUS
from ..utils.validation_utils import validate_smiles, validate_molecular_property, validate_molecule_flag_status


class MoleculeBase(BaseModel):
    """Base Pydantic model for molecule data with common fields"""
    smiles: str = Field(..., 
                        description="SMILES representation of the molecular structure",
                        example="CCO")
    flag_status: Optional[str] = Field(None, 
                                      description="Optional flag for priority review (e.g., 'HIGH', 'MEDIUM', 'LOW')",
                                      example="HIGH")
    properties: Optional[List[MoleculePropertyCreate]] = Field(default_factory=list,
                                                            description="List of molecular properties")

    @validator('smiles', pre=True)
    def validate_smiles(cls, v):
        """Validates that a SMILES string is properly formatted and represents a valid molecular structure"""
        if not v or not isinstance(v, str):
            raise ValueError("SMILES string cannot be empty")
        
        if not validate_smiles(v):
            raise ValueError(f"Invalid SMILES structure: {v}")
        
        return v

    @validator('properties', pre=True)
    def validate_properties(cls, v):
        """Validates that molecule properties are properly formatted and contain valid values"""
        if v is None:
            return []
        
        if not isinstance(v, list):
            raise ValueError("Properties must be a list")
        
        for prop in v:
            if not isinstance(prop, dict):
                raise ValueError("Each property must be a dictionary")
            
            if 'property_name' not in prop:
                raise ValueError("Each property must have a 'property_name' field")
            
            # Validate property name
            if not validate_molecular_property(prop['property_name'], prop.get('property_value')):
                raise ValueError(f"Invalid property: {prop['property_name']} with value {prop.get('property_value')}")
        
        return v

    @validator('flag_status', pre=True)
    def validate_flag_status(cls, v):
        """Validates that a flag status value is one of the allowed values"""
        if v is None:
            return None
        
        if not validate_molecule_flag_status(v):
            valid_statuses = list(MOLECULE_FLAG_STATUS.keys()) if MOLECULE_FLAG_STATUS else []
            raise ValueError(f"Invalid flag status: {v}. Must be one of {valid_statuses}")
        
        return v


class MoleculeCreate(MoleculeBase):
    """Pydantic model for creating a new molecule"""
    created_by: Optional[UUID] = Field(None, 
                                      description="ID of the user creating the molecule")


class MoleculeUpdate(BaseModel):
    """Pydantic model for updating an existing molecule"""
    flag_status: Optional[str] = Field(None, 
                                      description="Optional flag for priority review (e.g., 'HIGH', 'MEDIUM', 'LOW')",
                                      example="HIGH")
    properties: Optional[List[MoleculePropertyCreate]] = Field(None,
                                                           description="List of molecular properties to update")


class MoleculeInDBBase(MoleculeBase):
    """Base Pydantic model for molecule data as stored in the database"""
    id: UUID = Field(..., description="Unique identifier for the molecule")
    created_by: UUID = Field(..., description="ID of the user who created this molecule")
    created_at: datetime = Field(..., description="Timestamp when the molecule was created")
    updated_at: Optional[datetime] = Field(None, description="Timestamp when the molecule was last updated")

    class Config:
        from_attributes = True


class MoleculeResponse(MoleculeInDBBase):
    """Pydantic model for molecule data in API responses"""
    properties: List[MoleculePropertyResponse] = Field(default_factory=list,
                                                   description="List of molecular properties")
    creator: Optional[UserResponse] = Field(None, description="User who created this molecule")

    class Config:
        from_attributes = True


class MoleculeDetailResponse(MoleculeResponse):
    """Pydantic model for detailed molecule data in API responses"""
    image_url: Optional[str] = Field(None, 
                                    description="URL to a rendered image of the molecular structure",
                                    example="/api/v1/molecules/1/image")

    class Config:
        from_attributes = True


class MoleculeListResponse(BaseModel):
    """Pydantic model for paginated list of molecules in API responses"""
    items: List[MoleculeResponse] = Field(..., description="List of molecules")
    total: int = Field(..., description="Total number of molecules matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")


class MoleculeFilter(BaseModel):
    """Pydantic model for filtering molecules by various criteria"""
    smiles: Optional[str] = Field(None, description="Filter by SMILES string (partial match)")
    created_by: Optional[UUID] = Field(None, description="Filter by creator user ID")
    flag_status: Optional[str] = Field(None, description="Filter by flag status")
    properties: Optional[List[Dict]] = Field(None, description="Filter by property values")
    created_after: Optional[datetime] = Field(None, description="Filter by creation date (after)")
    created_before: Optional[datetime] = Field(None, description="Filter by creation date (before)")
    library_ids: Optional[List[UUID]] = Field(None, description="Filter by library membership")
    experiment_ids: Optional[List[UUID]] = Field(None, description="Filter by experiment membership")
    sort_by: Optional[str] = Field(None, description="Field to sort by")
    sort_desc: Optional[bool] = Field(False, description="Sort in descending order if true")