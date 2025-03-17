"""
Pydantic schema models for library data validation, serialization, and API responses.

This module defines the data models used for library-related operations in the
Molecular Data Management and CRO Integration Platform, with appropriate validation
rules, type annotations, and documentation.
"""

from pydantic import BaseModel, Field, validator  # version 2.0+
from typing import Optional, List, Dict, Any, Union  # standard library
from uuid import UUID  # standard library
from datetime import datetime  # standard library

from ..models.library import Library
from .user import UserResponse
from .molecule import MoleculeResponse
from ..utils.validation_utils import validate_string_length, validate_library_name, validate_description


class LibraryBase(BaseModel):
    """Base Pydantic model for library data with common fields."""
    name: str = Field(
        ..., 
        description="Name of the library",
        example="High Activity Compounds",
        min_length=1,
        max_length=100
    )
    description: Optional[str] = Field(
        None, 
        description="Detailed description of the library's purpose",
        example="Collection of compounds with high binding affinity to target protein A"
    )

    @validator('name')
    def validate_library_name(cls, v):
        """Validates that a library name is properly formatted and not empty."""
        if v is None or v.strip() == "":
            raise ValueError("Library name cannot be empty")
        
        if not validate_string_length(v, 1, 100):
            raise ValueError("Library name must be between 1 and 100 characters")
        
        if not validate_library_name(v):
            raise ValueError("Library name contains invalid characters")
        
        return v
    
    @validator('description')
    def validate_description(cls, v):
        """Validates that a library description is properly formatted."""
        if v is None:
            return None
        
        if not validate_description(v):
            raise ValueError("Library description contains invalid characters")
        
        return v
    
    class Config:
        from_attributes = True


class LibraryCreate(LibraryBase):
    """Pydantic model for creating a new library."""
    created_by: Optional[UUID] = Field(
        None, 
        description="ID of the user creating the library"
    )


class LibraryUpdate(BaseModel):
    """Pydantic model for updating an existing library."""
    name: Optional[str] = Field(
        None, 
        description="Updated name of the library",
        example="Optimized High Activity Compounds"
    )
    description: Optional[str] = Field(
        None, 
        description="Updated description of the library's purpose",
        example="Refined collection of high-affinity compounds"
    )

    @validator('name')
    def validate_library_name(cls, v):
        """Validates that a library name is properly formatted if provided."""
        if v is None:
            return None
        
        if v.strip() == "":
            raise ValueError("Library name cannot be empty")
        
        if not validate_string_length(v, 1, 100):
            raise ValueError("Library name must be between 1 and 100 characters")
        
        if not validate_library_name(v):
            raise ValueError("Library name contains invalid characters")
        
        return v
    
    @validator('description')
    def validate_description(cls, v):
        """Validates that a library description is properly formatted if provided."""
        if v is None:
            return None
        
        if not validate_description(v):
            raise ValueError("Library description contains invalid characters")
        
        return v
    
    class Config:
        from_attributes = True


class LibraryInDBBase(LibraryBase):
    """Base Pydantic model for library data as stored in the database."""
    id: UUID = Field(..., description="Unique identifier for the library")
    created_by: UUID = Field(..., description="ID of the user who created this library")
    created_at: datetime = Field(..., description="Timestamp when the library was created")
    updated_at: Optional[datetime] = Field(None, description="Timestamp when the library was last updated")

    class Config:
        from_attributes = True


class LibraryResponse(LibraryInDBBase):
    """Pydantic model for library data in API responses."""
    creator: Optional[UserResponse] = Field(
        None, 
        description="User who created this library"
    )
    molecule_count: int = Field(
        0, 
        description="Number of molecules in this library",
        ge=0
    )

    class Config:
        from_attributes = True


class LibraryDetailResponse(LibraryInDBBase):
    """Pydantic model for detailed library data in API responses, including molecules."""
    creator: Optional[UserResponse] = Field(
        None, 
        description="User who created this library"
    )
    molecules: List[MoleculeResponse] = Field(
        default_factory=list,
        description="List of molecules in this library"
    )

    class Config:
        from_attributes = True


class LibraryListResponse(BaseModel):
    """Pydantic model for paginated list of libraries in API responses."""
    items: List[LibraryResponse] = Field(..., description="List of libraries")
    total: int = Field(..., description="Total number of libraries matching the criteria")
    page: int = Field(..., description="Current page number")
    size: int = Field(..., description="Number of items per page")


class LibraryFilter(BaseModel):
    """Pydantic model for filtering libraries by various criteria."""
    name: Optional[str] = Field(
        None, 
        description="Filter by library name (partial match)",
        example="activity"
    )
    created_by: Optional[UUID] = Field(
        None, 
        description="Filter by creator user ID"
    )
    created_after: Optional[datetime] = Field(
        None, 
        description="Filter by creation date (after)"
    )
    created_before: Optional[datetime] = Field(
        None, 
        description="Filter by creation date (before)"
    )
    contains_molecule_id: Optional[UUID] = Field(
        None, 
        description="Filter libraries containing a specific molecule"
    )
    sort_by: Optional[str] = Field(
        None, 
        description="Field to sort by",
        example="name"
    )
    sort_desc: Optional[bool] = Field(
        False, 
        description="Sort in descending order if true"
    )


class LibraryMoleculeOperation(BaseModel):
    """Pydantic model for adding or removing molecules from a library."""
    library_id: UUID = Field(
        ..., 
        description="ID of the library to modify"
    )
    molecule_ids: List[UUID] = Field(
        ..., 
        description="List of molecule IDs to add or remove"
    )