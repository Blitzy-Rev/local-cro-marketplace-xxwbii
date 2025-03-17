"""
Defines Pydantic schema models for library-related API requests and responses in the
Molecular Data Management and CRO Integration Platform.

These schemas are used for data validation, serialization, and OpenAPI documentation
in the library endpoints.
"""

from pydantic import BaseModel, Field, validator  # version 2.0+
from typing import List, Dict, Optional, Any, Union  # standard library
from uuid import UUID  # standard library

from ...schemas.library import (
    LibraryBase, LibraryCreate, LibraryUpdate, LibraryInDBBase, 
    LibraryResponse, LibraryDetailResponse, LibraryListResponse,
    LibraryFilter, LibraryMoleculeOperation
)
from ...schemas.molecule import MoleculeResponse
from ...schemas.msg import Msg


def validate_operation(v: str, values: dict) -> str:
    """
    Validates that the operation is one of the allowed values for library molecule operations.
    
    Args:
        v: The operation string to validate.
        values: Dictionary of other values in the model (not used in this validator).
        
    Returns:
        The validated operation if valid.
        
    Raises:
        ValueError: If the operation is not valid.
    """
    if v is None or v.strip() == "":
        raise ValueError("Operation cannot be empty")
    
    allowed_operations = ["add", "remove"]
    if v.lower() not in allowed_operations:
        raise ValueError(f"Operation must be one of {allowed_operations}")
    
    return v.lower()


class LibraryRead(LibraryResponse):
    """Schema for library data in API responses, extending the base LibraryResponse."""
    pass


class LibraryDetailRead(LibraryDetailResponse):
    """
    Schema for detailed library data including molecules in API responses, 
    extending the base LibraryDetailResponse.
    """
    pass


class LibraryList(LibraryListResponse):
    """
    Schema for paginated list of libraries in API responses, 
    extending the base LibraryListResponse.
    """
    pass


class LibraryCreateRequest(LibraryCreate):
    """Schema for creating a new library in API requests."""
    pass


class LibraryUpdateRequest(LibraryUpdate):
    """Schema for updating a library in API requests."""
    pass


class LibraryMoleculeOperationRequest(BaseModel):
    """Schema for adding or removing molecules from a library in API requests."""
    molecule_ids: List[UUID] = Field(
        ...,
        description="List of molecule IDs to add or remove from the library"
    )
    operation: str = Field(
        ...,
        description="Operation to perform: 'add' or 'remove'",
        example="add"
    )
    
    @validator("operation")
    def validate_operation(cls, v):
        """Validates that the operation is one of the allowed values."""
        return validate_operation(v, {})
    
    @validator("molecule_ids")
    def validate_molecule_ids(cls, v):
        """Validates that the molecule_ids list is not empty."""
        if v is None or len(v) == 0:
            raise ValueError("Must provide at least one molecule ID")
        return v


class LibraryOperationResponse(BaseModel):
    """Schema for library operation response."""
    success: bool = Field(..., description="Whether the operation was successful")
    message: str = Field(..., description="Response message")
    processed_count: int = Field(..., description="Total number of molecules processed")
    success_count: int = Field(..., description="Number of successfully processed molecules")
    failure_count: int = Field(..., description="Number of failed molecule operations")
    failures: Optional[List[Dict[str, Any]]] = Field(
        None,
        description="Details of failures if any occurred"
    )