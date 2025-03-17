from enum import Enum  # standard library
from typing import List, Dict, Optional, Any, Union  # standard library
import re  # standard library

from pydantic import BaseModel, Field, validator, constr  # pydantic 2.0+

from ....schemas.msg import Msg  # Base message schema for simple responses
from ....constants import (
    MAX_CSV_SIZE_MB,
    CSV_CHUNK_SIZE,
    DEFAULT_MOLECULE_PROPERTIES,
)

# Define CSV processing status enum
CSV_PROCESSING_STATUS = Enum('CSV_PROCESSING_STATUS', [
    'PENDING',
    'PROCESSING',
    'COMPLETED',
    'FAILED',
    'CANCELLED'
])


def validate_file_id(v: str, values: dict) -> str:
    """Validates that the file_id is a non-empty string"""
    if not v or not isinstance(v, str) or v.strip() == "":
        raise ValueError("File ID is required and must be a non-empty string")
    return v


def validate_mapping(v: dict, values: dict) -> dict:
    """Validates that the mapping contains a SMILES column mapping"""
    if not v or not isinstance(v, dict):
        raise ValueError("Mapping must be a non-empty dictionary")
    
    # Check if any of the values is 'SMILES' (case-insensitive)
    has_smiles = any(val.upper() == 'SMILES' for val in v.values())
    
    if not has_smiles:
        raise ValueError("Mapping must include a column mapped to 'SMILES'")
    
    return v


def validate_property_name(v: str, values: dict) -> str:
    """Validates that a property name is a valid identifier"""
    if not v or not isinstance(v, str) or v.strip() == "":
        raise ValueError("Property name is required and must be a non-empty string")
    
    # Check if the property name is a valid Python identifier
    if not re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', v):
        raise ValueError("Property name must be a valid identifier (letters, numbers, underscore, starting with letter or underscore)")
    
    return v


class CSVUploadResponse(BaseModel):
    """Response schema for CSV file upload endpoint"""
    file_id: str = Field(..., description="Unique identifier for the uploaded CSV file")
    headers: List[str] = Field(..., description="List of CSV column headers detected in the file")
    row_count: int = Field(..., description="Number of data rows in the CSV file")
    message: str = Field("File uploaded successfully", description="Status message about the upload")


class PropertyMapping(BaseModel):
    """Schema for mapping a CSV column to a system property"""
    csv_column: str = Field(..., description="Name of the column in the CSV file")
    system_property: str = Field(..., description="Name of the corresponding system property")


class CSVMappingRequest(BaseModel):
    """Request schema for mapping CSV columns to system properties"""
    file_id: str = Field(..., description="Unique identifier for the uploaded CSV file")
    mapping: Dict[str, str] = Field(
        ..., 
        description="Mapping of CSV column names to system property names (key: CSV column, value: system property)"
    )
    
    # Validators
    @validator('file_id')
    def _validate_file_id(cls, v):
        return validate_file_id(v, {})
    
    @validator('mapping')
    def _validate_mapping(cls, v):
        return validate_mapping(v, {})


class CSVMappingResponse(BaseModel):
    """Response schema for CSV mapping endpoint"""
    job_id: str = Field(..., description="Unique identifier for the CSV processing job")
    status: str = Field(..., description="Initial status of the processing job")
    message: str = Field("CSV mapping submitted for processing", description="Status message")


class CSVProcessingStatus(BaseModel):
    """Schema for CSV processing status information"""
    status: str = Field(..., description="Current status of the CSV processing job")
    progress: Optional[int] = Field(None, description="Processing progress as a percentage (0-100)")
    errors: Optional[Dict[str, Any]] = Field(None, description="Error information if processing failed")


class CSVImportSummary(BaseModel):
    """Schema for CSV import summary information"""
    total_rows: int = Field(..., description="Total number of rows in the CSV file")
    processed_rows: int = Field(..., description="Number of rows that were processed")
    successful_rows: int = Field(..., description="Number of rows successfully imported")
    failed_rows: int = Field(..., description="Number of rows that failed to import")
    errors: List[Dict[str, Any]] = Field([], description="List of error details for failed rows")


class CSVStatusResponse(BaseModel):
    """Response schema for CSV processing status endpoint"""
    job_id: str = Field(..., description="Unique identifier for the CSV processing job")
    status: str = Field(..., description="Current status of the processing job")
    progress: Optional[int] = Field(None, description="Processing progress as a percentage (0-100)")
    summary: Optional[CSVImportSummary] = Field(None, description="Summary of import results if processing is complete")


class SystemProperty(BaseModel):
    """Schema for system property information"""
    name: str = Field(..., description="Name of the system property")
    description: str = Field(..., description="Description of the property")
    data_type: str = Field(..., description="Data type of the property (string, number, boolean, etc.)")
    required: Optional[bool] = Field(False, description="Whether this property is required")


class AvailablePropertiesResponse(BaseModel):
    """Response schema for available system properties endpoint"""
    properties: List[SystemProperty] = Field(..., description="List of available system properties")