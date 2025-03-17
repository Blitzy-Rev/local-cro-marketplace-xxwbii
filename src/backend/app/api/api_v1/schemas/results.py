"""
Pydantic schema models for result-related API requests and responses in the Molecular Data Management
and CRO Integration Platform.

These schemas are used for data validation, serialization, and OpenAPI documentation in the result endpoints.
"""

from pydantic import BaseModel, Field, validator  # pydantic version 2.0+
from typing import List, Dict, Optional, Any, Union  # standard library
from uuid import UUID  # standard library
from enum import Enum  # standard library

from ...schemas.result import (
    ResultBase, ResultCreate, ResultUpdate, ResultRead, ResultDetailedRead,
    ResultList, ResultFilter, ResultApproval, 
    ResultFileBase, ResultFileCreate, ResultFileRead,
    ResultDataBase, ResultDataCreate, ResultDataRead,
    RESULT_STATUS_VALUES
)
from ...schemas.msg import Msg

# Create an Enum for result status values for better API documentation
RESULT_STATUS = Enum('RESULT_STATUS', {status: status for status in RESULT_STATUS_VALUES})


def validate_file_list(v: List[Dict[str, Any]], values: dict) -> List[Dict[str, Any]]:
    """
    Validates that the file list is not empty when uploading results.
    
    Args:
        v: The file list to validate
        values: The values dict from pydantic
        
    Returns:
        The validated file list
        
    Raises:
        ValueError: If the file list is empty
    """
    if v is None or len(v) == 0:
        raise ValueError("At least one file must be provided when uploading results")
    return v


def validate_status(v: str, values: dict) -> str:
    """
    Validates that a result status is one of the allowed values.
    
    Args:
        v: The status value to validate
        values: The values dict from pydantic
        
    Returns:
        The validated status if valid
        
    Raises:
        ValueError: If the status is not valid
    """
    if v is None:
        return 'PENDING'
    
    if v not in RESULT_STATUS_VALUES:
        raise ValueError(f"Invalid status: {v}. Must be one of {RESULT_STATUS_VALUES}")
    
    return v


class ResultDetailRead(ResultDetailedRead):
    """
    Schema for detailed result data in API responses, extending the base ResultDetailedRead.
    
    This schema includes all fields from ResultDetailedRead and may add additional
    API-specific fields as needed.
    """
    pass


class ResultUpload(BaseModel):
    """
    Schema for uploading result files and structured data.
    
    This schema is used when a CRO uploads experimental results, including both
    file references and structured data points.
    """
    submission_id: UUID = Field(..., description="ID of the submission this result belongs to")
    status: Optional[str] = Field("UPLOADED", description="Status of the result (default: UPLOADED)")
    structured_data: List[Dict[str, Any]] = Field(default_factory=list, 
                                               description="Structured data points from the experiment")
    notes: Optional[str] = Field(None, description="Additional notes or comments about the result")
    
    @validator('submission_id')
    def validate_submission_id(cls, v):
        """Validates that the submission_id is provided"""
        if v is None:
            raise ValueError("Submission ID is required")
        return v
    
    @validator('status')
    def validate_status(cls, v):
        """Validates that the status is a valid result status"""
        return validate_status(v, {})


class ResultStatusUpdate(BaseModel):
    """
    Schema for updating the status of a result.
    
    This schema is used when updating just the status of a result, such as
    when approving or rejecting results.
    """
    status: str = Field(..., description="New status for the result")
    notes: Optional[str] = Field(None, description="Additional notes about the status change")
    
    @validator('status')
    def validate_status(cls, v):
        """Validates that the status is a valid result status"""
        return validate_status(v, {})


class ResultFileUpload(BaseModel):
    """
    Schema for uploading a result file.
    
    This schema captures metadata about a file being uploaded as part of
    experimental results.
    """
    file_name: str = Field(..., description="Original name of the uploaded file")
    file_type: str = Field(..., description="MIME type or file extension")
    description: Optional[str] = Field(None, description="Optional description of the file contents")


class ResultDataPoint(BaseModel):
    """
    Schema for a single data point in structured result data.
    
    This schema captures a single experimental measurement linking a molecule
    to a specific data value.
    """
    molecule_id: UUID = Field(..., description="ID of the molecule this data is associated with")
    data_name: str = Field(..., description="Name or identifier of the data point (e.g., 'IC50', 'Solubility')")
    data_value: float = Field(..., description="Numerical value of the measurement")
    data_unit: Optional[str] = Field(None, description="Unit of measurement (e.g., 'nM', 'mg/mL')")
    
    @validator('data_name')
    def validate_data_name(cls, v):
        """Validates that the data_name is not empty"""
        if v is None or v.strip() == "":
            raise ValueError("Data name cannot be empty")
        return v


class ResultAnalysis(BaseModel):
    """
    Schema for result data analysis.
    
    This schema is used for requesting or returning analysis of experimental
    results, including summary statistics and visualizations.
    """
    result_id: UUID = Field(..., description="ID of the result to analyze")
    summary_statistics: Dict[str, Dict[str, Any]] = Field(..., 
                                                     description="Summary statistics grouped by data types")
    molecule_data: Dict[str, List[Dict[str, Any]]] = Field(..., 
                                                      description="Molecule-specific data for each measurement")
    visualizations: Optional[Dict[str, Any]] = Field(None, 
                                                 description="Optional visualization configurations or data")


class ResultExport(BaseModel):
    """
    Schema for exporting result data.
    
    This schema is used when requesting an export of result data in a specific
    format, with options for included fields and files.
    """
    result_id: UUID = Field(..., description="ID of the result to export")
    export_format: str = Field(..., description="Format for export (csv, xlsx, json)")
    include_fields: Optional[List[str]] = Field(None, description="Specific fields to include in export")
    include_files: Optional[bool] = Field(True, description="Whether to include result files in export")
    
    @validator('export_format')
    def validate_export_format(cls, v):
        """Validates that the export_format is one of the allowed values"""
        allowed_formats = ["csv", "xlsx", "json"]
        if v not in allowed_formats:
            raise ValueError(f"Export format must be one of {allowed_formats}")
        return v


class ResultBatchApproval(BaseModel):
    """
    Schema for batch approval of multiple results.
    
    This schema is used when approving or rejecting multiple results at once
    for efficiency in workflows with many results.
    """
    result_ids: List[UUID] = Field(..., description="List of result IDs to approve or reject")
    approved: bool = Field(..., description="Whether to approve (True) or reject (False) the results")
    notes: Optional[str] = Field(None, description="Additional notes about the approval decision")
    
    @validator('result_ids')
    def validate_result_ids(cls, v):
        """Validates that the result_ids list is not empty"""
        if v is None or len(v) == 0:
            raise ValueError("At least one result ID must be provided")
        return v