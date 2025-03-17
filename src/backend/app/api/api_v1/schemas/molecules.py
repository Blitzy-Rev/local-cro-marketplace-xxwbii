"""
Pydantic schema models for molecule-related API requests and responses in the 
Molecular Data Management and CRO Integration Platform.

These schemas are used for data validation, serialization, and OpenAPI documentation
in the molecule endpoints.
"""

from enum import Enum  # standard library
from typing import List, Dict, Optional, Any, Union  # standard library
from uuid import UUID  # standard library

from pydantic import BaseModel, Field, validator  # pydantic 2.0+

# Import base schemas
from ....schemas.molecule import (
    MoleculeBase, MoleculeCreate, MoleculeUpdate, MoleculeResponse, 
    MoleculeListResponse, MoleculeFilter
)
from ....schemas.molecule import (
    MoleculePropertyCreate, MoleculePropertyUpdate, MoleculePropertyResponse
)
from ....schemas.msg import Msg


# Define operation types for bulk operations
BULK_OPERATION_TYPE = Enum('BULK_OPERATION_TYPE', ['FLAG', 'DELETE', 'ADD_TO_LIBRARY', 'ADD_TO_EXPERIMENT'])


def validate_molecule_ids(v: List[UUID], values: dict) -> List[UUID]:
    """
    Validates that the molecule_ids list is not empty.
    
    Args:
        v: List of molecule UUIDs
        values: Other validated values
        
    Returns:
        The validated molecule_ids list
        
    Raises:
        ValueError: If the list is empty
    """
    if not v or len(v) == 0:
        raise ValueError("molecule_ids cannot be empty")
    return v


def validate_operation_type(v: str, values: dict) -> str:
    """
    Validates that the operation_type is a valid bulk operation type.
    
    Args:
        v: Operation type string
        values: Other validated values
        
    Returns:
        The validated operation_type
        
    Raises:
        ValueError: If the operation type is invalid
    """
    if not v or v.strip() == "":
        raise ValueError("operation_type cannot be empty")
    
    try:
        BULK_OPERATION_TYPE[v]
    except KeyError:
        valid_operations = [op.name for op in BULK_OPERATION_TYPE]
        raise ValueError(f"Invalid operation_type: {v}. Must be one of {valid_operations}")
    
    return v


class MoleculeRead(MoleculeResponse):
    """
    Schema for molecule data in API responses, extending the base MoleculeResponse.
    
    This model is used for returning molecule data in API responses with all
    the necessary fields and relationships.
    """
    pass


class MoleculeList(MoleculeListResponse):
    """
    Schema for paginated list of molecules in API responses, extending the base MoleculeListResponse.
    
    This model provides a standardized format for paginated molecule data with total count
    and additional metadata.
    """
    pass


class PropertyFilter(BaseModel):
    """
    Schema for filtering molecules by property values.
    
    This model defines the structure for property-based filtering in molecule queries,
    supporting minimum and maximum value ranges or discrete value lists.
    """
    property_name: str = Field(..., description="Name of the property to filter by", example="LogP")
    min_value: Optional[float] = Field(None, description="Minimum value (inclusive)", example=1.5)
    max_value: Optional[float] = Field(None, description="Maximum value (inclusive)", example=4.5)
    values: Optional[List[str]] = Field(None, description="List of values to match (OR condition)", example=["aromatic", "aliphatic"])


class MoleculeFlagRequest(BaseModel):
    """
    Schema for flagging a molecule with a specific status.
    
    This model is used for updating the flag status of a molecule, supporting
    priorities like "important", "review", "rejected", or "none".
    """
    flag_status: str = Field(..., description="Flag status to apply", example="important")
    
    @validator('flag_status')
    def validate_flag_status(cls, v):
        """
        Validates that the flag_status is a valid value.
        
        Args:
            v: Flag status value
            
        Returns:
            The validated flag_status
            
        Raises:
            ValueError: If the flag status is invalid
        """
        if not v or v.strip() == "":
            raise ValueError("flag_status cannot be empty")
        
        allowed_values = ["important", "review", "rejected", "none"]
        if v.lower() not in allowed_values:
            raise ValueError(f"Invalid flag_status: {v}. Must be one of {allowed_values}")
        
        return v.lower()


class MoleculeBulkOperation(BaseModel):
    """
    Schema for bulk operations on multiple molecules.
    
    This model is used for performing operations on multiple molecules at once,
    such as flagging, deleting, or adding to libraries or experiments.
    """
    operation_type: str = Field(..., description="Type of operation to perform", example="FLAG")
    molecule_ids: List[UUID] = Field(..., description="List of molecule IDs to operate on")
    flag_status: Optional[str] = Field(None, description="Flag status for FLAG operations", example="important")
    library_id: Optional[UUID] = Field(None, description="Library ID for ADD_TO_LIBRARY operations")
    experiment_id: Optional[UUID] = Field(None, description="Experiment ID for ADD_TO_EXPERIMENT operations")
    
    @validator('molecule_ids')
    def _validate_molecule_ids(cls, v):
        """Validates the molecule_ids field"""
        return validate_molecule_ids(v, {})
    
    @validator('operation_type')
    def _validate_operation_type(cls, v):
        """Validates the operation_type field"""
        return validate_operation_type(v, {})
    
    @validator('flag_status', 'operation_type', always=True)
    def validate_flag_operation(cls, values):
        """
        Validates that flag_status is provided for FLAG operations.
        
        Args:
            values: Dict of validated values
            
        Returns:
            Validated values dictionary
            
        Raises:
            ValueError: If flag_status is missing for FLAG operations
        """
        flag_status = values.get('flag_status')
        operation_type = values.get('operation_type')
        
        if operation_type == 'FLAG' and not flag_status:
            raise ValueError("flag_status is required for FLAG operations")
        
        return values
    
    @validator('library_id', 'operation_type', always=True)
    def validate_library_operation(cls, values):
        """
        Validates that library_id is provided for ADD_TO_LIBRARY operations.
        
        Args:
            values: Dict of validated values
            
        Returns:
            Validated values dictionary
            
        Raises:
            ValueError: If library_id is missing for ADD_TO_LIBRARY operations
        """
        library_id = values.get('library_id')
        operation_type = values.get('operation_type')
        
        if operation_type == 'ADD_TO_LIBRARY' and not library_id:
            raise ValueError("library_id is required for ADD_TO_LIBRARY operations")
        
        return values
    
    @validator('experiment_id', 'operation_type', always=True)
    def validate_experiment_operation(cls, values):
        """
        Validates that experiment_id is provided for ADD_TO_EXPERIMENT operations.
        
        Args:
            values: Dict of validated values
            
        Returns:
            Validated values dictionary
            
        Raises:
            ValueError: If experiment_id is missing for ADD_TO_EXPERIMENT operations
        """
        experiment_id = values.get('experiment_id')
        operation_type = values.get('operation_type')
        
        if operation_type == 'ADD_TO_EXPERIMENT' and not experiment_id:
            raise ValueError("experiment_id is required for ADD_TO_EXPERIMENT operations")
        
        return values


class BulkOperationResponse(BaseModel):
    """
    Schema for bulk operation response.
    
    This model is used for returning the results of a bulk operation, including
    success status, counts, and details of any failures.
    """
    success: bool = Field(..., description="Whether the operation was successful", example=True)
    message: str = Field(..., description="Operation result message", example="Operation completed successfully")
    processed_count: int = Field(..., description="Number of molecules processed", example=10)
    success_count: int = Field(..., description="Number of molecules successfully processed", example=9)
    failure_count: int = Field(..., description="Number of molecules that failed to process", example=1)
    failures: Optional[List[Dict[str, Any]]] = Field(None, description="Details of failed operations")


class SimilaritySearchRequest(BaseModel):
    """
    Schema for similarity search request.
    
    This model is used for searching molecules by structural similarity, either
    by providing a SMILES string or referencing an existing molecule.
    """
    smiles: Optional[str] = Field(None, description="SMILES string to compare against", example="CCO")
    molecule_id: Optional[UUID] = Field(None, description="ID of molecule to compare against")
    threshold: float = Field(0.7, description="Similarity threshold (0-1)", example=0.7)
    limit: int = Field(50, description="Maximum number of results to return", example=50)
    
    @validator('smiles', 'molecule_id', always=True)
    def validate_search_parameters(cls, values):
        """
        Validates that either smiles or molecule_id is provided.
        
        Args:
            values: Dict of validated values
            
        Returns:
            Validated values dictionary
            
        Raises:
            ValueError: If neither smiles nor molecule_id is provided
        """
        smiles = values.get('smiles')
        molecule_id = values.get('molecule_id')
        
        if not smiles and not molecule_id:
            raise ValueError("Either smiles or molecule_id must be provided")
        
        return values
    
    @validator('threshold')
    def validate_threshold(cls, v):
        """
        Validates that threshold is between 0 and 1.
        
        Args:
            v: Threshold value
            
        Returns:
            Validated threshold
            
        Raises:
            ValueError: If threshold is not between 0 and 1
        """
        if not (0 <= v <= 1):
            raise ValueError("threshold must be between 0 and 1")
        
        return v
    
    @validator('limit')
    def validate_limit(cls, v):
        """
        Validates that limit is a positive integer.
        
        Args:
            v: Limit value
            
        Returns:
            Validated limit
            
        Raises:
            ValueError: If limit is not a positive integer
        """
        if v <= 0:
            raise ValueError("limit must be a positive integer")
        
        return v


class SubstructureSearchRequest(BaseModel):
    """
    Schema for substructure search request.
    
    This model is used for searching molecules containing a specific substructure
    pattern defined using SMARTS notation.
    """
    smarts: str = Field(..., description="SMARTS pattern to search for", example="[#6]-[#8]")
    limit: int = Field(50, description="Maximum number of results to return", example=50)
    
    @validator('smarts')
    def validate_smarts(cls, v):
        """
        Validates that smarts is a valid SMARTS pattern.
        
        Args:
            v: SMARTS pattern
            
        Returns:
            Validated SMARTS pattern
            
        Raises:
            ValueError: If SMARTS pattern is invalid
        """
        if not v or v.strip() == "":
            raise ValueError("smarts cannot be empty")
        
        # Note: In a real implementation, we would use RDKit to validate the SMARTS pattern
        # try:
        #     mol = Chem.MolFromSmarts(v)
        #     if mol is None:
        #         raise ValueError(f"Invalid SMARTS pattern: {v}")
        # except Exception as e:
        #     raise ValueError(f"Error parsing SMARTS pattern: {str(e)}")
        
        return v
    
    @validator('limit')
    def validate_limit(cls, v):
        """
        Validates that limit is a positive integer.
        
        Args:
            v: Limit value
            
        Returns:
            Validated limit
            
        Raises:
            ValueError: If limit is not a positive integer
        """
        if v <= 0:
            raise ValueError("limit must be a positive integer")
        
        return v