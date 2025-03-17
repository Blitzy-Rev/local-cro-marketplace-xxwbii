from typing import List, Optional, Dict, Any
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status, Path, Query
from sqlalchemy.orm import Session

from ..deps import get_current_pharma_user, get_current_admin_user, get_db_session
from ..schemas.libraries import (
    LibraryRead, 
    LibraryDetailRead, 
    LibraryList, 
    LibraryCreateRequest, 
    LibraryUpdateRequest, 
    LibraryMoleculeOperationRequest, 
    LibraryOperationResponse
)
from ...schemas.library import LibraryFilter
from ...schemas.msg import Msg
from ...services.library_service import library_service
from ...exceptions import NotFoundException, ValidationException, ConflictException, AuthorizationException

import logging

# Create router with prefix and tag
libraries_router = APIRouter(prefix="/libraries", tags=["libraries"])
logger = logging.getLogger(__name__)

@libraries_router.get("/", response_model=LibraryList)
def get_libraries(
    name_contains: Optional[str] = Query(None, description="Filter libraries by name containing this string"),
    skip: Optional[int] = Query(0, description="Number of libraries to skip"),
    limit: Optional[int] = Query(100, description="Maximum number of libraries to return"),
    sort_by: Optional[str] = Query(None, description="Field to sort by"),
    sort_desc: Optional[bool] = Query(False, description="Sort in descending order"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> LibraryList:
    """
    Get a paginated list of libraries with optional filtering.
    
    Args:
        name_contains: Filter libraries by name containing this string
        skip: Number of libraries to skip for pagination
        limit: Maximum number of libraries to return
        sort_by: Field to sort by (e.g., "name", "created_at")
        sort_desc: Sort in descending order if True
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        LibraryList: Paginated list of libraries
    """
    # Create filter object with the provided parameters
    filters = LibraryFilter(
        name=name_contains,
        created_by=current_user["user_id"],
        sort_by=sort_by,
        sort_desc=sort_desc
    )
    
    # Get libraries using the library service
    return library_service.get_libraries(
        filters=filters,
        skip=skip,
        limit=limit
    )

@libraries_router.post("/", response_model=LibraryRead, status_code=status.HTTP_201_CREATED)
def create_library(
    library_data: LibraryCreateRequest,
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> LibraryRead:
    """
    Create a new library.
    
    Args:
        library_data: Library creation data
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        LibraryRead: Created library
        
    Raises:
        HTTPException: If library creation fails
    """
    try:
        return library_service.create_library(
            library_data=library_data, 
            user_id=current_user["user_id"]
        )
    except ConflictException as e:
        logger.warning(
            f"Conflict while creating library: {str(e)}", 
            {"user_id": current_user["user_id"], "library_name": library_data.name}
        )
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=str(e)
        )
    except ValidationException as e:
        logger.warning(
            f"Validation error while creating library: {str(e)}",
            {"user_id": current_user["user_id"]}
        )
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )

@libraries_router.get("/{library_id}", response_model=LibraryDetailRead)
def get_library(
    library_id: UUID = Path(..., description="The ID of the library to get"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> LibraryDetailRead:
    """
    Get a library by ID.
    
    Args:
        library_id: The ID of the library to retrieve
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        LibraryDetailRead: Library with detailed information
        
    Raises:
        HTTPException: If library not found or user not authorized
    """
    try:
        return library_service.get_by_id(
            db=db, 
            id=library_id, 
            user_id=current_user["user_id"]
        )
    except NotFoundException as e:
        logger.warning(
            f"Library not found: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except AuthorizationException as e:
        logger.warning(
            f"Authorization error while getting library: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=str(e)
        )

@libraries_router.put("/{library_id}", response_model=LibraryRead)
def update_library(
    library_id: UUID = Path(..., description="The ID of the library to update"),
    library_data: LibraryUpdateRequest = ...,
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> LibraryRead:
    """
    Update a library by ID.
    
    Args:
        library_id: The ID of the library to update
        library_data: Updated library data
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        LibraryRead: Updated library
        
    Raises:
        HTTPException: If library not found, user not authorized, or validation error
    """
    try:
        return library_service.update_library(
            library_id=library_id,
            library_data=library_data,
            user_id=current_user["user_id"]
        )
    except NotFoundException as e:
        logger.warning(
            f"Library not found for update: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except ConflictException as e:
        logger.warning(
            f"Conflict while updating library: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_409_CONFLICT,
            detail=str(e)
        )
    except AuthorizationException as e:
        logger.warning(
            f"Authorization error while updating library: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=str(e)
        )
    except ValidationException as e:
        logger.warning(
            f"Validation error while updating library: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )

@libraries_router.delete("/{library_id}", response_model=Msg)
def delete_library(
    library_id: UUID = Path(..., description="The ID of the library to delete"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> Msg:
    """
    Delete a library by ID.
    
    Args:
        library_id: The ID of the library to delete
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        Msg: Success message
        
    Raises:
        HTTPException: If library not found or user not authorized
    """
    try:
        library_service.delete_library(
            library_id=library_id,
            user_id=current_user["user_id"]
        )
        return {"msg": "Library successfully deleted"}
    except NotFoundException as e:
        logger.warning(
            f"Library not found for deletion: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except AuthorizationException as e:
        logger.warning(
            f"Authorization error while deleting library: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=str(e)
        )

@libraries_router.get("/{library_id}/molecules")
def get_library_molecules(
    library_id: UUID = Path(..., description="The ID of the library"),
    skip: int = Query(0, description="Number of molecules to skip"),
    limit: int = Query(100, description="Maximum number of molecules to return"),
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> Dict[str, Any]:
    """
    Get molecules in a library with pagination.
    
    Args:
        library_id: The ID of the library to get molecules from
        skip: Number of molecules to skip for pagination
        limit: Maximum number of molecules to return
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        Dict containing molecules list and pagination info
        
    Raises:
        HTTPException: If library not found or user not authorized
    """
    try:
        return library_service.get_library_molecules(
            library_id=library_id,
            skip=skip,
            limit=limit
        )
    except NotFoundException as e:
        logger.warning(
            f"Library not found for getting molecules: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except AuthorizationException as e:
        logger.warning(
            f"Authorization error while getting library molecules: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=str(e)
        )

@libraries_router.post("/{library_id}/molecules", response_model=LibraryOperationResponse)
def library_molecule_operation(
    library_id: UUID = Path(..., description="The ID of the library"),
    operation_data: LibraryMoleculeOperationRequest = ...,
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    db: Session = Depends(get_db_session)
) -> LibraryOperationResponse:
    """
    Add or remove molecules from a library.
    
    Args:
        library_id: The ID of the library to modify
        operation_data: Operation details (add/remove molecules)
        current_user: Current authenticated user (from dependency)
        db: Database session (from dependency)
        
    Returns:
        LibraryOperationResponse: Result of the operation
        
    Raises:
        HTTPException: If library not found, user not authorized, or validation error
    """
    try:
        if operation_data.operation.lower() == "add":
            result = library_service.add_molecules_to_library(
                library_id=library_id,
                molecule_ids=operation_data.molecule_ids,
                user_id=current_user["user_id"]
            )
        elif operation_data.operation.lower() == "remove":
            result = library_service.remove_molecules_from_library(
                library_id=library_id,
                molecule_ids=operation_data.molecule_ids,
                user_id=current_user["user_id"]
            )
        else:
            raise ValidationException(
                f"Invalid operation: {operation_data.operation}. Must be 'add' or 'remove'",
                {"operation": operation_data.operation}
            )
        
        # Convert result to LibraryOperationResponse format
        return LibraryOperationResponse(
            success=True,
            message=f"Successfully {operation_data.operation}ed molecules to library",
            processed_count=len(operation_data.molecule_ids),
            success_count=result.get("success", 0),
            failure_count=len(result.get("failures", [])),
            failures=result.get("failures", [])
        )
    except NotFoundException as e:
        logger.warning(
            f"Library not found for molecule operation: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except ValidationException as e:
        logger.warning(
            f"Validation error in molecule operation: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
            detail=str(e)
        )
    except AuthorizationException as e:
        logger.warning(
            f"Authorization error in molecule operation: {str(e)}",
            {"user_id": current_user["user_id"], "library_id": str(library_id)}
        )
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=str(e)
        )