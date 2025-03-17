"""
Implements the API endpoints for molecule management in the Molecular Data Management and CRO Integration Platform.
This module provides RESTful endpoints for creating, retrieving, updating, and deleting molecules,
as well as specialized endpoints for molecule filtering, flagging, bulk operations, similarity search,
and substructure search.
"""

from typing import List, Dict, Optional, Any
from uuid import UUID

from fastapi import APIRouter, Depends, HTTPException, status, Response, Query  # FastAPI 0.95+
from sqlalchemy.orm import Session  # SQLAlchemy 2.0+

# Internal imports for dependencies, schemas, and services
from ..deps import get_db_session, get_current_user_from_token, get_current_pharma_user, get_current_admin_user
from ..schemas.molecules import (
    MoleculeRead,
    MoleculeList,
    PropertyFilter,
    MoleculeFlagRequest,
    MoleculeBulkOperation,
    BulkOperationResponse,
    SimilaritySearchRequest,
    SubstructureSearchRequest,
    BULK_OPERATION_TYPE
)
from '../../schemas/molecule' import MoleculeCreate, MoleculeUpdate, MoleculeFilter
from '../../schemas/msg' import Msg
from '../../services/molecule_service' import molecule_service
from '../../services/library_service' import library_service
from '../../services/experiment_service' import experiment_service
from '../../exceptions' import NotFoundException, ValidationException, BusinessLogicException

# Router instance for molecule endpoints
router = APIRouter(prefix='/molecules', tags=['molecules'])


@router.get('/{molecule_id}', response_model=MoleculeRead)
def get_molecule(molecule_id: UUID, current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Get a single molecule by ID

    Args:
        molecule_id: UUID of the molecule to retrieve
        current_user: Current user information (dependency)

    Returns:
        MoleculeRead: Molecule data
    """
    try:
        # Get molecule by ID using molecule_service.get_by_id
        molecule = molecule_service.get_molecule(molecule_id)
    except NotFoundException as e:
        # If molecule not found, raise HTTPException with 404 status
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the molecule data
    return molecule


@router.post('/', response_model=MoleculeRead, status_code=status.HTTP_201_CREATED)
def create_molecule(molecule_data: MoleculeCreate, current_user: Dict[str, Any] = Depends(get_current_pharma_user)):
    """
    Create a new molecule

    Args:
        molecule_data: Molecule data for creation
        current_user: Current user information (dependency)

    Returns:
        MoleculeRead: Created molecule data
    """
    try:
        # Extract user_id from current_user
        user_id = current_user.get("user_id")

        # Create molecule using molecule_service.create with molecule_data and user_id
        molecule = molecule_service.create_molecule(molecule_data, user_id)
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the created molecule data
    return molecule


@router.put('/{molecule_id}', response_model=MoleculeRead)
def update_molecule(molecule_id: UUID, molecule_data: MoleculeUpdate, current_user: Dict[str, Any] = Depends(get_current_pharma_user)):
    """
    Update an existing molecule

    Args:
        molecule_id: UUID of the molecule to update
        molecule_data: Molecule data for update
        current_user: Current user information (dependency)

    Returns:
        MoleculeRead: Updated molecule data
    """
    try:
        # Get molecule by ID using molecule_service.get_by_id
        molecule = molecule_service.get_molecule(molecule_id)
        if not molecule:
            # If molecule not found, raise HTTPException with 404 status
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")

        # Check if current user is the creator of the molecule or an admin
        user_id = current_user.get("user_id")
        if molecule.get("created_by") != user_id and current_user.get("role") != "admin":
            # If not authorized, raise HTTPException with 403 status
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized to update this molecule")

        # Update molecule using molecule_service.update with molecule_id and molecule_data
        updated_molecule = molecule_service.update_molecule(molecule_id, molecule_data)
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the updated molecule data
    return updated_molecule


@router.delete('/{molecule_id}', response_model=Msg)
def delete_molecule(molecule_id: UUID, current_user: Dict[str, Any] = Depends(get_current_admin_user)):
    """
    Delete a molecule

    Args:
        molecule_id: UUID of the molecule to delete
        current_user: Current user information (dependency)

    Returns:
        Msg: Success message
    """
    try:
        # Get molecule by ID using molecule_service.get_by_id
        molecule = molecule_service.get_molecule(molecule_id)
        if not molecule:
            # If molecule not found, raise HTTPException with 404 status
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")

        # Delete molecule using molecule_service.delete with molecule_id
        molecule_service.delete_molecule(molecule_id)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return success message
    return Msg(msg="Molecule deleted successfully")


@router.get('/', response_model=MoleculeList)
def get_molecules(
    smiles_pattern: Optional[str] = None,
    flag_status: Optional[str] = None,
    property_filters: Optional[List[PropertyFilter]] = None,
    library_id: Optional[UUID] = None,
    experiment_id: Optional[UUID] = None,
    skip: int = 0,
    limit: int = 100,
    sort_by: Optional[str] = None,
    sort_desc: bool = False,
    current_user: Dict[str, Any] = Depends(get_current_user_from_token)
):
    """
    Get molecules with filtering, sorting, and pagination

    Args:
        smiles_pattern: Optional SMILES pattern for filtering
        flag_status: Optional flag status for filtering
        property_filters: Optional list of property filters
        library_id: Optional library ID for filtering
        experiment_id: Optional experiment ID for filtering
        skip: Number of records to skip (for pagination)
        limit: Maximum number of records to return
        sort_by: Field to sort by
        sort_desc: Sort in descending order if True
        current_user: Current user information (dependency)

    Returns:
        MoleculeList: Paginated list of molecules
    """
    try:
        # Create MoleculeFilter object with all filter parameters
        molecule_filter = MoleculeFilter(
            smiles=smiles_pattern,
            flag_status=flag_status,
            # property_filters=property_filters, # TODO: Implement property filters
            # library_id=library_id, # TODO: Implement library ID filter
            # experiment_id=experiment_id, # TODO: Implement experiment ID filter
            # sort_by=sort_by, # TODO: Implement sorting
            # sort_desc=sort_desc # TODO: Implement sorting direction
        )

        # Get molecules using molecule_service.get_molecules with filter object
        molecules_data = molecule_service.get_molecules(filters=molecule_filter, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the molecules list with pagination info
    return MoleculeList(**molecules_data)


@router.get('/me', response_model=MoleculeList)
def get_my_molecules(
    skip: int = 0,
    limit: int = 100,
    current_user: Dict[str, Any] = Depends(get_current_user_from_token)
):
    """
    Get molecules created by the current user

    Args:
        skip: Number of records to skip (for pagination)
        limit: Maximum number of records to return
        current_user: Current user information (dependency)

    Returns:
        MoleculeList: Paginated list of user's molecules
    """
    try:
        # Extract user_id from current_user
        user_id = current_user.get("user_id")

        # Get user's molecules using molecule_service.get_by_user with user_id, skip, and limit
        molecules_data = molecule_service.get_molecules_by_user(user_id=user_id, skip=skip, limit=limit)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the molecules list with pagination info
    return MoleculeList(**molecules_data)


@router.put('/{molecule_id}/flag', response_model=MoleculeRead)
def flag_molecule(molecule_id: UUID, flag_data: MoleculeFlagRequest, current_user: Dict[str, Any] = Depends(get_current_pharma_user)):
    """
    Update the flag status of a molecule

    Args:
        molecule_id: UUID of the molecule to update
        flag_data: Flag data for update
        current_user: Current user information (dependency)

    Returns:
        MoleculeRead: Updated molecule data
    """
    try:
        # Get molecule by ID using molecule_service.get_by_id
        molecule = molecule_service.get_molecule(molecule_id)
        if not molecule:
            # If molecule not found, raise HTTPException with 404 status
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")

        # Check if current user is the creator of the molecule or an admin
        user_id = current_user.get("user_id")
        if molecule.get("created_by") != user_id and current_user.get("role") != "admin":
            # If not authorized, raise HTTPException with 403 status
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized to update this molecule")

        # Update molecule flag status using molecule_service.update_flag_status with molecule_id and flag_status
        updated_molecule = molecule_service.update_molecule_flag(molecule_id=molecule_id, flag_status=flag_data.flag_status)
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return the updated molecule data
    return updated_molecule


@router.post('/bulk', response_model=BulkOperationResponse)
def bulk_operation(operation_data: MoleculeBulkOperation, current_user: Dict[str, Any] = Depends(get_current_pharma_user), db: Session = Depends(get_db_session)):
    """
    Perform bulk operations on multiple molecules

    Args:
        operation_data: Bulk operation data
        current_user: Current user information (dependency)
        db: Database session (dependency)

    Returns:
        BulkOperationResponse: Result of the bulk operation
    """
    try:
        # Extract user_id from current_user
        user_id = current_user.get("user_id")

        # Extract operation_type and molecule_ids from operation_data
        operation_type = operation_data.operation_type
        molecule_ids = operation_data.molecule_ids

        if operation_type == BULK_OPERATION_TYPE.FLAG.name:
            # Update flag status for all molecules
            success_count = 0
            failure_count = 0
            failures = []
            for molecule_id in molecule_ids:
                try:
                    molecule_service.update_molecule_flag(molecule_id=molecule_id, flag_status=operation_data.flag_status)
                    success_count += 1
                except Exception as e:
                    failure_count += 1
                    failures.append({"molecule_id": str(molecule_id), "error": str(e)})

            return BulkOperationResponse(
                success=True,
                message="Flag operation completed",
                processed_count=len(molecule_ids),
                success_count=success_count,
                failure_count=failure_count,
                failures=failures,
            )
        elif operation_type == BULK_OPERATION_TYPE.DELETE.name:
            # Delete all molecules (admin only)
            if current_user.get("role") != "admin":
                raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Only admin users can perform delete operation")

            success_count = 0
            failure_count = 0
            failures = []
            for molecule_id in molecule_ids:
                try:
                    molecule_service.delete_molecule(molecule_id=molecule_id)
                    success_count += 1
                except Exception as e:
                    failure_count += 1
                    failures.append({"molecule_id": str(molecule_id), "error": str(e)})

            return BulkOperationResponse(
                success=True,
                message="Delete operation completed",
                processed_count=len(molecule_ids),
                success_count=success_count,
                failure_count=failure_count,
                failures=failures,
            )
        elif operation_type == BULK_OPERATION_TYPE.ADD_TO_LIBRARY.name:
            # Add molecules to library
            library_id = operation_data.library_id
            result = library_service.add_molecules_to_library(library_id=library_id, molecule_ids=molecule_ids, user_id=user_id)

            return BulkOperationResponse(
                success=True,
                message="Add to library operation completed",
                processed_count=len(molecule_ids),
                success_count=result["success"],
                failure_count=len(result["failures"]),
                failures=result["failures"],
            )
        elif operation_type == BULK_OPERATION_TYPE.ADD_TO_EXPERIMENT.name:
            # Add molecules to experiment
            experiment_id = operation_data.experiment_id
            result = experiment_service.add_molecules_to_experiment(experiment_id=experiment_id, molecule_ids=molecule_ids, user_id=user_id)

            return BulkOperationResponse(
                success=True,
                message="Add to experiment operation completed",
                processed_count=len(molecule_ids),
                success_count=result["success"],
                failure_count=len(result["failures"]),
                failures=result["failures"],
            )
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid operation type")
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))


@router.get('/properties/{property_name}/range')
def get_property_ranges(property_name: str, current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Get min and max values for a specific property

    Args:
        property_name: Name of the property to get range for
        current_user: Current user information (dependency)

    Returns:
        Dict[str, float]: Min and max values for the property
    """
    try:
        # Get property ranges using molecule_service.get_property_ranges with property_name
        property_ranges = molecule_service.get_property_ranges(property_name=property_name)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return dictionary with min_value and max_value
    return property_ranges


@router.get('/properties/available')
def get_available_properties(current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Get list of all available property names

    Args:
        current_user: Current user information (dependency)

    Returns:
        List[str]: List of property names
    """
    try:
        # Get available properties using molecule_service.get_available_properties
        properties = molecule_service.get_available_properties()
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return list of property names
    return properties


@router.post('/search/similarity')
def search_by_similarity(search_params: SimilaritySearchRequest, current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Search for molecules similar to a reference molecule

    Args:
        search_params: Similarity search parameters
        current_user: Current user information (dependency)

    Returns:
        Dict[str, Any]: Search results with similarity scores
    """
    try:
        # Extract search parameters (smiles, molecule_id, threshold, limit)
        smiles = search_params.smiles
        molecule_id = search_params.molecule_id
        threshold = search_params.threshold
        limit = search_params.limit

        # If smiles is provided, use it directly
        if smiles:
            reference_smiles = smiles
        elif molecule_id:
            # If molecule_id is provided, get the molecule's SMILES
            molecule = molecule_service.get_molecule(molecule_id)
            if not molecule:
                raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")
            reference_smiles = molecule.get("smiles")
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Either smiles or molecule_id must be provided")

        # Perform similarity search using molecule_service.search_by_similarity
        search_results = molecule_service.search_by_similarity(reference_smiles=reference_smiles, threshold=threshold, limit=limit)
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return search results with molecules and similarity scores
    return search_results


@router.post('/search/substructure')
def search_by_substructure(search_params: SubstructureSearchRequest, current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Search for molecules containing a specific substructure

    Args:
        search_params: Substructure search parameters
        current_user: Current user information (dependency)

    Returns:
        Dict[str, Any]: Search results with matching molecules
    """
    try:
        # Extract search parameters (smarts, limit)
        smarts = search_params.smarts
        limit = search_params.limit

        # Perform substructure search using molecule_service.search_by_substructure
        search_results = molecule_service.search_by_substructure(substructure_smarts=smarts, limit=limit)
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

    # Return search results with matching molecules
    return search_results


@router.get('/{molecule_id}/image')
def get_molecule_image(molecule_id: UUID, width: int = 300, height: int = 200, format: str = 'svg', current_user: Dict[str, Any] = Depends(get_current_user_from_token)):
    """
    Generate an image of a molecule structure

    Args:
        molecule_id: UUID of the molecule
        width: Width of the image
        height: Height of the image
        format: Format of the image (svg or png)
        current_user: Current user information (dependency)

    Returns:
        Response: Image response with appropriate content type
    """
    try:
        # Get molecule by ID using molecule_service.get_by_id
        molecule = molecule_service.get_molecule(molecule_id)
        if not molecule:
            # If molecule not found, raise HTTPException with 404 status
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")

        # Generate molecule image using molecule_service.generate_image with molecule's SMILES
        image_data = molecule_service.generate_image(smiles=molecule.get("smiles"), width=width, height=height)

        # Set appropriate content type based on format (image/svg+xml for SVG, image/png for PNG)
        content_type = "image/svg+xml" if format.lower() == "svg" else "image/png"

        # Return Response with image data and content type
        return Response(content=image_data, media_type=content_type)
    except NotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=e.message)
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))