from typing import List, Dict, Any, Optional
from uuid import UUID
from fastapi import APIRouter, Depends, HTTPException, status, Path, Query, Body
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+

from ..deps import get_current_pharma_user, get_current_admin_user, get_db_session
from ..schemas.experiments import (
    ExperimentCreateRequest,
    ExperimentUpdateRequest,
    ExperimentResponse,
    ExperimentDetailResponse,
    ExperimentListResponse,
    ExperimentFilterParams,
    ExperimentMoleculeOperationRequest,
    ExperimentStatusUpdateRequest,
    ExperimentTypeResponse
)
from ...services import experiment_service
from ...exceptions import ValidationException, ResourceNotFoundException

experiments_router = APIRouter(prefix='/experiments', tags=['experiments'])

@experiments_router.get('/', response_model=ExperimentListResponse)
def get_experiments(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    filters: ExperimentFilterParams = Depends()
) -> ExperimentListResponse:
    """
    Get a paginated list of experiments with optional filtering
    """
    user_id = current_user.get("user_id")
    filters.created_by = user_id
    experiment_list_response = experiment_service.get_experiments(filters)
    return experiment_list_response

@experiments_router.get('/{experiment_id}', response_model=ExperimentDetailResponse)
def get_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to retrieve')
) -> ExperimentDetailResponse:
    """
    Get detailed information about a specific experiment
    """
    try:
        experiment_detail = experiment_service.get_experiment_by_id(experiment_id, include_molecules=True)
        if not experiment_detail:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_id = current_user.get("user_id")
        user_role = current_user.get("role")

        if experiment_detail["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        return experiment_detail
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.post('/', response_model=ExperimentResponse, status_code=status.HTTP_201_CREATED)
def create_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_data: ExperimentCreateRequest = Body(...)
) -> ExperimentResponse:
    """
    Create a new experiment
    """
    user_id = current_user.get("user_id")
    try:
        created_experiment = experiment_service.create_experiment(experiment_data, user_id)
        return created_experiment
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.put('/{experiment_id}', response_model=ExperimentResponse)
def update_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to update'),
    experiment_data: ExperimentUpdateRequest = Body(...)
) -> ExperimentResponse:
    """
    Update an existing experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        updated_experiment = experiment_service.update_experiment(experiment_id, experiment_data, user_id)
        return updated_experiment
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.delete('/{experiment_id}', response_model=Dict[str, Any])
def delete_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to delete')
) -> Dict[str, Any]:
    """
    Delete an experiment (only allowed for DRAFT or CANCELLED status)
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        experiment_service.delete_experiment(experiment_id, user_id)
        return {"message": "Experiment deleted successfully"}
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.put('/{experiment_id}/status', response_model=ExperimentResponse)
def update_experiment_status(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to update'),
    status_data: ExperimentStatusUpdateRequest = Body(...)
) -> ExperimentResponse:
    """
    Update the status of an experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        updated_experiment = experiment_service.update_experiment_status(experiment_id, status_data.status, user_id)
        return updated_experiment
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.post('/{experiment_id}/queue', response_model=ExperimentResponse)
def queue_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to queue')
) -> ExperimentResponse:
    """
    Change experiment status to QUEUED for submission to CRO
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        queued_experiment = experiment_service.queue_experiment(experiment_id, user_id)
        return queued_experiment
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.get('/{experiment_id}/prepare-submission', response_model=Dict[str, Any])
def prepare_for_submission(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment to prepare')
) -> Dict[str, Any]:
    """
    Prepare an experiment for submission to CRO by validating requirements
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        prepared_experiment = experiment_service.prepare_for_submission(experiment_id, user_id)
        return prepared_experiment
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.get('/{experiment_id}/molecules', response_model=Dict[str, Any])
def get_experiment_molecules(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment'),
    skip: int = Query(0, ge=0, description='Skip N records'),
    limit: int = Query(100, ge=1, le=1000, description='Limit to N records')
) -> Dict[str, Any]:
    """
    Get molecules associated with an experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        molecules = experiment_service.get_experiment_molecules(experiment_id, skip, limit)
        return molecules
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.post('/{experiment_id}/molecules', response_model=Dict[str, Any])
def add_molecules_to_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment'),
    molecule_data: ExperimentMoleculeOperationRequest = Body(...)
) -> Dict[str, Any]:
    """
    Add molecules to an experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        if molecule_data.operation != "add":
            raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail="Invalid operation")

        result = experiment_service.add_molecules_to_experiment(experiment_id, molecule_data.molecule_ids, user_id)
        return result
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.delete('/{experiment_id}/molecules', response_model=Dict[str, Any])
def remove_molecules_from_experiment(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment'),
    molecule_data: ExperimentMoleculeOperationRequest = Body(...)
) -> Dict[str, Any]:
    """
    Remove molecules from an experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        if molecule_data.operation != "remove":
            raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail="Invalid operation")

        result = experiment_service.remove_molecules_from_experiment(experiment_id, molecule_data.molecule_ids, user_id)
        return result
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.get('/{experiment_id}/parameters', response_model=Dict[str, str])
def get_experiment_parameters(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment')
) -> Dict[str, str]:
    """
    Get parameters for a specific experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        parameters = experiment_service.get_experiment_parameters(experiment_id)
        return parameters
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.put('/{experiment_id}/parameters/{parameter_name}', response_model=Dict[str, Any])
def set_experiment_parameter(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user),
    experiment_id: UUID = Path(..., description='The ID of the experiment'),
    parameter_name: str = Path(..., description='The name of the parameter'),
    parameter_data: Dict[str, str] = Body(...)
) -> Dict[str, Any]:
    """
    Set a parameter value for an experiment
    """
    user_id = current_user.get("user_id")
    try:
        existing_experiment = experiment_service.get_experiment_by_id(experiment_id)
        if not existing_experiment:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Experiment not found")

        user_role = current_user.get("role")
        if existing_experiment["created_by"] != user_id and user_role != "admin":
            raise HTTPException(status_code=status.HTTP_403_FORBIDDEN, detail="Not authorized")

        parameter_value = parameter_data.get("parameter_value")
        updated_parameter = experiment_service.set_experiment_parameter(experiment_id, parameter_name, parameter_value)
        return updated_parameter
    except ResourceNotFoundException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except ValidationException as e:
        raise HTTPException(status_code=status.HTTP_422_UNPROCESSABLE_ENTITY, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))

@experiments_router.get('/types', response_model=List[ExperimentTypeResponse])
def get_experiment_types(
    current_user: Dict[str, Any] = Depends(get_current_pharma_user)
) -> List[ExperimentTypeResponse]:
    """
    Get all available experiment types
    """
    experiment_types = experiment_service.get_experiment_types()
    return experiment_types

def handle_validation_exception(exc: ValidationException):
    """
    Utility function to handle ValidationException and convert to HTTPException
    """
    raise HTTPException(status_code=422, detail=exc.message)

def handle_resource_not_found_exception(exc: ResourceNotFoundException):
    """
    Utility function to handle ResourceNotFoundException and convert to HTTPException
    """
    raise HTTPException(status_code=404, detail=exc.message)