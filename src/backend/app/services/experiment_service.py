"""
Service layer for experiment management in the Molecular Data Management and CRO Integration Platform.
This service provides high-level business logic for experiment operations, including creation,
retrieval, status management, molecule association, and preparation for CRO submission.
It acts as an intermediary between API endpoints and the data access layer, implementing
core experiment workflow functionality.
"""

import logging  # standard library
from typing import List, Dict, Optional, Any, Tuple, Union  # standard library
from uuid import UUID  # standard library

from ..crud import crud_experiment as experiment  # Internal import: Data access layer for experiment operations
from ..models.experiment import Experiment, ExperimentStatus  # Internal import: Database model for experiments
from ..schemas.experiment import ExperimentCreate, ExperimentUpdate, ExperimentFilter  # Internal import: Schemas for experiment data validation
from ..db.session import get_db  # Internal import: Database session context manager
from . import molecule_service  # Internal import: Service for molecule operations
from .notification_service import send_experiment_status_change  # Internal import: Send notifications for experiment status changes
from ..exceptions import ValidationException, ResourceNotFoundException  # Internal import: Exception classes

# Set up logger
logger = logging.getLogger(__name__)


class ExperimentService:
    """
    Service class for experiment management operations
    """

    def __init__(self):
        """
        Initializes the ExperimentService
        """
        self.logger = logger

    def get_experiment_by_id(self, experiment_id: UUID, include_molecules: bool = False) -> Optional[Dict[str, Any]]:
        """
        Retrieves an experiment by its ID

        Args:
            experiment_id: ID of the experiment
            include_molecules: Whether to include molecule details

        Returns:
            Optional[Dict[str, Any]]: Experiment data if found, None otherwise
        """
        return get_experiment_by_id(experiment_id, include_molecules)

    def get_experiment_by_name(self, name: str, user_id: int) -> Optional[Dict[str, Any]]:
        """
        Retrieves an experiment by its name

        Args:
            name: Name of the experiment
            user_id: ID of the user

        Returns:
            Optional[Dict[str, Any]]: Experiment data if found, None otherwise
        """
        return get_experiment_by_name(name, user_id)

    def create_experiment(self, experiment_data: ExperimentCreate, user_id: int) -> Dict[str, Any]:
        """
        Creates a new experiment

        Args:
            experiment_data: Data for the new experiment
            user_id: ID of the user creating the experiment

        Returns:
            Dict[str, Any]: Created experiment data
        """
        return create_experiment(experiment_data, user_id)

    def update_experiment(self, experiment_id: UUID, experiment_data: ExperimentUpdate, user_id: int) -> Dict[str, Any]:
        """
        Updates an existing experiment

        Args:
            experiment_id: ID of the experiment to update
            experiment_data: Updated experiment data
            user_id: ID of the user making the update

        Returns:
            Dict[str, Any]: Updated experiment data
        """
        return update_experiment(experiment_id, experiment_data, user_id)

    def delete_experiment(self, experiment_id: UUID, user_id: int) -> Dict[str, Any]:
        """
        Deletes an experiment

        Args:
            experiment_id: ID of the experiment to delete
            user_id: ID of the user making the request

        Returns:
            Dict[str, Any]: Result of deletion operation
        """
        return delete_experiment(experiment_id, user_id)

    def update_experiment_status(self, experiment_id: UUID, status: str, user_id: int) -> Dict[str, Any]:
        """
        Updates the status of an experiment

        Args:
            experiment_id: ID of the experiment to update
            status: New status for the experiment
            user_id: ID of the user making the update

        Returns:
            Dict[str, Any]: Updated experiment data
        """
        return update_experiment_status(experiment_id, status, user_id)

    def add_molecules_to_experiment(self, experiment_id: UUID, molecule_ids: List[UUID], user_id: int) -> Dict[str, Any]:
        """
        Adds molecules to an experiment

        Args:
            experiment_id: ID of the experiment
            molecule_ids: List of molecule IDs to add
            user_id: ID of the user making the request

        Returns:
            Dict[str, Any]: Result with success count and failures
        """
        return add_molecules_to_experiment(experiment_id, molecule_ids, user_id)

    def remove_molecules_from_experiment(self, experiment_id: UUID, molecule_ids: List[UUID], user_id: int) -> Dict[str, Any]:
        """
        Removes molecules from an experiment

        Args:
            experiment_id: ID of the experiment
            molecule_ids: List of molecule IDs to remove
            user_id: ID of the user making the request

        Returns:
            Dict[str, Any]: Result with success count and failures
        """
        return remove_molecules_from_experiment(experiment_id, molecule_ids, user_id)

    def get_experiment_molecules(self, experiment_id: UUID, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves molecules associated with an experiment

        Args:
            experiment_id: ID of the experiment
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return

        Returns:
            Dict[str, Any]: Dictionary with molecules list and pagination info
        """
        return get_experiment_molecules(experiment_id, skip, limit)

    def get_experiments(self, filters: ExperimentFilter, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves experiments based on filter criteria

        Args:
            filters: Filter criteria for experiments
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return

        Returns:
            Dict[str, Any]: Dictionary with experiments list and pagination info
        """
        return get_experiments(filters, skip, limit)

    def get_experiments_by_user(self, user_id: int, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves experiments created by a specific user

        Args:
            user_id: ID of the user
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return

        Returns:
            Dict[str, Any]: Dictionary with experiments list and pagination info
        """
        return get_experiments_by_user(user_id, skip, limit)

    def get_experiments_by_status(self, status: str, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves experiments with a specific status

        Args:
            status: Status to filter by
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return

        Returns:
            Dict[str, Any]: Dictionary with experiments list and pagination info
        """
        return get_experiments_by_status(status, skip, limit)

    def queue_experiment(self, experiment_id: UUID, user_id: int) -> Dict[str, Any]:
        """
        Changes experiment status from DRAFT to QUEUED

        Args:
            experiment_id: ID of the experiment to queue
            user_id: ID of the user making the request

        Returns:
            Dict[str, Any]: Updated experiment data
        """
        return queue_experiment(experiment_id, user_id)

    def prepare_for_submission(self, experiment_id: UUID, user_id: int) -> Dict[str, Any]:
        """
        Validates that an experiment is ready for submission to a CRO

        Args:
            experiment_id: ID of the experiment to prepare
            user_id: ID of the user making the request

        Returns:
            Dict[str, Any]: Experiment data if ready for submission
        """
        return prepare_for_submission(experiment_id, user_id)

    def validate_status_transition(self, current_status: str, new_status: str) -> bool:
        """
        Validates if a status transition is allowed

        Args:
            current_status: Current status of the experiment
            new_status: New status to transition to

        Returns:
            bool: True if transition is allowed, False otherwise
        """
        return validate_status_transition(current_status, new_status)


experiment_service = ExperimentService()


def get_experiment_by_id(experiment_id: UUID, include_molecules: bool = False) -> Optional[Dict[str, Any]]:
    """
    Retrieves an experiment by its ID with parameters and molecules
    """
    logger.debug(f"Getting experiment by ID: {experiment_id}")
    with get_db() as db:
        if include_molecules:
            db_experiment = experiment.get_with_molecules(db, experiment_id)
        else:
            db_experiment = experiment.get_with_parameters(db, experiment_id)

        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            return None

        experiment_data = format_experiment_data(db_experiment, include_molecules)
        return experiment_data


def get_experiment_by_name(name: str, user_id: int) -> Optional[Dict[str, Any]]:
    """
    Retrieves an experiment by its name and creator ID
    """
    logger.debug(f"Getting experiment by name: {name} and user_id: {user_id}")
    with get_db() as db:
        db_experiment = experiment.get_by_name(db, name, user_id)
        if not db_experiment:
            logger.warning(f"Experiment with name {name} and user_id {user_id} not found")
            return None

        experiment_data = format_experiment_data(db_experiment)
        return experiment_data


def create_experiment(experiment_data: ExperimentCreate, user_id: int) -> Dict[str, Any]:
    """
    Creates a new experiment with parameters and optional molecules
    """
    logger.info(f"Creating new experiment: {experiment_data.name} for user {user_id}")
    with get_db() as db:
        # Check if experiment with same name already exists for this user
        existing_experiment = experiment.get_by_name(db, experiment_data.name, user_id)
        if existing_experiment:
            logger.warning(f"Experiment with name {experiment_data.name} already exists for user {user_id}")
            raise ValidationException(
                f"Experiment with name '{experiment_data.name}' already exists",
                {"name": experiment_data.name}
            )

        # Set default status to DRAFT if not provided
        if not experiment_data.status:
            experiment_data.status = ExperimentStatus.DRAFT.value

        db_experiment = experiment.create(db, obj_in=experiment_data, user_id=user_id)

        experiment_data = format_experiment_data(db_experiment)
        logger.info(f"Created experiment with ID {db_experiment.id}")
        return experiment_data


def update_experiment(experiment_id: UUID, experiment_data: ExperimentUpdate, user_id: int) -> Dict[str, Any]:
    """
    Updates an existing experiment's data, parameters, and molecules
    """
    logger.info(f"Updating experiment with ID {experiment_id} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment
        db_experiment = experiment.get_with_parameters(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # If experiment status is not DRAFT, only allow status updates
        if db_experiment.status != ExperimentStatus.DRAFT and experiment_data.status is None:
            logger.warning(f"Experiment {experiment_id} is not in DRAFT status, only status updates allowed")
            raise ValidationException(
                f"Experiment is not in DRAFT status, only status updates allowed",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        db_experiment = experiment.update(db, db_obj=db_experiment, obj_in=experiment_data)

        experiment_data = format_experiment_data(db_experiment)
        logger.info(f"Updated experiment with ID {experiment_id}")
        return experiment_data


def delete_experiment(experiment_id: UUID, user_id: int) -> Dict[str, Any]:
    """
    Deletes an experiment if it's in DRAFT status
    """
    logger.info(f"Deleting experiment with ID {experiment_id} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment
        db_experiment = experiment.get(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to delete experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to delete this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Verify experiment status is DRAFT (only draft experiments can be deleted)
        if db_experiment.status != ExperimentStatus.DRAFT:
            logger.warning(f"Experiment {experiment_id} is not in DRAFT status, cannot be deleted")
            raise ValidationException(
                f"Experiment is not in DRAFT status, cannot be deleted",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        experiment.remove(db, id=experiment_id)
        logger.info(f"Deleted experiment with ID {experiment_id}")
        return {"message": "Experiment deleted successfully"}


def update_experiment_status(experiment_id: UUID, status: str, user_id: int) -> Dict[str, Any]:
    """
    Updates the status of an experiment and sends notifications
    """
    logger.info(f"Updating status of experiment with ID {experiment_id} to {status} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment
        db_experiment = experiment.get(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Validate status transition is allowed based on current status
        is_valid_transition = validate_status_transition(db_experiment.status, status)
        if not is_valid_transition:
            logger.warning(f"Invalid status transition from {db_experiment.status} to {status}")
            raise ValidationException(
                f"Invalid status transition from {db_experiment.status} to {status}",
                {"experiment_id": experiment_id, "current_status": db_experiment.status, "new_status": status}
            )

        db_experiment = experiment.update_status(db, experiment_id, status)

        # Send experiment status change notification
        send_experiment_status_change(user_id, experiment_id, status, db_experiment.name)

        experiment_data = format_experiment_data(db_experiment)
        logger.info(f"Updated status of experiment with ID {experiment_id} to {status}")
        return experiment_data


def add_molecules_to_experiment(experiment_id: UUID, molecule_ids: List[UUID], user_id: int) -> Dict[str, Any]:
    """
    Adds molecules to an experiment
    """
    logger.info(f"Adding molecules to experiment with ID {experiment_id} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment
        db_experiment = experiment.get(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Verify experiment status is DRAFT or QUEUED (can only add molecules in these states)
        if db_experiment.status not in [ExperimentStatus.DRAFT, ExperimentStatus.QUEUED]:
            logger.warning(f"Experiment {experiment_id} is not in DRAFT or QUEUED status, cannot add molecules")
            raise ValidationException(
                f"Experiment is not in DRAFT or QUEUED status, cannot add molecules",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        result = experiment.add_molecules(db, experiment_id, molecule_ids)
        logger.info(f"Added molecules to experiment with ID {experiment_id}")
        return result


def remove_molecules_from_experiment(experiment_id: UUID, molecule_ids: List[UUID], user_id: int) -> Dict[str, Any]:
    """
    Removes molecules from an experiment
    """
    logger.info(f"Removing molecules from experiment with ID {experiment_id} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment
        db_experiment = experiment.get(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Verify experiment status is DRAFT or QUEUED (can only remove molecules in these states)
        if db_experiment.status not in [ExperimentStatus.DRAFT, ExperimentStatus.QUEUED]:
            logger.warning(f"Experiment {experiment_id} is not in DRAFT or QUEUED status, cannot remove molecules")
            raise ValidationException(
                f"Experiment is not in DRAFT or QUEUED status, cannot remove molecules",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        result = experiment.remove_molecules(db, experiment_id, molecule_ids)
        logger.info(f"Removed molecules from experiment with ID {experiment_id}")
        return result


def get_experiment_molecules(experiment_id: UUID, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves molecules associated with an experiment
    """
    logger.debug(f"Getting molecules for experiment with ID {experiment_id} (skip={skip}, limit={limit})")
    with get_db() as db:
        db_experiment = experiment.get(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        molecules = experiment.get_molecules(db, experiment_id, skip, limit)
        total = experiment.count_molecules(db, experiment_id)

        result = {
            "items": molecules,
            "total": total,
            "page": (skip // limit) + 1,
            "size": limit
        }
        logger.debug(f"Returning {len(molecules)} molecules for experiment {experiment_id}")
        return result


def get_experiments(filters: ExperimentFilter, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves experiments based on filter criteria
    """
    logger.debug(f"Getting experiments with filters: {filters} (skip={skip}, limit={limit})")
    with get_db() as db:
        # Convert filter parameters to dictionary if needed
        filter_params = filters.dict(exclude_unset=True) if filters else {}

        experiments, total = experiment.filter_experiments(db, filter_params, skip, limit)

        # Format each experiment as a dictionary with related data
        formatted_experiments = [format_experiment_data(exp) for exp in experiments]

        result = {
            "items": formatted_experiments,
            "total": total,
            "page": (skip // limit) + 1,
            "size": limit
        }
        logger.debug(f"Returning {len(experiments)} experiments")
        return result


def get_experiments_by_user(user_id: int, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves experiments created by a specific user
    """
    logger.debug(f"Getting experiments for user {user_id} (skip={skip}, limit={limit})")
    with get_db() as db:
        experiments = experiment.get_by_user(db, user_id, skip, limit)
        total = experiment.count_by_user(db, user_id)

        formatted_experiments = [format_experiment_data(exp) for exp in experiments]

        result = {
            "items": formatted_experiments,
            "total": total,
            "page": (skip // limit) + 1,
            "size": limit
        }
        logger.debug(f"Returning {len(experiments)} experiments for user {user_id}")
        return result


def get_experiments_by_status(status: str, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves experiments with a specific status
    """
    logger.debug(f"Getting experiments with status {status} (skip={skip}, limit={limit})")
    with get_db() as db:
        experiments = experiment.get_by_status(db, status, skip, limit)
        total = experiment.count_by_status(db, status)

        formatted_experiments = [format_experiment_data(exp) for exp in experiments]

        result = {
            "items": formatted_experiments,
            "total": total,
            "page": (skip // limit) + 1,
            "size": limit
        }
        logger.debug(f"Returning {len(experiments)} experiments with status {status}")
        return result


def queue_experiment(experiment_id: UUID, user_id: int) -> Dict[str, Any]:
    """
    Changes experiment status from DRAFT to QUEUED
    """
    logger.info(f"Queuing experiment with ID {experiment_id} for user {user_id}")
    with get_db() as db:
        # Get the existing experiment with molecules
        db_experiment = experiment.get_with_molecules(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Verify experiment status is DRAFT
        if db_experiment.status != ExperimentStatus.DRAFT:
            logger.warning(f"Experiment {experiment_id} is not in DRAFT status, cannot be queued")
            raise ValidationException(
                f"Experiment is not in DRAFT status, cannot be queued",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        # Verify experiment has at least one molecule associated
        if not db_experiment.molecules:
            logger.warning(f"Experiment {experiment_id} has no molecules, cannot be queued")
            raise ValidationException(
                f"Experiment has no molecules, cannot be queued",
                {"experiment_id": experiment_id}
            )

        # Update status to QUEUED
        db_experiment = experiment.update_status(db, experiment_id, ExperimentStatus.QUEUED)

        # Send experiment status change notification
        send_experiment_status_change(user_id, experiment_id, ExperimentStatus.QUEUED, db_experiment.name)

        experiment_data = format_experiment_data(db_experiment)
        logger.info(f"Queued experiment with ID {experiment_id}")
        return experiment_data


def prepare_for_submission(experiment_id: UUID, user_id: int) -> Dict[str, Any]:
    """
    Validates that an experiment is ready for submission to a CRO
    """
    logger.info(f"Preparing experiment with ID {experiment_id} for submission for user {user_id}")
    with get_db() as db:
        # Get the existing experiment with molecules
        db_experiment = experiment.get_with_molecules(db, experiment_id)
        if not db_experiment:
            logger.warning(f"Experiment with ID {experiment_id} not found")
            raise ResourceNotFoundException(
                f"Experiment with ID {experiment_id} not found",
                {"experiment_id": experiment_id}
            )

        # Verify experiment belongs to the specified user
        if db_experiment.created_by != user_id:
            logger.warning(f"User {user_id} does not have permission to update experiment {experiment_id}")
            raise AuthorizationException(
                f"User does not have permission to update this experiment",
                {"experiment_id": experiment_id, "user_id": user_id}
            )

        # Verify experiment status is QUEUED (only queued experiments can be submitted)
        if db_experiment.status != ExperimentStatus.QUEUED:
            logger.warning(f"Experiment {experiment_id} is not in QUEUED status, cannot be submitted")
            raise ValidationException(
                f"Experiment is not in QUEUED status, cannot be submitted",
                {"experiment_id": experiment_id, "status": db_experiment.status}
            )

        # Verify experiment has at least one molecule associated
        if not db_experiment.molecules:
            logger.warning(f"Experiment {experiment_id} has no molecules, cannot be submitted")
            raise ValidationException(
                f"Experiment has no molecules, cannot be submitted",
                {"experiment_id": experiment_id}
            )

        # Verify experiment has required parameters set
        if not db_experiment.parameters:
            logger.warning(f"Experiment {experiment_id} has no parameters, cannot be submitted")
            raise ValidationException(
                f"Experiment has no parameters, cannot be submitted",
                {"experiment_id": experiment_id}
            )

        experiment_data = format_experiment_data(db_experiment)
        logger.info(f"Experiment with ID {experiment_id} is ready for submission")
        return experiment_data


def validate_status_transition(current_status: str, new_status: str) -> bool:
    """
    Validates if a status transition is allowed based on current status
    """
    logger.debug(f"Validating status transition from {current_status} to {new_status}")
    # Define allowed transitions
    allowed_transitions = {
        ExperimentStatus.DRAFT: [ExperimentStatus.QUEUED, ExperimentStatus.CANCELLED],
        ExperimentStatus.QUEUED: [ExperimentStatus.SUBMITTED, ExperimentStatus.DRAFT, ExperimentStatus.CANCELLED],
        ExperimentStatus.SUBMITTED: [ExperimentStatus.QUOTE_PENDING, ExperimentStatus.REJECTED, ExperimentStatus.CANCELLED],
        ExperimentStatus.QUOTE_PENDING: [ExperimentStatus.IN_PROGRESS, ExperimentStatus.QUOTE_REJECTED, ExperimentStatus.CANCELLED],
        ExperimentStatus.QUOTE_REJECTED: [ExperimentStatus.QUOTE_PENDING, ExperimentStatus.CANCELLED],
        ExperimentStatus.IN_PROGRESS: [ExperimentStatus.RESULTS_PENDING, ExperimentStatus.CANCELLED],
        ExperimentStatus.RESULTS_PENDING: [ExperimentStatus.RESULTS_AVAILABLE, ExperimentStatus.RESULTS_REJECTED, ExperimentStatus.CANCELLED],
        ExperimentStatus.RESULTS_AVAILABLE: [ExperimentStatus.COMPLETED, ExperimentStatus.RESULTS_REJECTED, ExperimentStatus.CANCELLED],
        ExperimentStatus.RESULTS_REJECTED: [ExperimentStatus.RESULTS_PENDING, ExperimentStatus.CANCELLED],
        ExperimentStatus.COMPLETED: [],
        ExperimentStatus.CANCELLED: []
    }

    # Check if new_status is in allowed transitions for current_status
    if current_status in allowed_transitions and new_status in allowed_transitions[current_status]:
        logger.debug(f"Status transition from {current_status} to {new_status} is valid")
        return True
    else:
        logger.warning(f"Status transition from {current_status} to {new_status} is invalid")
        return False


def format_experiment_data(experiment_obj: Experiment, include_molecules: bool = False) -> Dict[str, Any]:
    """
    Formats experiment data as a dictionary with related information
    """
    logger.debug(f"Formatting experiment data for experiment {experiment_obj.id}")
    # Create base dictionary with experiment properties
    experiment_data = {
        "id": experiment_obj.id,
        "name": experiment_obj.name,
        "type_id": experiment_obj.type_id,
        "status": experiment_obj.status,
        "created_by": experiment_obj.created_by,
        "created_at": experiment_obj.created_at,
        "updated_at": experiment_obj.updated_at,
        "description": experiment_obj.description,
    }

    # Add experiment type information
    experiment_data["experiment_type"] = {
        "id": experiment_obj.experiment_type.id,
        "name": experiment_obj.experiment_type.name,
        "description": experiment_obj.experiment_type.description,
        "category": experiment_obj.experiment_type.category
    }

    # Add creator information
    experiment_data["creator"] = {
        "id": experiment_obj.created_by_user.id,
        "email": experiment_obj.created_by_user.email,
        "role": experiment_obj.created_by_user.role
    }

    # Add parameters as list of dictionaries
    experiment_data["parameters"] = [
        {"parameter_name": param.parameter_name, "parameter_value": param.parameter_value}
        for param in experiment_obj.parameters
    ]

    # If include_molecules is True, add molecules data
    if include_molecules:
        experiment_data["molecules"] = [
            {"id": mol.id, "smiles": mol.smiles}
            for mol in experiment_obj.molecules
        ]

    # Add molecule count
    experiment_data["molecule_count"] = len(experiment_obj.molecules)

    return experiment_data