"""
Service layer implementation for CRO submission management in the Molecular Data Management and CRO Integration Platform.
This service provides high-level business logic for submission operations, including creation, status management, quote handling, and integration with experiments and results. It acts as an intermediary between API endpoints and the data access layer, implementing core submission workflow functionality.
"""

import logging  # standard library
from typing import List, Dict, Optional, Any, Tuple, Union  # standard library
from uuid import UUID  # standard library

from ..crud import crud_submission as submission  # Internal import: Data access layer for submission operations
from ..models.submission import Submission  # Internal import: Database model for submissions
from ..constants import SubmissionStatus  # Internal import: Enumeration of submission status values
from ..schemas.submission import SubmissionCreate, SubmissionUpdate, SubmissionFilter, QuoteProvide, QuoteResponse  # Internal import: Schemas for submission data validation
from ..db.session import get_db  # Internal import: Database session context manager
from .experiment_service import experiment_service  # Internal import: Service for experiment operations
from .notification_service import send_submission_created, send_quote_provided  # Internal import: Send notifications for submission events
from ..exceptions import ValidationException, ResourceNotFoundException, SubmissionException  # Internal import: Exception classes

# Set up logger
logger = logging.getLogger(__name__)


class SubmissionService:
    """
    Service class for submission management operations
    """

    def __init__(self):
        """
        Initializes the SubmissionService
        """
        self.logger = logger

    def get_submission_by_id(self, submission_id: UUID, include_results: bool = False) -> Optional[Dict[str, Any]]:
        """
        Retrieves a submission by its ID with all details and related data
        """
        with get_db() as db:
            db_submission = submission.get_with_all_relations(db, submission_id) if include_results else submission.get_with_details(db, submission_id)
            if not db_submission:
                return None
            return self.format_submission_data(db_submission, include_results)

    def create_submission(self, submission_data: SubmissionCreate, user_id: int) -> Dict[str, Any]:
        """
        Creates a new submission for an experiment to a CRO
        """
        with get_db() as db:
            # Call experiment_service.prepare_for_submission to validate experiment is ready
            try:
                experiment_service.get_experiment_by_id(submission_data.experiment_id)
            except ResourceNotFoundException as e:
                raise ValidationException(f"Experiment not found: {submission_data.experiment_id}", details={"experiment_id": submission_data.experiment_id}) from e

            # Verify CRO user exists and has CRO role - omitted for local deployment
            db_submission = submission.create(db, obj_in=submission_data, user_id=user_id)

            # Update experiment status to SUBMITTED
            experiment_service.update_experiment_status(submission_data.experiment_id, SubmissionStatus.SUBMITTED.name, user_id)

            # Get CRO user details for notification - omitted for local deployment
            # Call send_submission_created to notify CRO of new submission
            send_submission_created(db_submission.cro_id, db_submission.id, db_submission.experiment.name, db_submission.experiment.created_by_user.email)

            return self.format_submission_data(db_submission)

    def update_submission(self, submission_id: UUID, submission_data: SubmissionUpdate, user_id: int) -> Dict[str, Any]:
        """
        Updates an existing submission's data
        """
        with get_db() as db:
            # Get the existing submission
            db_submission = submission.get_with_details(db, submission_id)
            if not db_submission:
                raise ResourceNotFoundException(f"Submission not found: {submission_id}", details={"submission_id": submission_id})

            # Verify user has permission to update this submission (creator or assigned CRO) - omitted for local deployment

            # Update the submission
            db_submission = submission.update(db, db_obj=db_submission, obj_in=submission_data)

            # If status changed, update experiment status accordingly
            if submission_data.status and db_submission.experiment_id:
                experiment_service.update_experiment_status(db_submission.experiment_id, db_submission.status, user_id)

            return self.format_submission_data(db_submission)

    def cancel_submission(self, submission_id: UUID, user_id: int) -> Dict[str, Any]:
        """
        Cancels a submission if it's in an appropriate state
        """
        with get_db() as db:
            # Get the existing submission
            db_submission = submission.get_with_experiment(db, submission_id)
            if not db_submission:
                raise ResourceNotFoundException(f"Submission not found: {submission_id}", details={"submission_id": submission_id})

            # Verify user has permission to cancel this submission (creator or assigned CRO) - omitted for local deployment

            # Verify submission status allows cancellation (not COMPLETED or already CANCELLED)
            if db_submission.status in [SubmissionStatus.COMPLETED, SubmissionStatus.CANCELLED]:
                raise ValidationException(f"Submission cannot be cancelled in its current state: {db_submission.status}", details={"submission_id": submission_id, "status": db_submission.status})

            # Update status to CANCELLED
            db_submission = submission.update_status(db, submission_id, SubmissionStatus.CANCELLED.name)

            # Update experiment status to reflect cancellation
            experiment_service.update_experiment_status(db_submission.experiment_id, SubmissionStatus.CANCELLED.name, user_id)

            return {"message": "Submission cancelled successfully"}

    def update_submission_status(self, submission_id: UUID, status: str, user_id: int) -> Dict[str, Any]:
        """
        Updates the status of a submission
        """
        with get_db() as db:
            # Get the existing submission
            db_submission = submission.get_with_experiment(db, submission_id)
            if not db_submission:
                raise ResourceNotFoundException(f"Submission not found: {submission_id}", details={"submission_id": submission_id})

            # Verify user has permission to update this submission (creator or assigned CRO) - omitted for local deployment

            # Validate status transition is allowed based on current status
            self.validate_status_transition(db_submission.status, status)

            # Update the submission status
            db_submission = submission.update_status(db, submission_id, status)

            # Update experiment status to reflect submission status change
            experiment_service.update_experiment_status(db_submission.experiment_id, status, user_id)

            return self.format_submission_data(db_submission)

    def provide_quote(self, submission_id: UUID, quote_data: QuoteProvide, cro_user_id: int) -> Dict[str, Any]:
        """
        Allows a CRO to provide a quote for a submission
        """
        with get_db() as db:
            # Get the existing submission
            db_submission = submission.get_with_experiment(db, submission_id)
            if not db_submission:
                raise ResourceNotFoundException(f"Submission not found: {submission_id}", details={"submission_id": submission_id})

            # Verify submission is assigned to the specified CRO user - omitted for local deployment

            # Verify submission status is PENDING (awaiting quote)
            if db_submission.status != SubmissionStatus.PENDING.name:
                raise ValidationException(f"Submission is not awaiting a quote: {db_submission.status}", details={"submission_id": submission_id, "status": db_submission.status})

            # Update the submission with the quote
            db_submission = submission.provide_quote(db, submission_id, quote_data.price, quote_data.turnaround_days, notes=quote_data.notes)

            # Get experiment creator details for notification - omitted for local deployment
            # Call send_quote_provided to notify pharma user of quote
            send_quote_provided(db_submission.experiment.created_by, db_submission.id, db_submission.experiment.name, db_submission.cro.email, quote_data.price)

            return self.format_submission_data(db_submission)

    def respond_to_quote(self, submission_id: UUID, response_data: QuoteResponse, user_id: int) -> Dict[str, Any]:
        """
        Allows a pharma user to approve or reject a quote
        """
        with get_db() as db:
            # Get the existing submission
            db_submission = submission.get_with_experiment(db, submission_id)
            if not db_submission:
                raise ResourceNotFoundException(f"Submission not found: {submission_id}", details={"submission_id": submission_id})

            # Verify user is the creator of the associated experiment - omitted for local deployment

            # Verify submission status is QUOTE_PROVIDED (awaiting response)
            if db_submission.status != SubmissionStatus.QUOTE_PROVIDED.name:
                raise ValidationException(f"Submission is not awaiting a response to the quote: {db_submission.status}", details={"submission_id": submission_id, "status": db_submission.status})

            # Update the submission with the response
            db_submission = submission.respond_to_quote(db, submission_id, response_data.approved, notes=response_data.notes)

            # Get CRO user details for notification - omitted for local deployment
            # Send appropriate notification based on approval decision - omitted for local deployment

            return self.format_submission_data(db_submission)

    def get_submissions(self, filters: SubmissionFilter, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves submissions based on filter criteria with pagination
        """
        with get_db() as db:
            # Convert filter parameters to dictionary if needed
            filter_params = filters.dict(exclude_unset=True) if filters else {}

            # Call submission.filter_submissions with filter parameters
            submissions, total = submission.filter_submissions(db, filter_params, skip, limit)

            # Format each submission as a dictionary with related data
            formatted_submissions = [self.format_submission_data(sub) for sub in submissions]

            # Return dictionary with items, total count, page number, and page size
            return {
                "items": formatted_submissions,
                "total": total,
                "page": (skip // limit) + 1,
                "size": limit
            }

    def get_submissions_by_experiment(self, experiment_id: UUID, user_id: int, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves submissions for a specific experiment with pagination
        """
        with get_db() as db:
            # Verify user has access to the experiment - omitted for local deployment

            # Call submission.get_by_experiment to get submissions for the experiment
            submissions = submission.get_by_experiment(db, experiment_id, skip, limit)

            # Calculate total count of experiment's submissions
            total = submission.count_by_experiment(db, experiment_id)

            # Format each submission as a dictionary
            formatted_submissions = [self.format_submission_data(sub) for sub in submissions]

            # Return dictionary with items, total count, page number, and page size
            return {
                "items": formatted_submissions,
                "total": total,
                "page": (skip // limit) + 1,
                "size": limit
            }

    def get_submissions_by_cro(self, cro_user_id: int, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves submissions assigned to a specific CRO with pagination
        """
        with get_db() as db:
            # Verify user has CRO role - omitted for local deployment

            # Call submission.get_by_cro to get submissions assigned to the CRO
            submissions = submission.get_by_cro(db, cro_user_id, skip, limit)

            # Calculate total count of CRO's submissions
            total = submission.count_by_cro(db, cro_user_id)

            # Format each submission as a dictionary
            formatted_submissions = [self.format_submission_data(sub) for sub in submissions]

            # Return dictionary with items, total count, page number, and page size
            return {
                "items": formatted_submissions,
                "total": total,
                "page": (skip // limit) + 1,
                "size": limit
            }

    def get_submissions_by_status(self, status: str, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
        """
        Retrieves submissions with a specific status with pagination
        """
        with get_db() as db:
            # Call submission.get_by_status to get submissions with the specified status
            submissions = submission.get_by_status(db, status, skip, limit)

            # Calculate total count of submissions with this status
            total = submission.count_by_status(db, status)

            # Format each submission as a dictionary
            formatted_submissions = [self.format_submission_data(sub) for sub in submissions]

            # Return dictionary with items, total count, page number, and page size
            return {
                "items": formatted_submissions,
                "total": total,
                "page": (skip // limit) + 1,
                "size": limit
            }

    def validate_status_transition(self, current_status: str, new_status: str) -> bool:
        """
        Validates if a status transition is allowed based on current status
        """
        # Define allowed transitions
        allowed_transitions = {
            SubmissionStatus.PENDING.name: [SubmissionStatus.QUOTE_PROVIDED.name, SubmissionStatus.REJECTED.name, SubmissionStatus.CANCELLED.name],
            SubmissionStatus.QUOTE_PROVIDED.name: [SubmissionStatus.APPROVED.name, SubmissionStatus.QUOTE_REJECTED.name, SubmissionStatus.CANCELLED.name],
            SubmissionStatus.APPROVED.name: [SubmissionStatus.IN_PROGRESS.name, SubmissionStatus.CANCELLED.name],
            SubmissionStatus.IN_PROGRESS.name: [SubmissionStatus.COMPLETED.name, SubmissionStatus.CANCELLED.name],
            SubmissionStatus.COMPLETED.name: [],
            SubmissionStatus.REJECTED.name: [],
            SubmissionStatus.QUOTE_REJECTED.name: [],
            SubmissionStatus.CANCELLED.name: []
        }

        # Check if new_status is in allowed transitions for current_status
        if current_status in allowed_transitions and new_status in allowed_transitions[current_status]:
            return True
        else:
            raise ValidationException(f"Invalid status transition from {current_status} to {new_status}", details={"current_status": current_status, "new_status": new_status})

    def format_submission_data(self, submission_obj: Submission, include_results: bool = False) -> Dict[str, Any]:
        """
        Formats submission data as a dictionary with related information
        """
        # Create base dictionary with submission properties
        submission_data = {
            "id": submission_obj.id,
            "experiment_id": submission_obj.experiment_id,
            "cro_id": submission_obj.cro_id,
            "status": submission_obj.status,
            "submitted_at": submission_obj.submitted_at,
            "updated_at": submission_obj.updated_at,
            "price": submission_obj.price,
            "turnaround_days": submission_obj.turnaround_days,
            "notes": submission_obj.notes
        }

        # Add experiment information
        submission_data["experiment"] = {
            "id": submission_obj.experiment.id,
            "name": submission_obj.experiment.name,
            "type": submission_obj.experiment.type_id
        }

        # Add CRO information
        submission_data["cro"] = {
            "id": submission_obj.cro.id,
            "email": submission_obj.cro.email
        }

        # Add details as list of dictionaries
        submission_data["details"] = [{"detail_name": detail.detail_name, "detail_value": detail.detail_value} for detail in submission_obj.details]

        # If include_results is True, add results data
        if include_results:
            submission_data["results"] = [{"id": result.id, "status": result.status} for result in submission_obj.results]

        return submission_data


submission_service = SubmissionService()