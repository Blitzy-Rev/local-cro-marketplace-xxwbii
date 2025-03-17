from typing import Dict, List, Optional, Any, Tuple, Union
from datetime import datetime
from uuid import UUID, uuid4
import logging
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select, update, delete, and_, or_, func, desc, asc

from .base import CRUDBase
from ..models.submission import Submission, SubmissionStatus
from ..models.submission_detail import SubmissionDetail
from ..models.experiment import Experiment, ExperimentStatus
from ..schemas.submission import SubmissionCreate, SubmissionUpdate

# Set up logging
logger = logging.getLogger(__name__)


class CRUDSubmission(CRUDBase[Submission, SubmissionCreate, SubmissionUpdate]):
    """
    CRUD operations for submission management with extended functionality for details, 
    status tracking, and filtering.
    
    This class extends the base CRUD operations with submission-specific functions
    for handling the complex workflows involved in submitting experiments to CROs,
    including detail management, status transitions, and relationship handling.
    """
    
    def __init__(self):
        """Initialize the CRUD object with Submission model."""
        super().__init__(Submission)
    
    def create(self, db: Session, obj_in: Union[SubmissionCreate, Dict[str, Any]]) -> Submission:
        """
        Create a new submission with details.
        
        Args:
            db: Database session
            obj_in: Submission data to create
        
        Returns:
            Created submission instance
        """
        # Extract details if present
        if isinstance(obj_in, dict):
            obj_data = obj_in.copy()
            details = obj_data.pop("details", [])
        else:
            obj_data = obj_in.dict()
            details = obj_data.pop("details", [])
        
        # Set defaults
        if "submitted_at" not in obj_data:
            obj_data["submitted_at"] = datetime.utcnow()
        if "status" not in obj_data or not obj_data["status"]:
            obj_data["status"] = SubmissionStatus.PENDING.name
        
        # Create submission object
        db_obj = Submission(**obj_data)
        db.add(db_obj)
        db.flush()  # Flush to get ID while still in transaction
        
        # Create detail records if provided
        if details:
            for detail in details:
                if isinstance(detail, dict):
                    detail_data = detail
                else:
                    detail_data = detail.dict()
                
                detail_obj = SubmissionDetail(
                    submission_id=db_obj.id,
                    detail_name=detail_data["detail_name"],
                    detail_value=detail_data["detail_value"]
                )
                db.add(detail_obj)
        
        # Update experiment status to SUBMITTED
        if db_obj.experiment_id:
            stmt = select(Experiment).where(Experiment.id == db_obj.experiment_id)
            experiment = db.execute(stmt).scalar_one_or_none()
            if experiment:
                experiment.status = ExperimentStatus.SUBMITTED.name
        
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def update(self, db: Session, db_obj: Submission, obj_in: Union[SubmissionUpdate, Dict[str, Any]]) -> Submission:
        """
        Update a submission and its details.
        
        Args:
            db: Database session
            db_obj: Existing submission object to update
            obj_in: Updated submission data
        
        Returns:
            Updated submission instance
        """
        # Extract details if present
        if isinstance(obj_in, dict):
            obj_data = obj_in.copy()
            details = obj_data.pop("details", None)
        else:
            obj_data = obj_in.dict(exclude_unset=True)
            details = obj_data.pop("details", None) if "details" in obj_in.__dict__ else None
        
        # Update submission object
        for field in obj_data:
            setattr(db_obj, field, obj_data[field])
        
        # Update details if provided
        if details is not None:
            # Remove existing details
            stmt = delete(SubmissionDetail).where(SubmissionDetail.submission_id == db_obj.id)
            db.execute(stmt)
            
            # Add new details
            for detail in details:
                if isinstance(detail, dict):
                    detail_data = detail
                else:
                    detail_data = detail.dict()
                
                detail_obj = SubmissionDetail(
                    submission_id=db_obj.id,
                    detail_name=detail_data["detail_name"],
                    detail_value=detail_data["detail_value"]
                )
                db.add(detail_obj)
        
        # Update timestamp
        db_obj.updated_at = datetime.utcnow()
        
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def get_with_details(self, db: Session, submission_id: UUID) -> Optional[Submission]:
        """
        Get a submission with its details.
        
        Args:
            db: Database session
            submission_id: ID of the submission to retrieve
            
        Returns:
            Submission with details if found, None otherwise
        """
        stmt = (
            select(Submission)
            .where(Submission.id == submission_id)
        )
        submission = db.execute(stmt).scalar_one_or_none()
        
        # Access details to ensure they're loaded
        if submission:
            _ = submission.details
            
        return submission
    
    def get_with_experiment(self, db: Session, submission_id: UUID) -> Optional[Submission]:
        """
        Get a submission with its associated experiment.
        
        Args:
            db: Database session
            submission_id: ID of the submission to retrieve
            
        Returns:
            Submission with experiment if found, None otherwise
        """
        stmt = (
            select(Submission)
            .where(Submission.id == submission_id)
        )
        submission = db.execute(stmt).scalar_one_or_none()
        
        # Access experiment to ensure it's loaded
        if submission:
            _ = submission.experiment
            
        return submission
    
    def get_with_cro(self, db: Session, submission_id: UUID) -> Optional[Submission]:
        """
        Get a submission with its associated CRO user.
        
        Args:
            db: Database session
            submission_id: ID of the submission to retrieve
            
        Returns:
            Submission with CRO user if found, None otherwise
        """
        stmt = (
            select(Submission)
            .where(Submission.id == submission_id)
        )
        submission = db.execute(stmt).scalar_one_or_none()
        
        # Access CRO to ensure it's loaded
        if submission:
            _ = submission.cro
            
        return submission
    
    def get_with_results(self, db: Session, submission_id: UUID) -> Optional[Submission]:
        """
        Get a submission with its associated results.
        
        Args:
            db: Database session
            submission_id: ID of the submission to retrieve
            
        Returns:
            Submission with results if found, None otherwise
        """
        stmt = (
            select(Submission)
            .where(Submission.id == submission_id)
        )
        submission = db.execute(stmt).scalar_one_or_none()
        
        # Access results to ensure they're loaded
        if submission:
            _ = submission.results
            
        return submission
    
    def get_with_all_relations(self, db: Session, submission_id: UUID) -> Optional[Submission]:
        """
        Get a submission with all its related data.
        
        Args:
            db: Database session
            submission_id: ID of the submission to retrieve
            
        Returns:
            Submission with all relations if found, None otherwise
        """
        stmt = (
            select(Submission)
            .where(Submission.id == submission_id)
        )
        submission = db.execute(stmt).scalar_one_or_none()
        
        # Access all related objects to ensure they're loaded
        if submission:
            _ = submission.experiment
            _ = submission.cro
            _ = submission.details
            _ = submission.results
            
        return submission
    
    def get_by_status(self, db: Session, status: str, skip: int = 0, limit: int = 100) -> List[Submission]:
        """
        Get submissions with a specific status.
        
        Args:
            db: Database session
            status: Status to filter by
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of submissions with the specified status
        """
        stmt = (
            select(Submission)
            .where(Submission.status == status)
            .offset(skip)
            .limit(limit)
        )
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_by_status(self, db: Session, status: str) -> int:
        """
        Count submissions with a specific status.
        
        Args:
            db: Database session
            status: Status to count
            
        Returns:
            Count of submissions with the specified status
        """
        stmt = (
            select(func.count())
            .select_from(Submission)
            .where(Submission.status == status)
        )
        return db.execute(stmt).scalar_one()
    
    def update_status(self, db: Session, submission_id: UUID, status: str) -> Optional[Submission]:
        """
        Update the status of a submission.
        
        Args:
            db: Database session
            submission_id: ID of the submission to update
            status: New status value
            
        Returns:
            Updated submission if found, None otherwise
        """
        submission = self.get(db, submission_id)
        if not submission:
            return None
        
        submission.status = status
        submission.updated_at = datetime.utcnow()
        
        # Update experiment status based on submission status
        if submission.experiment_id:
            stmt = select(Experiment).where(Experiment.id == submission.experiment_id)
            experiment = db.execute(stmt).scalar_one_or_none()
            if experiment:
                # Map submission status to experiment status
                if status == SubmissionStatus.QUOTE_PROVIDED.name:
                    experiment.status = ExperimentStatus.QUOTE_PENDING.name
                elif status == SubmissionStatus.APPROVED.name:
                    experiment.status = ExperimentStatus.IN_PROGRESS.name
                elif status == SubmissionStatus.QUOTE_REJECTED.name:
                    experiment.status = ExperimentStatus.QUOTE_REJECTED.name
                elif status == SubmissionStatus.COMPLETED.name:
                    experiment.status = ExperimentStatus.RESULTS_AVAILABLE.name
                elif status == SubmissionStatus.CANCELLED.name:
                    experiment.status = ExperimentStatus.CANCELLED.name
        
        db.add(submission)
        db.commit()
        db.refresh(submission)
        return submission
    
    def get_by_experiment(self, db: Session, experiment_id: UUID, skip: int = 0, limit: int = 100) -> List[Submission]:
        """
        Get submissions for a specific experiment.
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of submissions for the experiment
        """
        stmt = (
            select(Submission)
            .where(Submission.experiment_id == experiment_id)
            .offset(skip)
            .limit(limit)
        )
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_by_experiment(self, db: Session, experiment_id: UUID) -> int:
        """
        Count submissions for a specific experiment.
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            
        Returns:
            Count of submissions for the experiment
        """
        stmt = (
            select(func.count())
            .select_from(Submission)
            .where(Submission.experiment_id == experiment_id)
        )
        return db.execute(stmt).scalar_one()
    
    def get_by_cro(self, db: Session, cro_id: UUID, skip: int = 0, limit: int = 100) -> List[Submission]:
        """
        Get submissions assigned to a specific CRO.
        
        Args:
            db: Database session
            cro_id: ID of the CRO user
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of submissions assigned to the CRO
        """
        stmt = (
            select(Submission)
            .where(Submission.cro_id == cro_id)
            .offset(skip)
            .limit(limit)
        )
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_by_cro(self, db: Session, cro_id: UUID) -> int:
        """
        Count submissions assigned to a specific CRO.
        
        Args:
            db: Database session
            cro_id: ID of the CRO user
            
        Returns:
            Count of submissions assigned to the CRO
        """
        stmt = (
            select(func.count())
            .select_from(Submission)
            .where(Submission.cro_id == cro_id)
        )
        return db.execute(stmt).scalar_one()
    
    def add_detail(self, db: Session, submission_id: UUID, detail_name: str, detail_value: str) -> SubmissionDetail:
        """
        Add a detail to a submission.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            detail_name: Name of the detail
            detail_value: Value of the detail
            
        Returns:
            Created detail object
        """
        # Check if detail already exists
        stmt = (
            select(SubmissionDetail)
            .where(
                SubmissionDetail.submission_id == submission_id,
                SubmissionDetail.detail_name == detail_name
            )
        )
        existing_detail = db.execute(stmt).scalar_one_or_none()
        
        if existing_detail:
            # Update existing detail
            existing_detail.detail_value = detail_value
            db.add(existing_detail)
            db.commit()
            db.refresh(existing_detail)
            return existing_detail
        else:
            # Create new detail
            detail_obj = SubmissionDetail(
                submission_id=submission_id,
                detail_name=detail_name,
                detail_value=detail_value
            )
            db.add(detail_obj)
            db.commit()
            db.refresh(detail_obj)
            return detail_obj
    
    def get_detail(self, db: Session, submission_id: UUID, detail_name: str) -> Optional[str]:
        """
        Get a specific detail from a submission.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            detail_name: Name of the detail to retrieve
            
        Returns:
            Detail value if found, None otherwise
        """
        stmt = (
            select(SubmissionDetail.detail_value)
            .where(
                SubmissionDetail.submission_id == submission_id,
                SubmissionDetail.detail_name == detail_name
            )
        )
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def get_details(self, db: Session, submission_id: UUID) -> Dict[str, str]:
        """
        Get all details for a submission.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            
        Returns:
            Dictionary of detail name-value pairs
        """
        stmt = (
            select(SubmissionDetail.detail_name, SubmissionDetail.detail_value)
            .where(SubmissionDetail.submission_id == submission_id)
        )
        results = db.execute(stmt).all()
        return {detail.detail_name: detail.detail_value for detail in results}
    
    def remove_detail(self, db: Session, submission_id: UUID, detail_name: str) -> bool:
        """
        Remove a detail from a submission.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            detail_name: Name of the detail to remove
            
        Returns:
            True if detail was removed, False if not found
        """
        stmt = (
            select(SubmissionDetail)
            .where(
                SubmissionDetail.submission_id == submission_id,
                SubmissionDetail.detail_name == detail_name
            )
        )
        detail = db.execute(stmt).scalar_one_or_none()
        
        if not detail:
            return False
        
        db.delete(detail)
        db.commit()
        return True
    
    def provide_quote(
        self, db: Session, submission_id: UUID, price: float, turnaround_time: int, 
        currency: Optional[str] = "USD", notes: Optional[str] = None
    ) -> Optional[Submission]:
        """
        Update a submission with quote information.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            price: Quoted price value
            turnaround_time: Estimated turnaround time in days
            currency: Currency code (default: USD)
            notes: Additional notes about the quote
            
        Returns:
            Updated submission if found, None otherwise
        """
        submission = self.get(db, submission_id)
        if not submission:
            return None
        
        # Update submission status and fields
        submission.status = SubmissionStatus.QUOTE_PROVIDED.name
        submission.price = price
        submission.turnaround_days = turnaround_time
        submission.updated_at = datetime.utcnow()
        
        # Add or update detail records
        self.add_detail(db, submission_id, "price", str(price))
        self.add_detail(db, submission_id, "currency", currency)
        self.add_detail(db, submission_id, "turnaround_days", str(turnaround_time))
        if notes:
            self.add_detail(db, submission_id, "quote_notes", notes)
        
        # Update experiment status
        stmt = select(Experiment).where(Experiment.id == submission.experiment_id)
        experiment = db.execute(stmt).scalar_one_or_none()
        if experiment:
            experiment.status = ExperimentStatus.QUOTE_PENDING.name
        
        db.add(submission)
        db.commit()
        db.refresh(submission)
        return submission
    
    def respond_to_quote(
        self, db: Session, submission_id: UUID, approved: bool, notes: Optional[str] = None
    ) -> Optional[Submission]:
        """
        Respond to a quote with approval or rejection.
        
        Args:
            db: Database session
            submission_id: ID of the submission
            approved: Whether the quote is approved
            notes: Additional notes about the decision
            
        Returns:
            Updated submission if found, None otherwise
        """
        submission = self.get(db, submission_id)
        if not submission:
            return None
        
        if approved:
            submission.status = SubmissionStatus.APPROVED.name
            # Update experiment status
            stmt = select(Experiment).where(Experiment.id == submission.experiment_id)
            experiment = db.execute(stmt).scalar_one_or_none()
            if experiment:
                experiment.status = ExperimentStatus.IN_PROGRESS.name
        else:
            submission.status = SubmissionStatus.QUOTE_REJECTED.name
            # Update experiment status
            stmt = select(Experiment).where(Experiment.id == submission.experiment_id)
            experiment = db.execute(stmt).scalar_one_or_none()
            if experiment:
                experiment.status = ExperimentStatus.QUOTE_REJECTED.name
        
        submission.updated_at = datetime.utcnow()
        
        # Add response notes if provided
        if notes:
            self.add_detail(db, submission_id, "response_notes", notes)
        
        db.add(submission)
        db.commit()
        db.refresh(submission)
        return submission
    
    def filter_submissions(
        self, db: Session, filters: Dict[str, Any], skip: int = 0, limit: int = 100
    ) -> Tuple[List[Submission], int]:
        """
        Filter submissions based on various criteria.
        
        Args:
            db: Database session
            filters: Dictionary of filter parameters
            skip: Number of records to skip (pagination)
            limit: Maximum number of records to return
            
        Returns:
            Tuple of (list of filtered submissions, total count)
        """
        query = select(Submission)
        
        # Apply filters
        if "experiment_id" in filters and filters["experiment_id"]:
            query = query.where(Submission.experiment_id == filters["experiment_id"])
        
        if "cro_id" in filters and filters["cro_id"]:
            query = query.where(Submission.cro_id == filters["cro_id"])
        
        if "status" in filters and filters["status"]:
            query = query.where(Submission.status == filters["status"])
        
        if "submitted_after" in filters and filters["submitted_after"]:
            query = query.where(Submission.submitted_at >= filters["submitted_after"])
        
        if "submitted_before" in filters and filters["submitted_before"]:
            query = query.where(Submission.submitted_at <= filters["submitted_before"])
        
        # Get total count before pagination
        count_query = select(func.count()).select_from(query.subquery())
        total = db.execute(count_query).scalar_one()
        
        # Apply sorting
        if "sort_by" in filters and filters["sort_by"]:
            sort_field = filters["sort_by"]
            if hasattr(Submission, sort_field):
                sort_column = getattr(Submission, sort_field)
                if "sort_desc" in filters and filters["sort_desc"]:
                    query = query.order_by(desc(sort_column))
                else:
                    query = query.order_by(asc(sort_column))
            else:
                # Default to submitted_at if invalid field
                query = query.order_by(desc(Submission.submitted_at))
        else:
            # Default sort
            query = query.order_by(desc(Submission.submitted_at))
        
        # Apply pagination
        query = query.offset(skip).limit(limit)
        
        # Execute query
        results = db.execute(query).scalars().all()
        
        return list(results), total


# Create singleton instance
submission = CRUDSubmission()