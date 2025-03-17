from typing import List, Dict, Any, Optional, Union, Tuple
from uuid import UUID
from datetime import datetime
from sqlalchemy.orm import Session
from sqlalchemy import select, join, and_, or_

from .base import CRUDBase
from ..models.result import Result, ResultStatus
from ..models.result_file import ResultFile
from ..models.result_data import ResultData
from ..schemas.result import ResultCreate, ResultUpdate


class CRUDResult(CRUDBase[Result, ResultCreate, ResultUpdate]):
    """
    CRUD operations for the Result model with specialized methods for result management
    
    This class extends the base CRUD operations with additional methods for managing
    experimental results, their associated files, and structured data. It provides
    comprehensive data access capabilities for the CRO result management workflow.
    """
    
    def __init__(self):
        """Initialize the CRUD object with the Result model"""
        super().__init__(Result)
    
    def get_with_details(self, db: Session, result_id: UUID) -> Optional[Result]:
        """
        Get a result by ID with all related details (submission, files, data)
        
        Args:
            db: Database session
            result_id: ID of the result to retrieve
            
        Returns:
            Result with loaded relationships if found, None otherwise
        """
        query = (
            select(Result)
            .where(Result.id == result_id)
            .join(Result.submission)
            .join(Result.files, isouter=True)
            .join(Result.data_points, isouter=True)
        )
        
        result = db.execute(query).unique().scalar_one_or_none()
        return result
    
    def get_with_files(self, db: Session, result_id: UUID) -> Optional[Result]:
        """
        Get a result by ID with its associated files
        
        Args:
            db: Database session
            result_id: ID of the result to retrieve
            
        Returns:
            Result with loaded file relationships if found, None otherwise
        """
        query = (
            select(Result)
            .where(Result.id == result_id)
            .join(Result.files, isouter=True)
        )
        
        result = db.execute(query).unique().scalar_one_or_none()
        return result
    
    def get_with_data(self, db: Session, result_id: UUID) -> Optional[Result]:
        """
        Get a result by ID with its associated structured data
        
        Args:
            db: Database session
            result_id: ID of the result to retrieve
            
        Returns:
            Result with loaded data relationships if found, None otherwise
        """
        query = (
            select(Result)
            .where(Result.id == result_id)
            .join(Result.data_points, isouter=True)
        )
        
        result = db.execute(query).unique().scalar_one_or_none()
        return result
    
    def get_by_submission(self, db: Session, submission_id: UUID, skip: int = 0, limit: int = 100) -> List[Result]:
        """
        Get results for a specific submission with pagination
        
        Args:
            db: Database session
            submission_id: ID of the submission to get results for
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of results for the submission
        """
        query = (
            select(Result)
            .where(Result.submission_id == submission_id)
            .offset(skip)
            .limit(limit)
        )
        
        results = db.execute(query).scalars().all()
        return list(results)
    
    def get_by_status(self, db: Session, status: str, skip: int = 0, limit: int = 100) -> List[Result]:
        """
        Get results with a specific status with pagination
        
        Args:
            db: Database session
            status: Status to filter by (one of ResultStatus enum values)
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of results with the specified status
        """
        query = (
            select(Result)
            .where(Result.status == status)
            .offset(skip)
            .limit(limit)
        )
        
        results = db.execute(query).scalars().all()
        return list(results)
    
    def update_status(
        self, db: Session, result_id: UUID, status: str, approved_at: Optional[datetime] = None
    ) -> Optional[Result]:
        """
        Update the status of a result
        
        Args:
            db: Database session
            result_id: ID of the result to update
            status: New status value (one of ResultStatus enum values)
            approved_at: Timestamp for approval (if status is APPROVED)
            
        Returns:
            Updated result if found, None otherwise
        """
        result = self.get(db, result_id)
        if not result:
            return None
        
        # Update the status
        result.status = status
        
        # If the status is being set to APPROVED and no approved_at time is provided,
        # set it to the current time
        if status == ResultStatus.APPROVED and approved_at is None:
            result.approved_at = datetime.utcnow()
        elif approved_at is not None:
            result.approved_at = approved_at
        
        db.add(result)
        db.commit()
        db.refresh(result)
        return result
    
    def add_file(
        self, db: Session, result_id: UUID, file_name: str, file_path: str, file_size: int, file_type: str
    ) -> ResultFile:
        """
        Add a file to a result
        
        Args:
            db: Database session
            result_id: ID of the result this file belongs to
            file_name: Original name of the uploaded file
            file_path: Path to the file in the object storage system
            file_size: Size of the file in bytes
            file_type: MIME type or file extension
            
        Returns:
            Created result file instance
        """
        file = ResultFile(
            result_id=result_id,
            file_name=file_name,
            file_path=file_path,
            file_size=file_size,
            file_type=file_type,
            uploaded_at=datetime.utcnow()
        )
        
        db.add(file)
        db.commit()
        db.refresh(file)
        return file
    
    def get_file(self, db: Session, file_id: UUID) -> Optional[ResultFile]:
        """
        Get a result file by ID
        
        Args:
            db: Database session
            file_id: ID of the file to retrieve
            
        Returns:
            Result file if found, None otherwise
        """
        query = select(ResultFile).where(ResultFile.id == file_id)
        file = db.execute(query).scalar_one_or_none()
        return file
    
    def delete_file(self, db: Session, file_id: UUID) -> bool:
        """
        Delete a result file
        
        Args:
            db: Database session
            file_id: ID of the file to delete
            
        Returns:
            True if file was deleted, False otherwise
        """
        file = self.get_file(db, file_id)
        if not file:
            return False
        
        db.delete(file)
        db.commit()
        return True
    
    def add_data(
        self, db: Session, result_id: UUID, molecule_id: UUID, data_name: str, 
        data_value: float, data_unit: Optional[str] = None
    ) -> ResultData:
        """
        Add structured data to a result for a specific molecule
        
        Args:
            db: Database session
            result_id: ID of the result this data belongs to
            molecule_id: ID of the molecule this data is associated with
            data_name: Name or identifier of the data point (e.g., "IC50", "Solubility")
            data_value: Numerical value of the measurement
            data_unit: Unit of measurement (e.g., "nM", "mg/mL")
            
        Returns:
            Created result data instance
        """
        data = ResultData(
            result_id=result_id,
            molecule_id=molecule_id,
            data_name=data_name,
            data_value=data_value,
            data_unit=data_unit
        )
        
        db.add(data)
        db.commit()
        db.refresh(data)
        return data
    
    def get_data(self, db: Session, result_id: UUID) -> List[ResultData]:
        """
        Get all structured data for a result
        
        Args:
            db: Database session
            result_id: ID of the result to get data for
            
        Returns:
            List of structured data records for the result
        """
        query = (
            select(ResultData)
            .where(ResultData.result_id == result_id)
            .join(ResultData.molecule)
        )
        
        data = db.execute(query).scalars().all()
        return list(data)
    
    def get_data_by_molecule(self, db: Session, result_id: UUID, molecule_id: UUID) -> List[ResultData]:
        """
        Get structured data for a result filtered by molecule
        
        Args:
            db: Database session
            result_id: ID of the result to get data for
            molecule_id: ID of the molecule to filter by
            
        Returns:
            List of structured data records for the result and molecule
        """
        query = (
            select(ResultData)
            .where(and_(
                ResultData.result_id == result_id,
                ResultData.molecule_id == molecule_id
            ))
            .join(ResultData.molecule)
        )
        
        data = db.execute(query).scalars().all()
        return list(data)
    
    def filter_results(self, db: Session, filters: Dict[str, Any]) -> Tuple[List[Result], int]:
        """
        Filter results based on various criteria with pagination
        
        Args:
            db: Database session
            filters: Dictionary of filter criteria, which may include:
                - submission_id: Filter by submission ID
                - experiment_id: Filter by experiment ID via submission relationship
                - status: Filter by result status
                - uploaded_after: Filter by upload date greater than or equal to
                - uploaded_before: Filter by upload date less than or equal to
                - user_id: Filter by submission's experiment creator
                - skip: Number of records to skip (pagination)
                - limit: Maximum number of records to return
                - sort_by: Field to sort by
                - sort_desc: Sort in descending order if true
            
        Returns:
            Tuple containing (list of filtered results, total count)
        """
        # Build base query
        query = select(Result)
        
        # Apply filters
        filter_conditions = []
        
        if filters.get("submission_id"):
            filter_conditions.append(Result.submission_id == filters["submission_id"])
        
        if filters.get("experiment_id"):
            # This requires joining on the submission relationship
            query = query.join(Result.submission)
            filter_conditions.append(Result.submission.has(experiment_id=filters["experiment_id"]))
        
        if filters.get("status"):
            filter_conditions.append(Result.status == filters["status"])
        
        if filters.get("uploaded_after"):
            filter_conditions.append(Result.uploaded_at >= filters["uploaded_after"])
        
        if filters.get("uploaded_before"):
            filter_conditions.append(Result.uploaded_at <= filters["uploaded_before"])
        
        if filters.get("user_id"):
            # This requires joining through submission to experiment to user
            query = query.join(Result.submission).join(Submission.experiment)
            filter_conditions.append(Experiment.created_by == filters["user_id"])
        
        # Apply all filters
        if filter_conditions:
            query = query.where(and_(*filter_conditions))
        
        # Get total count
        count_query = query.with_only_columns(func.count())
        total = db.execute(count_query).scalar_one()
        
        # Apply sorting
        sort_by = filters.get("sort_by", "uploaded_at")
        sort_desc = filters.get("sort_desc", True)
        
        if hasattr(Result, sort_by):
            sort_field = getattr(Result, sort_by)
            if sort_desc:
                query = query.order_by(sort_field.desc())
            else:
                query = query.order_by(sort_field)
        
        # Apply pagination
        skip = filters.get("skip", 0)
        limit = filters.get("limit", 100)
        query = query.offset(skip).limit(limit)
        
        # Execute query
        results = db.execute(query).scalars().all()
        
        return list(results), total


# Create singleton instance
result = CRUDResult()