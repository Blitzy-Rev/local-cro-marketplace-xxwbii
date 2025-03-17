from typing import Dict, List, Optional, Any, Tuple, Union
from uuid import UUID, uuid4
from datetime import datetime
import logging
from sqlalchemy import select, update, delete, and_, or_, func, desc, asc
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+

from .base import CRUDBase
from ..models.experiment import Experiment, ExperimentStatus
from ..models.experiment_parameter import ExperimentParameter
from ..models.experiment_molecule import ExperimentMolecule
from ..models.molecule import Molecule
from ..schemas.experiment import ExperimentCreate, ExperimentUpdate

# Set up logger
logger = logging.getLogger(__name__)


class CRUDExperiment(CRUDBase[Experiment, ExperimentCreate, ExperimentUpdate]):
    """
    CRUD operations for experiment management with extended functionality for parameters,
    molecules, and filtering.
    """
    
    def __init__(self):
        """Initialize the CRUD object with Experiment model"""
        super().__init__(Experiment)
    
    def get_by_name(self, db: Session, name: str, user_id: int) -> Optional[Experiment]:
        """
        Get an experiment by name and creator ID
        
        Args:
            db: Database session
            name: Experiment name to search for
            user_id: ID of the creator user
            
        Returns:
            Experiment if found, None otherwise
        """
        stmt = select(self.model).where(
            and_(
                self.model.name == name,
                self.model.created_by == user_id
            )
        )
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def create(
        self, db: Session, obj_in: Union[ExperimentCreate, Dict[str, Any]], user_id: int
    ) -> Experiment:
        """
        Create a new experiment with parameters
        
        Args:
            db: Database session
            obj_in: Input data (schema or dictionary)
            user_id: ID of the creator user
            
        Returns:
            Created experiment instance
        """
        # Extract parameters if they exist
        if isinstance(obj_in, dict):
            parameters = obj_in.pop("parameters", [])
            obj_data = obj_in
        else:
            parameters = obj_in.parameters
            obj_data = obj_in.dict(exclude={"parameters", "molecule_ids"})
        
        # Create experiment object
        obj_data["created_by"] = user_id
        db_obj = self.model(**obj_data)
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        
        # Create parameter records if any
        if parameters:
            for param in parameters:
                if isinstance(param, dict):
                    param_data = param
                else:
                    param_data = param.dict()
                
                param_obj = ExperimentParameter(
                    experiment_id=db_obj.id,
                    parameter_name=param_data["parameter_name"],
                    parameter_value=param_data["parameter_value"]
                )
                db.add(param_obj)
        
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def update(
        self, db: Session, db_obj: Experiment, obj_in: Union[ExperimentUpdate, Dict[str, Any]]
    ) -> Experiment:
        """
        Update an experiment and its parameters
        
        Args:
            db: Database session
            db_obj: Existing experiment object
            obj_in: Update data (schema or dictionary)
            
        Returns:
            Updated experiment instance
        """
        # Extract parameters if they exist
        if isinstance(obj_in, dict):
            parameters = obj_in.pop("parameters", None)
            update_data = obj_in
        else:
            parameters = obj_in.parameters if hasattr(obj_in, "parameters") else None
            update_data = obj_in.dict(exclude_unset=True, exclude={"parameters", "molecule_ids"})
        
        # Update the experiment object with base fields
        for field in update_data:
            if hasattr(db_obj, field):
                setattr(db_obj, field, update_data[field])
        
        # Update parameters if provided
        if parameters is not None:
            # First remove existing parameters
            db.execute(delete(ExperimentParameter).where(
                ExperimentParameter.experiment_id == db_obj.id
            ))
            
            # Then add new parameters
            for param in parameters:
                if isinstance(param, dict):
                    param_data = param
                else:
                    param_data = param.dict()
                
                param_obj = ExperimentParameter(
                    experiment_id=db_obj.id,
                    parameter_name=param_data["parameter_name"],
                    parameter_value=param_data["parameter_value"]
                )
                db.add(param_obj)
        
        # Set updated timestamp
        db_obj.updated_at = datetime.utcnow()
        
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def get_with_parameters(self, db: Session, experiment_id: UUID) -> Optional[Experiment]:
        """
        Get an experiment with its parameters
        
        Args:
            db: Database session
            experiment_id: ID of the experiment to retrieve
            
        Returns:
            Experiment with parameters if found, None otherwise
        """
        stmt = select(self.model).where(self.model.id == experiment_id)
        result = db.execute(stmt).scalar_one_or_none()
        
        if result:
            # Ensure parameters are loaded
            param_stmt = select(ExperimentParameter).where(
                ExperimentParameter.experiment_id == experiment_id
            )
            result.parameters = db.execute(param_stmt).scalars().all()
        
        return result
    
    def get_with_molecules(self, db: Session, experiment_id: UUID) -> Optional[Experiment]:
        """
        Get an experiment with its associated molecules
        
        Args:
            db: Database session
            experiment_id: ID of the experiment to retrieve
            
        Returns:
            Experiment with molecules if found, None otherwise
        """
        stmt = select(self.model).where(self.model.id == experiment_id)
        result = db.execute(stmt).scalar_one_or_none()
        
        if result:
            # Ensure molecules are loaded
            molecules_stmt = select(Molecule).join(
                ExperimentMolecule, 
                and_(
                    ExperimentMolecule.molecule_id == Molecule.id,
                    ExperimentMolecule.experiment_id == experiment_id
                )
            )
            result.molecules = db.execute(molecules_stmt).scalars().all()
        
        return result
    
    def get_by_status(
        self, db: Session, status: str, skip: int = 0, limit: int = 100
    ) -> List[Experiment]:
        """
        Get experiments with a specific status
        
        Args:
            db: Database session
            status: Status to filter by
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of experiments with the specified status
        """
        stmt = select(self.model).where(
            self.model.status == status
        ).offset(skip).limit(limit)
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_by_status(self, db: Session, status: str) -> int:
        """
        Count experiments with a specific status
        
        Args:
            db: Database session
            status: Status to count
            
        Returns:
            Count of experiments with the specified status
        """
        stmt = select(func.count()).select_from(self.model).where(
            self.model.status == status
        )
        result = db.execute(stmt).scalar_one()
        return result
    
    def update_status(
        self, db: Session, experiment_id: UUID, status: str
    ) -> Optional[Experiment]:
        """
        Update the status of an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment to update
            status: New status value
            
        Returns:
            Updated experiment if found, None otherwise
        """
        db_obj = self.get(db, experiment_id)
        if not db_obj:
            return None
        
        db_obj.status = status
        db_obj.updated_at = datetime.utcnow()
        
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def get_by_user(
        self, db: Session, user_id: int, skip: int = 0, limit: int = 100
    ) -> List[Experiment]:
        """
        Get experiments created by a specific user
        
        Args:
            db: Database session
            user_id: ID of the creator user
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of experiments created by the user
        """
        stmt = select(self.model).where(
            self.model.created_by == user_id
        ).offset(skip).limit(limit)
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_by_user(self, db: Session, user_id: int) -> int:
        """
        Count experiments created by a specific user
        
        Args:
            db: Database session
            user_id: ID of the creator user
            
        Returns:
            Count of experiments created by the user
        """
        stmt = select(func.count()).select_from(self.model).where(
            self.model.created_by == user_id
        )
        result = db.execute(stmt).scalar_one()
        return result
    
    def add_molecules(
        self, db: Session, experiment_id: UUID, molecule_ids: List[UUID]
    ) -> Dict[str, Any]:
        """
        Add molecules to an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            molecule_ids: List of molecule IDs to add
            
        Returns:
            Result with success count and failures
        """
        # Check if experiment exists
        experiment = self.get(db, experiment_id)
        if not experiment:
            return {
                "success": False,
                "error": "Experiment not found",
                "experiment_id": experiment_id,
                "added": 0,
                "failed": len(molecule_ids)
            }
        
        success_count = 0
        failed_count = 0
        failed_ids = []
        
        for molecule_id in molecule_ids:
            # Check if molecule exists
            molecule_stmt = select(Molecule).where(Molecule.id == molecule_id)
            molecule = db.execute(molecule_stmt).scalar_one_or_none()
            
            if not molecule:
                failed_count += 1
                failed_ids.append(str(molecule_id))
                continue
            
            # Check if association already exists
            assoc_stmt = select(ExperimentMolecule).where(
                and_(
                    ExperimentMolecule.experiment_id == experiment_id,
                    ExperimentMolecule.molecule_id == molecule_id
                )
            )
            existing_assoc = db.execute(assoc_stmt).scalar_one_or_none()
            
            if existing_assoc:
                # Already associated, count as success but do nothing
                success_count += 1
                continue
            
            # Create new association
            new_assoc = ExperimentMolecule(
                experiment_id=experiment_id,
                molecule_id=molecule_id
            )
            db.add(new_assoc)
            success_count += 1
        
        db.commit()
        
        return {
            "success": True,
            "experiment_id": str(experiment_id),
            "added": success_count,
            "failed": failed_count,
            "failed_ids": failed_ids
        }
    
    def remove_molecules(
        self, db: Session, experiment_id: UUID, molecule_ids: List[UUID]
    ) -> Dict[str, Any]:
        """
        Remove molecules from an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            molecule_ids: List of molecule IDs to remove
            
        Returns:
            Result with success count and failures
        """
        # Check if experiment exists
        experiment = self.get(db, experiment_id)
        if not experiment:
            return {
                "success": False,
                "error": "Experiment not found",
                "experiment_id": experiment_id,
                "removed": 0,
                "failed": len(molecule_ids)
            }
        
        success_count = 0
        failed_count = 0
        failed_ids = []
        
        for molecule_id in molecule_ids:
            # Delete association if exists
            delete_stmt = delete(ExperimentMolecule).where(
                and_(
                    ExperimentMolecule.experiment_id == experiment_id,
                    ExperimentMolecule.molecule_id == molecule_id
                )
            )
            result = db.execute(delete_stmt)
            
            if result.rowcount > 0:
                success_count += 1
            else:
                failed_count += 1
                failed_ids.append(str(molecule_id))
        
        db.commit()
        
        return {
            "success": True,
            "experiment_id": str(experiment_id),
            "removed": success_count,
            "failed": failed_count,
            "failed_ids": failed_ids
        }
    
    def get_molecules(
        self, db: Session, experiment_id: UUID, skip: int = 0, limit: int = 100
    ) -> List[Molecule]:
        """
        Get molecules associated with an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of molecules in the experiment
        """
        stmt = select(Molecule).join(
            ExperimentMolecule, 
            and_(
                ExperimentMolecule.molecule_id == Molecule.id,
                ExperimentMolecule.experiment_id == experiment_id
            )
        ).offset(skip).limit(limit)
        
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def count_molecules(self, db: Session, experiment_id: UUID) -> int:
        """
        Count molecules associated with an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            
        Returns:
            Count of molecules in the experiment
        """
        stmt = select(func.count()).select_from(Molecule).join(
            ExperimentMolecule, 
            and_(
                ExperimentMolecule.molecule_id == Molecule.id,
                ExperimentMolecule.experiment_id == experiment_id
            )
        )
        
        result = db.execute(stmt).scalar_one()
        return result
    
    def is_molecule_in_experiment(
        self, db: Session, experiment_id: UUID, molecule_id: UUID
    ) -> bool:
        """
        Check if a molecule is associated with an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            molecule_id: ID of the molecule
            
        Returns:
            True if the molecule is in the experiment, False otherwise
        """
        stmt = select(ExperimentMolecule).where(
            and_(
                ExperimentMolecule.experiment_id == experiment_id,
                ExperimentMolecule.molecule_id == molecule_id
            )
        )
        
        result = db.execute(stmt).scalar_one_or_none()
        return result is not None
    
    def get_by_type(
        self, db: Session, type_id: UUID, skip: int = 0, limit: int = 100
    ) -> List[Experiment]:
        """
        Get experiments of a specific type
        
        Args:
            db: Database session
            type_id: ID of the experiment type
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of experiments with the specified type
        """
        stmt = select(self.model).where(
            self.model.type_id == type_id
        ).offset(skip).limit(limit)
        
        results = db.execute(stmt).scalars().all()
        return list(results)
    
    def get_parameter_value(
        self, db: Session, experiment_id: UUID, parameter_name: str
    ) -> Optional[str]:
        """
        Get the value of a specific parameter for an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            parameter_name: Name of the parameter to retrieve
            
        Returns:
            Parameter value if found, None otherwise
        """
        stmt = select(ExperimentParameter.parameter_value).where(
            and_(
                ExperimentParameter.experiment_id == experiment_id,
                ExperimentParameter.parameter_name == parameter_name
            )
        )
        
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def set_parameter_value(
        self, db: Session, experiment_id: UUID, parameter_name: str, parameter_value: str
    ) -> ExperimentParameter:
        """
        Set the value of a specific parameter for an experiment
        
        Args:
            db: Database session
            experiment_id: ID of the experiment
            parameter_name: Name of the parameter
            parameter_value: Value to set
            
        Returns:
            Updated or created parameter
        """
        # Check if parameter exists
        stmt = select(ExperimentParameter).where(
            and_(
                ExperimentParameter.experiment_id == experiment_id,
                ExperimentParameter.parameter_name == parameter_name
            )
        )
        
        existing_param = db.execute(stmt).scalar_one_or_none()
        
        if existing_param:
            # Update existing parameter
            existing_param.parameter_value = parameter_value
            db.add(existing_param)
            db.commit()
            db.refresh(existing_param)
            return existing_param
        else:
            # Create new parameter
            new_param = ExperimentParameter(
                experiment_id=experiment_id,
                parameter_name=parameter_name,
                parameter_value=parameter_value
            )
            db.add(new_param)
            db.commit()
            db.refresh(new_param)
            return new_param
    
    def filter_experiments(
        self, db: Session, filters: Dict[str, Any], skip: int = 0, limit: int = 100
    ) -> Tuple[List[Experiment], int]:
        """
        Filter experiments based on various criteria
        
        Args:
            db: Database session
            filters: Dictionary of filter criteria
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            Tuple of (list of filtered experiments, total count)
        """
        # Start with base query
        query = select(self.model)
        
        # Apply filters
        if filters.get("name"):
            query = query.where(self.model.name.ilike(f"%{filters['name']}%"))
        
        if filters.get("status"):
            query = query.where(self.model.status == filters["status"])
        
        if filters.get("type_id"):
            query = query.where(self.model.type_id == filters["type_id"])
        
        if filters.get("created_by"):
            query = query.where(self.model.created_by == filters["created_by"])
        
        if filters.get("molecule_id"):
            # Filter by associated molecule
            query = query.join(
                ExperimentMolecule,
                ExperimentMolecule.experiment_id == self.model.id
            ).where(ExperimentMolecule.molecule_id == filters["molecule_id"])
        
        # Get total count before pagination
        count_query = select(func.count()).select_from(query.subquery())
        total = db.execute(count_query).scalar_one()
        
        # Apply sorting
        if filters.get("sort_by"):
            sort_column = getattr(self.model, filters["sort_by"], self.model.created_at)
            if filters.get("sort_desc", False):
                query = query.order_by(desc(sort_column))
            else:
                query = query.order_by(asc(sort_column))
        else:
            # Default sort by created_at descending
            query = query.order_by(desc(self.model.created_at))
        
        # Apply pagination
        query = query.offset(skip).limit(limit)
        
        # Execute query
        results = db.execute(query).scalars().all()
        
        return list(results), total


# Create a singleton instance
experiment = CRUDExperiment()