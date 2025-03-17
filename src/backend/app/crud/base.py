from typing import TypeVar, Generic, Type, Any, Dict, List, Optional, Union
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select, update, delete, func

from ..db.base_class import Base

# Define type variables for generic typing
ModelType = TypeVar('ModelType', bound=Base)
CreateSchemaType = TypeVar('CreateSchemaType')
UpdateSchemaType = TypeVar('UpdateSchemaType')


class CRUDBase(Generic[ModelType, CreateSchemaType, UpdateSchemaType]):
    """
    Base class for CRUD operations on SQLAlchemy models.
    
    This class provides standard create, read, update, and delete operations
    for SQLAlchemy models. It is designed to be extended by model-specific
    CRUD classes to provide consistent data access patterns across the application.
    
    Type Parameters:
        ModelType: The SQLAlchemy model class
        CreateSchemaType: The schema type used for model creation
        UpdateSchemaType: The schema type used for model updates
    """
    
    def __init__(self, model: Type[ModelType]):
        """
        Initialize a CRUD instance with a specific model class.
        
        Args:
            model: The SQLAlchemy model class this CRUD instance will operate on
        """
        self.model = model
    
    def get(self, db: Session, id: Any) -> Optional[ModelType]:
        """
        Get a single record by ID.
        
        Args:
            db: Database session
            id: Primary key value of the record to retrieve
            
        Returns:
            The model instance if found, None otherwise
        """
        stmt = select(self.model).where(self.model.id == id)
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def get_multi(self, db: Session, skip: int = 0, limit: int = 100) -> List[ModelType]:
        """
        Get multiple records with pagination.
        
        Args:
            db: Database session
            skip: Number of records to skip (for pagination)
            limit: Maximum number of records to return
            
        Returns:
            List of model instances
        """
        stmt = select(self.model).offset(skip).limit(limit)
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def create(self, db: Session, obj_in: Union[CreateSchemaType, Dict[str, Any]]) -> ModelType:
        """
        Create a new record.
        
        Args:
            db: Database session
            obj_in: Input data (schema instance or dictionary)
            
        Returns:
            The created model instance
        """
        if isinstance(obj_in, dict):
            obj_data = obj_in
        else:
            obj_data = obj_in.dict()
        
        db_obj = self.model(**obj_data)
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def update(
        self, db: Session, db_obj: ModelType, obj_in: Union[UpdateSchemaType, Dict[str, Any]]
    ) -> ModelType:
        """
        Update an existing record.
        
        Args:
            db: Database session
            db_obj: Database model instance to update
            obj_in: Update data (schema instance or dictionary)
            
        Returns:
            The updated model instance
        """
        obj_data = {
            key: val for key, val in db_obj.__dict__.items() 
            if not key.startswith("_")
        }
        
        if isinstance(obj_in, dict):
            update_data = obj_in
        else:
            update_data = obj_in.dict(exclude_unset=True)
        
        for field in update_data:
            if field in obj_data:
                setattr(db_obj, field, update_data[field])
        
        db.add(db_obj)
        db.commit()
        db.refresh(db_obj)
        return db_obj
    
    def remove(self, db: Session, id: Any) -> Optional[ModelType]:
        """
        Delete a record.
        
        Args:
            db: Database session
            id: Primary key value of the record to delete
            
        Returns:
            The deleted model instance if found, None otherwise
        """
        obj = self.get(db, id)
        if obj is None:
            return None
        
        db.delete(obj)
        db.commit()
        return obj
    
    def count(self, db: Session) -> int:
        """
        Count total number of records.
        
        Args:
            db: Database session
            
        Returns:
            Total count of records
        """
        stmt = select(func.count()).select_from(self.model)
        result = db.execute(stmt).scalar_one()
        return result