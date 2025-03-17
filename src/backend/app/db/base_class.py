from typing import Any
from sqlalchemy import Column, Integer
from sqlalchemy.orm import as_declarative, declared_attr, declarative_base  # sqlalchemy.orm version 2.0+


@as_declarative()
class Base:
    """
    SQLAlchemy declarative base class for all database models in the application.
    
    This base class provides common functionality and properties that will be inherited
    by all model classes, ensuring consistent structure and behavior across the database schema.
    
    Features:
    - Automatic table name generation from class name
    - Common id primary key field for all tables
    """
    
    # All tables will have an 'id' primary key column
    id = Column(Integer, primary_key=True, index=True)
    
    @declared_attr
    @classmethod
    def __tablename__(cls) -> str:
        """
        Automatically generates the table name based on the class name.
        
        This method is called during the SQLAlchemy model class creation process
        to determine the table name in the database.
        
        Returns:
            str: The lowercase version of the class name to be used as table name
        """
        return cls.__name__.lower()