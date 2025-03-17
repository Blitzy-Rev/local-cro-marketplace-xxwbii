from sqlalchemy import Column, ForeignKey, String, Integer, DateTime  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+
from datetime import datetime  # standard library
from uuid import UUID  # standard library

from ..db.base_class import Base
from .result import Result


class ResultFile(Base):
    """
    SQLAlchemy model representing a file associated with an experimental result.
    
    This model stores metadata about files uploaded as part of experimental results,
    including file name, path, size, type, and upload timestamp. Files are stored in
    the MinIO object storage system, with this model tracking the metadata and
    maintaining the relationship to the parent Result entity.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        result_id (int): Foreign key referencing the result this file belongs to
        file_name (str): Original name of the uploaded file
        file_path (str): Path to the file in the object storage system
        file_size (int): Size of the file in bytes
        file_type (str): MIME type or file extension
        uploaded_at (datetime): When the file was uploaded
        
        result (Result): Relationship to the associated result
    """
    
    # Foreign keys
    result_id = Column(ForeignKey('results.id'), nullable=False, index=True)
    
    # File metadata
    file_name = Column(String(255), nullable=False)
    file_path = Column(String(1024), nullable=False, unique=True)
    file_size = Column(Integer, nullable=False)
    file_type = Column(String(100), nullable=False)
    uploaded_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    
    # Relationships
    result = relationship("Result", back_populates="files")