from sqlalchemy import Column, ForeignKey, String, DateTime, Enum, Text  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship, backref  # sqlalchemy.orm version 2.0+
from datetime import datetime  # standard library
from typing import List, Optional, Union
from uuid import UUID  # standard library

from ..db.base_class import Base
from ..constants import ResultStatus
from .submission import Submission


class Result(Base):
    """
    SQLAlchemy model representing experimental results uploaded by CRO users.
    
    This model stores information about experimental results including their status,
    submission relationship, upload/approval timestamps, and maintains relationships
    with result files and structured data points.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        submission_id (int): Foreign key referencing the submission this result belongs to
        status (ResultStatus): Current status of the result (PENDING, UPLOADED, APPROVED, REJECTED)
        uploaded_at (datetime): When the result was uploaded
        approved_at (datetime): When the result was approved (if applicable)
        notes (str): Additional notes or comments about the result
        
        submission (Submission): Relationship to the associated submission
        files (list): One-to-many relationship with result files
        data_points (list): One-to-many relationship with structured result data points
    """
    
    # Foreign keys
    submission_id = Column(ForeignKey('submissions.id'), nullable=False, index=True)
    
    # Status and timestamps
    status = Column(Enum(ResultStatus), default=ResultStatus.PENDING, nullable=False, index=True)
    uploaded_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    approved_at = Column(DateTime, nullable=True)
    notes = Column(Text, nullable=True)
    
    # Relationships
    submission = relationship("Submission", back_populates="results")
    
    # One-to-many relationships with result files and data points
    # These will be defined in their respective models with back_populates
    files = relationship("ResultFile", back_populates="result", cascade="all, delete-orphan")
    data_points = relationship("ResultData", back_populates="result", cascade="all, delete-orphan")