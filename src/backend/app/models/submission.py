from sqlalchemy import Column, String, DateTime, ForeignKey, Enum, Integer, Float, Text  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship, backref  # sqlalchemy.orm version 2.0+
from datetime import datetime  # standard library
from uuid import UUID  # standard library

from ..db.base_class import Base
from ..constants import SubmissionStatus
from .user import User
from .experiment import Experiment


class Submission(Base):
    """
    SQLAlchemy model representing a submission of an experiment to a CRO for testing.
    
    This model tracks submissions of experiments to Contract Research Organizations (CROs),
    including status tracking, pricing information, and relationships to experiments,
    CRO users, submission details, and results.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        experiment_id (int): Foreign key referencing the experiment being submitted
        cro_id (int): Foreign key referencing the CRO user assigned to the submission
        status (SubmissionStatus): Current status of the submission
        submitted_at (datetime): When the submission was created
        updated_at (datetime): When the submission was last updated
        price (float): Price quoted for the experiment
        turnaround_days (int): Estimated turnaround time in days
        notes (str): Additional notes about the submission
        
        experiment (Experiment): Relationship to the experiment being submitted
        cro (User): Relationship to the CRO user assigned to the submission
        details (list): Relationship to submission details
        results (list): Relationship to experimental results
    """
    
    # Foreign keys
    experiment_id = Column(Integer, ForeignKey('experiments.id'), nullable=False)
    cro_id = Column(Integer, ForeignKey('users.id'), nullable=False)
    
    # Status and timestamps
    status = Column(Enum(SubmissionStatus), default=SubmissionStatus.PENDING, nullable=False)
    submitted_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    
    # Pricing and turnaround information
    price = Column(Float, nullable=True)
    turnaround_days = Column(Integer, nullable=True)
    notes = Column(Text, nullable=True)
    
    # Relationships
    experiment = relationship("Experiment", back_populates="submissions")
    cro = relationship("User", foreign_keys=[cro_id], back_populates="cro_submissions")
    
    # One-to-many relationships with submission details and results
    details = relationship("SubmissionDetail", back_populates="submission", cascade="all, delete-orphan")
    results = relationship("Result", back_populates="submission", cascade="all, delete-orphan")