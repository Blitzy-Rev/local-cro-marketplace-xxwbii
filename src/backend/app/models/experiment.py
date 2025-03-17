from sqlalchemy import Column, String, DateTime, ForeignKey, Enum, Integer, Text  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship, backref  # sqlalchemy.orm version 2.0+
from datetime import datetime  # standard library
from uuid import UUID  # standard library

from ..db.base_class import Base
from ..constants import ExperimentStatus
from .user import User
from .experiment_type import ExperimentType


class Experiment(Base):
    """
    SQLAlchemy model representing an experiment configuration for molecular testing.
    
    This model tracks experiments throughout their lifecycle from creation to results
    and manages relationships with users, molecules, parameters, and submissions.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        name (str): Name of the experiment
        type_id (int): Foreign key referencing the experiment type
        status (ExperimentStatus): Current status of the experiment
        created_by (int): Foreign key referencing the user who created the experiment
        created_at (datetime): When the experiment was created
        updated_at (datetime): When the experiment was last updated
        description (str): Detailed description of the experiment
        
        created_by_user (User): Relationship to the user who created the experiment
        experiment_type (ExperimentType): Relationship to the experiment type
        parameters (list): Relationship to experiment parameters
        molecules (list): Relationship to molecules in this experiment
        molecule_associations (list): Relationship to the association table
        submissions (list): Relationship to CRO submissions for this experiment
    """
    
    # Basic properties
    name = Column(String(255), nullable=False, index=True)
    type_id = Column(Integer, ForeignKey('experiment_types.id'), nullable=False)
    status = Column(Enum(ExperimentStatus), default=ExperimentStatus.DRAFT, nullable=False)
    created_by = Column(Integer, ForeignKey('users.id'), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    description = Column(Text, nullable=True)
    
    # Relationships
    created_by_user = relationship("User", back_populates="experiments", foreign_keys=[created_by])
    experiment_type = relationship("ExperimentType", back_populates="experiments")
    
    # One-to-many relationship with experiment parameters
    parameters = relationship("ExperimentParameter", back_populates="experiment", cascade="all, delete-orphan")
    
    # Many-to-many relationship with molecules
    molecules = relationship(
        "Molecule",
        secondary="experiment_molecules",
        back_populates="experiments"
    )
    
    # Relationship to the association table for experiment-molecule relationships
    molecule_associations = relationship("ExperimentMolecule", back_populates="experiment", cascade="all, delete-orphan")
    
    # One-to-many relationship with submissions
    submissions = relationship("Submission", back_populates="experiment", cascade="all, delete-orphan")