from sqlalchemy import Column, String  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+
from uuid import UUID  # standard library

from ..db.base_class import Base


class ExperimentType(Base):
    """
    SQLAlchemy model representing a type of experiment that can be conducted by CROs.
    
    This model serves as a catalog of available experimental assessments, including binding assays,
    ADME panels, toxicity assessments, and other experimental procedures supported by the platform.
    
    Attributes:
        id (int): Primary key for the experiment type (inherited from Base)
        name (str): Name of the experiment type
        description (str): Detailed description of what the experiment measures or assesses
        category (str): Category of the experiment (e.g., 'ADME', 'Binding', 'Toxicity')
        experiments (relationship): Relationship to experiments of this type
    """
    
    name = Column(String(100), nullable=False, index=True, unique=True)
    description = Column(String(500), nullable=True)
    category = Column(String(100), nullable=False, index=True)
    
    # Define relationship to experiments (one-to-many)
    experiments = relationship("Experiment", back_populates="experiment_type")