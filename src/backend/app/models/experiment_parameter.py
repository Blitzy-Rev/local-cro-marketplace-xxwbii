from sqlalchemy import Column, String, ForeignKey, Integer, Index  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+

from ..db.base_class import Base
from .experiment import Experiment  # Imported for type annotations

class ExperimentParameter(Base):
    """
    SQLAlchemy model representing a parameter for an experiment, stored as a key-value pair.
    
    This model allows for flexible configuration of experiment settings without requiring
    schema changes for new parameter types. Parameters are stored as name-value pairs
    associated with specific experiments.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        experiment_id (int): Foreign key to the experiment table
        parameter_name (str): Name of the parameter
        parameter_value (str): Value of the parameter
        experiment (Experiment): Relationship to the parent experiment
    """
    
    experiment_id = Column(Integer, ForeignKey("experiments.id"), nullable=False, index=True)
    parameter_name = Column(String(255), nullable=False)
    parameter_value = Column(String(500), nullable=True)
    
    # Create a composite index on experiment_id and parameter_name for efficient lookups
    __table_args__ = (
        Index("ix_experiment_parameter_exp_name", "experiment_id", "parameter_name", unique=True),
    )
    
    # Relationship back to the experiment
    experiment = relationship("Experiment", back_populates="parameters")