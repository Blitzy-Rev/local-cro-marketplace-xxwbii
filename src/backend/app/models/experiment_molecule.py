from sqlalchemy import Column, ForeignKey, DateTime, Index  # sqlalchemy version 2.0+
from sqlalchemy.orm import relationship  # sqlalchemy.orm version 2.0+
from datetime import datetime  # standard library

from ..db.base_class import Base


class ExperimentMolecule(Base):
    """
    SQLAlchemy model representing the many-to-many relationship between experiments and molecules.
    
    This junction table associates molecules with experiments, allowing molecules to be grouped
    into experiments for testing by CROs. It includes metadata such as when the molecule was
    added to the experiment.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        experiment_id (int): Foreign key to the experiment
        molecule_id (int): Foreign key to the molecule
        added_at (datetime): Timestamp when the molecule was added to the experiment
        experiment (Experiment): Relationship to the experiment
        molecule (Molecule): Relationship to the molecule
    """
    # Override table name to match what's expected in relationship definitions
    __tablename__ = "experiment_molecules"
    
    experiment_id = Column(ForeignKey('experiments.id'), nullable=False)
    molecule_id = Column(ForeignKey('molecules.id'), nullable=False)
    added_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    
    # Define relationships
    experiment = relationship("Experiment", back_populates="molecule_associations")
    molecule = relationship("Molecule", back_populates="experiment_associations")
    
    # Add unique constraint to ensure a molecule is only added once to an experiment
    __table_args__ = (
        Index('ix_experiment_molecules_exp_mol', 'experiment_id', 'molecule_id', unique=True),
    )