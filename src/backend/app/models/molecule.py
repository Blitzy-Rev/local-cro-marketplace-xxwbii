from sqlalchemy import Column, String, DateTime, ForeignKey, Text, Index, Integer
from sqlalchemy.orm import relationship, backref
from datetime import datetime
from uuid import UUID  # For type hints in docstrings

from ..db.base_class import Base


class Molecule(Base):
    """
    SQLAlchemy model representing a molecular structure with its SMILES representation and associated metadata.
    
    This model serves as a central entity in the system's data model, enabling the organization,
    filtering, and experimental testing of chemical compounds.
    
    Attributes:
        id (int): Primary key (inherited from Base)
        smiles (str): SMILES representation of the molecular structure
        created_by (int): Foreign key to the user who created this molecule
        created_at (datetime): Timestamp when the molecule was created
        updated_at (datetime): Timestamp when the molecule was last updated
        flag_status (str): Optional flag for priority review (e.g., "HIGH", "MEDIUM", "LOW")
        
        created_by_user (User): Relationship to the user who created this molecule
        properties (List[MoleculeProperty]): Relationship to molecular properties
        libraries (List[Library]): Relationship to libraries containing this molecule
        library_associations (List[LibraryMolecule]): Relationship to library junction table
        experiments (List[Experiment]): Relationship to experiments using this molecule
        experiment_associations (List[ExperimentMolecule]): Relationship to experiment junction table
        result_data (List[ResultData]): Relationship to experimental results for this molecule
    """
    # Column definitions
    smiles = Column(String(4000), unique=True, nullable=False, index=True)
    created_by = Column(Integer, ForeignKey('users.id'), nullable=False)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=True)
    flag_status = Column(String(50), nullable=True, index=True)
    
    # Relationships
    created_by_user = relationship("User", back_populates="molecules")
    properties = relationship("MoleculeProperty", back_populates="molecule", cascade="all, delete-orphan")
    
    # Library relationships
    libraries = relationship("Library", secondary="library_molecule", back_populates="molecules")
    library_associations = relationship("LibraryMolecule", back_populates="molecule", cascade="all, delete-orphan")
    
    # Experiment relationships
    experiments = relationship("Experiment", secondary="experiment_molecule", back_populates="molecules")
    experiment_associations = relationship("ExperimentMolecule", back_populates="molecule", cascade="all, delete-orphan")
    
    # Result relationships
    result_data = relationship("ResultData", back_populates="molecule", cascade="all, delete-orphan")
    
    # Indexes for performance optimization
    __table_args__ = (
        Index('ix_molecule_smiles_flag', 'smiles', 'flag_status'),
    )