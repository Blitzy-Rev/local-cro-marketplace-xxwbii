from sqlalchemy import Column, ForeignKey, DateTime, Index
from sqlalchemy.orm import relationship
from datetime import datetime  # standard library

from ..db.base_class import Base
from .library import Library
from .molecule import Molecule


class LibraryMolecule(Base):
    """
    SQLAlchemy model representing the many-to-many relationship between libraries and molecules.
    
    This model serves as a junction table that connects libraries and molecules, allowing 
    molecules to be organized into user-defined collections. It includes metadata such as 
    when a molecule was added to a library.
    
    Attributes:
        library_id (int): Foreign key reference to the library
        molecule_id (int): Foreign key reference to the molecule
        added_at (datetime): When the molecule was added to the library
        library (Library): Relationship to the library
        molecule (Molecule): Relationship to the molecule
    """
    
    __tablename__ = "library_molecule"
    
    # Override the 'id' column from Base
    id = None
    
    # Foreign keys that form the composite primary key
    library_id = Column(ForeignKey('libraries.id'), primary_key=True, nullable=False)
    molecule_id = Column(ForeignKey('molecules.id'), primary_key=True, nullable=False)
    
    # Metadata
    added_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    
    # Relationships
    library = relationship("Library", back_populates="library_molecules")
    molecule = relationship("Molecule", back_populates="library_associations")
    
    # Define indexes for performance optimization
    __table_args__ = (
        Index('ix_library_molecule_library_id', 'library_id'),
        Index('ix_library_molecule_molecule_id', 'molecule_id'),
    )