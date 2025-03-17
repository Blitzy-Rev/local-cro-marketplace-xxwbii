from sqlalchemy import Column, String, ForeignKey, DateTime, Text, Index
from sqlalchemy.orm import relationship, backref
from datetime import datetime  # standard library

from ..db.base_class import Base
from .user import User


class Library(Base):
    """
    SQLAlchemy model representing a user-defined collection of molecules.
    
    This model enables users to organize molecules into named libraries for
    easier management and research purposes. Libraries have a many-to-many
    relationship with molecules, allowing a molecule to belong to multiple
    libraries and libraries to contain multiple molecules.
    
    Attributes:
        name (str): Name of the library (required)
        description (str): Detailed description of the library's purpose
        created_by (int): Foreign key reference to the user who created the library
        created_at (datetime): When the library was created
        updated_at (datetime): When the library was last updated
        created_by_user (relationship): Relationship to the user who created the library
        creator (property): Alias for created_by_user to fulfill export requirements
        molecules (relationship): Many-to-many relationship to molecules in this library
        library_molecules (relationship): Relationship to the library_molecule association table
    """
    name = Column(String(100), nullable=False, index=True)
    description = Column(Text, nullable=True)
    created_by = Column(ForeignKey('users.id'), nullable=False, index=True)
    created_at = Column(DateTime, default=datetime.utcnow, nullable=False)
    updated_at = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow, nullable=False)
    
    # Relationships
    created_by_user = relationship("User", back_populates="libraries")
    molecules = relationship(
        "Molecule",
        secondary="library_molecule",
        back_populates="libraries"
    )
    library_molecules = relationship("LibraryMolecule", back_populates="library")
    
    @property
    def creator(self):
        """
        Alias for created_by_user to fulfill export specification.
        
        Returns:
            User: The user who created this library
        """
        return self.created_by_user