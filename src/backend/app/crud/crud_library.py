from typing import List, Dict, Any, Optional, Tuple, Union
from uuid import UUID
from datetime import datetime
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select, update, delete, join, func, and_, or_

from .base import CRUDBase
from ..models.library import Library
from ..models.library_molecule import LibraryMolecule
from ..models.molecule import Molecule
from ..schemas.library import LibraryCreate, LibraryUpdate


class CRUDLibrary(CRUDBase[Library, LibraryCreate, LibraryUpdate]):
    """
    CRUD operations for Library model with specialized methods for library management.
    
    This class extends the base CRUD operations to provide specialized functionality
    for managing libraries, including library-molecule associations, filtering by
    various criteria, and organizing molecules into custom collections.
    """
    
    def __init__(self):
        """Initialize the CRUD object with the Library model"""
        super().__init__(Library)
    
    def get_by_name(self, db: Session, name: str, user_id: UUID) -> Optional[Library]:
        """
        Get a library by name and creator ID.
        
        Args:
            db: Database session
            name: Name of the library to retrieve
            user_id: ID of the user who created the library
            
        Returns:
            Library instance if found, None otherwise
        """
        stmt = select(Library).where(
            and_(
                Library.name == name,
                Library.created_by == user_id
            )
        )
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def get_with_molecules(
        self, db: Session, library_id: UUID, skip: int = 0, limit: int = 100
    ) -> Optional[Tuple[Library, List[Any]]]:
        """
        Get a library with its associated molecules.
        
        Args:
            db: Database session
            library_id: ID of the library to retrieve
            skip: Number of molecules to skip (for pagination)
            limit: Maximum number of molecules to return
            
        Returns:
            Tuple containing the library and its molecules if found, None otherwise
        """
        # Get the library
        library = self.get(db, library_id)
        if not library:
            return None
        
        # Get the molecules in the library
        molecules = self.get_molecules(db, library_id, skip, limit)
        
        return (library, molecules)
    
    def get_by_user(
        self, db: Session, user_id: UUID, skip: int = 0, limit: int = 100
    ) -> List[Library]:
        """
        Get libraries created by a specific user.
        
        Args:
            db: Database session
            user_id: ID of the user who created the libraries
            skip: Number of libraries to skip (for pagination)
            limit: Maximum number of libraries to return
            
        Returns:
            List of libraries created by the user
        """
        stmt = select(Library).where(
            Library.created_by == user_id
        ).offset(skip).limit(limit)
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def count_by_user(self, db: Session, user_id: UUID) -> int:
        """
        Count libraries created by a specific user.
        
        Args:
            db: Database session
            user_id: ID of the user who created the libraries
            
        Returns:
            Count of libraries created by the user
        """
        stmt = select(func.count()).select_from(Library).where(
            Library.created_by == user_id
        )
        result = db.execute(stmt).scalar_one()
        return result
    
    def add_molecules(
        self, db: Session, library_id: UUID, molecule_ids: List[UUID]
    ) -> Dict[str, Any]:
        """
        Add molecules to a library.
        
        Args:
            db: Database session
            library_id: ID of the library to add molecules to
            molecule_ids: List of molecule IDs to add to the library
            
        Returns:
            Dictionary with success count and failures
        """
        # Initialize counters and result
        success_count = 0
        failures = []
        
        # Check if library exists
        library = self.get(db, library_id)
        if not library:
            return {
                "success": success_count,
                "failures": [{"error": "Library not found", "library_id": str(library_id)}]
            }
        
        # Add molecules one by one to handle failures gracefully
        for molecule_id in molecule_ids:
            try:
                # Check if molecule exists
                stmt = select(Molecule).where(Molecule.id == molecule_id)
                molecule = db.execute(stmt).scalar_one_or_none()
                if not molecule:
                    failures.append({
                        "molecule_id": str(molecule_id),
                        "error": "Molecule not found"
                    })
                    continue
                
                # Check if molecule is already in the library
                if self.is_molecule_in_library(db, library_id, molecule_id):
                    failures.append({
                        "molecule_id": str(molecule_id),
                        "error": "Molecule already in library"
                    })
                    continue
                
                # Add molecule to library
                library_molecule = LibraryMolecule(
                    library_id=library_id,
                    molecule_id=molecule_id,
                    added_at=datetime.utcnow()
                )
                db.add(library_molecule)
                success_count += 1
            
            except Exception as e:
                failures.append({
                    "molecule_id": str(molecule_id),
                    "error": str(e)
                })
        
        # Commit transaction if any molecules were added successfully
        if success_count > 0:
            db.commit()
        
        return {
            "success": success_count,
            "failures": failures
        }
    
    def remove_molecules(
        self, db: Session, library_id: UUID, molecule_ids: List[UUID]
    ) -> Dict[str, Any]:
        """
        Remove molecules from a library.
        
        Args:
            db: Database session
            library_id: ID of the library to remove molecules from
            molecule_ids: List of molecule IDs to remove from the library
            
        Returns:
            Dictionary with success count and failures
        """
        failures = []
        
        # Check if library exists
        library = self.get(db, library_id)
        if not library:
            return {
                "success": 0,
                "failures": [{"error": "Library not found", "library_id": str(library_id)}]
            }
        
        try:
            # Create delete statement for multiple molecules at once
            stmt = delete(LibraryMolecule).where(
                and_(
                    LibraryMolecule.library_id == library_id,
                    LibraryMolecule.molecule_id.in_(molecule_ids)
                )
            )
            result = db.execute(stmt)
            success_count = result.rowcount
            
            # Commit the transaction
            db.commit()
            
        except Exception as e:
            db.rollback()
            return {
                "success": 0,
                "failures": [{"error": str(e)}]
            }
        
        # Check if any molecules were not found in the library
        if success_count < len(molecule_ids):
            failures.append({
                "error": "Some molecules were not found in the library",
                "count": len(molecule_ids) - success_count
            })
        
        return {
            "success": success_count,
            "failures": failures
        }
    
    def get_molecules(
        self, db: Session, library_id: UUID, skip: int = 0, limit: int = 100
    ) -> List[Any]:
        """
        Get molecules in a library with pagination.
        
        Args:
            db: Database session
            library_id: ID of the library to get molecules from
            skip: Number of molecules to skip (for pagination)
            limit: Maximum number of molecules to return
            
        Returns:
            List of molecules in the library
        """
        # Join Molecule and LibraryMolecule tables to get molecules in the library
        stmt = select(Molecule).join(
            LibraryMolecule, 
            Molecule.id == LibraryMolecule.molecule_id
        ).where(
            LibraryMolecule.library_id == library_id
        ).offset(skip).limit(limit)
        
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def count_molecules(self, db: Session, library_id: UUID) -> int:
        """
        Count molecules in a library.
        
        Args:
            db: Database session
            library_id: ID of the library to count molecules in
            
        Returns:
            Count of molecules in the library
        """
        stmt = select(func.count()).select_from(LibraryMolecule).where(
            LibraryMolecule.library_id == library_id
        )
        result = db.execute(stmt).scalar_one()
        return result
    
    def is_molecule_in_library(self, db: Session, library_id: UUID, molecule_id: UUID) -> bool:
        """
        Check if a molecule is in a library.
        
        Args:
            db: Database session
            library_id: ID of the library to check
            molecule_id: ID of the molecule to check
            
        Returns:
            True if the molecule is in the library, False otherwise
        """
        stmt = select(LibraryMolecule).where(
            and_(
                LibraryMolecule.library_id == library_id,
                LibraryMolecule.molecule_id == molecule_id
            )
        )
        result = db.execute(stmt).first()
        return result is not None
    
    def filter_libraries(
        self, 
        db: Session, 
        name_contains: Optional[str] = None,
        created_by: Optional[UUID] = None,
        skip: int = 0, 
        limit: int = 100,
        sort_by: Optional[str] = None,
        sort_desc: bool = False
    ) -> Tuple[List[Library], int]:
        """
        Filter libraries based on criteria.
        
        Args:
            db: Database session
            name_contains: Filter libraries by name (case-insensitive partial match)
            created_by: Filter libraries by creator user ID
            skip: Number of libraries to skip (for pagination)
            limit: Maximum number of libraries to return
            sort_by: Field to sort by (name, created_at)
            sort_desc: Sort in descending order if True
            
        Returns:
            Tuple containing filtered libraries and total count
        """
        # Base query
        query = select(Library)
        count_query = select(func.count()).select_from(Library)
        
        # Apply filters
        filters = []
        
        if name_contains:
            filters.append(Library.name.ilike(f"%{name_contains}%"))
        
        if created_by:
            filters.append(Library.created_by == created_by)
        
        if filters:
            query = query.where(and_(*filters))
            count_query = count_query.where(and_(*filters))
        
        # Apply sorting
        if sort_by:
            if sort_by == "name":
                order_by = Library.name.desc() if sort_desc else Library.name
            elif sort_by == "created_at":
                order_by = Library.created_at.desc() if sort_desc else Library.created_at
            else:
                # Default to sort by name
                order_by = Library.name.desc() if sort_desc else Library.name
            
            query = query.order_by(order_by)
        
        # Apply pagination
        query = query.offset(skip).limit(limit)
        
        # Execute queries
        libraries = list(db.execute(query).scalars().all())
        total = db.execute(count_query).scalar_one()
        
        return libraries, total


# Singleton instance for application-wide use
library = CRUDLibrary()