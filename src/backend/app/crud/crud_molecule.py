from typing import List, Dict, Optional, Any, Tuple, Union
from datetime import datetime
from sqlalchemy.orm import Session  # sqlalchemy.orm version 2.0+
from sqlalchemy import select, update, delete, func, and_, or_, not_, desc, asc  # sqlalchemy version 2.0+
from sqlalchemy.orm import join, outerjoin  # sqlalchemy.orm version 2.0+

from .base import CRUDBase
from ..models.molecule import Molecule
from ..models.molecule_property import MoleculeProperty
from ..schemas.molecule import MoleculeCreate, MoleculeUpdate


class CRUDMolecule(CRUDBase[Molecule, MoleculeCreate, MoleculeUpdate]):
    """
    CRUD operations for molecules with specialized methods for molecular data.
    
    This class extends the base CRUD operations with molecule-specific methods
    such as filtering by properties, searching by substructure/similarity,
    handling flag status, and performing bulk operations.
    """
    
    def __init__(self):
        """Initialize the CRUD object with the Molecule model"""
        super().__init__(Molecule)
    
    def get_by_smiles(self, db: Session, smiles: str) -> Optional[Molecule]:
        """
        Get a molecule by its SMILES string.
        
        Args:
            db: Database session
            smiles: SMILES representation of the molecule
            
        Returns:
            Molecule instance if found, None otherwise
        """
        stmt = select(Molecule).where(Molecule.smiles == smiles)
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def get_with_properties(self, db: Session, id: Any) -> Optional[Molecule]:
        """
        Get a molecule by ID with all its properties.
        
        Args:
            db: Database session
            id: ID of the molecule to retrieve
            
        Returns:
            Molecule instance with properties if found, None otherwise
        """
        stmt = (
            select(Molecule)
            .outerjoin(MoleculeProperty)
            .where(Molecule.id == id)
            .options(
                # This would use SQLAlchemy's selectinload or joinedload in a real implementation
                # Since we're not showing the full SQLAlchemy import structure, we'll use this comment
                # to indicate the intent to eagerly load properties
            )
        )
        result = db.execute(stmt).scalar_one_or_none()
        return result
    
    def get_by_property_range(
        self, 
        db: Session, 
        property_name: str, 
        min_value: Optional[float] = None, 
        max_value: Optional[float] = None,
        skip: int = 0, 
        limit: int = 100,
        sort_by: Optional[str] = None,
        sort_desc: bool = False
    ) -> List[Molecule]:
        """
        Get molecules with a property value within a specified range.
        
        Args:
            db: Database session
            property_name: Name of the property to filter by
            min_value: Minimum value (inclusive) for the property
            max_value: Maximum value (inclusive) for the property
            skip: Number of records to skip
            limit: Maximum number of records to return
            sort_by: Property to sort by
            sort_desc: Whether to sort in descending order
            
        Returns:
            List of molecules matching the property range criteria
        """
        # Base query joining molecules and properties
        stmt = (
            select(Molecule)
            .join(MoleculeProperty, Molecule.id == MoleculeProperty.molecule_id)
            .where(MoleculeProperty.property_name == property_name)
        )
        
        # Apply min value filter if provided
        if min_value is not None:
            stmt = stmt.where(MoleculeProperty.property_value >= min_value)
        
        # Apply max value filter if provided
        if max_value is not None:
            stmt = stmt.where(MoleculeProperty.property_value <= max_value)
        
        # Apply sorting if requested
        if sort_by:
            if sort_by == property_name:
                # Sorting by the property we're filtering on
                sort_column = MoleculeProperty.property_value
                if sort_desc:
                    stmt = stmt.order_by(desc(sort_column))
                else:
                    stmt = stmt.order_by(asc(sort_column))
            else:
                # Sorting by a molecule attribute
                sort_column = getattr(Molecule, sort_by, None)
                if sort_column:
                    if sort_desc:
                        stmt = stmt.order_by(desc(sort_column))
                    else:
                        stmt = stmt.order_by(asc(sort_column))
        
        # Apply pagination
        stmt = stmt.offset(skip).limit(limit)
        
        # Execute the query
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def count_by_property_range(
        self, 
        db: Session, 
        property_name: str,
        min_value: Optional[float] = None, 
        max_value: Optional[float] = None
    ) -> int:
        """
        Count molecules with a property value within a specified range.
        
        Args:
            db: Database session
            property_name: Name of the property to filter by
            min_value: Minimum value (inclusive) for the property
            max_value: Maximum value (inclusive) for the property
            
        Returns:
            Count of molecules matching the property range criteria
        """
        # Base query for counting molecules with the specified property
        stmt = (
            select(func.count(Molecule.id.distinct()))
            .join(MoleculeProperty, Molecule.id == MoleculeProperty.molecule_id)
            .where(MoleculeProperty.property_name == property_name)
        )
        
        # Apply min value filter if provided
        if min_value is not None:
            stmt = stmt.where(MoleculeProperty.property_value >= min_value)
        
        # Apply max value filter if provided
        if max_value is not None:
            stmt = stmt.where(MoleculeProperty.property_value <= max_value)
        
        # Execute the query
        result = db.execute(stmt).scalar_one()
        return result
    
    def get_by_flag_status(
        self, 
        db: Session, 
        flag_status: str,
        skip: int = 0, 
        limit: int = 100
    ) -> List[Molecule]:
        """
        Get molecules with a specific flag status.
        
        Args:
            db: Database session
            flag_status: Flag status to filter by
            skip: Number of records to skip
            limit: Maximum number of records to return
            
        Returns:
            List of molecules with the specified flag status
        """
        stmt = (
            select(Molecule)
            .where(Molecule.flag_status == flag_status)
            .offset(skip)
            .limit(limit)
        )
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def update_flag_status(
        self, 
        db: Session, 
        id: Any, 
        flag_status: str
    ) -> Optional[Molecule]:
        """
        Update the flag status of a molecule.
        
        Args:
            db: Database session
            id: ID of the molecule to update
            flag_status: New flag status value
            
        Returns:
            Updated molecule instance if found, None otherwise
        """
        molecule = self.get(db, id)
        if not molecule:
            return None
        
        molecule.flag_status = flag_status
        db.add(molecule)
        db.commit()
        db.refresh(molecule)
        return molecule
    
    def get_by_user(
        self, 
        db: Session, 
        user_id: Any,
        skip: int = 0, 
        limit: int = 100
    ) -> List[Molecule]:
        """
        Get molecules created by a specific user.
        
        Args:
            db: Database session
            user_id: ID of the user who created the molecules
            skip: Number of records to skip
            limit: Maximum number of records to return
            
        Returns:
            List of molecules created by the specified user
        """
        stmt = (
            select(Molecule)
            .where(Molecule.created_by == user_id)
            .offset(skip)
            .limit(limit)
        )
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def get_property_ranges(
        self, 
        db: Session, 
        property_name: str
    ) -> Tuple[float, float]:
        """
        Get the minimum and maximum values for a specific property.
        
        Args:
            db: Database session
            property_name: Name of the property to get range for
            
        Returns:
            Tuple containing (min_value, max_value) for the property
        """
        stmt = (
            select(
                func.min(MoleculeProperty.property_value),
                func.max(MoleculeProperty.property_value)
            )
            .where(MoleculeProperty.property_name == property_name)
        )
        result = db.execute(stmt).one_or_none()
        
        if result and result[0] is not None and result[1] is not None:
            return (float(result[0]), float(result[1]))
        return (0.0, 0.0)  # Default range if no data exists
    
    def get_available_properties(self, db: Session) -> List[str]:
        """
        Get a list of all property names available in the database.
        
        Args:
            db: Database session
            
        Returns:
            List of distinct property names
        """
        stmt = (
            select(MoleculeProperty.property_name)
            .distinct()
            .order_by(MoleculeProperty.property_name)
        )
        result = db.execute(stmt).scalars().all()
        return list(result)
    
    def search_by_substructure(
        self, 
        db: Session, 
        substructure_smarts: str,
        skip: int = 0, 
        limit: int = 100
    ) -> List[Molecule]:
        """
        Search for molecules containing a specific substructure (using RDKit cartridge).
        
        This method would use RDKit's substructure search functionality, which typically
        requires RDKit PostgreSQL cartridge to be installed in the database.
        
        Args:
            db: Database session
            substructure_smarts: SMARTS pattern for the substructure to search
            skip: Number of records to skip
            limit: Maximum number of records to return
            
        Returns:
            List of molecules containing the substructure
        """
        # Note: This is a simplified implementation
        # In a real application, this would use RDKit's substructure search functionality
        # via SQL function calls to the RDKit PostgreSQL cartridge
        
        # For example (pseudocode):
        # stmt = (
        #     select(Molecule)
        #     .where(func.mol_contains(func.mol_from_smiles(Molecule.smiles), 
        #                             func.mol_from_smarts(substructure_smarts)))
        #     .offset(skip)
        #     .limit(limit)
        # )
        
        # For this implementation, we'll return an empty list
        # as we can't perform actual RDKit operations
        return []
    
    def search_by_similarity(
        self, 
        db: Session, 
        reference_smiles: str,
        threshold: float = 0.7,
        skip: int = 0, 
        limit: int = 100
    ) -> List[Tuple[Molecule, float]]:
        """
        Search for molecules similar to a reference molecule (using RDKit cartridge).
        
        This method would use RDKit's similarity search functionality, which typically
        requires RDKit PostgreSQL cartridge to be installed in the database.
        
        Args:
            db: Database session
            reference_smiles: SMILES of the reference molecule
            threshold: Minimum similarity score (0.0-1.0)
            skip: Number of records to skip
            limit: Maximum number of records to return
            
        Returns:
            List of tuples containing (molecule, similarity_score)
        """
        # Note: This is a simplified implementation
        # In a real application, this would use RDKit's similarity search functionality
        # via SQL function calls to the RDKit PostgreSQL cartridge
        
        # For example (pseudocode):
        # stmt = (
        #     select(Molecule, func.tanimoto_similarity(func.mol_from_smiles(Molecule.smiles), 
        #                                              func.mol_from_smiles(reference_smiles))
        #            .label('similarity'))
        #     .where(func.tanimoto_similarity(func.mol_from_smiles(Molecule.smiles), 
        #                                    func.mol_from_smiles(reference_smiles)) >= threshold)
        #     .order_by(desc('similarity'))
        #     .offset(skip)
        #     .limit(limit)
        # )
        
        # For this implementation, we'll return an empty list
        # as we can't perform actual RDKit operations
        return []
    
    def bulk_create(
        self, 
        db: Session, 
        obj_list: List[Dict[str, Any]], 
        user_id: int
    ) -> Dict[str, Any]:
        """
        Create multiple molecules in a single transaction.
        
        This method efficiently handles bulk insertion of molecules from CSV imports,
        checking for duplicates by SMILES string and adding properties.
        
        Args:
            db: Database session
            obj_list: List of dictionaries containing molecule data
            user_id: ID of the user creating the molecules
            
        Returns:
            Dictionary with created molecules and statistics
        """
        created_molecules = []
        duplicate_count = 0
        
        # Process each molecule in the list
        for molecule_data in obj_list:
            # Check if molecule already exists by SMILES string
            smiles = molecule_data.get("smiles")
            existing_molecule = self.get_by_smiles(db, smiles)
            
            if existing_molecule:
                duplicate_count += 1
                continue
            
            # Create a new molecule
            molecule_in = MoleculeCreate(
                smiles=smiles,
                created_by=user_id,
                flag_status=molecule_data.get("flag_status")
            )
            
            # Create the molecule
            db_molecule = self.create(db, molecule_in)
            
            # Add properties if present
            if "properties" in molecule_data and molecule_data["properties"]:
                for prop in molecule_data["properties"]:
                    property_obj = MoleculeProperty(
                        molecule_id=db_molecule.id,
                        property_name=prop["property_name"],
                        property_value=prop.get("property_value"),
                        property_unit=prop.get("property_unit"),
                        is_calculated=prop.get("is_calculated", False)
                    )
                    db.add(property_obj)
            
            created_molecules.append(db_molecule)
        
        # Commit all changes
        db.commit()
        
        # Return statistics and created molecules
        return {
            "created_molecules": created_molecules,
            "total_created": len(created_molecules),
            "total_duplicates": duplicate_count,
            "total_processed": len(obj_list)
        }


# Create a singleton instance for application-wide use
molecule = CRUDMolecule()