from typing import List, Dict, Optional, Union, Any, Tuple
from uuid import UUID
from datetime import datetime

from ..crud.crud_library import library
from ..crud.crud_molecule import molecule
from ..models.library import Library
from ..schemas.library import LibraryCreate, LibraryUpdate, LibraryFilter
from ..db.session import get_db
from ..exceptions import LibraryServiceException
from ..logging_config import logger


def get_library(library_id: Any) -> Optional[Dict[str, Any]]:
    """
    Retrieves a library by ID.
    
    Args:
        library_id: Unique identifier of the library
        
    Returns:
        Library data if found, None otherwise
    """
    with get_db() as db:
        db_library = library.get(db, library_id)
        if not db_library:
            return None
        
        # Convert library to dictionary
        return {
            "id": db_library.id,
            "name": db_library.name,
            "description": db_library.description,
            "created_by": db_library.created_by,
            "created_at": db_library.created_at,
            "updated_at": db_library.updated_at
        }


def get_library_by_name(name: str, user_id: UUID) -> Optional[Dict[str, Any]]:
    """
    Retrieves a library by name and creator ID.
    
    Args:
        name: Name of the library
        user_id: ID of the user who created the library
        
    Returns:
        Library data if found, None otherwise
    """
    with get_db() as db:
        db_library = library.get_by_name(db, name, user_id)
        if not db_library:
            return None
        
        # Convert library to dictionary
        return {
            "id": db_library.id,
            "name": db_library.name,
            "description": db_library.description,
            "created_by": db_library.created_by,
            "created_at": db_library.created_at,
            "updated_at": db_library.updated_at
        }


def create_library(library_data: LibraryCreate, user_id: UUID) -> Dict[str, Any]:
    """
    Creates a new library.
    
    Args:
        library_data: Library creation data
        user_id: ID of the user creating the library
        
    Returns:
        Created library data
        
    Raises:
        LibraryServiceException: If a library with the same name already exists for this user
    """
    with get_db() as db:
        # Check if library with same name already exists for this user
        existing_library = library.get_by_name(db, library_data.name, user_id)
        if existing_library:
            raise LibraryServiceException(
                f"Library with name '{library_data.name}' already exists for this user",
                {"name": library_data.name, "user_id": str(user_id)}
            )
        
        # Create library object with user ID
        lib_in = LibraryCreate(
            name=library_data.name,
            description=library_data.description,
            created_by=user_id
        )
        
        # Create library in database
        db_library = library.create(db, lib_in)
        
        # Return created library data
        return {
            "id": db_library.id,
            "name": db_library.name,
            "description": db_library.description,
            "created_by": db_library.created_by,
            "created_at": db_library.created_at,
            "updated_at": db_library.updated_at
        }


def update_library(library_id: Any, library_data: LibraryUpdate, user_id: UUID) -> Dict[str, Any]:
    """
    Updates an existing library.
    
    Args:
        library_id: ID of the library to update
        library_data: Updated library data
        user_id: ID of the user performing the update
        
    Returns:
        Updated library data
        
    Raises:
        LibraryServiceException: If library not found or user is not authorized
    """
    with get_db() as db:
        # Get the library
        db_library = library.get(db, library_id)
        if not db_library:
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        # Check if user is authorized to update this library
        if db_library.created_by != user_id:
            raise LibraryServiceException(
                "User not authorized to update this library",
                {"library_id": str(library_id), "user_id": str(user_id)}
            )
        
        # If name is being updated, check if new name already exists for this user
        if library_data.name and library_data.name != db_library.name:
            existing_library = library.get_by_name(db, library_data.name, user_id)
            if existing_library and existing_library.id != library_id:
                raise LibraryServiceException(
                    f"Library with name '{library_data.name}' already exists for this user",
                    {"name": library_data.name, "user_id": str(user_id)}
                )
        
        # Update library
        updated_library = library.update(db, db_library, library_data)
        
        # Return updated library data
        return {
            "id": updated_library.id,
            "name": updated_library.name,
            "description": updated_library.description,
            "created_by": updated_library.created_by,
            "created_at": updated_library.created_at,
            "updated_at": updated_library.updated_at
        }


def delete_library(library_id: Any, user_id: UUID) -> bool:
    """
    Deletes a library by ID.
    
    Args:
        library_id: ID of the library to delete
        user_id: ID of the user performing the delete
        
    Returns:
        True if library was deleted successfully
        
    Raises:
        LibraryServiceException: If library not found or user is not authorized
    """
    with get_db() as db:
        # Get the library
        db_library = library.get(db, library_id)
        if not db_library:
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        # Check if user is authorized to delete this library
        if db_library.created_by != user_id:
            raise LibraryServiceException(
                "User not authorized to delete this library",
                {"library_id": str(library_id), "user_id": str(user_id)}
            )
        
        # Delete the library
        library.remove(db, library_id)
        return True


def get_libraries(filters: Optional[LibraryFilter] = None, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves libraries with optional filtering and pagination.
    
    Args:
        filters: Optional filter criteria
        skip: Number of records to skip
        limit: Maximum number of records to return
        
    Returns:
        Dictionary with libraries list, total count, and pagination info
    """
    with get_db() as db:
        # Initialize filter parameters
        name_contains = None
        created_by = None
        contains_molecule_id = None
        sort_by = None
        sort_desc = False
        
        if filters:
            name_contains = filters.name
            created_by = filters.created_by
            contains_molecule_id = filters.contains_molecule_id
            sort_by = filters.sort_by
            sort_desc = filters.sort_desc or False
        
        # If filtering by molecule, check if molecule exists
        if contains_molecule_id:
            mol = molecule.get(db, contains_molecule_id)
            if not mol:
                return {
                    "items": [],
                    "total": 0,
                    "skip": skip,
                    "limit": limit
                }
        
        # Get libraries with filtering
        libraries_list, total = library.filter_libraries(
            db,
            name_contains=name_contains,
            created_by=created_by,
            skip=skip,
            limit=limit,
            sort_by=sort_by,
            sort_desc=sort_desc
        )
        
        # Add molecule count to each library
        result_libraries = []
        for lib in libraries_list:
            lib_dict = {
                "id": lib.id,
                "name": lib.name,
                "description": lib.description,
                "created_by": lib.created_by,
                "created_at": lib.created_at,
                "updated_at": lib.updated_at,
                "molecule_count": library.count_molecules(db, lib.id)
            }
            result_libraries.append(lib_dict)
        
        return {
            "items": result_libraries,
            "total": total,
            "skip": skip,
            "limit": limit
        }


def get_user_libraries(user_id: UUID, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves libraries created by a specific user.
    
    Args:
        user_id: ID of the user who created the libraries
        skip: Number of records to skip
        limit: Maximum number of records to return
        
    Returns:
        Dictionary with user's libraries and count
    """
    with get_db() as db:
        # Get libraries created by user
        user_libraries = library.get_by_user(db, user_id, skip, limit)
        total = library.count_by_user(db, user_id)
        
        # Add molecule count to each library
        result_libraries = []
        for lib in user_libraries:
            lib_dict = {
                "id": lib.id,
                "name": lib.name,
                "description": lib.description,
                "created_by": lib.created_by,
                "created_at": lib.created_at,
                "updated_at": lib.updated_at,
                "molecule_count": library.count_molecules(db, lib.id)
            }
            result_libraries.append(lib_dict)
        
        return {
            "items": result_libraries,
            "total": total,
            "skip": skip,
            "limit": limit
        }


def get_library_with_molecules(library_id: Any, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves a library with its associated molecules.
    
    Args:
        library_id: ID of the library to retrieve
        skip: Number of molecules to skip
        limit: Maximum number of molecules to return
        
    Returns:
        Library data with molecules if found
        
    Raises:
        LibraryServiceException: If library not found
    """
    with get_db() as db:
        # Get library with molecules
        result = library.get_with_molecules(db, library_id, skip, limit)
        if not result:
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        lib, molecules_list = result
        
        # Convert library to dictionary
        lib_dict = {
            "id": lib.id,
            "name": lib.name,
            "description": lib.description,
            "created_by": lib.created_by,
            "created_at": lib.created_at,
            "updated_at": lib.updated_at
        }
        
        # Convert molecules to list of dictionaries
        molecules_data = []
        for mol in molecules_list:
            mol_dict = {
                "id": mol.id,
                "smiles": mol.smiles,
                "created_by": mol.created_by,
                "created_at": mol.created_at,
                "updated_at": mol.updated_at,
                "flag_status": mol.flag_status
            }
            molecules_data.append(mol_dict)
        
        # Get total molecule count
        total_molecules = library.count_molecules(db, library_id)
        
        return {
            "library": lib_dict,
            "molecules": molecules_data,
            "total_molecules": total_molecules,
            "skip": skip,
            "limit": limit
        }


def add_molecules_to_library(library_id: Any, molecule_ids: List[Any], user_id: UUID) -> Dict[str, Any]:
    """
    Adds molecules to a library.
    
    Args:
        library_id: ID of the library to add molecules to
        molecule_ids: List of molecule IDs to add
        user_id: ID of the user performing the operation
        
    Returns:
        Result of the operation with success count and failures
        
    Raises:
        LibraryServiceException: If library not found or user is not authorized
    """
    with get_db() as db:
        # Get the library
        db_library = library.get(db, library_id)
        if not db_library:
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        # Check if user is authorized to modify this library
        if db_library.created_by != user_id:
            raise LibraryServiceException(
                "User not authorized to modify this library",
                {"library_id": str(library_id), "user_id": str(user_id)}
            )
        
        # Validate that all molecule IDs exist
        invalid_molecules = []
        for mol_id in molecule_ids:
            if not molecule.get(db, mol_id):
                invalid_molecules.append(str(mol_id))
        
        if invalid_molecules:
            raise LibraryServiceException(
                "Some molecules do not exist",
                {"molecule_ids": invalid_molecules}
            )
        
        # Add molecules to library
        result = library.add_molecules(db, library_id, molecule_ids)
        
        return result


def remove_molecules_from_library(library_id: Any, molecule_ids: List[Any], user_id: UUID) -> Dict[str, Any]:
    """
    Removes molecules from a library.
    
    Args:
        library_id: ID of the library to remove molecules from
        molecule_ids: List of molecule IDs to remove
        user_id: ID of the user performing the operation
        
    Returns:
        Result of the operation with success count and failures
        
    Raises:
        LibraryServiceException: If library not found or user is not authorized
    """
    with get_db() as db:
        # Get the library
        db_library = library.get(db, library_id)
        if not db_library:
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        # Check if user is authorized to modify this library
        if db_library.created_by != user_id:
            raise LibraryServiceException(
                "User not authorized to modify this library",
                {"library_id": str(library_id), "user_id": str(user_id)}
            )
        
        # Remove molecules from library
        result = library.remove_molecules(db, library_id, molecule_ids)
        
        return result


def check_molecule_in_library(library_id: Any, molecule_id: Any) -> bool:
    """
    Checks if a molecule is in a library.
    
    Args:
        library_id: ID of the library to check
        molecule_id: ID of the molecule to check
        
    Returns:
        True if the molecule is in the library, False otherwise
    """
    with get_db() as db:
        return library.is_molecule_in_library(db, library_id, molecule_id)


def get_library_molecules(library_id: Any, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Gets molecules in a library with pagination.
    
    Args:
        library_id: ID of the library to get molecules from
        skip: Number of molecules to skip
        limit: Maximum number of molecules to return
        
    Returns:
        Dictionary with molecules list and total count
        
    Raises:
        LibraryServiceException: If library not found
    """
    with get_db() as db:
        # Check if library exists
        if not library.get(db, library_id):
            raise LibraryServiceException(
                f"Library not found: {library_id}",
                {"library_id": str(library_id)}
            )
        
        # Get molecules in library
        molecules_list = library.get_molecules(db, library_id, skip, limit)
        total = library.count_molecules(db, library_id)
        
        # Convert molecules to list of dictionaries
        molecules_data = []
        for mol in molecules_list:
            mol_dict = {
                "id": mol.id,
                "smiles": mol.smiles,
                "created_by": mol.created_by,
                "created_at": mol.created_at,
                "updated_at": mol.updated_at,
                "flag_status": mol.flag_status
            }
            molecules_data.append(mol_dict)
        
        return {
            "items": molecules_data,
            "total": total,
            "skip": skip,
            "limit": limit
        }