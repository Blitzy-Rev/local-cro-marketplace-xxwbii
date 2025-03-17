"""
Molecule Service Module

This module provides the service layer for molecule management operations in the 
Molecular Data Management and CRO Integration Platform. It handles business logic 
for molecule creation, retrieval, filtering, organization, and property calculation,
serving as an intermediary between API endpoints and the data access layer.
"""

from typing import List, Dict, Optional, Union, Any, Tuple
import io
from io import BytesIO
import uuid
from uuid import UUID

from rdkit import Chem  # version 2023.03+
from rdkit.Chem import Draw  # version 2023.03+

from ..crud.crud_molecule import molecule
from ..models.molecule import Molecule
from ..models.molecule_property import MoleculeProperty
from ..schemas.molecule import MoleculeCreate, MoleculeUpdate, MoleculeFilter
from ..molecular.processor import (
    process_molecule, 
    process_molecules_batch, 
    process_molecules_stream, 
    standardize_molecules,
    remove_duplicates,
    check_drug_likeness
)
from .file_storage_service import upload_molecule_image, get_molecule_image_url
from ..db.session import get_db
from ..exceptions import MoleculeServiceException
from ..logging_config import logger
from ..constants import BATCH_SIZE


def get_molecule(molecule_id: Any) -> Optional[Dict[str, Any]]:
    """
    Retrieves a molecule by ID with its properties.
    
    Args:
        molecule_id: The ID of the molecule to retrieve
        
    Returns:
        Optional[Dict[str, Any]]: Molecule data with properties if found, None otherwise
    """
    with get_db() as db:
        db_molecule = molecule.get_with_properties(db, molecule_id)
        if not db_molecule:
            return None
        
        # Convert molecule model to dictionary with properties
        result = {
            "id": db_molecule.id,
            "smiles": db_molecule.smiles,
            "flag_status": db_molecule.flag_status,
            "created_at": db_molecule.created_at,
            "updated_at": db_molecule.updated_at,
            "created_by": db_molecule.created_by,
            "properties": []
        }
        
        # Add properties if available
        if hasattr(db_molecule, 'properties') and db_molecule.properties:
            result["properties"] = [
                {
                    "property_name": prop.property_name,
                    "property_value": prop.property_value,
                    "property_unit": prop.property_unit,
                    "is_calculated": prop.is_calculated
                }
                for prop in db_molecule.properties
            ]
            
        return result


def get_molecule_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    """
    Retrieves a molecule by its SMILES string.
    
    Args:
        smiles: SMILES representation of the molecule
        
    Returns:
        Optional[Dict[str, Any]]: Molecule data if found, None otherwise
    """
    with get_db() as db:
        db_molecule = molecule.get_by_smiles(db, smiles)
        if not db_molecule:
            return None
        
        # Convert molecule model to dictionary
        result = {
            "id": db_molecule.id,
            "smiles": db_molecule.smiles,
            "flag_status": db_molecule.flag_status,
            "created_at": db_molecule.created_at,
            "updated_at": db_molecule.updated_at,
            "created_by": db_molecule.created_by
        }
        
        return result


def create_molecule(molecule_data: MoleculeCreate, user_id: UUID) -> Dict[str, Any]:
    """
    Creates a new molecule with properties.
    
    Args:
        molecule_data: Data for the new molecule
        user_id: ID of the user creating the molecule
        
    Returns:
        Dict[str, Any]: Created molecule data
        
    Raises:
        MoleculeServiceException: If molecule already exists or validation fails
    """
    with get_db() as db:
        # Check if molecule already exists by SMILES
        existing_molecule = molecule.get_by_smiles(db, molecule_data.smiles)
        if existing_molecule:
            raise MoleculeServiceException(
                f"Molecule with SMILES '{molecule_data.smiles}' already exists",
                {"smiles": molecule_data.smiles}
            )
        
        # Process molecule to validate structure and calculate properties
        try:
            processed_data = process_molecule(molecule_data.smiles)
            if not processed_data.get('validation', {}).get('structure_valid', False):
                raise MoleculeServiceException(
                    f"Invalid molecular structure: {molecule_data.smiles}",
                    {"smiles": molecule_data.smiles}
                )
        except Exception as e:
            logger.error(f"Error processing molecule: {str(e)}")
            raise MoleculeServiceException(
                f"Error processing molecule: {str(e)}",
                {"smiles": molecule_data.smiles}
            )
        
        # Create molecule object
        molecule_in = MoleculeCreate(
            smiles=molecule_data.smiles,
            flag_status=molecule_data.flag_status,
            created_by=user_id
        )
        
        # Add to database
        db_molecule = molecule.create(db, obj_in=molecule_in)
        
        # Add properties if provided
        if molecule_data.properties:
            for prop in molecule_data.properties:
                property_obj = MoleculeProperty(
                    molecule_id=db_molecule.id,
                    property_name=prop.property_name,
                    property_value=prop.property_value,
                    property_unit=prop.property_unit,
                    is_calculated=prop.is_calculated
                )
                db.add(property_obj)
                
            db.commit()
            db.refresh(db_molecule)
        
        # Convert created molecule to dictionary
        result = {
            "id": db_molecule.id,
            "smiles": db_molecule.smiles,
            "flag_status": db_molecule.flag_status,
            "created_at": db_molecule.created_at,
            "updated_at": db_molecule.updated_at,
            "created_by": db_molecule.created_by,
            "properties": []
        }
        
        # Add properties if available
        if hasattr(db_molecule, 'properties') and db_molecule.properties:
            result["properties"] = [
                {
                    "property_name": prop.property_name,
                    "property_value": prop.property_value,
                    "property_unit": prop.property_unit,
                    "is_calculated": prop.is_calculated
                }
                for prop in db_molecule.properties
            ]
            
        return result


def update_molecule(molecule_id: Any, molecule_data: MoleculeUpdate) -> Dict[str, Any]:
    """
    Updates an existing molecule's properties or flag status.
    
    Args:
        molecule_id: ID of the molecule to update
        molecule_data: Updated molecule data
        
    Returns:
        Dict[str, Any]: Updated molecule data
        
    Raises:
        MoleculeServiceException: If molecule not found
    """
    with get_db() as db:
        # Get existing molecule
        db_molecule = molecule.get_with_properties(db, molecule_id)
        if not db_molecule:
            raise MoleculeServiceException(
                f"Molecule with ID {molecule_id} not found",
                {"molecule_id": molecule_id}
            )
        
        # Update dictionary with new data
        update_dict = {}
        
        # Update flag status if provided
        if molecule_data.flag_status is not None:
            update_dict["flag_status"] = molecule_data.flag_status
        
        # Apply updates if any
        if update_dict:
            db_molecule = molecule.update(db, db_obj=db_molecule, obj_in=update_dict)
        
        # Update properties if provided
        if molecule_data.properties:
            # Create a map of existing properties for efficient lookup
            existing_properties = {
                prop.property_name: prop for prop in db_molecule.properties
            }
            
            for prop in molecule_data.properties:
                if prop.property_name in existing_properties:
                    # Update existing property
                    existing_prop = existing_properties[prop.property_name]
                    existing_prop.property_value = prop.property_value
                    existing_prop.property_unit = prop.property_unit
                    existing_prop.is_calculated = prop.is_calculated
                    db.add(existing_prop)
                else:
                    # Add new property
                    property_obj = MoleculeProperty(
                        molecule_id=db_molecule.id,
                        property_name=prop.property_name,
                        property_value=prop.property_value,
                        property_unit=prop.property_unit,
                        is_calculated=prop.is_calculated
                    )
                    db.add(property_obj)
            
            db.commit()
            db.refresh(db_molecule)
        
        # Convert updated molecule to dictionary
        result = {
            "id": db_molecule.id,
            "smiles": db_molecule.smiles,
            "flag_status": db_molecule.flag_status,
            "created_at": db_molecule.created_at,
            "updated_at": db_molecule.updated_at,
            "created_by": db_molecule.created_by,
            "properties": []
        }
        
        # Add properties if available
        if hasattr(db_molecule, 'properties') and db_molecule.properties:
            result["properties"] = [
                {
                    "property_name": prop.property_name,
                    "property_value": prop.property_value,
                    "property_unit": prop.property_unit,
                    "is_calculated": prop.is_calculated
                }
                for prop in db_molecule.properties
            ]
            
        return result


def delete_molecule(molecule_id: Any) -> bool:
    """
    Deletes a molecule by ID.
    
    Args:
        molecule_id: ID of the molecule to delete
        
    Returns:
        bool: True if molecule was deleted successfully
    """
    with get_db() as db:
        deleted = molecule.remove(db, molecule_id)
        return deleted is not None


def get_molecules(filters: Optional[MoleculeFilter] = None, skip: int = 0, limit: int = 100) -> Dict[str, Any]:
    """
    Retrieves molecules with optional filtering and pagination.
    
    Args:
        filters: Filter criteria for molecules
        skip: Number of records to skip (for pagination)
        limit: Maximum number of records to return
        
    Returns:
        Dict[str, Any]: Dictionary with molecules list, total count, and pagination info
    """
    with get_db() as db:
        query_filters = {}
        total = 0
        
        if filters:
            # Apply SMILES filter if provided
            if filters.smiles:
                # This would be a partial match in a real implementation
                # For this service layer, we'll use exact match
                smiles_molecule = molecule.get_by_smiles(db, filters.smiles)
                if smiles_molecule:
                    return {
                        "items": [get_molecule(smiles_molecule.id)],
                        "total": 1,
                        "skip": 0,
                        "limit": 1
                    }
                else:
                    return {
                        "items": [],
                        "total": 0,
                        "skip": skip,
                        "limit": limit
                    }
            
            # Apply user filter if provided
            if filters.created_by:
                user_molecules = molecule.get_by_user(db, filters.created_by, skip, limit)
                total = len(user_molecules)  # In a real implementation, would use a count query
                result = {
                    "items": [get_molecule(mol.id) for mol in user_molecules],
                    "total": total,
                    "skip": skip,
                    "limit": limit
                }
                return result
            
            # Apply flag status filter if provided
            if filters.flag_status:
                flag_molecules = molecule.get_by_flag_status(db, filters.flag_status, skip, limit)
                total = len(flag_molecules)  # In a real implementation, would use a count query
                result = {
                    "items": [get_molecule(mol.id) for mol in flag_molecules],
                    "total": total,
                    "skip": skip,
                    "limit": limit
                }
                return result
            
            # Apply property filters if provided
            if filters.properties:
                for prop_filter in filters.properties:
                    property_name = prop_filter.get("name")
                    min_value = prop_filter.get("min_value")
                    max_value = prop_filter.get("max_value")
                    
                    if property_name:
                        # Get molecules by property range
                        property_molecules = molecule.get_by_property_range(
                            db, 
                            property_name, 
                            min_value, 
                            max_value, 
                            skip, 
                            limit,
                            filters.sort_by,
                            filters.sort_desc
                        )
                        
                        # Get count for the filter
                        total = molecule.count_by_property_range(
                            db, 
                            property_name, 
                            min_value, 
                            max_value
                        )
                        
                        result = {
                            "items": [get_molecule(mol.id) for mol in property_molecules],
                            "total": total,
                            "skip": skip,
                            "limit": limit
                        }
                        return result
        
        # If no specific filters applied, get all molecules with pagination
        all_molecules = molecule.get_multi(db, skip=skip, limit=limit)
        total = molecule.count(db)
        
        result = {
            "items": [get_molecule(mol.id) for mol in all_molecules],
            "total": total,
            "skip": skip,
            "limit": limit
        }
        
        return result


def update_molecule_flag(molecule_id: Any, flag_status: str) -> Dict[str, Any]:
    """
    Updates the flag status of a molecule.
    
    Args:
        molecule_id: ID of the molecule to update
        flag_status: New flag status value
        
    Returns:
        Dict[str, Any]: Updated molecule data
        
    Raises:
        MoleculeServiceException: If molecule not found
    """
    with get_db() as db:
        db_molecule = molecule.update_flag_status(db, molecule_id, flag_status)
        
        if not db_molecule:
            raise MoleculeServiceException(
                f"Molecule with ID {molecule_id} not found",
                {"molecule_id": molecule_id}
            )
        
        # Convert updated molecule to dictionary
        result = {
            "id": db_molecule.id,
            "smiles": db_molecule.smiles,
            "flag_status": db_molecule.flag_status,
            "created_at": db_molecule.created_at,
            "updated_at": db_molecule.updated_at,
            "created_by": db_molecule.created_by
        }
        
        return result


def get_property_ranges(property_name: str) -> Dict[str, float]:
    """
    Gets the minimum and maximum values for a specific property.
    
    Args:
        property_name: Name of the property to get range for
        
    Returns:
        Dict[str, float]: Dictionary with min and max values for the property
    """
    with get_db() as db:
        min_val, max_val = molecule.get_property_ranges(db, property_name)
        return {
            "min": min_val,
            "max": max_val
        }


def get_available_properties() -> List[str]:
    """
    Gets a list of all property names available in the database.
    
    Returns:
        List[str]: List of distinct property names
    """
    with get_db() as db:
        properties = molecule.get_available_properties(db)
        return properties


def bulk_create_molecules(
    molecules_data: List[Dict[str, Any]], 
    user_id: UUID, 
    process_properties: bool = True,
    batch_size: int = None
) -> Dict[str, Any]:
    """
    Creates multiple molecules in a single transaction.
    
    Args:
        molecules_data: List of molecule data dictionaries
        user_id: ID of the user creating the molecules
        process_properties: Whether to validate and calculate properties
        batch_size: Size of each processing batch
        
    Returns:
        Dict[str, Any]: Dictionary with created molecules and statistics
    """
    # Set default batch size if not provided
    if batch_size is None:
        batch_size = BATCH_SIZE
    
    # Extract SMILES strings from data
    smiles_list = [mol_data.get("smiles") for mol_data in molecules_data if mol_data.get("smiles")]
    
    # Process molecules if requested
    if process_properties:
        logger.info(f"Processing {len(smiles_list)} molecules in batches of {batch_size}")
        
        # Process molecules in stream for validation and property calculation
        processed_molecules = process_molecules_stream(smiles_list, batch_size=batch_size)
        
        # Standardize SMILES to canonical form
        standardized_molecules = standardize_molecules(processed_molecules)
        
        # Remove duplicates
        unique_molecules, duplicates = remove_duplicates(standardized_molecules)
        
        # Prepare molecules for database insertion
        molecules_for_db = []
        for mol in unique_molecules:
            if mol.get('validation', {}).get('structure_valid', False):
                properties = []
                
                # Convert properties to required format
                if 'properties' in mol and mol['properties']:
                    for prop_name, prop_value in mol['properties'].items():
                        properties.append({
                            "property_name": prop_name,
                            "property_value": prop_value,
                            "property_unit": None,  # Unit would be determined by property type
                            "is_calculated": True   # These properties were calculated
                        })
                
                molecules_for_db.append({
                    "smiles": mol["smiles"],
                    "created_by": user_id,
                    "flag_status": None,
                    "properties": properties
                })
    else:
        # Use the raw molecule data without processing
        molecules_for_db = []
        for mol_data in molecules_data:
            if "smiles" in mol_data:
                properties = []
                
                # Convert properties to required format if present
                if "properties" in mol_data and mol_data["properties"]:
                    for prop in mol_data["properties"]:
                        if isinstance(prop, dict) and "property_name" in prop:
                            properties.append({
                                "property_name": prop["property_name"],
                                "property_value": prop.get("property_value"),
                                "property_unit": prop.get("property_unit"),
                                "is_calculated": prop.get("is_calculated", False)
                            })
                
                molecules_for_db.append({
                    "smiles": mol_data["smiles"],
                    "created_by": user_id,
                    "flag_status": mol_data.get("flag_status"),
                    "properties": properties
                })
    
    # Create molecules in database
    with get_db() as db:
        logger.info(f"Creating {len(molecules_for_db)} molecules in database")
        result = molecule.bulk_create(db, molecules_for_db, user_id)
        
        return {
            "created_molecules": result["created_molecules"],
            "total_created": result["total_created"],
            "total_duplicates": result["total_duplicates"],
            "total_processed": result["total_processed"],
            "total_failed": len(smiles_list) - result["total_created"] - result["total_duplicates"]
        }


def generate_molecule_image(smiles: str, width: int = 300, height: int = 200, store_image: bool = True) -> Dict[str, Any]:
    """
    Generates and stores a 2D image of a molecule structure.
    
    Args:
        smiles: SMILES representation of the molecule
        width: Image width in pixels
        height: Image height in pixels
        store_image: Whether to store the image in file storage
        
    Returns:
        Dict[str, Any]: Dictionary with image URL or data
        
    Raises:
        MoleculeServiceException: If molecule structure is invalid
    """
    try:
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise MoleculeServiceException(
                f"Invalid SMILES structure: {smiles}",
                {"smiles": smiles}
            )
        
        # Generate 2D coordinates
        mol = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol, Chem.AllChem.ETKDGV3())
        mol = Chem.RemoveHs(mol)
        
        # Generate molecule image
        img = Draw.MolToImage(mol, size=(width, height))
        
        # Convert image to bytes
        img_bytes = BytesIO()
        img.save(img_bytes, format='PNG')
        img_bytes.seek(0)
        
        if store_image:
            # Generate unique filename
            filename = f"{uuid.uuid4()}.png"
            
            # Upload image to storage
            object_name = upload_molecule_image(img_bytes, filename)
            
            # Get URL for the image
            image_url = get_molecule_image_url(object_name)
            
            return {
                "image_url": image_url,
                "width": width,
                "height": height,
                "smiles": smiles
            }
        else:
            # Return image data directly
            import base64
            img_base64 = base64.b64encode(img_bytes.getvalue()).decode('utf-8')
            
            return {
                "image_data": f"data:image/png;base64,{img_base64}",
                "width": width,
                "height": height,
                "smiles": smiles
            }
            
    except Exception as e:
        logger.error(f"Error generating molecule image: {str(e)}")
        raise MoleculeServiceException(
            f"Error generating molecule image: {str(e)}",
            {"smiles": smiles}
        )


def check_molecule_drug_likeness(smiles: str) -> Dict[str, Any]:
    """
    Checks if a molecule follows Lipinski's Rule of Five for drug-likeness.
    
    Args:
        smiles: SMILES representation of the molecule
        
    Returns:
        Dict[str, Any]: Dictionary with drug-likeness assessment
        
    Raises:
        MoleculeServiceException: If molecule structure is invalid
    """
    try:
        # Process molecule to validate structure
        processed_data = process_molecule(smiles)
        if not processed_data.get('validation', {}).get('structure_valid', False):
            raise MoleculeServiceException(
                f"Invalid molecular structure: {smiles}",
                {"smiles": smiles}
            )
        
        # Convert SMILES to RDKit molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise MoleculeServiceException(
                f"Invalid SMILES structure: {smiles}",
                {"smiles": smiles}
            )
        
        # Check drug-likeness
        lipinski_results = check_drug_likeness(mol)
        
        return {
            "smiles": smiles,
            "molecular_weight": lipinski_results["molecular_weight"],
            "logp": lipinski_results["logp"],
            "h_bond_donors": lipinski_results["h_bond_donors"],
            "h_bond_acceptors": lipinski_results["h_bond_acceptors"],
            "violations": lipinski_results["violations"],
            "is_drug_like": lipinski_results["is_drug_like"],
            "details": {
                "mw_pass": lipinski_results["mw_pass"],
                "logp_pass": lipinski_results["logp_pass"],
                "donors_pass": lipinski_results["donors_pass"],
                "acceptors_pass": lipinski_results["acceptors_pass"]
            }
        }
        
    except MoleculeServiceException:
        # Re-raise MoleculeServiceException
        raise
    except Exception as e:
        logger.error(f"Error checking drug-likeness: {str(e)}")
        raise MoleculeServiceException(
            f"Error checking drug-likeness: {str(e)}",
            {"smiles": smiles}
        )