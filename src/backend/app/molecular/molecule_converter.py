"""
Molecule Converter Module

This module provides utilities for converting between different molecular representations
in the Molecular Data Management and CRO Integration Platform. It includes functions for:

- Converting SMILES strings to RDKit molecule objects and vice versa
- Generating molecular structure images in various formats
- Converting molecules to different chemical identifiers (InChI, InChIKey)
- Generating molecular fingerprints for similarity calculations
- Converting between various molecular file formats

All functions include robust error handling with detailed error messages to help
diagnose issues with molecular processing.
"""

from typing import Optional, Union, Dict, List, Tuple, Any
import io
import base64
import logging

# RDKit version 2023.03+
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolDescriptors

from ..exceptions import MolecularProcessingException

# Set up module logger
logger = logging.getLogger(__name__)

# Dictionary mapping fingerprint types to their respective RDKit functions
FINGERPRINT_TYPES = {
    'morgan': rdMolDescriptors.GetMorganFingerprintAsBitVect,
    'maccs': rdMolDescriptors.GetMACCSKeysFingerprint,
    'topological': rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect,
    'pattern': rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect
}


def smiles_to_mol(smiles: str) -> Optional[Chem.rdchem.Mol]:
    """
    Converts a SMILES string to an RDKit molecule object.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        RDKit molecule object or None if conversion fails
    """
    if not smiles:
        return None
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol
    except Exception as e:
        logger.error(f"Error converting SMILES to molecule: {str(e)}")
        return None


def smiles_to_mol_with_exception(smiles: str) -> Chem.rdchem.Mol:
    """
    Converts a SMILES string to an RDKit molecule object with exception handling.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        RDKit molecule object
        
    Raises:
        MolecularProcessingException: If SMILES conversion fails
    """
    if not smiles:
        raise MolecularProcessingException(
            "Cannot convert empty SMILES string to molecule",
            details={"smiles": smiles}
        )
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise MolecularProcessingException(
            "Invalid SMILES string",
            details={"smiles": smiles}
        )
    
    return mol


def mol_to_smiles(mol: Chem.rdchem.Mol, canonical: bool = True) -> str:
    """
    Converts an RDKit molecule object to a SMILES string.
    
    Args:
        mol: RDKit molecule object
        canonical: Whether to return the canonical SMILES (default: True)
        
    Returns:
        SMILES string representation of the molecule
        
    Raises:
        MolecularProcessingException: If molecule conversion fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot convert None to SMILES",
            details={"mol": "None"}
        )
    
    try:
        smiles = Chem.MolToSmiles(mol, canonical=canonical)
        return smiles
    except Exception as e:
        raise MolecularProcessingException(
            f"Error converting molecule to SMILES: {str(e)}",
            details={"error": str(e)}
        )


def canonicalize_smiles(smiles: str) -> str:
    """
    Converts a SMILES string to its canonical form.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Canonical SMILES string
        
    Raises:
        MolecularProcessingException: If SMILES conversion fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_smiles(mol, canonical=True)


def mol_to_image(
    mol: Chem.rdchem.Mol,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_atoms: List[int] = None,
    highlight_bonds: List[int] = None
) -> bytes:
    """
    Converts an RDKit molecule to an image.
    
    Args:
        mol: RDKit molecule object
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_atoms: List of atom indices to highlight (default: None)
        highlight_bonds: List of bond indices to highlight (default: None)
        
    Returns:
        Image data in the specified format
        
    Raises:
        MolecularProcessingException: If image generation fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot generate image from None molecule",
            details={"mol": "None"}
        )
    
    try:
        if highlight_atoms or highlight_bonds:
            # Create highlighting
            highlights = {}
            if highlight_atoms:
                highlights['highlightAtoms'] = highlight_atoms
            if highlight_bonds:
                highlights['highlightBonds'] = highlight_bonds
            
            drawer = Draw.MolDraw2DCairo(width, height)
            drawer.DrawMolecule(mol, **highlights)
            drawer.FinishDrawing()
            img_data = drawer.GetDrawingText()
            return img_data
        else:
            # Standard drawing without highlighting
            if format.lower() == 'svg':
                drawer = Draw.MolDraw2DSVG(width, height)
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                return drawer.GetDrawingText().encode('utf-8')
            else:  # PNG or other formats
                img = Draw.MolToImage(mol, size=(width, height))
                img_bytesio = io.BytesIO()
                img.save(img_bytesio, format=format.upper())
                return img_bytesio.getvalue()
    except Exception as e:
        raise MolecularProcessingException(
            f"Error generating molecule image: {str(e)}",
            details={"error": str(e)}
        )


def smiles_to_image(
    smiles: str,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_atoms: List[int] = None,
    highlight_bonds: List[int] = None
) -> bytes:
    """
    Converts a SMILES string to a molecular structure image.
    
    Args:
        smiles: SMILES string representation of a molecule
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_atoms: List of atom indices to highlight (default: None)
        highlight_bonds: List of bond indices to highlight (default: None)
        
    Returns:
        Image data in the specified format
        
    Raises:
        MolecularProcessingException: If SMILES conversion or image generation fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_image(mol, width, height, format, highlight_atoms, highlight_bonds)


def image_to_base64(image_data: bytes, format: str = 'png') -> str:
    """
    Converts image bytes to a base64-encoded string.
    
    Args:
        image_data: Image data as bytes
        format: Image format (default: 'png')
        
    Returns:
        Base64-encoded image string with MIME type prefix
    """
    base64_data = base64.b64encode(image_data).decode('utf-8')
    
    # Determine MIME type based on format
    mime_type = 'image/png'
    if format.lower() == 'svg':
        mime_type = 'image/svg+xml'
    elif format.lower() == 'jpg' or format.lower() == 'jpeg':
        mime_type = 'image/jpeg'
    
    return f"data:{mime_type};base64,{base64_data}"


def smiles_to_base64_image(
    smiles: str,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_atoms: List[int] = None,
    highlight_bonds: List[int] = None
) -> str:
    """
    Converts a SMILES string to a base64-encoded image.
    
    Args:
        smiles: SMILES string representation of a molecule
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_atoms: List of atom indices to highlight (default: None)
        highlight_bonds: List of bond indices to highlight (default: None)
        
    Returns:
        Base64-encoded image string with MIME type prefix
        
    Raises:
        MolecularProcessingException: If SMILES conversion or image generation fails
    """
    image_data = smiles_to_image(smiles, width, height, format, highlight_atoms, highlight_bonds)
    return image_to_base64(image_data, format)


def mol_to_fingerprint(
    mol: Chem.rdchem.Mol,
    fingerprint_type: str = 'morgan',
    params: Dict[str, Any] = None
) -> Any:
    """
    Generates a molecular fingerprint from an RDKit molecule.
    
    Args:
        mol: RDKit molecule object
        fingerprint_type: Type of fingerprint to generate, one of: 
                          'morgan', 'maccs', 'topological', 'pattern' (default: 'morgan')
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        Molecular fingerprint object
        
    Raises:
        MolecularProcessingException: If fingerprint generation fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot generate fingerprint from None molecule",
            details={"mol": "None"}
        )
    
    if fingerprint_type not in FINGERPRINT_TYPES:
        raise MolecularProcessingException(
            f"Unsupported fingerprint type: {fingerprint_type}",
            details={"supported_types": list(FINGERPRINT_TYPES.keys())}
        )
    
    try:
        fp_function = FINGERPRINT_TYPES[fingerprint_type]
        
        # Set default parameters if none provided
        if params is None:
            params = {}
            if fingerprint_type == 'morgan':
                params = {'radius': 2, 'nBits': 2048}
            elif fingerprint_type == 'topological' or fingerprint_type == 'pattern':
                params = {'nBits': 2048}
        
        # Generate fingerprint with appropriate parameters
        if fingerprint_type == 'morgan':
            fp = fp_function(mol, radius=params.get('radius', 2), nBits=params.get('nBits', 2048))
        elif fingerprint_type == 'maccs':
            fp = fp_function(mol)
        else:  # topological or pattern
            fp = fp_function(mol, nBits=params.get('nBits', 2048))
        
        return fp
    except Exception as e:
        raise MolecularProcessingException(
            f"Error generating molecular fingerprint: {str(e)}",
            details={"error": str(e), "fingerprint_type": fingerprint_type}
        )


def smiles_to_fingerprint(
    smiles: str,
    fingerprint_type: str = 'morgan',
    params: Dict[str, Any] = None
) -> Any:
    """
    Generates a molecular fingerprint from a SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        fingerprint_type: Type of fingerprint to generate, one of: 
                          'morgan', 'maccs', 'topological', 'pattern' (default: 'morgan')
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        Molecular fingerprint object
        
    Raises:
        MolecularProcessingException: If SMILES conversion or fingerprint generation fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_fingerprint(mol, fingerprint_type, params)


def mol_to_inchi(mol: Chem.rdchem.Mol) -> str:
    """
    Converts an RDKit molecule to an InChI string.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        InChI string representation of the molecule
        
    Raises:
        MolecularProcessingException: If InChI conversion fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot convert None to InChI",
            details={"mol": "None"}
        )
    
    try:
        inchi = Chem.MolToInchi(mol)
        if not inchi:
            raise MolecularProcessingException(
                "Failed to convert molecule to InChI",
                details={"mol": Chem.MolToSmiles(mol) if mol else "None"}
            )
        return inchi
    except Exception as e:
        raise MolecularProcessingException(
            f"Error converting molecule to InChI: {str(e)}",
            details={"error": str(e)}
        )


def smiles_to_inchi(smiles: str) -> str:
    """
    Converts a SMILES string to an InChI string.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        InChI string representation of the molecule
        
    Raises:
        MolecularProcessingException: If SMILES conversion or InChI conversion fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_inchi(mol)


def mol_to_inchikey(mol: Chem.rdchem.Mol) -> str:
    """
    Converts an RDKit molecule to an InChIKey.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        InChIKey representation of the molecule
        
    Raises:
        MolecularProcessingException: If InChIKey conversion fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot convert None to InChIKey",
            details={"mol": "None"}
        )
    
    try:
        inchi = mol_to_inchi(mol)
        inchikey = Chem.InchiToInchiKey(inchi)
        if not inchikey:
            raise MolecularProcessingException(
                "Failed to convert InChI to InChIKey",
                details={"inchi": inchi}
            )
        return inchikey
    except Exception as e:
        raise MolecularProcessingException(
            f"Error converting molecule to InChIKey: {str(e)}",
            details={"error": str(e)}
        )


def smiles_to_inchikey(smiles: str) -> str:
    """
    Converts a SMILES string to an InChIKey.
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        InChIKey representation of the molecule
        
    Raises:
        MolecularProcessingException: If SMILES conversion or InChIKey conversion fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_inchikey(mol)


def mol_to_molblock(mol: Chem.rdchem.Mol) -> str:
    """
    Converts an RDKit molecule to a molblock string (MDL MOL format).
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Molblock string representation of the molecule
        
    Raises:
        MolecularProcessingException: If molblock conversion fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot convert None to molblock",
            details={"mol": "None"}
        )
    
    try:
        molblock = Chem.MolToMolBlock(mol)
        if not molblock:
            raise MolecularProcessingException(
                "Failed to convert molecule to molblock",
                details={"mol": Chem.MolToSmiles(mol) if mol else "None"}
            )
        return molblock
    except Exception as e:
        raise MolecularProcessingException(
            f"Error converting molecule to molblock: {str(e)}",
            details={"error": str(e)}
        )


def smiles_to_molblock(smiles: str) -> str:
    """
    Converts a SMILES string to a molblock string (MDL MOL format).
    
    Args:
        smiles: SMILES string representation of a molecule
        
    Returns:
        Molblock string representation of the molecule
        
    Raises:
        MolecularProcessingException: If SMILES conversion or molblock conversion fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return mol_to_molblock(mol)