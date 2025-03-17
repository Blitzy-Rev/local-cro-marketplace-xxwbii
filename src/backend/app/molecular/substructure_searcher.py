"""
Substructure Searcher Module

This module provides functionality for searching molecules for substructures, finding
substructure matches, and visualizing substructure matches within molecules. It supports
both individual molecules and collections of molecules, and offers both functional
and object-oriented interfaces for substructure searching.

Key functionalities:
- Searching for the presence of substructures using SMARTS patterns
- Finding all matches of a substructure within a molecule
- Filtering collections of molecules based on substructure presence
- Generating visualizations with highlighted substructure matches
- Finding common substructures among a set of molecules

All functions include comprehensive error handling to provide clear feedback on
any issues during substructure operations.
"""

from typing import List, Dict, Optional, Union, Any, Tuple, Set
import logging

# RDKit version 2023.03+
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, MCS

from ..exceptions import MolecularProcessingException
from .molecule_converter import (
    smiles_to_mol,
    smiles_to_mol_with_exception,
    mol_to_image,
    image_to_base64
)

# Set up module logger
logger = logging.getLogger(__name__)

# Default color for substructure highlighting
DEFAULT_HIGHLIGHT_COLOR = (0.7, 1.0, 0.7)  # Light green


def search_substructure(mol: Chem.rdchem.Mol, smarts_pattern: str) -> bool:
    """
    Searches for a substructure in a molecule using a SMARTS pattern.
    
    Args:
        mol: RDKit molecule object
        smarts_pattern: SMARTS pattern representing the substructure to search for
        
    Returns:
        True if the molecule contains the substructure, False otherwise
        
    Raises:
        MolecularProcessingException: If the molecule is None or SMARTS pattern is invalid
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot search for substructure in None molecule"
        )
    
    if not smarts_pattern:
        raise MolecularProcessingException(
            "SMARTS pattern cannot be empty"
        )
    
    query = Chem.MolFromSmarts(smarts_pattern)
    if query is None:
        raise MolecularProcessingException(
            f"Invalid SMARTS pattern: {smarts_pattern}"
        )
    
    return mol.HasSubstructMatch(query)


def search_substructure_smiles(smiles: str, smarts_pattern: str) -> bool:
    """
    Searches for a substructure in a molecule specified by a SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        smarts_pattern: SMARTS pattern representing the substructure to search for
        
    Returns:
        True if the molecule contains the substructure, False otherwise
        
    Raises:
        MolecularProcessingException: If the SMILES string is invalid or SMARTS pattern is invalid
    """
    mol = smiles_to_mol_with_exception(smiles)
    return search_substructure(mol, smarts_pattern)


def find_substructure_matches(
    mol: Chem.rdchem.Mol, 
    smarts_pattern: str, 
    unique_matches: bool = True
) -> List[Tuple[int, ...]]:
    """
    Finds all matches of a substructure in a molecule.
    
    Args:
        mol: RDKit molecule object
        smarts_pattern: SMARTS pattern representing the substructure to search for
        unique_matches: Whether to return only unique matches (default: True)
        
    Returns:
        List of tuples containing atom indices for each match
        
    Raises:
        MolecularProcessingException: If the molecule is None or SMARTS pattern is invalid
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot find substructure matches in None molecule"
        )
    
    if not smarts_pattern:
        raise MolecularProcessingException(
            "SMARTS pattern cannot be empty"
        )
    
    query = Chem.MolFromSmarts(smarts_pattern)
    if query is None:
        raise MolecularProcessingException(
            f"Invalid SMARTS pattern: {smarts_pattern}"
        )
    
    return mol.GetSubstructMatches(query, uniquify=unique_matches)


def filter_molecules_with_substructure(
    molecules: List[Chem.rdchem.Mol], 
    smarts_pattern: str
) -> List[Chem.rdchem.Mol]:
    """
    Filters a list of molecules to include only those with a specific substructure.
    
    Args:
        molecules: List of RDKit molecule objects
        smarts_pattern: SMARTS pattern representing the substructure to filter by
        
    Returns:
        Filtered list of molecules containing the substructure
        
    Raises:
        MolecularProcessingException: If the SMARTS pattern is invalid
    """
    if not molecules:
        return []
    
    if not smarts_pattern:
        raise MolecularProcessingException(
            "SMARTS pattern cannot be empty"
        )
    
    query = Chem.MolFromSmarts(smarts_pattern)
    if query is None:
        raise MolecularProcessingException(
            f"Invalid SMARTS pattern: {smarts_pattern}"
        )
    
    return [mol for mol in molecules if mol is not None and mol.HasSubstructMatch(query)]


def filter_smiles_with_substructure(
    smiles_list: List[str], 
    smarts_pattern: str
) -> List[str]:
    """
    Filters a list of SMILES strings to include only those with a specific substructure.
    
    Args:
        smiles_list: List of SMILES strings
        smarts_pattern: SMARTS pattern representing the substructure to filter by
        
    Returns:
        Filtered list of SMILES strings containing the substructure
        
    Raises:
        MolecularProcessingException: If the SMARTS pattern is invalid
    """
    if not smiles_list:
        return []
    
    # Convert SMILES strings to molecules, ignoring invalid ones
    molecules = []
    valid_smiles = []
    
    for i, smiles in enumerate(smiles_list):
        try:
            mol = smiles_to_mol(smiles)
            if mol is not None:
                molecules.append(mol)
                valid_smiles.append(smiles)
        except Exception as e:
            logger.warning(f"Skipping invalid SMILES at index {i}: {smiles}. Error: {str(e)}")
    
    # Filter molecules
    filtered_molecules = filter_molecules_with_substructure(molecules, smarts_pattern)
    
    # Map filtered molecules back to their SMILES strings
    indices = []
    for mol in filtered_molecules:
        for i, orig_mol in enumerate(molecules):
            if mol is orig_mol:
                indices.append(i)
                break
    
    return [valid_smiles[i] for i in indices]


def generate_substructure_highlight_image(
    mol: Chem.rdchem.Mol,
    smarts_pattern: str,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_color: tuple = DEFAULT_HIGHLIGHT_COLOR
) -> bytes:
    """
    Generates an image with highlighted substructure matches.
    
    Args:
        mol: RDKit molecule object
        smarts_pattern: SMARTS pattern representing the substructure to highlight
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_color: RGB tuple for highlight color (default: light green)
        
    Returns:
        Image data with highlighted substructure
        
    Raises:
        MolecularProcessingException: If the molecule is None, SMARTS pattern is invalid,
                                     or image generation fails
    """
    # Find all substructure matches
    matches = find_substructure_matches(mol, smarts_pattern)
    
    # Extract all atoms involved in matches
    highlight_atoms = set()
    for match in matches:
        highlight_atoms.update(match)
    
    # Generate image with highlighted atoms
    return mol_to_image(
        mol, 
        width=width, 
        height=height, 
        format=format, 
        highlight_atoms=list(highlight_atoms)
    )


def generate_substructure_highlight_base64(
    mol: Chem.rdchem.Mol,
    smarts_pattern: str,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_color: tuple = DEFAULT_HIGHLIGHT_COLOR
) -> str:
    """
    Generates a base64-encoded image with highlighted substructure matches.
    
    Args:
        mol: RDKit molecule object
        smarts_pattern: SMARTS pattern representing the substructure to highlight
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_color: RGB tuple for highlight color (default: light green)
        
    Returns:
        Base64-encoded image string with highlighted substructure
        
    Raises:
        MolecularProcessingException: If the molecule is None, SMARTS pattern is invalid,
                                     or image generation fails
    """
    image_data = generate_substructure_highlight_image(
        mol, 
        smarts_pattern, 
        width, 
        height, 
        format, 
        highlight_color
    )
    
    return image_to_base64(image_data, format)


def generate_smiles_highlight_base64(
    smiles: str,
    smarts_pattern: str,
    width: int = 300,
    height: int = 200,
    format: str = 'png',
    highlight_color: tuple = DEFAULT_HIGHLIGHT_COLOR
) -> str:
    """
    Generates a base64-encoded image with highlighted substructure matches for a SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        smarts_pattern: SMARTS pattern representing the substructure to highlight
        width: Image width in pixels (default: 300)
        height: Image height in pixels (default: 200)
        format: Image format, 'png' or 'svg' (default: 'png')
        highlight_color: RGB tuple for highlight color (default: light green)
        
    Returns:
        Base64-encoded image string with highlighted substructure
        
    Raises:
        MolecularProcessingException: If the SMILES string is invalid, SMARTS pattern is invalid,
                                     or image generation fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    
    return generate_substructure_highlight_base64(
        mol, 
        smarts_pattern, 
        width, 
        height, 
        format, 
        highlight_color
    )


def find_common_substructures(
    molecules: List[Chem.rdchem.Mol],
    threshold: float = 0.7,
    maximize_bonds: bool = True
) -> List[Dict[str, Any]]:
    """
    Identifies common substructures in a set of molecules.
    
    Args:
        molecules: List of RDKit molecule objects
        threshold: Minimum fraction of molecules that must contain the substructure (default: 0.7)
        maximize_bonds: Whether to maximize the number of bonds (default: True)
        
    Returns:
        List of dictionaries containing SMARTS patterns and match information
        
    Raises:
        MolecularProcessingException: If fewer than 2 valid molecules are provided
    """
    # Filter out None molecules
    valid_molecules = [mol for mol in molecules if mol is not None]
    
    if len(valid_molecules) < 2:
        raise MolecularProcessingException(
            "At least 2 valid molecules are required to find common substructures"
        )
    
    try:
        # Create MCS finder
        mcs_finder = MCS.FindMCS(
            valid_molecules,
            threshold=threshold,
            maximizeBonds=maximize_bonds,
            completeRingsOnly=True,
            matchValences=True
        )
        
        # Get SMARTS pattern for MCS
        mcs_smarts = mcs_finder.smartsString
        
        if not mcs_smarts:
            return []
        
        # Find matches in each molecule
        results = []
        for i, mol in enumerate(valid_molecules):
            matches = find_substructure_matches(mol, mcs_smarts)
            if matches:
                results.append({
                    "molecule_index": i,
                    "smarts": mcs_smarts,
                    "num_atoms": mcs_finder.numAtoms,
                    "num_bonds": mcs_finder.numBonds,
                    "matches": matches,
                    "match_count": len(matches)
                })
        
        return results
    except Exception as e:
        raise MolecularProcessingException(
            f"Error finding common substructures: {str(e)}",
            details={"error": str(e)}
        )


def find_common_substructures_from_smiles(
    smiles_list: List[str],
    threshold: float = 0.7,
    maximize_bonds: bool = True
) -> List[Dict[str, Any]]:
    """
    Identifies common substructures in a set of molecules specified by SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        threshold: Minimum fraction of molecules that must contain the substructure (default: 0.7)
        maximize_bonds: Whether to maximize the number of bonds (default: True)
        
    Returns:
        List of dictionaries containing SMARTS patterns and match information
        
    Raises:
        MolecularProcessingException: If fewer than 2 valid molecules are provided
    """
    # Convert SMILES strings to molecules, ignoring invalid ones
    molecules = []
    
    for smiles in smiles_list:
        try:
            mol = smiles_to_mol(smiles)
            if mol is not None:
                molecules.append(mol)
        except Exception as e:
            logger.warning(f"Skipping invalid SMILES: {smiles}. Error: {str(e)}")
    
    return find_common_substructures(molecules, threshold, maximize_bonds)


class SubstructureSearcher:
    """
    Class for performing substructure searches and highlighting.
    
    This class provides an object-oriented interface for working with molecular
    substructures, including searching for substructures, finding all matches,
    filtering collections of molecules, and generating visualizations.
    """
    
    def __init__(self, smarts_pattern: str = None, highlight_color: tuple = DEFAULT_HIGHLIGHT_COLOR):
        """
        Initializes the SubstructureSearcher with a SMARTS pattern.
        
        Args:
            smarts_pattern: SMARTS pattern representing the substructure to search for (default: None)
            highlight_color: RGB tuple for highlight color (default: light green)
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is invalid
        """
        self._smarts_pattern = None
        self._highlight_color = highlight_color
        
        # Initialize logger
        self.logger = logging.getLogger(__name__)
        
        # Validate and set SMARTS pattern if provided
        if smarts_pattern:
            self.set_smarts_pattern(smarts_pattern)
    
    def set_smarts_pattern(self, smarts_pattern: str) -> None:
        """
        Sets the SMARTS pattern to use for substructure searches.
        
        Args:
            smarts_pattern: SMARTS pattern representing the substructure to search for
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is invalid
        """
        if not smarts_pattern:
            raise MolecularProcessingException(
                "SMARTS pattern cannot be empty"
            )
        
        query = Chem.MolFromSmarts(smarts_pattern)
        if query is None:
            raise MolecularProcessingException(
                f"Invalid SMARTS pattern: {smarts_pattern}"
            )
        
        self._smarts_pattern = smarts_pattern
    
    def set_highlight_color(self, highlight_color: tuple) -> None:
        """
        Sets the color to use for highlighting substructure matches.
        
        Args:
            highlight_color: RGB tuple for highlight color (values from 0 to 1)
            
        Raises:
            MolecularProcessingException: If the highlight color is invalid
        """
        if not isinstance(highlight_color, tuple) or len(highlight_color) != 3:
            raise MolecularProcessingException(
                "Highlight color must be a tuple of 3 values"
            )
        
        for value in highlight_color:
            if not isinstance(value, (int, float)) or value < 0 or value > 1:
                raise MolecularProcessingException(
                    "Highlight color values must be between 0 and 1"
                )
        
        self._highlight_color = highlight_color
    
    def search(self, molecule: Union[str, Chem.rdchem.Mol]) -> bool:
        """
        Searches for the configured substructure in a molecule.
        
        Args:
            molecule: SMILES string or RDKit molecule object
            
        Returns:
            True if the molecule contains the substructure, False otherwise
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is not set,
                                         molecule is invalid, or search fails
        """
        if self._smarts_pattern is None:
            raise MolecularProcessingException(
                "SMARTS pattern has not been set"
            )
        
        # Convert SMILES to molecule if needed
        mol = molecule
        if isinstance(molecule, str):
            mol = smiles_to_mol_with_exception(molecule)
        
        return search_substructure(mol, self._smarts_pattern)
    
    def find_matches(
        self, 
        molecule: Union[str, Chem.rdchem.Mol], 
        unique_matches: bool = True
    ) -> List[Tuple[int, ...]]:
        """
        Finds all matches of the configured substructure in a molecule.
        
        Args:
            molecule: SMILES string or RDKit molecule object
            unique_matches: Whether to return only unique matches (default: True)
            
        Returns:
            List of tuples containing atom indices for each match
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is not set,
                                         molecule is invalid, or search fails
        """
        if self._smarts_pattern is None:
            raise MolecularProcessingException(
                "SMARTS pattern has not been set"
            )
        
        # Convert SMILES to molecule if needed
        mol = molecule
        if isinstance(molecule, str):
            mol = smiles_to_mol_with_exception(molecule)
        
        return find_substructure_matches(mol, self._smarts_pattern, unique_matches)
    
    def filter_molecules(
        self, 
        molecules: List[Union[str, Chem.rdchem.Mol]]
    ) -> List[Union[str, Chem.rdchem.Mol]]:
        """
        Filters a list of molecules to include only those with the configured substructure.
        
        Args:
            molecules: List of SMILES strings or RDKit molecule objects
            
        Returns:
            Filtered list of molecules containing the substructure
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is not set or filtering fails
        """
        if self._smarts_pattern is None:
            raise MolecularProcessingException(
                "SMARTS pattern has not been set"
            )
        
        if not molecules:
            return []
        
        # Determine input type (SMILES or RDKit molecules)
        if all(isinstance(m, str) for m in molecules):
            return filter_smiles_with_substructure(molecules, self._smarts_pattern)
        elif all(isinstance(m, Chem.rdchem.Mol) or m is None for m in molecules):
            return filter_molecules_with_substructure(molecules, self._smarts_pattern)
        else:
            raise MolecularProcessingException(
                "All molecules must be of the same type (either all SMILES strings or all RDKit molecules)"
            )
    
    def generate_highlight_image(
        self,
        molecule: Union[str, Chem.rdchem.Mol],
        width: int = 300,
        height: int = 200,
        format: str = 'png'
    ) -> bytes:
        """
        Generates an image with highlighted substructure matches.
        
        Args:
            molecule: SMILES string or RDKit molecule object
            width: Image width in pixels (default: 300)
            height: Image height in pixels (default: 200)
            format: Image format, 'png' or 'svg' (default: 'png')
            
        Returns:
            Image data with highlighted substructure
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is not set,
                                         molecule is invalid, or image generation fails
        """
        if self._smarts_pattern is None:
            raise MolecularProcessingException(
                "SMARTS pattern has not been set"
            )
        
        # Convert SMILES to molecule if needed
        mol = molecule
        if isinstance(molecule, str):
            mol = smiles_to_mol_with_exception(molecule)
        
        return generate_substructure_highlight_image(
            mol, 
            self._smarts_pattern, 
            width, 
            height, 
            format, 
            self._highlight_color
        )
    
    def generate_highlight_base64(
        self,
        molecule: Union[str, Chem.rdchem.Mol],
        width: int = 300,
        height: int = 200,
        format: str = 'png'
    ) -> str:
        """
        Generates a base64-encoded image with highlighted substructure matches.
        
        Args:
            molecule: SMILES string or RDKit molecule object
            width: Image width in pixels (default: 300)
            height: Image height in pixels (default: 200)
            format: Image format, 'png' or 'svg' (default: 'png')
            
        Returns:
            Base64-encoded image string with highlighted substructure
            
        Raises:
            MolecularProcessingException: If the SMARTS pattern is not set,
                                         molecule is invalid, or image generation fails
        """
        if self._smarts_pattern is None:
            raise MolecularProcessingException(
                "SMARTS pattern has not been set"
            )
        
        # Convert SMILES to molecule if needed
        mol = molecule
        if isinstance(molecule, str):
            mol = smiles_to_mol_with_exception(molecule)
        
        return generate_substructure_highlight_base64(
            mol, 
            self._smarts_pattern, 
            width, 
            height, 
            format, 
            self._highlight_color
        )