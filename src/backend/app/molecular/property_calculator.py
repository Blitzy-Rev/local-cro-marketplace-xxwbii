"""
Property Calculator Module

This module provides functionality for calculating molecular properties in the
Molecular Data Management and CRO Integration Platform. It includes functions for
computing various physicochemical properties of molecules, including molecular weight,
LogP, hydrogen bond donors/acceptors, and other drug-like properties.

The module supports both individual molecule processing and batch calculations with
efficient parallel execution for performance optimization.
"""

from typing import List, Dict, Optional, Union, Any, Callable, Tuple
import concurrent.futures
import logging

# RDKit version 2023.03+
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED

from ..exceptions import MolecularProcessingException
from .molecule_converter import smiles_to_mol, smiles_to_mol_with_exception
from .validator import PROPERTY_RANGES, validate_property_value

# Set up module logger
logger = logging.getLogger(__name__)

# Default number of workers for parallel processing
DEFAULT_NUM_WORKERS = 4

# Dictionary mapping property names to their respective RDKit calculation functions
PROPERTY_CALCULATORS = {
    'molecular_weight': Descriptors.MolWt,
    'logp': Descriptors.MolLogP,
    'h_bond_donors': Descriptors.NumHDonors,
    'h_bond_acceptors': Descriptors.NumHAcceptors,
    'rotatable_bonds': Descriptors.NumRotatableBonds,
    'polar_surface_area': Descriptors.TPSA,
    'heavy_atom_count': Descriptors.HeavyAtomCount,
    'ring_count': Descriptors.RingCount,
    'aromatic_rings': Lipinski.NumAromaticRings,
    'qed': QED.qed,
    'fraction_csp3': Descriptors.FractionCSP3,
    'num_stereocenters': Chem.FindMolChiralCenters
}


def calculate_property(mol: Chem.rdchem.Mol, property_name: str) -> float:
    """
    Calculates a single molecular property for a given molecule.
    
    Args:
        mol: RDKit molecule object
        property_name: Name of the property to calculate
        
    Returns:
        float: Calculated property value
        
    Raises:
        MolecularProcessingException: If property calculation fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot calculate property for None molecule",
            details={"property_name": property_name}
        )
    
    if property_name not in PROPERTY_CALCULATORS:
        raise MolecularProcessingException(
            f"Unknown property: {property_name}",
            details={"available_properties": list(PROPERTY_CALCULATORS.keys())}
        )
    
    calculator = PROPERTY_CALCULATORS[property_name]
    
    try:
        # Special handling for num_stereocenters which returns a tuple list
        if property_name == 'num_stereocenters':
            stereocenters = calculator(mol)
            return float(len(stereocenters))
        
        # Calculate the property
        value = calculator(mol)
        return float(value)
    except Exception as e:
        raise MolecularProcessingException(
            f"Error calculating property {property_name}: {str(e)}",
            details={"property_name": property_name, "error": str(e)}
        )


def calculate_properties(mol: Chem.rdchem.Mol, property_names: List[str] = None) -> Dict[str, float]:
    """
    Calculates multiple properties for a given molecule.
    
    Args:
        mol: RDKit molecule object
        property_names: List of property names to calculate (default: all available properties)
        
    Returns:
        Dict[str, float]: Dictionary of property names and their calculated values
        
    Raises:
        MolecularProcessingException: If molecule is None or property calculation fails
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot calculate properties for None molecule"
        )
    
    # If no property names specified, calculate all available properties
    if property_names is None:
        property_names = list(PROPERTY_CALCULATORS.keys())
    
    results = {}
    
    for prop_name in property_names:
        try:
            value = calculate_property(mol, prop_name)
            results[prop_name] = value
        except MolecularProcessingException as e:
            logger.warning(f"Failed to calculate property {prop_name}: {e.message}")
            # Skip this property but continue with others
            continue
    
    return results


def calculate_properties_from_smiles(smiles: str, property_names: List[str] = None) -> Dict[str, float]:
    """
    Calculates properties for a molecule specified by SMILES string.
    
    Args:
        smiles: SMILES string representation of a molecule
        property_names: List of property names to calculate (default: all available properties)
        
    Returns:
        Dict[str, float]: Dictionary of property names and their calculated values
        
    Raises:
        MolecularProcessingException: If SMILES conversion fails or property calculation fails
    """
    mol = smiles_to_mol_with_exception(smiles)
    return calculate_properties(mol, property_names)


def batch_calculate_properties(
    molecules: List[Chem.rdchem.Mol],
    property_names: List[str] = None,
    num_workers: int = None
) -> List[Dict[str, float]]:
    """
    Calculates properties for a batch of molecules in parallel.
    
    Args:
        molecules: List of RDKit molecule objects
        property_names: List of property names to calculate (default: all available properties)
        num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
        
    Returns:
        List[Dict[str, float]]: List of dictionaries containing calculated properties for each molecule
    """
    if num_workers is None:
        num_workers = DEFAULT_NUM_WORKERS
    
    if not molecules:
        return []
    
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit property calculation tasks
        future_to_index = {}
        for i, mol in enumerate(molecules):
            if mol is not None:
                future = executor.submit(calculate_properties, mol, property_names)
                future_to_index[future] = i
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_index):
            idx = future_to_index[future]
            try:
                prop_dict = future.result()
                results.append((idx, prop_dict))
            except Exception as e:
                logger.error(f"Error calculating properties for molecule at index {idx}: {str(e)}")
                # Add empty result for failed calculation
                results.append((idx, {}))
    
    # Sort results by original index and return just the property dictionaries
    results.sort(key=lambda x: x[0])
    return [r[1] for r in results]


def batch_calculate_properties_from_smiles(
    smiles_list: List[str],
    property_names: List[str] = None,
    num_workers: int = None
) -> List[Dict[str, float]]:
    """
    Calculates properties for a batch of molecules specified by SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        property_names: List of property names to calculate (default: all available properties)
        num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
        
    Returns:
        List[Dict[str, float]]: List of dictionaries containing calculated properties for each molecule
    """
    if num_workers is None:
        num_workers = DEFAULT_NUM_WORKERS
    
    if not smiles_list:
        return []
    
    results = []
    
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit property calculation tasks
        future_to_index = {}
        for i, smiles in enumerate(smiles_list):
            if smiles:
                future = executor.submit(calculate_properties_from_smiles, smiles, property_names)
                future_to_index[future] = i
        
        # Collect results as they complete
        for future in concurrent.futures.as_completed(future_to_index):
            idx = future_to_index[future]
            try:
                prop_dict = future.result()
                results.append((idx, prop_dict))
            except Exception as e:
                logger.error(f"Error calculating properties for SMILES at index {idx}: {str(e)}")
                # Add empty result for failed calculation
                results.append((idx, {}))
    
    # Sort results by original index and return just the property dictionaries
    results.sort(key=lambda x: x[0])
    return [r[1] for r in results]


def check_lipinski_rule_of_five(mol: Chem.rdchem.Mol) -> Dict[str, Any]:
    """
    Checks if a molecule follows Lipinski's Rule of Five for drug-likeness.
    
    Lipinski's rules:
    - Molecular weight < 500 Da
    - LogP < 5
    - H-bond donors <= 5
    - H-bond acceptors <= 10
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dict[str, Any]: Dictionary with Lipinski properties and rule compliance
        
    Raises:
        MolecularProcessingException: If molecule is None
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot check Lipinski rules for None molecule"
        )
    
    # Calculate properties
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    h_donors = Descriptors.NumHDonors(mol)
    h_acceptors = Descriptors.NumHAcceptors(mol)
    
    # Check rules
    mw_pass = mw <= 500
    logp_pass = logp <= 5
    donors_pass = h_donors <= 5
    acceptors_pass = h_acceptors <= 10
    
    # Count violations
    violations = sum(not x for x in [mw_pass, logp_pass, donors_pass, acceptors_pass])
    
    # Generally, up to 1 violation is still considered drug-like
    is_drug_like = violations <= 1
    
    return {
        'molecular_weight': mw,
        'logp': logp,
        'h_bond_donors': h_donors,
        'h_bond_acceptors': h_acceptors,
        'mw_pass': mw_pass,
        'logp_pass': logp_pass,
        'donors_pass': donors_pass,
        'acceptors_pass': acceptors_pass,
        'violations': violations,
        'is_drug_like': is_drug_like
    }


def check_veber_rules(mol: Chem.rdchem.Mol) -> Dict[str, Any]:
    """
    Checks if a molecule follows Veber's rules for oral bioavailability.
    
    Veber's rules:
    - Rotatable bonds <= 10
    - Polar surface area <= 140 Å²
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        Dict[str, Any]: Dictionary with Veber properties and rule compliance
        
    Raises:
        MolecularProcessingException: If molecule is None
    """
    if mol is None:
        raise MolecularProcessingException(
            "Cannot check Veber rules for None molecule"
        )
    
    # Calculate properties
    rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    polar_surface_area = Descriptors.TPSA(mol)
    
    # Check rules
    rotatable_pass = rotatable_bonds <= 10
    psa_pass = polar_surface_area <= 140
    
    # Both rules must be satisfied
    passes_veber_rules = rotatable_pass and psa_pass
    
    return {
        'rotatable_bonds': rotatable_bonds,
        'polar_surface_area': polar_surface_area,
        'rotatable_pass': rotatable_pass,
        'psa_pass': psa_pass,
        'passes_veber_rules': passes_veber_rules
    }


def register_property_calculator(property_name: str, calculator_function: Callable) -> None:
    """
    Registers a new property calculator function.
    
    Args:
        property_name: Name of the property
        calculator_function: Function that calculates the property from an RDKit molecule
        
    Returns:
        None
    """
    PROPERTY_CALCULATORS[property_name] = calculator_function
    logger.info(f"Registered new property calculator: {property_name}")


def get_available_properties() -> List[str]:
    """
    Returns a list of all available property calculators.
    
    Returns:
        List[str]: List of available property names
    """
    return list(PROPERTY_CALCULATORS.keys())


class PropertyCalculator:
    """
    Class for calculating and managing molecular properties with caching capabilities.
    """
    
    def __init__(self):
        """
        Initialize the PropertyCalculator with default calculators and property ranges.
        """
        self._calculators = dict(PROPERTY_CALCULATORS)
        self._property_ranges = dict(PROPERTY_RANGES)
        self._cache = {}
    
    def register_calculator(self, property_name: str, calculator_function: Callable) -> None:
        """
        Register a new property calculator function.
        
        Args:
            property_name: Name of the property
            calculator_function: Function that calculates the property from an RDKit molecule
            
        Returns:
            None
        """
        self._calculators[property_name] = calculator_function
    
    def set_property_range(self, property_name: str, min_value: float, max_value: float) -> None:
        """
        Set the valid range for a property.
        
        Args:
            property_name: Name of the property
            min_value: Minimum valid value
            max_value: Maximum valid value
            
        Returns:
            None
        """
        self._property_ranges[property_name] = {'min': min_value, 'max': max_value}
    
    def calculate(self, mol: Chem.rdchem.Mol, property_name: str) -> float:
        """
        Calculate a single property for a molecule.
        
        Args:
            mol: RDKit molecule object
            property_name: Name of the property to calculate
            
        Returns:
            float: Calculated property value
            
        Raises:
            MolecularProcessingException: If property calculation fails
        """
        if mol is None:
            raise MolecularProcessingException(
                "Cannot calculate property for None molecule",
                details={"property_name": property_name}
            )
        
        # Generate a unique identifier for the molecule (using InChI Key)
        try:
            mol_id = Chem.MolToInchiKey(mol)
        except:
            # If InChI Key generation fails, use a fallback approach
            mol_id = f"MOL_{id(mol)}"
        
        # Check cache
        if mol_id in self._cache and property_name in self._cache[mol_id]:
            return self._cache[mol_id][property_name]
        
        if property_name not in self._calculators:
            raise MolecularProcessingException(
                f"Unknown property: {property_name}",
                details={"available_properties": list(self._calculators.keys())}
            )
        
        calculator = self._calculators[property_name]
        
        try:
            # Special handling for num_stereocenters which returns a tuple list
            if property_name == 'num_stereocenters':
                stereocenters = calculator(mol)
                value = float(len(stereocenters))
            else:
                # Calculate the property
                value = calculator(mol)
                value = float(value)
            
            # Cache the result
            if mol_id not in self._cache:
                self._cache[mol_id] = {}
            self._cache[mol_id][property_name] = value
            
            return value
        except Exception as e:
            raise MolecularProcessingException(
                f"Error calculating property {property_name}: {str(e)}",
                details={"property_name": property_name, "error": str(e)}
            )
    
    def calculate_all(self, mol: Chem.rdchem.Mol, property_names: List[str] = None) -> Dict[str, float]:
        """
        Calculate all or specified properties for a molecule.
        
        Args:
            mol: RDKit molecule object
            property_names: List of property names to calculate (default: all available properties)
            
        Returns:
            Dict[str, float]: Dictionary of property names and their calculated values
            
        Raises:
            MolecularProcessingException: If molecule is None or property calculation fails
        """
        if mol is None:
            raise MolecularProcessingException(
                "Cannot calculate properties for None molecule"
            )
        
        # If no property names specified, calculate all available properties
        if property_names is None:
            property_names = list(self._calculators.keys())
        
        results = {}
        
        for prop_name in property_names:
            try:
                value = self.calculate(mol, prop_name)
                results[prop_name] = value
            except MolecularProcessingException as e:
                logger.warning(f"Failed to calculate property {prop_name}: {e.message}")
                # Skip this property but continue with others
                continue
        
        return results
    
    def calculate_from_smiles(self, smiles: str, property_names: List[str] = None) -> Dict[str, float]:
        """
        Calculate properties for a molecule specified by SMILES string.
        
        Args:
            smiles: SMILES string representation of a molecule
            property_names: List of property names to calculate (default: all available properties)
            
        Returns:
            Dict[str, float]: Dictionary of property names and their calculated values
            
        Raises:
            MolecularProcessingException: If SMILES conversion fails or property calculation fails
        """
        mol = smiles_to_mol_with_exception(smiles)
        return self.calculate_all(mol, property_names)
    
    def batch_calculate(
        self,
        molecules: List[Chem.rdchem.Mol],
        property_names: List[str] = None,
        num_workers: int = None
    ) -> List[Dict[str, float]]:
        """
        Calculate properties for a batch of molecules in parallel.
        
        Args:
            molecules: List of RDKit molecule objects
            property_names: List of property names to calculate (default: all available properties)
            num_workers: Number of parallel workers to use (default: DEFAULT_NUM_WORKERS)
            
        Returns:
            List[Dict[str, float]]: List of dictionaries containing calculated properties for each molecule
        """
        if num_workers is None:
            num_workers = DEFAULT_NUM_WORKERS
        
        if not molecules:
            return []
        
        results = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
            # Submit property calculation tasks
            future_to_index = {}
            for i, mol in enumerate(molecules):
                if mol is not None:
                    future = executor.submit(self.calculate_all, mol, property_names)
                    future_to_index[future] = i
            
            # Collect results as they complete
            for future in concurrent.futures.as_completed(future_to_index):
                idx = future_to_index[future]
                try:
                    prop_dict = future.result()
                    results.append((idx, prop_dict))
                except Exception as e:
                    logger.error(f"Error calculating properties for molecule at index {idx}: {str(e)}")
                    # Add empty result for failed calculation
                    results.append((idx, {}))
        
        # Sort results by original index and return just the property dictionaries
        results.sort(key=lambda x: x[0])
        return [r[1] for r in results]
    
    def validate_property(self, property_name: str, value: float) -> bool:
        """
        Validate if a property value is within the expected range.
        
        Args:
            property_name: Name of the property to validate
            value: Value to validate
            
        Returns:
            bool: True if value is within range, False otherwise
        """
        return validate_property_value(property_name, value)
    
    def clear_cache(self) -> None:
        """
        Clear the property calculation cache.
        
        Returns:
            None
        """
        self._cache = {}