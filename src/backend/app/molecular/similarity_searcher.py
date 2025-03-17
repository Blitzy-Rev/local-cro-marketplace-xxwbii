"""
Molecular Similarity Searcher Module

This module provides functionality for molecular similarity searching in the
Molecular Data Management and CRO Integration Platform. It includes functions for:

- Calculating similarity between molecules using various fingerprint types
- Finding similar molecules in a dataset based on a query molecule
- Clustering molecules based on structural similarity
- Selecting diverse subsets of molecules for analysis

The module supports different fingerprint types including Morgan fingerprints,
MACCS keys, topological fingerprints, and pattern fingerprints.
"""

from typing import List, Dict, Optional, Union, Any, Tuple, Set
import numpy as np
import logging

# RDKit imports (version 2023.03+)
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import rdMolDescriptors

# Internal imports
from ..exceptions import MolecularProcessingException
from .molecule_converter import (
    smiles_to_mol,
    smiles_to_mol_with_exception,
    mol_to_fingerprint,
    FINGERPRINT_TYPES
)

# Set up module logger
logger = logging.getLogger(__name__)

# Default values for similarity searching
DEFAULT_SIMILARITY_THRESHOLD = 0.7
DEFAULT_FINGERPRINT_TYPE = "morgan"
DEFAULT_MORGAN_RADIUS = 2
DEFAULT_MORGAN_BITS = 2048


def calculate_similarity(
    mol1: Chem.rdchem.Mol,
    mol2: Chem.rdchem.Mol,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> float:
    """
    Calculates the Tanimoto similarity between two molecules.
    
    Args:
        mol1: First RDKit molecule
        mol2: Second RDKit molecule
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        Similarity score between 0 and 1
        
    Raises:
        MolecularProcessingException: If either molecule is None or fingerprint generation fails
    """
    if mol1 is None or mol2 is None:
        raise MolecularProcessingException(
            "Cannot calculate similarity between None molecules",
            details={"mol1": "None" if mol1 is None else "Valid", 
                     "mol2": "None" if mol2 is None else "Valid"}
        )
    
    # Generate fingerprints for both molecules
    fp1 = mol_to_fingerprint(mol1, fingerprint_type, params)
    fp2 = mol_to_fingerprint(mol2, fingerprint_type, params)
    
    # Calculate Tanimoto similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity


def calculate_similarity_from_smiles(
    smiles1: str,
    smiles2: str,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> float:
    """
    Calculates the Tanimoto similarity between two molecules specified by SMILES strings.
    
    Args:
        smiles1: SMILES string of first molecule
        smiles2: SMILES string of second molecule
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        Similarity score between 0 and 1
        
    Raises:
        MolecularProcessingException: If SMILES conversion fails
    """
    # Convert SMILES to molecules
    mol1 = smiles_to_mol_with_exception(smiles1)
    mol2 = smiles_to_mol_with_exception(smiles2)
    
    return calculate_similarity(mol1, mol2, fingerprint_type, params)


def calculate_similarity_matrix(
    molecules: List[Chem.rdchem.Mol],
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> np.ndarray:
    """
    Calculates a pairwise similarity matrix for a list of molecules.
    
    Args:
        molecules: List of RDKit molecules
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        2D numpy array of pairwise similarity scores
        
    Raises:
        MolecularProcessingException: If fingerprint generation fails
    """
    if not molecules:
        return np.array([])
    
    n_mols = len(molecules)
    
    # Generate fingerprints for all molecules
    fingerprints = [mol_to_fingerprint(mol, fingerprint_type, params) for mol in molecules]
    
    # Initialize similarity matrix
    similarity_matrix = np.zeros((n_mols, n_mols))
    
    # Calculate pairwise similarities
    for i in range(n_mols):
        # Diagonal element (self-similarity) is always 1.0
        similarity_matrix[i, i] = 1.0
        
        # Calculate similarity with other molecules
        for j in range(i + 1, n_mols):
            sim = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
            similarity_matrix[i, j] = sim
            similarity_matrix[j, i] = sim  # Similarity matrix is symmetric
    
    return similarity_matrix


def find_similar_molecules(
    query_mol: Chem.rdchem.Mol,
    target_mols: List[Chem.rdchem.Mol],
    threshold: float = DEFAULT_SIMILARITY_THRESHOLD,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> List[Tuple[Chem.rdchem.Mol, float]]:
    """
    Finds molecules similar to a query molecule based on a similarity threshold.
    
    Args:
        query_mol: Query molecule
        target_mols: List of target molecules to search
        threshold: Minimum similarity threshold (default: DEFAULT_SIMILARITY_THRESHOLD)
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        List of (molecule, similarity) tuples for molecules above the threshold,
        sorted by decreasing similarity
        
    Raises:
        MolecularProcessingException: If query_mol is None or fingerprint generation fails
    """
    if query_mol is None:
        raise MolecularProcessingException(
            "Query molecule cannot be None",
            details={"query_mol": "None"}
        )
    
    if not target_mols:
        return []
    
    # Generate fingerprint for query molecule
    query_fp = mol_to_fingerprint(query_mol, fingerprint_type, params)
    
    # Calculate similarity for each target molecule
    results = []
    for target_mol in target_mols:
        if target_mol is not None:
            # Generate fingerprint for target molecule
            target_fp = mol_to_fingerprint(target_mol, fingerprint_type, params)
            
            # Calculate Tanimoto similarity
            similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
            
            # Add to results if above threshold
            if similarity >= threshold:
                results.append((target_mol, similarity))
    
    # Sort results by decreasing similarity
    results.sort(key=lambda x: x[1], reverse=True)
    
    return results


def find_similar_molecules_from_smiles(
    query_smiles: str,
    target_mols: List[Union[str, Chem.rdchem.Mol]],
    threshold: float = DEFAULT_SIMILARITY_THRESHOLD,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> List[Tuple[Union[str, Chem.rdchem.Mol], float]]:
    """
    Finds molecules similar to a query SMILES string based on a similarity threshold.
    
    Args:
        query_smiles: SMILES string of query molecule
        target_mols: List of target molecules (as SMILES strings or RDKit molecules)
        threshold: Minimum similarity threshold (default: DEFAULT_SIMILARITY_THRESHOLD)
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        List of (molecule, similarity) tuples for molecules above the threshold,
        sorted by decreasing similarity
        
    Raises:
        MolecularProcessingException: If SMILES conversion fails
    """
    # Convert query SMILES to molecule
    query_mol = smiles_to_mol_with_exception(query_smiles)
    
    # Convert any SMILES strings in target_mols to molecules
    converted_targets = []
    original_to_converted = {}  # Maps original index to converted index
    
    for idx, mol in enumerate(target_mols):
        if isinstance(mol, str):
            # It's a SMILES string, convert to molecule
            rdkit_mol = smiles_to_mol(mol)
            if rdkit_mol:  # Only add valid molecules
                converted_targets.append(rdkit_mol)
                original_to_converted[idx] = len(converted_targets) - 1
        else:
            # It's already a molecule
            if mol:  # Only add non-None molecules
                converted_targets.append(mol)
                original_to_converted[idx] = len(converted_targets) - 1
    
    # Find similar molecules
    similar_mols = find_similar_molecules(
        query_mol, converted_targets, threshold, fingerprint_type, params
    )
    
    # Map back to original molecules
    results = []
    converted_to_original = {v: k for k, v in original_to_converted.items()}
    
    for i, (mol, similarity) in enumerate(similar_mols):
        # Find the original molecule corresponding to this similar molecule
        for j, conv_mol in enumerate(converted_targets):
            if mol is conv_mol:  # Identity comparison for RDKit molecules
                orig_idx = converted_to_original.get(j)
                if orig_idx is not None:
                    results.append((target_mols[orig_idx], similarity))
                break
    
    return results


def cluster_molecules(
    molecules: List[Chem.rdchem.Mol],
    algorithm: str = "butina",
    threshold: float = DEFAULT_SIMILARITY_THRESHOLD,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> List[List[int]]:
    """
    Clusters molecules based on structural similarity using a specified algorithm.
    
    Args:
        molecules: List of RDKit molecules
        algorithm: Clustering algorithm to use, either "butina" or "hierarchical" (default: "butina")
        threshold: Similarity threshold for clustering (default: DEFAULT_SIMILARITY_THRESHOLD)
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        List of clusters, where each cluster is a list of indices into the original molecules list
        
    Raises:
        MolecularProcessingException: If an unsupported algorithm is specified or clustering fails
    """
    if not molecules:
        return []
    
    # Calculate similarity matrix
    similarity_matrix = calculate_similarity_matrix(molecules, fingerprint_type, params)
    
    # Perform clustering based on specified algorithm
    if algorithm.lower() == "butina":
        # Butina clustering (based on RDKit's Butina algorithm)
        # First, convert similarity matrix to distance matrix (distance = 1 - similarity)
        distance_matrix = 1 - similarity_matrix
        
        # Flatten the distance matrix to a list of distances
        n_mols = len(molecules)
        distances = []
        for i in range(n_mols):
            for j in range(i + 1, n_mols):
                distances.append((i, j, distance_matrix[i, j]))
        
        # Sort distances by increasing distance
        distances.sort(key=lambda x: x[2])
        
        # Perform clustering
        # Initialize all molecules as singletons
        clusters = [{i} for i in range(n_mols)]
        for i, j, dist in distances:
            if dist <= 1 - threshold:  # dist <= cutoff means sim >= threshold
                # Find clusters containing molecules i and j
                cluster_i = None
                cluster_j = None
                for k, cluster in enumerate(clusters):
                    if i in cluster:
                        cluster_i = k
                    if j in cluster:
                        cluster_j = k
                    if cluster_i is not None and cluster_j is not None:
                        break
                
                if cluster_i != cluster_j:  # Only merge if they're in different clusters
                    # Merge clusters
                    clusters[cluster_i].update(clusters[cluster_j])
                    # Remove the now-empty cluster
                    clusters.pop(cluster_j)
        
        # Convert sets to lists
        return [sorted(list(cluster)) for cluster in clusters]
    
    elif algorithm.lower() == "hierarchical":
        try:
            # Hierarchical clustering using scipy
            from scipy.cluster.hierarchy import linkage, fcluster
            
            # Convert similarity matrix to distance matrix
            distance_matrix = 1 - similarity_matrix
            
            # We only need the upper triangle of the distance matrix
            n_mols = len(molecules)
            condensed_dist = []
            for i in range(n_mols):
                for j in range(i + 1, n_mols):
                    condensed_dist.append(distance_matrix[i, j])
            
            # Perform hierarchical clustering
            Z = linkage(condensed_dist, method='average')
            
            # Cut the dendrogram at the specified threshold
            clusters_raw = fcluster(Z, 1 - threshold, criterion='distance')
            
            # Organize molecules into clusters
            cluster_dict = {}
            for i, cluster_id in enumerate(clusters_raw):
                if cluster_id not in cluster_dict:
                    cluster_dict[cluster_id] = []
                cluster_dict[cluster_id].append(i)
            
            return sorted([sorted(cluster) for cluster in cluster_dict.values()])
        
        except ImportError:
            raise MolecularProcessingException(
                "Hierarchical clustering requires scipy. Please install scipy.",
                details={"algorithm": algorithm}
            )
    
    else:
        raise MolecularProcessingException(
            f"Unsupported clustering algorithm: {algorithm}",
            details={"supported_algorithms": ["butina", "hierarchical"]}
        )


def select_diverse_subset(
    molecules: List[Chem.rdchem.Mol],
    n_select: int,
    fingerprint_type: str = "morgan",
    params: Dict[str, Any] = None
) -> List[int]:
    """
    Selects a diverse subset of molecules using the MaxMin algorithm.
    
    Args:
        molecules: List of RDKit molecules
        n_select: Number of molecules to select
        fingerprint_type: Type of fingerprint to use (default: "morgan")
        params: Additional parameters for fingerprint generation (default: None)
        
    Returns:
        Indices of selected diverse molecules
        
    Raises:
        MolecularProcessingException: If fingerprint generation fails
    """
    if not molecules:
        return []
    
    n_mols = len(molecules)
    if n_select >= n_mols:
        # If we want to select as many or more molecules than we have,
        # just return all indices
        return list(range(n_mols))
    
    # Generate fingerprints for all molecules
    fingerprints = [mol_to_fingerprint(mol, fingerprint_type, params) for mol in molecules]
    
    # Initialize with a random molecule
    import random
    selected = [random.randint(0, n_mols - 1)]
    
    # Iteratively add the molecule that is most dissimilar to the already selected molecules
    while len(selected) < n_select:
        max_min_dist = -1
        max_idx = -1
        
        # For each unselected molecule
        for i in range(n_mols):
            if i in selected:
                continue
            
            # Find minimum distance (maximum similarity) to any selected molecule
            min_dist = float('inf')
            for j in selected:
                similarity = DataStructs.TanimotoSimilarity(fingerprints[i], fingerprints[j])
                dist = 1 - similarity
                min_dist = min(min_dist, dist)
            
            # Update max_min_dist if this molecule is more dissimilar
            if min_dist > max_min_dist:
                max_min_dist = min_dist
                max_idx = i
        
        # Add the most dissimilar molecule to the selected list
        selected.append(max_idx)
    
    return sorted(selected)


class SimilaritySearcher:
    """
    Class for performing molecular similarity searches and clustering.
    
    This class provides methods for calculating similarity between molecules,
    finding similar molecules in a dataset, clustering molecules based on
    structural similarity, and selecting diverse subsets of molecules.
    """
    
    def __init__(
        self,
        fingerprint_type: str = DEFAULT_FINGERPRINT_TYPE,
        fingerprint_params: Dict[str, Any] = None,
        similarity_threshold: float = DEFAULT_SIMILARITY_THRESHOLD
    ):
        """
        Initializes the SimilaritySearcher with configuration parameters.
        
        Args:
            fingerprint_type: Type of fingerprint to use (default: DEFAULT_FINGERPRINT_TYPE)
            fingerprint_params: Additional parameters for fingerprint generation (default: None)
            similarity_threshold: Default similarity threshold (default: DEFAULT_SIMILARITY_THRESHOLD)
            
        Raises:
            MolecularProcessingException: If an unsupported fingerprint type is specified
        """
        if fingerprint_type not in FINGERPRINT_TYPES:
            raise MolecularProcessingException(
                f"Unsupported fingerprint type: {fingerprint_type}",
                details={"supported_types": list(FINGERPRINT_TYPES.keys())}
            )
        
        self._fingerprint_type = fingerprint_type
        
        # Set default fingerprint parameters if none provided
        if fingerprint_params is None:
            fingerprint_params = {}
            if fingerprint_type == 'morgan':
                fingerprint_params = {'radius': DEFAULT_MORGAN_RADIUS, 'nBits': DEFAULT_MORGAN_BITS}
            elif fingerprint_type in ['topological', 'pattern']:
                fingerprint_params = {'nBits': DEFAULT_MORGAN_BITS}
        
        self._fingerprint_params = fingerprint_params
        self._similarity_threshold = similarity_threshold
        
        # Initialize logger
        self.logger = logging.getLogger(f"{__name__}.{self.__class__.__name__}")
        self.logger.debug(
            f"Initialized SimilaritySearcher with fingerprint_type={fingerprint_type}, "
            f"threshold={similarity_threshold}"
        )
    
    def set_fingerprint_type(self, fingerprint_type: str, params: Dict[str, Any] = None) -> None:
        """
        Sets the fingerprint type to use for similarity calculations.
        
        Args:
            fingerprint_type: Type of fingerprint to use
            params: Additional parameters for fingerprint generation (default: None)
            
        Raises:
            MolecularProcessingException: If an unsupported fingerprint type is specified
        """
        if fingerprint_type not in FINGERPRINT_TYPES:
            raise MolecularProcessingException(
                f"Unsupported fingerprint type: {fingerprint_type}",
                details={"supported_types": list(FINGERPRINT_TYPES.keys())}
            )
        
        self._fingerprint_type = fingerprint_type
        
        # Set default fingerprint parameters if none provided
        if params is None:
            params = {}
            if fingerprint_type == 'morgan':
                params = {'radius': DEFAULT_MORGAN_RADIUS, 'nBits': DEFAULT_MORGAN_BITS}
            elif fingerprint_type in ['topological', 'pattern']:
                params = {'nBits': DEFAULT_MORGAN_BITS}
        
        self._fingerprint_params = params
        self.logger.debug(f"Set fingerprint type to {fingerprint_type} with params {params}")
    
    def set_similarity_threshold(self, threshold: float) -> None:
        """
        Sets the similarity threshold for finding similar molecules.
        
        Args:
            threshold: Similarity threshold (between 0 and 1)
            
        Raises:
            MolecularProcessingException: If threshold is not between 0 and 1
        """
        if threshold < 0 or threshold > 1:
            raise MolecularProcessingException(
                f"Similarity threshold must be between 0 and 1, got {threshold}",
                details={"threshold": threshold}
            )
        
        self._similarity_threshold = threshold
        self.logger.debug(f"Set similarity threshold to {threshold}")
    
    def calculate_similarity(self, mol1: Chem.rdchem.Mol, mol2: Chem.rdchem.Mol) -> float:
        """
        Calculates similarity between two molecules.
        
        Args:
            mol1: First RDKit molecule
            mol2: Second RDKit molecule
            
        Returns:
            Similarity score between 0 and 1
            
        Raises:
            MolecularProcessingException: If either molecule is None
        """
        return calculate_similarity(mol1, mol2, self._fingerprint_type, self._fingerprint_params)
    
    def calculate_similarity_from_smiles(self, smiles1: str, smiles2: str) -> float:
        """
        Calculates similarity between two SMILES strings.
        
        Args:
            smiles1: SMILES string of first molecule
            smiles2: SMILES string of second molecule
            
        Returns:
            Similarity score between 0 and 1
            
        Raises:
            MolecularProcessingException: If SMILES conversion fails
        """
        return calculate_similarity_from_smiles(
            smiles1, smiles2, self._fingerprint_type, self._fingerprint_params
        )
    
    def find_similar(
        self,
        query: Union[str, Chem.rdchem.Mol],
        targets: List[Union[str, Chem.rdchem.Mol]],
        threshold: float = None
    ) -> List[Tuple[Union[str, Chem.rdchem.Mol], float]]:
        """
        Finds molecules similar to a query molecule.
        
        Args:
            query: Query molecule (as SMILES string or RDKit molecule)
            targets: List of target molecules (as SMILES strings or RDKit molecules)
            threshold: Minimum similarity threshold (default: self._similarity_threshold)
            
        Returns:
            List of (molecule, similarity) tuples for molecules above the threshold,
            sorted by decreasing similarity
            
        Raises:
            MolecularProcessingException: If query conversion fails
        """
        if threshold is None:
            threshold = self._similarity_threshold
        
        # Process based on query type
        if isinstance(query, str):
            # If query is a SMILES string, use the SMILES-specific function
            return find_similar_molecules_from_smiles(
                query, targets, threshold, self._fingerprint_type, self._fingerprint_params
            )
        else:
            # Process RDKit molecule query
            if query is None:
                raise MolecularProcessingException(
                    "Query molecule cannot be None",
                    details={"query": "None"}
                )
            
            # Convert targets to RDKit molecules for processing
            rdkit_molecules = []
            original_indices = {}  # Maps RDKit molecule (by id) to original index
            
            for idx, target in enumerate(targets):
                if isinstance(target, str):
                    # It's a SMILES string, convert to RDKit molecule
                    rdkit_mol = smiles_to_mol(target)
                    if rdkit_mol:  # Only add valid molecules
                        rdkit_molecules.append(rdkit_mol)
                        original_indices[id(rdkit_mol)] = idx
                else:
                    # It's already a RDKit molecule
                    if target:  # Only add non-None molecules
                        rdkit_molecules.append(target)
                        original_indices[id(target)] = idx
            
            # Find similar molecules using RDKit molecules
            similar_mol_tuples = find_similar_molecules(
                query, rdkit_molecules, threshold, self._fingerprint_type, self._fingerprint_params
            )
            
            # Convert results back to original format
            results = []
            for mol, similarity in similar_mol_tuples:
                # Get original index
                orig_idx = original_indices.get(id(mol))
                if orig_idx is not None:
                    results.append((targets[orig_idx], similarity))
            
            return results
    
    def calculate_similarity_matrix(
        self,
        molecules: List[Union[str, Chem.rdchem.Mol]]
    ) -> np.ndarray:
        """
        Calculates a pairwise similarity matrix for a list of molecules.
        
        Args:
            molecules: List of molecules (as SMILES strings or RDKit molecules)
            
        Returns:
            2D numpy array of pairwise similarity scores
            
        Raises:
            MolecularProcessingException: If molecule conversion fails
        """
        # Convert molecules to RDKit molecules if they're SMILES strings
        converted_mols = []
        
        for mol in molecules:
            if isinstance(mol, str):
                rdkit_mol = smiles_to_mol_with_exception(mol)
                converted_mols.append(rdkit_mol)
            else:
                if mol is None:
                    raise MolecularProcessingException(
                        "Cannot calculate similarity matrix with None molecules",
                        details={"molecules": "Contains None"}
                    )
                converted_mols.append(mol)
        
        return calculate_similarity_matrix(
            converted_mols, self._fingerprint_type, self._fingerprint_params
        )
    
    def cluster(
        self,
        molecules: List[Union[str, Chem.rdchem.Mol]],
        algorithm: str = "butina",
        threshold: float = None
    ) -> List[List[int]]:
        """
        Clusters molecules based on structural similarity.
        
        Args:
            molecules: List of molecules (as SMILES strings or RDKit molecules)
            algorithm: Clustering algorithm to use (default: "butina")
            threshold: Similarity threshold for clustering (default: self._similarity_threshold)
            
        Returns:
            List of clusters, where each cluster is a list of indices
            
        Raises:
            MolecularProcessingException: If molecule conversion fails or clustering fails
        """
        if threshold is None:
            threshold = self._similarity_threshold
        
        # Convert molecules to RDKit molecules if they're SMILES strings
        converted_mols = []
        
        for mol in molecules:
            if isinstance(mol, str):
                rdkit_mol = smiles_to_mol_with_exception(mol)
                converted_mols.append(rdkit_mol)
            else:
                if mol is None:
                    raise MolecularProcessingException(
                        "Cannot cluster None molecules",
                        details={"molecules": "Contains None"}
                    )
                converted_mols.append(mol)
        
        return cluster_molecules(
            converted_mols, algorithm, threshold, self._fingerprint_type, self._fingerprint_params
        )
    
    def select_diverse(
        self,
        molecules: List[Union[str, Chem.rdchem.Mol]],
        n_select: int
    ) -> List[int]:
        """
        Selects a diverse subset of molecules.
        
        Args:
            molecules: List of molecules (as SMILES strings or RDKit molecules)
            n_select: Number of molecules to select
            
        Returns:
            Indices of selected diverse molecules
            
        Raises:
            MolecularProcessingException: If molecule conversion fails
        """
        # Convert molecules to RDKit molecules if they're SMILES strings
        converted_mols = []
        
        for mol in molecules:
            if isinstance(mol, str):
                rdkit_mol = smiles_to_mol_with_exception(mol)
                converted_mols.append(rdkit_mol)
            else:
                if mol is None:
                    raise MolecularProcessingException(
                        "Cannot select diverse subset with None molecules",
                        details={"molecules": "Contains None"}
                    )
                converted_mols.append(mol)
        
        return select_diverse_subset(
            converted_mols, n_select, self._fingerprint_type, self._fingerprint_params
        )