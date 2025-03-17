"""
Unit Tests for Similarity Searcher Module

This module contains tests for the similarity searcher functionality,
which provides methods for calculating molecular similarity, finding similar
molecules, clustering molecules based on structural similarity, and selecting
diverse subsets of molecules.
"""

import pytest
import numpy as np
from typing import List, Dict, Any
from unittest.mock import patch

from rdkit import Chem
from rdkit import DataStructs

from ../../app/molecular/similarity_searcher import (
    calculate_similarity,
    calculate_similarity_from_smiles,
    calculate_similarity_matrix,
    find_similar_molecules,
    find_similar_molecules_from_smiles,
    cluster_molecules,
    select_diverse_subset,
    SimilaritySearcher,
    DEFAULT_SIMILARITY_THRESHOLD,
    DEFAULT_FINGERPRINT_TYPE,
)
from ../../app/molecular/molecule_converter import FINGERPRINT_TYPES
from ../../app/exceptions import MolecularProcessingException
from ../../app/molecular/molecule_converter import smiles_to_mol

# Test data
TEST_SMILES = ['CCO', 'CCCCO', 'c1ccccc1', 'CC(=O)O', 'CCN']
TEST_MOLS = [Chem.MolFromSmiles(s) for s in TEST_SMILES]


def test_calculate_similarity():
    """Tests the calculate_similarity function with valid molecules."""
    # Create two RDKit molecules
    mol1 = Chem.MolFromSmiles("CCO")  # Ethanol
    mol2 = Chem.MolFromSmiles("CCCCO")  # 1-Butanol
    
    # Calculate similarity with default fingerprint
    similarity = calculate_similarity(mol1, mol2)
    
    # Check that the result is a float between 0 and 1
    assert isinstance(similarity, float)
    assert 0 <= similarity <= 1
    
    # Test with different fingerprint types
    morgan_sim = calculate_similarity(mol1, mol2, fingerprint_type="morgan")
    maccs_sim = calculate_similarity(mol1, mol2, fingerprint_type="maccs")
    topo_sim = calculate_similarity(mol1, mol2, fingerprint_type="topological")
    
    # Check that different fingerprint types give different similarity values
    assert morgan_sim != maccs_sim or morgan_sim != topo_sim or maccs_sim != topo_sim


def test_calculate_similarity_invalid_molecules():
    """Tests the calculate_similarity function with invalid molecules."""
    mol = Chem.MolFromSmiles("CCO")  # Valid molecule
    
    # Try with None as first molecule
    with pytest.raises(MolecularProcessingException):
        calculate_similarity(None, mol)
    
    # Try with None as second molecule
    with pytest.raises(MolecularProcessingException):
        calculate_similarity(mol, None)
    
    # Try with invalid fingerprint type
    with pytest.raises(MolecularProcessingException):
        calculate_similarity(mol, mol, fingerprint_type="invalid_type")


def test_calculate_similarity_from_smiles():
    """Tests the calculate_similarity_from_smiles function."""
    # Define SMILES strings
    smiles1 = "CCO"  # Ethanol
    smiles2 = "CCCCO"  # 1-Butanol
    
    # Calculate similarity with default fingerprint
    similarity = calculate_similarity_from_smiles(smiles1, smiles2)
    
    # Check that the result is a float between 0 and 1
    assert isinstance(similarity, float)
    assert 0 <= similarity <= 1
    
    # Test with different fingerprint types
    morgan_sim = calculate_similarity_from_smiles(smiles1, smiles2, fingerprint_type="morgan")
    maccs_sim = calculate_similarity_from_smiles(smiles1, smiles2, fingerprint_type="maccs")
    topo_sim = calculate_similarity_from_smiles(smiles1, smiles2, fingerprint_type="topological")
    
    # Check that different fingerprint types give different similarity values
    assert morgan_sim != maccs_sim or morgan_sim != topo_sim or maccs_sim != topo_sim


def test_calculate_similarity_from_smiles_invalid():
    """Tests the calculate_similarity_from_smiles function with invalid SMILES."""
    valid_smiles = "CCO"  # Valid SMILES
    
    # Try with invalid first SMILES
    with pytest.raises(MolecularProcessingException):
        calculate_similarity_from_smiles("invalid_smiles", valid_smiles)
    
    # Try with invalid second SMILES
    with pytest.raises(MolecularProcessingException):
        calculate_similarity_from_smiles(valid_smiles, "invalid_smiles")
    
    # Try with empty SMILES strings
    with pytest.raises(MolecularProcessingException):
        calculate_similarity_from_smiles("", valid_smiles)
    
    with pytest.raises(MolecularProcessingException):
        calculate_similarity_from_smiles(valid_smiles, "")


def test_calculate_similarity_matrix():
    """Tests the calculate_similarity_matrix function."""
    # Test with the predefined molecule list
    matrix = calculate_similarity_matrix(TEST_MOLS)
    
    # Check that the result is a numpy array
    assert isinstance(matrix, np.ndarray)
    
    # Check dimensions of the matrix
    assert matrix.shape == (len(TEST_MOLS), len(TEST_MOLS))
    
    # Check that diagonal values are 1.0 (self-similarity)
    for i in range(len(TEST_MOLS)):
        assert matrix[i, i] == 1.0
    
    # Check that the matrix is symmetric
    for i in range(len(TEST_MOLS)):
        for j in range(i + 1, len(TEST_MOLS)):
            assert matrix[i, j] == matrix[j, i]
    
    # Test with different fingerprint types
    matrix_maccs = calculate_similarity_matrix(TEST_MOLS, fingerprint_type="maccs")
    # Check that the dimensions are still correct
    assert matrix_maccs.shape == (len(TEST_MOLS), len(TEST_MOLS))


def test_calculate_similarity_matrix_empty():
    """Tests the calculate_similarity_matrix function with an empty list."""
    # Calculate similarity matrix for an empty list
    matrix = calculate_similarity_matrix([])
    
    # Check that the result is an empty numpy array
    assert isinstance(matrix, np.ndarray)
    assert matrix.size == 0


def test_find_similar_molecules():
    """Tests the find_similar_molecules function."""
    # Select a query molecule and target molecules
    query_mol = TEST_MOLS[0]  # Ethanol
    target_mols = TEST_MOLS[1:]  # All other test molecules
    
    # Find similar molecules with default threshold
    similar_mols = find_similar_molecules(query_mol, target_mols)
    
    # Check the result structure
    assert isinstance(similar_mols, list)
    if similar_mols:
        assert isinstance(similar_mols[0], tuple)
        assert len(similar_mols[0]) == 2
        mol, similarity = similar_mols[0]
        assert isinstance(mol, Chem.rdchem.Mol)
        assert isinstance(similarity, float)
        assert 0 <= similarity <= 1
    
    # Check that all similarities are above the threshold
    for mol, similarity in similar_mols:
        assert similarity >= DEFAULT_SIMILARITY_THRESHOLD
    
    # Check that results are sorted by similarity (descending)
    for i in range(1, len(similar_mols)):
        assert similar_mols[i-1][1] >= similar_mols[i][1]
    
    # Test with a higher threshold
    high_threshold = 0.9
    similar_mols_high = find_similar_molecules(query_mol, target_mols, threshold=high_threshold)
    
    # Higher threshold should give fewer results
    assert len(similar_mols_high) <= len(similar_mols)
    
    # Test with different fingerprint type
    similar_mols_maccs = find_similar_molecules(query_mol, target_mols, fingerprint_type="maccs")
    # Just verify the result type is correct
    assert isinstance(similar_mols_maccs, list)


def test_find_similar_molecules_invalid():
    """Tests the find_similar_molecules function with invalid inputs."""
    # Valid molecule for testing
    valid_mol = TEST_MOLS[0]
    
    # Try with None as query molecule
    with pytest.raises(MolecularProcessingException):
        find_similar_molecules(None, TEST_MOLS)
    
    # Try with empty target list
    result = find_similar_molecules(valid_mol, [])
    assert result == []
    
    # Try with invalid threshold values
    # Note: The function doesn't explicitly validate threshold range, so this would only
    # fail if the implementation changes to add this validation


def test_find_similar_molecules_from_smiles():
    """Tests the find_similar_molecules_from_smiles function."""
    # Select a query SMILES and target SMILES
    query_smiles = TEST_SMILES[0]  # Ethanol
    target_smiles = TEST_SMILES[1:]  # All other test SMILES
    
    # Find similar molecules with default threshold
    similar_mols = find_similar_molecules_from_smiles(query_smiles, target_smiles)
    
    # Check the result structure
    assert isinstance(similar_mols, list)
    if similar_mols:
        assert isinstance(similar_mols[0], tuple)
        assert len(similar_mols[0]) == 2
        smiles, similarity = similar_mols[0]
        assert isinstance(smiles, str)
        assert isinstance(similarity, float)
        assert 0 <= similarity <= 1
    
    # Check that all similarities are above the threshold
    for mol, similarity in similar_mols:
        assert similarity >= DEFAULT_SIMILARITY_THRESHOLD
    
    # Check that results are sorted by similarity (descending)
    for i in range(1, len(similar_mols)):
        assert similar_mols[i-1][1] >= similar_mols[i][1]
    
    # Test with a mix of SMILES and RDKit molecules
    mixed_targets = [TEST_SMILES[1], TEST_MOLS[2]]
    similar_mixed = find_similar_molecules_from_smiles(query_smiles, mixed_targets)
    
    # Verify the function handles mixed input types
    assert isinstance(similar_mixed, list)


def test_find_similar_molecules_from_smiles_invalid():
    """Tests the find_similar_molecules_from_smiles function with invalid inputs."""
    # Valid SMILES for testing
    valid_smiles = TEST_SMILES[0]
    
    # Try with invalid query SMILES
    with pytest.raises(MolecularProcessingException):
        find_similar_molecules_from_smiles("invalid_smiles", TEST_SMILES)
    
    # Try with empty target list
    result = find_similar_molecules_from_smiles(valid_smiles, [])
    assert result == []
    
    # Try with some invalid SMILES in target list
    mixed_targets = [valid_smiles, "invalid_smiles", TEST_SMILES[1]]
    similar_mols = find_similar_molecules_from_smiles(valid_smiles, mixed_targets)
    
    # Should skip invalid SMILES and process valid ones
    assert isinstance(similar_mols, list)
    # Maximum 2 results possible (skipping the invalid one)
    assert len(similar_mols) <= 2


def test_cluster_molecules():
    """Tests the cluster_molecules function."""
    # Cluster the test molecules using default algorithm (Butina)
    clusters = cluster_molecules(TEST_MOLS)
    
    # Check the result structure
    assert isinstance(clusters, list)
    if clusters:
        assert isinstance(clusters[0], list)
        assert all(isinstance(idx, int) for cluster in clusters for idx in cluster)
    
    # Check that all molecule indices are included in clusters
    all_indices = [idx for cluster in clusters for idx in cluster]
    assert sorted(all_indices) == list(range(len(TEST_MOLS)))
    
    # Test with different threshold values
    high_threshold_clusters = cluster_molecules(TEST_MOLS, threshold=0.9)
    low_threshold_clusters = cluster_molecules(TEST_MOLS, threshold=0.1)
    
    # Higher threshold should give more clusters with fewer molecules each
    assert len(high_threshold_clusters) >= len(low_threshold_clusters)
    
    # Test with hierarchical algorithm
    hier_clusters = cluster_molecules(TEST_MOLS, algorithm="hierarchical")
    
    # Verify that different algorithms produce different clustering results
    assert isinstance(hier_clusters, list)
    # Note: The clusterings might coincidentally be the same, but they're generally different


def test_cluster_molecules_empty():
    """Tests the cluster_molecules function with an empty list."""
    # Cluster an empty list of molecules
    clusters = cluster_molecules([])
    
    # Check that the result is an empty list
    assert clusters == []


def test_cluster_molecules_invalid():
    """Tests the cluster_molecules function with invalid inputs."""
    # Try with invalid algorithm name
    with pytest.raises(MolecularProcessingException):
        cluster_molecules(TEST_MOLS, algorithm="invalid_algorithm")
    
    # Try with some invalid molecules in the list
    molecules_with_none = TEST_MOLS + [None]
    with pytest.raises(MolecularProcessingException):
        cluster_molecules(molecules_with_none)


def test_select_diverse_subset():
    """Tests the select_diverse_subset function."""
    # Select a diverse subset of molecules
    n_select = 3
    subset_indices = select_diverse_subset(TEST_MOLS, n_select)
    
    # Check the result structure
    assert isinstance(subset_indices, list)
    assert len(subset_indices) == n_select
    assert all(isinstance(idx, int) for idx in subset_indices)
    assert all(0 <= idx < len(TEST_MOLS) for idx in subset_indices)
    
    # Check that selected molecules are diverse
    subset_mols = [TEST_MOLS[idx] for idx in subset_indices]
    
    # Calculate pairwise similarities between selected molecules
    similarity_sum = 0
    count = 0
    for i in range(n_select):
        for j in range(i + 1, n_select):
            similarity = calculate_similarity(subset_mols[i], subset_mols[j])
            similarity_sum += similarity
            count += 1
    
    # Average similarity should be relatively low for diverse subset
    if count > 0:  # Avoid division by zero
        avg_similarity = similarity_sum / count
        assert avg_similarity < 0.8  # This threshold might need adjustment
    
    # Test with n_select greater than number of molecules
    large_n = len(TEST_MOLS) + 5
    all_indices = select_diverse_subset(TEST_MOLS, large_n)
    
    # Should return all molecules
    assert len(all_indices) == len(TEST_MOLS)
    assert sorted(all_indices) == list(range(len(TEST_MOLS)))


def test_select_diverse_subset_empty():
    """Tests the select_diverse_subset function with an empty list."""
    # Select diverse subset from an empty list
    subset = select_diverse_subset([], 5)
    
    # Check that the result is an empty list
    assert subset == []


def test_similarity_searcher_class():
    """Tests the SimilaritySearcher class functionality."""
    # Create a SimilaritySearcher instance
    searcher = SimilaritySearcher()
    
    # Test the calculate_similarity method
    mol1 = TEST_MOLS[0]
    mol2 = TEST_MOLS[1]
    similarity = searcher.calculate_similarity(mol1, mol2)
    assert isinstance(similarity, float)
    assert 0 <= similarity <= 1
    
    # Test the calculate_similarity_from_smiles method
    smiles1 = TEST_SMILES[0]
    smiles2 = TEST_SMILES[1]
    similarity = searcher.calculate_similarity_from_smiles(smiles1, smiles2)
    assert isinstance(similarity, float)
    assert 0 <= similarity <= 1
    
    # Test the find_similar method
    similar_mols = searcher.find_similar(mol1, TEST_MOLS)
    assert isinstance(similar_mols, list)
    
    # Test the calculate_similarity_matrix method
    matrix = searcher.calculate_similarity_matrix(TEST_MOLS)
    assert isinstance(matrix, np.ndarray)
    
    # Test the cluster method
    clusters = searcher.cluster(TEST_MOLS)
    assert isinstance(clusters, list)
    
    # Test the select_diverse method
    subset = searcher.select_diverse(TEST_MOLS, 3)
    assert isinstance(subset, list)
    assert len(subset) == 3


def test_similarity_searcher_set_fingerprint_type():
    """Tests the set_fingerprint_type method of SimilaritySearcher."""
    # Create a SimilaritySearcher instance
    searcher = SimilaritySearcher()
    
    # Test setting different fingerprint types
    for fp_type in FINGERPRINT_TYPES.keys():
        searcher.set_fingerprint_type(fp_type)
        # Calculate similarity to verify the fingerprint type is used
        similarity = searcher.calculate_similarity(TEST_MOLS[0], TEST_MOLS[1])
        assert isinstance(similarity, float)
    
    # Test with invalid fingerprint type
    with pytest.raises(MolecularProcessingException):
        searcher.set_fingerprint_type("invalid_type")
    
    # Test with custom fingerprint parameters
    custom_params = {"radius": 3, "nBits": 1024}
    searcher.set_fingerprint_type("morgan", custom_params)
    
    # Calculate similarity to verify the parameters are used
    similarity1 = searcher.calculate_similarity(TEST_MOLS[0], TEST_MOLS[1])
    
    # Change parameters and recalculate
    different_params = {"radius": 2, "nBits": 2048}
    searcher.set_fingerprint_type("morgan", different_params)
    similarity2 = searcher.calculate_similarity(TEST_MOLS[0], TEST_MOLS[1])
    
    # Different parameters should (potentially) give different results
    # Note: This might not always be true, but is a reasonable expectation
    assert similarity1 != similarity2  # This might occasionally fail if parameters don't affect these specific molecules


def test_similarity_searcher_set_similarity_threshold():
    """Tests the set_similarity_threshold method of SimilaritySearcher."""
    # Create a SimilaritySearcher instance
    searcher = SimilaritySearcher()
    
    # Test setting different thresholds
    test_thresholds = [0.1, 0.5, 0.9]
    for threshold in test_thresholds:
        searcher.set_similarity_threshold(threshold)
        assert searcher._similarity_threshold == threshold
    
    # Test with invalid threshold values
    with pytest.raises(MolecularProcessingException):
        searcher.set_similarity_threshold(1.5)  # Greater than 1
    
    with pytest.raises(MolecularProcessingException):
        searcher.set_similarity_threshold(-0.2)  # Negative
    
    # Test the effect of threshold on find_similar results
    searcher.set_similarity_threshold(0.1)  # Low threshold
    similar_low = searcher.find_similar(TEST_MOLS[0], TEST_MOLS)
    
    searcher.set_similarity_threshold(0.9)  # High threshold
    similar_high = searcher.find_similar(TEST_MOLS[0], TEST_MOLS)
    
    # Higher threshold should result in fewer similar molecules
    assert len(similar_high) <= len(similar_low)


@pytest.mark.slow
def test_performance_large_dataset():
    """Tests the performance of similarity functions with larger datasets."""
    # This test is marked as slow and might be skipped in regular test runs
    
    # Generate a larger set of test molecules (100+)
    import random
    large_mols = []
    for _ in range(100):  # Create 100 molecules
        # Pick a random molecule from our test set
        base_mol = random.choice(TEST_MOLS)
        # We'll just use the same molecule multiple times for this test
        large_mols.append(base_mol)
    
    # Measure time to calculate similarity matrix
    import time
    
    start_time = time.time()
    _ = calculate_similarity_matrix(large_mols)
    matrix_time = time.time() - start_time
    
    # Measure time to cluster molecules
    start_time = time.time()
    _ = cluster_molecules(large_mols)
    cluster_time = time.time() - start_time
    
    # Measure time to select diverse subset
    start_time = time.time()
    _ = select_diverse_subset(large_mols, 10)
    diverse_time = time.time() - start_time
    
    # Assert that performance meets requirements
    assert matrix_time < 10.0, f"Similarity matrix calculation took {matrix_time:.2f} seconds"
    assert cluster_time < 10.0, f"Clustering took {cluster_time:.2f} seconds"
    assert diverse_time < 10.0, f"Diverse subset selection took {diverse_time:.2f} seconds"