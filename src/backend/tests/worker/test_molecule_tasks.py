"""Unit tests for the Celery tasks defined in the molecule_tasks module.
This file tests the asynchronous processing of molecular data, including
property calculation, batch processing, similarity calculations, and other
molecular operations performed in background tasks.
"""

# Standard library imports
from typing import List, Dict, Any

# External package imports
import pytest  # pytest 7.0+
from unittest import mock  # standard library
import numpy as np  # numpy 1.24+

# Internal imports
from app.worker.tasks.molecule_tasks import (
    calculate_molecule_properties,
    batch_process_molecules,
    stream_process_molecules,
    calculate_batch_properties,
    calculate_similarity,
    cluster_by_similarity,
    select_diverse_molecules,
    standardize_molecule_set,
    check_molecules_drug_likeness,
    bulk_update_molecule_properties,
    MolecularProcessingException,
)

# Define global test data
VALID_SMILES = ['CCO', 'CCCCO', 'c1ccccc1', 'CC(=O)O', 'CCN']
INVALID_SMILES = ['XX', 'C1C', 'invalid']
TEST_PROPERTIES = ['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors']


@pytest.mark.parametrize('smiles', VALID_SMILES)
def test_calculate_molecule_properties_success(smiles: str) -> None:
    """Tests successful property calculation for a single molecule"""
    with mock.patch('app.worker.tasks.molecule_tasks.process_molecule') as mock_process_molecule:
        mock_process_molecule.return_value = {'smiles': smiles, 'properties': {'molecular_weight': 46.07}}
        result = calculate_molecule_properties(smiles)
        mock_process_molecule.assert_called_once_with(smiles, None)
        assert 'molecular_weight' in result['properties']


@pytest.mark.parametrize('smiles', INVALID_SMILES)
def test_calculate_molecule_properties_failure(smiles: str) -> None:
    """Tests property calculation failure for an invalid molecule"""
    with mock.patch('app.worker.tasks.molecule_tasks.process_molecule') as mock_process_molecule:
        mock_process_molecule.side_effect = MolecularProcessingException("Invalid SMILES")
        with pytest.raises(MolecularProcessingException):
            calculate_molecule_properties(smiles)
        mock_process_molecule.assert_called_once_with(smiles, None)


def test_batch_process_molecules_success() -> None:
    """Tests successful batch processing of multiple molecules"""
    smiles_list = VALID_SMILES
    with mock.patch('app.worker.tasks.molecule_tasks.process_molecules_batch') as mock_process_molecules_batch:
        mock_process_molecules_batch.return_value = [{'smiles': smiles, 'properties': {'molecular_weight': 46.07}} for smiles in smiles_list]
        result = batch_process_molecules(smiles_list)
        mock_process_molecules_batch.assert_called_once_with(smiles_list, None, None)
        assert len(result['results']) == len(smiles_list)
        assert result['num_molecules'] == len(smiles_list)


def test_batch_process_molecules_partial_failure() -> None:
    """Tests batch processing with some failures"""
    smiles_list = VALID_SMILES + INVALID_SMILES
    with mock.patch('app.worker.tasks.molecule_tasks.process_molecules_batch') as mock_process_molecules_batch:
        valid_results = [{'smiles': smiles, 'properties': {'molecular_weight': 46.07}} for smiles in VALID_SMILES]
        invalid_results = [{'smiles': smiles, 'error': 'Invalid SMILES'} for smiles in INVALID_SMILES]
        mock_process_molecules_batch.return_value = valid_results + invalid_results
        result = batch_process_molecules(smiles_list)
        mock_process_molecules_batch.assert_called_once_with(smiles_list, None, None)
        assert len(result['results']) == len(smiles_list)
        assert result['num_molecules'] == len(smiles_list)


def test_stream_process_molecules_success() -> None:
    """Tests successful stream processing of large molecule datasets"""
    smiles_list = VALID_SMILES * 100
    with mock.patch('app.worker.tasks.molecule_tasks.process_molecules_stream') as mock_process_molecules_stream:
        mock_process_molecules_stream.return_value = [{'smiles': smiles, 'properties': {'molecular_weight': 46.07}} for smiles in smiles_list]
        result = stream_process_molecules(smiles_list)
        mock_process_molecules_stream.assert_called_once_with(smiles_list, None, None, None)
        assert len(result['results']) == len(smiles_list)
        assert result['num_molecules'] == len(smiles_list)


def test_calculate_batch_properties_success() -> None:
    """Tests successful property calculation for multiple molecules"""
    smiles_list = VALID_SMILES
    with mock.patch('app.worker.tasks.molecule_tasks.batch_calculate_properties') as mock_batch_calculate_properties:
        mock_batch_calculate_properties.return_value = [{'molecular_weight': 46.07} for _ in smiles_list]
        result = calculate_batch_properties(smiles_list, TEST_PROPERTIES)
        mock_batch_calculate_properties.assert_called_once_with(smiles_list, TEST_PROPERTIES, None)
        assert len(result['results']) == len(smiles_list)


def test_calculate_similarity_success() -> None:
    """Tests successful similarity matrix calculation"""
    smiles_list = VALID_SMILES
    similarity_matrix = np.random.rand(len(smiles_list), len(smiles_list))
    with mock.patch('app.worker.tasks.molecule_tasks.calculate_similarity_matrix') as mock_calculate_similarity_matrix:
        mock_calculate_similarity_matrix.return_value = similarity_matrix
        result = calculate_similarity(smiles_list)
        mock_calculate_similarity_matrix.assert_called_once_with(smiles_list, 'morgan')
        assert 'similarity_matrix' in result
        assert len(result['similarity_matrix']) == len(smiles_list)


def test_cluster_by_similarity_success() -> None:
    """Tests successful molecule clustering based on similarity"""
    smiles_list = VALID_SMILES
    cluster_assignments = [0, 1, 0, 2, 1]
    with mock.patch('app.worker.tasks.molecule_tasks.cluster_molecules') as mock_cluster_molecules:
        mock_cluster_molecules.return_value = cluster_assignments
        result = cluster_by_similarity(smiles_list, threshold=0.7)
        mock_cluster_molecules.assert_called_once_with(smiles_list, 0.7, 'morgan')
        assert result['cluster_assignments'] == cluster_assignments


def test_select_diverse_molecules_success() -> None:
    """Tests successful selection of diverse molecule subset"""
    smiles_list = VALID_SMILES
    subset_indices = [0, 2, 4]
    with mock.patch('app.worker.tasks.molecule_tasks.select_diverse_subset') as mock_select_diverse_subset:
        mock_select_diverse_subset.return_value = subset_indices
        result = select_diverse_molecules(smiles_list, subset_size=3)
        mock_select_diverse_subset.assert_called_once_with(smiles_list, 3, 'morgan')
        assert len(result['selected_molecules']) == 3
        assert result['subset_size'] == 3


def test_standardize_molecule_set_success() -> None:
    """Tests successful standardization of molecules to canonical form"""
    smiles_list = VALID_SMILES
    standardized_smiles = ['[H]C([H])([H])O', '[H]C([H])([H])C([H])([H])C([H])([H])O', 'c1([H])c([H])c([H])c([H])c([H])c1([H])', '[H]C([H])([H])C(=O)O', '[H]C([H])([H])N']
    with mock.patch('app.worker.tasks.molecule_tasks.standardize_molecules') as mock_standardize_molecules:
        mock_standardize_molecules.return_value = standardized_smiles
        result = standardize_molecule_set(smiles_list)
        mock_standardize_molecules.assert_called_once_with(smiles_list)
        assert result['standardized_molecules'] == standardized_smiles


def test_check_molecules_drug_likeness_success() -> None:
    """Tests successful drug-likeness checking for multiple molecules"""
    smiles_list = VALID_SMILES
    with mock.patch('app.worker.tasks.molecule_tasks.check_drug_likeness') as mock_check_drug_likeness:
        mock_check_drug_likeness.return_value = {'is_drug_like': True}
        result = check_molecules_drug_likeness(smiles_list)
        assert len(result['results']) == len(smiles_list)
        for molecule_result in result['results']:
            assert 'drug_likeness' in molecule_result


def test_bulk_update_molecule_properties_success() -> None:
    """Tests successful bulk update of molecule properties in database"""
    molecules_data = [{'id': 1, 'flag_status': 'HIGH'}, {'id': 2, 'flag_status': 'LOW'}]
    with mock.patch('app.worker.tasks.molecule_tasks.molecule_service.update_molecule') as mock_update_molecule:
        mock_update_molecule.return_value = {'id': 1, 'flag_status': 'HIGH'}
        result = bulk_update_molecule_properties(molecules_data, user_id=1)
        assert mock_update_molecule.call_count == len(molecules_data)
        assert result['success_count'] == len(molecules_data)
        assert result['failure_count'] == 0


def test_bulk_update_molecule_properties_partial_failure() -> None:
    """Tests bulk update with some failures"""
    molecules_data = [{'id': 1, 'flag_status': 'HIGH'}, {'id': 2, 'flag_status': 'INVALID'}]
    with mock.patch('app.worker.tasks.molecule_tasks.molecule_service.update_molecule') as mock_update_molecule:
        mock_update_molecule.side_effect = [
            {'id': 1, 'flag_status': 'HIGH'},
            Exception('Update failed')
        ]
        result = bulk_update_molecule_properties(molecules_data, user_id=1)
        assert mock_update_molecule.call_count == len(molecules_data)
        assert result['success_count'] == 1
        assert result['failure_count'] == 1


def test_task_retry_on_exception() -> None:
    """Tests that tasks retry on exception"""
    @app.task(name='test_task', bind=True, max_retries=3)
    def test_task(self):
        raise Exception("Test exception")

    with mock.patch.object(test_task, 'retry') as mock_retry:
        with pytest.raises(Exception, match="Test exception"):
            test_task()
        assert mock_retry.called