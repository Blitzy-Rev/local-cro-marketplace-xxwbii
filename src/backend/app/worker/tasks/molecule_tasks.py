# src/backend/app/worker/tasks/molecule_tasks.py
"""
Celery task module for asynchronous molecular data processing in the Molecular Data Management and CRO Integration Platform.
Implements background tasks for batch processing of molecules, property calculation, similarity searching, and other computationally intensive molecular operations.
"""

# Standard library imports
import typing
from typing import Dict, List, Any, Optional, Union

# External package imports
import celery  # celery 5.2+

# Internal imports
from app.worker.celery_app import app  # Import Celery application instance for task registration
from app.logging_config import logger  # Application logger for task logging
from app.molecular.processor import process_molecule  # Process a single molecule with validation and property calculation
from app.molecular.processor import process_molecules_batch  # Process a batch of molecules in parallel
from app.molecular.processor import process_molecules_stream  # Process large datasets of molecules efficiently in batches
from app.molecular.processor import calculate_similarity_matrix  # Calculate pairwise similarity matrix for a list of molecules
from app.molecular.processor import cluster_molecules  # Cluster molecules based on structural similarity
from app.molecular.processor import select_diverse_subset  # Select diverse subset of molecules
from app.molecular.processor import standardize_molecules  # Standardize molecules to canonical form
from app.molecular.processor import check_drug_likeness  # Check drug-likeness using Lipinski's Rule of Five
from app.molecular.property_calculator import calculate_properties  # Calculate multiple molecular properties
from app.molecular.property_calculator import batch_calculate_properties  # Calculate properties for multiple molecules in parallel
from app.molecular.property_calculator import stream_calculate_properties  # Calculate properties for large datasets in batches
from app.services.molecule_service import MoleculeService  # Service for molecule operations
from app.exceptions import MolecularProcessingException  # Exception class for molecular processing errors
from app.constants import BATCH_SIZE  # Default batch size for processing molecules

# Initialize global variables
logger = typing.cast(logging.Logger, logging.getLogger(__name__))
molecule_service = typing.cast(MoleculeService, MoleculeService())

@app.task(name='molecule_tasks.calculate_molecule_properties', max_retries=3, acks_late=True)
def calculate_molecule_properties(smiles: str, properties: List[str] = None) -> Dict[str, Any]:
    """
    Celery task for calculating properties for a single molecule
    
    Args:
        smiles (str): SMILES string representation of the molecule
        properties (List[str], optional): List of property names to calculate. Defaults to None.
    
    Returns:
        Dict[str, Any]: Dictionary with SMILES and calculated properties
    """
    logger.info(f"Starting calculate_molecule_properties task for SMILES: {smiles}")
    try:
        # Try to process the molecule using process_molecule
        result = process_molecule(smiles, properties)
        return result
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error calculating properties for SMILES {smiles}: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise calculate_molecule_properties.retry(exc=e)

@app.task(name='molecule_tasks.batch_process_molecules', max_retries=2, acks_late=True)
def batch_process_molecules(smiles_list: List[str], properties: List[str] = None, num_workers: int = None) -> Dict[str, Any]:
    """
    Celery task for processing a batch of molecules in parallel
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        properties (List[str], optional): List of property names to calculate. Defaults to None.
        num_workers (int, optional): Number of workers for parallel processing. Defaults to None.
    
    Returns:
        Dict[str, Any]: Processing results with statistics
    """
    logger.info(f"Starting batch_process_molecules task for {len(smiles_list)} molecules")
    try:
        # Try to process the batch using process_molecules_batch
        results = process_molecules_batch(smiles_list, properties, num_workers)
        return {"results": results, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error processing batch of molecules: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise batch_process_molecules.retry(exc=e)

@app.task(name='molecule_tasks.stream_process_molecules', max_retries=2, acks_late=True, time_limit=3600)
def stream_process_molecules(smiles_list: List[str], properties: List[str] = None, batch_size: int = None, num_workers: int = None) -> Dict[str, Any]:
    """
    Celery task for processing large datasets of molecules in batches
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        properties (List[str], optional): List of property names to calculate. Defaults to None.
        batch_size (int, optional): Batch size for processing. Defaults to None.
        num_workers (int, optional): Number of workers for parallel processing. Defaults to None.
    
    Returns:
        Dict[str, Any]: Processing results with statistics
    """
    logger.info(f"Starting stream_process_molecules task for {len(smiles_list)} molecules")
    try:
        # If batch_size is None, use BATCH_SIZE from constants
        if batch_size is None:
            batch_size = BATCH_SIZE
        # Try to process the molecules using process_molecules_stream
        results = process_molecules_stream(smiles_list, properties, batch_size, num_workers)
        return {"results": results, "num_molecules": len(smiles_list), "batch_size": batch_size}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error processing stream of molecules: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise stream_process_molecules.retry(exc=e)

@app.task(name='molecule_tasks.calculate_batch_properties', max_retries=2, acks_late=True)
def calculate_batch_properties(smiles_list: List[str], properties: List[str] = None, num_workers: int = None) -> Dict[str, Any]:
    """
    Celery task for calculating properties for a batch of molecules
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        properties (List[str], optional): List of property names to calculate. Defaults to None.
        num_workers (int, optional): Number of workers for parallel processing. Defaults to None.
    
    Returns:
        Dict[str, Any]: Dictionary with molecules and their properties
    """
    logger.info(f"Starting calculate_batch_properties task for {len(smiles_list)} molecules")
    try:
        # Try to calculate properties using batch_calculate_properties
        results = batch_calculate_properties(smiles_list, properties, num_workers)
        return {"results": results, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error calculating batch properties: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise calculate_batch_properties.retry(exc=e)

@app.task(name='molecule_tasks.calculate_similarity', max_retries=2, acks_late=True)
def calculate_similarity(smiles_list: List[str], fingerprint_type: str = 'morgan') -> Dict[str, Any]:
    """
    Celery task for calculating similarity matrix for a list of molecules
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        fingerprint_type (str, optional): Fingerprint type for similarity calculation. Defaults to 'morgan'.
    
    Returns:
        Dict[str, Any]: Dictionary with similarity matrix
    """
    logger.info(f"Starting calculate_similarity task for {len(smiles_list)} molecules")
    try:
        # Try to calculate similarity matrix using calculate_similarity_matrix
        similarity_matrix = calculate_similarity_matrix(smiles_list, fingerprint_type)
        # Convert numpy array to list for JSON serialization
        similarity_list = similarity_matrix.tolist()
        return {"similarity_matrix": similarity_list, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error calculating similarity matrix: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise calculate_similarity.retry(exc=e)

@app.task(name='molecule_tasks.cluster_by_similarity', max_retries=2, acks_late=True)
def cluster_by_similarity(smiles_list: List[str], threshold: float = 0.7, fingerprint_type: str = 'morgan') -> Dict[str, Any]:
    """
    Celery task for clustering molecules based on structural similarity
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        threshold (float, optional): Similarity threshold for clustering. Defaults to 0.7.
        fingerprint_type (str, optional): Fingerprint type for similarity calculation. Defaults to 'morgan'.
    
    Returns:
        Dict[str, Any]: Dictionary with cluster assignments
    """
    logger.info(f"Starting cluster_by_similarity task for {len(smiles_list)} molecules with threshold {threshold}")
    try:
        # Try to cluster molecules using cluster_molecules
        cluster_assignments = cluster_molecules(smiles_list, threshold, fingerprint_type)
        return {"cluster_assignments": cluster_assignments, "num_molecules": len(smiles_list), "threshold": threshold}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error clustering molecules by similarity: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise cluster_by_similarity.retry(exc=e)

@app.task(name='molecule_tasks.select_diverse_molecules', max_retries=2, acks_late=True)
def select_diverse_molecules(smiles_list: List[str], subset_size: int, fingerprint_type: str = 'morgan') -> Dict[str, Any]:
    """
    Celery task for selecting a diverse subset of molecules
    
    Args:
        smiles_list (List[str]): List of SMILES strings
        subset_size (int): Size of the diverse subset to select
        fingerprint_type (str, optional): Fingerprint type for diversity selection. Defaults to 'morgan'.
    
    Returns:
        Dict[str, Any]: Dictionary with selected diverse molecules
    """
    logger.info(f"Starting select_diverse_molecules task for {len(smiles_list)} molecules, subset size {subset_size}")
    try:
        # Try to select diverse subset using select_diverse_subset
        selected_indices = select_diverse_subset(smiles_list, subset_size, fingerprint_type)
        # Extract selected molecules using indices
        selected_molecules = [smiles_list[i] for i in selected_indices]
        return {"selected_molecules": selected_molecules, "subset_size": subset_size, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error selecting diverse molecules: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise select_diverse_molecules.retry(exc=e)

@app.task(name='molecule_tasks.standardize_molecule_set', max_retries=2, acks_late=True)
def standardize_molecule_set(smiles_list: List[str]) -> Dict[str, Any]:
    """
    Celery task for standardizing a set of molecules to canonical form
    
    Args:
        smiles_list (List[str]): List of SMILES strings
    
    Returns:
        Dict[str, Any]: Dictionary with standardized molecules
    """
    logger.info(f"Starting standardize_molecule_set task for {len(smiles_list)} molecules")
    try:
        # Try to standardize molecules using standardize_molecules
        standardized_molecules = standardize_molecules(smiles_list)
        return {"standardized_molecules": standardized_molecules, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error standardizing molecules: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise standardize_molecule_set.retry(exc=e)

@app.task(name='molecule_tasks.check_molecules_drug_likeness', max_retries=2, acks_late=True)
def check_molecules_drug_likeness(smiles_list: List[str]) -> Dict[str, Any]:
    """
    Celery task for checking drug-likeness of multiple molecules
    
    Args:
        smiles_list (List[str]): List of SMILES strings
    
    Returns:
        Dict[str, Any]: Dictionary with drug-likeness assessment for each molecule
    """
    logger.info(f"Starting check_molecules_drug_likeness task for {len(smiles_list)} molecules")
    results = []
    try:
        # For each SMILES string in the input list:
        for smiles in smiles_list:
            try:
                # Try to check drug-likeness using check_drug_likeness
                drug_likeness = check_drug_likeness(smiles)
                results.append({"smiles": smiles, "drug_likeness": drug_likeness})
            except Exception as e:
                logger.error(f"Error checking drug-likeness for {smiles}: {str(e)}")
                results.append({"smiles": smiles, "error": str(e)})
        return {"results": results, "num_molecules": len(smiles_list)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error checking drug-likeness for molecules: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise check_molecules_drug_likeness.retry(exc=e)

@app.task(name='molecule_tasks.bulk_update_molecule_properties', max_retries=2, acks_late=True)
def bulk_update_molecule_properties(molecules_data: List[Dict[str, Any]], user_id: int) -> Dict[str, Any]:
    """
    Celery task for updating properties of multiple molecules in the database
    
    Args:
        molecules_data (List[Dict[str, Any]]): List of molecule data dictionaries
        user_id (int): ID of the user performing the update
    
    Returns:
        Dict[str, Any]: Dictionary with update statistics
    """
    logger.info(f"Starting bulk_update_molecule_properties task for {len(molecules_data)} molecules")
    success_count = 0
    failure_count = 0
    try:
        # For each molecule in molecules_data:
        for molecule_data in molecules_data:
            try:
                # Try to update molecule using molecule_service.update_molecule
                molecule_id = molecule_data.get("id")
                if not molecule_id:
                    logger.warning(f"Molecule ID missing, skipping update")
                    failure_count += 1
                    continue
                
                # Create a MoleculeUpdate object from the molecule data
                molecule_update = MoleculeUpdate(**molecule_data)
                
                molecule_service.update_molecule(molecule_id, molecule_update)
                success_count += 1
            except Exception as e:
                logger.error(f"Error updating molecule {molecule_data.get('id')}: {str(e)}")
                failure_count += 1
        return {"success_count": success_count, "failure_count": failure_count, "num_molecules": len(molecules_data)}
    except Exception as e:
        # If exception occurs, log error details
        logger.error(f"Error updating molecules in bulk: {str(e)}")
        # Retry the task if appropriate (up to max_retries)
        raise bulk_update_molecule_properties.retry(exc=e)