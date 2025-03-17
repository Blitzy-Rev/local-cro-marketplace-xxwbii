"""
Unit tests for the molecular processor module which handles batch processing of molecules,
validation, property calculation, filtering, and sorting.
"""

import pytest
from typing import List, Dict, Any
import time
from unittest.mock import Mock, patch
from rdkit import Chem

from ../../app.molecular.processor import (
    process_molecule,
    process_molecules_batch,
    process_molecules_stream,
    filter_molecules,
    sort_molecules,
    deduplicate_molecules,
    enrich_molecules,
    validate_and_standardize_molecules,
    generate_processing_summary,
    MoleculeProcessor,
    MolecularProcessingException
)
from ../../app.molecular.validator import PROPERTY_RANGES
from ../../app.molecular.property_calculator import PROPERTY_CALCULATORS


# Valid SMILES for testing
VALID_SMILES = [
    "CCO",              # Ethanol
    "CC(=O)O",          # Acetic acid
    "c1ccccc1",         # Benzene
    "C1=CC=CC=C1",      # Benzene (alternate)
    "CC1=CC=CC=C1",     # Toluene
    "COC",              # Dimethyl ether
    "CCC",              # Propane
    "CCCC",             # Butane
    "C1CCCCC1",         # Cyclohexane
    "OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO"  # Glucose
]

# Invalid SMILES for testing
INVALID_SMILES = [
    "",                 # Empty string
    "invalid",          # Not a SMILES string
    "C1CC",             # Incomplete ring
    "C%C",              # Invalid character
    "C(C)(C)(C)(C)",    # Too many bonds
    "CC(CC",            # Unclosed parenthesis
]


def test_process_molecule_valid():
    """Tests processing a valid molecule with default properties."""
    # Define a valid SMILES string
    smiles = "CCO"  # Ethanol
    
    # Call process_molecule with the SMILES string
    result = process_molecule(smiles)
    
    # Assert that the result contains the SMILES string
    assert result["smiles"] == smiles
    
    # Assert that the structure validation result is True
    assert result["validation"]["structure_valid"] is True
    
    # Properties are not calculated unless requested
    assert "properties" not in result


def test_process_molecule_invalid():
    """Tests processing an invalid molecule."""
    # Define an invalid SMILES string
    smiles = "C1CC"  # Incomplete ring
    
    # Call process_molecule with the invalid SMILES string
    result = process_molecule(smiles)
    
    # Assert that the result contains the SMILES string
    assert result["smiles"] == smiles
    
    # Assert that the structure validation result is False
    assert result["validation"]["structure_valid"] is False
    
    # Assert that the result does not contain calculated properties
    assert "properties" not in result


def test_process_molecule_with_properties():
    """Tests processing a molecule with specific properties."""
    # Define a valid SMILES string
    smiles = "CCO"  # Ethanol
    
    # Define a list of specific properties to calculate
    properties = ["molecular_weight", "logp"]
    
    # Call process_molecule with the SMILES string and properties list
    result = process_molecule(smiles, properties)
    
    # Assert that the result contains the SMILES string
    assert result["smiles"] == smiles
    
    # Assert that the validation result is True
    assert result["validation"]["structure_valid"] is True
    
    # Assert that the result contains properties
    assert "properties" in result
    
    # Assert that the result contains the specified properties
    for prop in properties:
        assert prop in result["properties"]
    
    # Assert that the result does not contain other properties
    for prop in result["properties"]:
        assert prop in properties


def test_process_molecules_batch():
    """Tests batch processing of multiple molecules."""
    # Define a list of valid SMILES strings
    smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
    
    # Call process_molecules_batch with the SMILES list
    results = process_molecules_batch(smiles_list)
    
    # Assert that the result list has the same length as the input list
    assert len(results) == len(smiles_list)
    
    # Assert that each result contains the corresponding SMILES string
    for i, result in enumerate(results):
        assert result["smiles"] == smiles_list[i]
    
    # Assert that valid molecules have validation results
    for result in results:
        assert "validation" in result
        assert result["validation"]["structure_valid"] is True


def test_process_molecules_batch_empty():
    """Tests batch processing with an empty list."""
    # Call process_molecules_batch with an empty list
    results = process_molecules_batch([])
    
    # Assert that the result is an empty list
    assert results == []


def test_process_molecules_batch_with_properties():
    """Tests batch processing with specific properties."""
    # Define a list of valid SMILES strings
    smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
    
    # Define a list of specific properties to calculate
    properties = ["molecular_weight", "logp"]
    
    # Call process_molecules_batch with the SMILES list and properties list
    results = process_molecules_batch(smiles_list, properties)
    
    # Assert that the result list has the same length as the input list
    assert len(results) == len(smiles_list)
    
    # Assert that each result contains only the specified properties
    for result in results:
        assert "properties" in result
        for prop in properties:
            assert prop in result["properties"]
    
    # Assert that results do not contain other properties
    for result in results:
        if "properties" in result:
            for prop in result["properties"]:
                assert prop in properties


def test_process_molecules_stream():
    """Tests stream processing of a large number of molecules."""
    # Define a list of valid SMILES strings (larger dataset)
    smiles_list = VALID_SMILES * 10  # Multiply the list to get more molecules
    
    # Call process_molecules_stream with the SMILES list
    results = process_molecules_stream(smiles_list, batch_size=20)
    
    # Assert that the result list has the same length as the input list
    assert len(results) == len(smiles_list)
    
    # Assert that each result contains a SMILES string and validation
    for result in results:
        assert "smiles" in result
        assert "validation" in result
        assert result["validation"]["structure_valid"] is True


def test_process_molecules_stream_performance():
    """Tests the performance of stream processing for large datasets."""
    # Generate a list of 1,000 valid SMILES strings
    smiles_list = VALID_SMILES * 100  # 100 copies of our 10 valid SMILES gives 1,000
    
    # Measure the time taken to process the stream
    start_time = time.time()
    results = process_molecules_stream(smiles_list, batch_size=250)
    end_time = time.time()
    
    # Calculate the processing time
    processing_time = end_time - start_time
    
    # Calculate the estimated time for 10,000 molecules (linear scaling assumption)
    estimated_time_for_10k = processing_time * 10  # Scale up by factor of 10
    
    # Assert that the estimated time is less than 30 seconds
    assert estimated_time_for_10k < 30, f"Estimated time for 10,000 molecules: {estimated_time_for_10k:.2f}s > 30s"
    
    # Also verify that all molecules were processed correctly
    assert len(results) == len(smiles_list)


def test_filter_molecules():
    """Tests filtering molecules based on property values."""
    # Define a list of processed molecules with various properties
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07, "logp": -0.14}
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 58.12, "logp": 1.0}
        },
        {
            "smiles": "c1ccccc1",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 78.11, "logp": 2.0}
        }
    ]
    
    # Define property filters with min/max ranges
    property_filters = {
        "molecular_weight": {"min": 50, "max": 80},
        "logp": {"min": 0, "max": 2.5}
    }
    
    # Call filter_molecules with the molecules and property filters
    filtered_molecules = filter_molecules(molecules, property_filters)
    
    # Assert that the filtered list contains only molecules within the specified ranges
    assert len(filtered_molecules) == 2
    for molecule in filtered_molecules:
        assert molecule["properties"]["molecular_weight"] >= 50
        assert molecule["properties"]["molecular_weight"] <= 80
        assert molecule["properties"]["logp"] >= 0
        assert molecule["properties"]["logp"] <= 2.5


def test_filter_molecules_with_substructure():
    """Tests filtering molecules based on substructure."""
    # Define a list of processed molecules with various structures
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True}
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True}
        },
        {
            "smiles": "c1ccccc1",
            "validation": {"structure_valid": True}
        },
        {
            "smiles": "CC(=O)O",
            "validation": {"structure_valid": True}
        }
    ]
    
    # Define a substructure SMARTS pattern - matches molecules containing oxygen
    substructure = "O"
    
    # Call filter_molecules with the molecules and substructure
    filtered_molecules = filter_molecules(molecules, substructure=substructure)
    
    # Assert that the filtered list contains only molecules with the substructure
    assert len(filtered_molecules) == 2
    assert any(mol["smiles"] == "CCO" for mol in filtered_molecules)
    assert any(mol["smiles"] == "CC(=O)O" for mol in filtered_molecules)


def test_sort_molecules():
    """Tests sorting molecules based on property values."""
    # Define a list of processed molecules with various properties
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07, "logp": -0.14}
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 58.12, "logp": 1.0}
        },
        {
            "smiles": "c1ccccc1",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 78.11, "logp": 2.0}
        }
    ]
    
    # Call sort_molecules with the molecules and a property name
    sorted_molecules = sort_molecules(molecules, "molecular_weight")
    
    # Assert that the sorted list is in ascending order by the specified property
    assert sorted_molecules[0]["properties"]["molecular_weight"] < sorted_molecules[1]["properties"]["molecular_weight"]
    assert sorted_molecules[1]["properties"]["molecular_weight"] < sorted_molecules[2]["properties"]["molecular_weight"]
    
    # Call sort_molecules with ascending=False
    sorted_molecules_desc = sort_molecules(molecules, "molecular_weight", ascending=False)
    
    # Assert that the sorted list is in descending order by the specified property
    assert sorted_molecules_desc[0]["properties"]["molecular_weight"] > sorted_molecules_desc[1]["properties"]["molecular_weight"]
    assert sorted_molecules_desc[1]["properties"]["molecular_weight"] > sorted_molecules_desc[2]["properties"]["molecular_weight"]


def test_deduplicate_molecules():
    """Tests removing duplicate molecules based on SMILES."""
    # Define a list of processed molecules with some duplicate SMILES
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07}
        },
        {
            "smiles": "CCO",  # Duplicate
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07, "logp": -0.14}  # More properties
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 58.12}
        },
        {
            "smiles": "c1ccccc1",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 78.11}
        }
    ]
    
    # Call deduplicate_molecules with the molecules list
    deduplicated_molecules = deduplicate_molecules(molecules)
    
    # Assert that the deduplicated list contains only unique SMILES
    assert len(deduplicated_molecules) == 3
    
    # Check that duplicate with more properties is preserved (CCO with logp)
    for molecule in deduplicated_molecules:
        if molecule["smiles"] == "CCO":
            assert "logp" in molecule["properties"]


def test_enrich_molecules():
    """Tests enriching molecules with additional properties."""
    # Define a list of processed molecules with some properties
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07}
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 58.12}
        }
    ]
    
    # Define a list of additional properties to calculate
    additional_properties = ["logp", "h_bond_donors"]
    
    # Call enrich_molecules with the molecules and additional properties
    enriched_molecules = enrich_molecules(molecules, additional_properties)
    
    # Assert that the enriched molecules contain the additional properties
    for molecule in enriched_molecules:
        assert "properties" in molecule
        for prop in additional_properties:
            assert prop in molecule["properties"]


def test_validate_and_standardize_molecules():
    """Tests validating and standardizing a list of molecules."""
    # Define a list of molecules with both valid and invalid SMILES
    molecules = [
        {"smiles": "CCO"},                 # Valid
        {"smiles": "C1CC"},                # Invalid
        {"smiles": "c1ccccc1"},            # Valid
        {"smiles": "invalid_smiles"},      # Invalid
        {"smiles": "CC(=O)O"}              # Valid
    ]
    
    # Call validate_and_standardize_molecules with the molecules list
    valid_molecules, invalid_molecules = validate_and_standardize_molecules(molecules)
    
    # Assert that the valid molecules list contains only valid molecules
    assert len(valid_molecules) == 3
    for molecule in valid_molecules:
        assert molecule["validation"]["structure_valid"] is True
    
    # Assert that the invalid molecules list contains only invalid molecules
    assert len(invalid_molecules) == 2
    for molecule in invalid_molecules:
        assert molecule["validation"]["structure_valid"] is False
    
    # Assert that the valid molecules have standardized SMILES
    for molecule in valid_molecules:
        assert "smiles" in molecule


def test_generate_processing_summary():
    """Tests generating a summary of molecule processing results."""
    # Define a list of processed molecules with various properties and validation results
    processed_molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07, "logp": -0.14}
        },
        {
            "smiles": "C1CC",
            "validation": {"structure_valid": False},
            "error": "Invalid structure"
        },
        {
            "smiles": "c1ccccc1",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 78.11, "logp": 2.0}
        },
        {
            "smiles": "invalid",
            "validation": {"structure_valid": False},
            "error": "Invalid SMILES"
        },
        {
            "smiles": "CC(=O)O",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 60.05, "logp": -0.17}
        }
    ]
    
    # Call generate_processing_summary with the processed molecules
    summary = generate_processing_summary(processed_molecules)
    
    # Assert that the summary contains the correct total count
    assert summary["total_molecules"] == 5
    
    # Assert that the summary contains the correct valid/invalid counts
    assert summary["valid_molecules"] == 3
    assert summary["invalid_molecules"] == 2
    
    # Assert that the summary contains property statistics
    assert "property_statistics" in summary
    assert "molecular_weight" in summary["property_statistics"]
    assert "logp" in summary["property_statistics"]
    
    # Check specific statistics
    mw_stats = summary["property_statistics"]["molecular_weight"]
    assert mw_stats["min"] == pytest.approx(46.07)
    assert mw_stats["max"] == pytest.approx(78.11)
    assert mw_stats["average"] == pytest.approx((46.07 + 78.11 + 60.05) / 3)
    assert mw_stats["count"] == 3


def test_molecule_processor_class():
    """Tests the MoleculeProcessor class functionality."""
    # Create a MoleculeProcessor instance
    processor = MoleculeProcessor()
    
    # Test the process method with a valid SMILES
    result = processor.process("CCO")
    assert result["smiles"] == "CCO"
    assert result["validation"]["structure_valid"] is True
    
    # Test the process_batch method with a list of SMILES
    smiles_list = ["CCO", "CC(=O)O", "c1ccccc1"]
    batch_results = processor.process_batch(smiles_list)
    assert len(batch_results) == len(smiles_list)
    
    # Test the process_stream method with a larger list
    stream_results = processor.process_stream(VALID_SMILES * 5)
    assert len(stream_results) == len(VALID_SMILES) * 5
    
    # Test the filter method with property filters
    molecules = [
        {
            "smiles": "CCO",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 46.07, "logp": -0.14}
        },
        {
            "smiles": "CCCC",
            "validation": {"structure_valid": True},
            "properties": {"molecular_weight": 58.12, "logp": 1.0}
        }
    ]
    property_filters = {"molecular_weight": {"min": 50}}
    filtered = processor.filter(molecules, property_filters)
    assert len(filtered) == 1
    assert filtered[0]["smiles"] == "CCCC"
    
    # Test the sort method with a property name
    sorted_molecules = processor.sort(molecules, "logp")
    assert sorted_molecules[0]["smiles"] == "CCO"
    assert sorted_molecules[1]["smiles"] == "CCCC"
    
    # Test the deduplicate method with duplicate molecules
    duplicate_molecules = [
        {"smiles": "CCO", "validation": {"structure_valid": True}},
        {"smiles": "CCO", "validation": {"structure_valid": True}}
    ]
    deduplicated = processor.deduplicate(duplicate_molecules)
    assert len(deduplicated) == 1
    
    # Test the enrich method with additional properties
    enriched = processor.enrich(molecules, ["h_bond_donors"])
    assert "h_bond_donors" in enriched[0]["properties"]
    
    # Test the validate_and_standardize method with mixed molecules
    mixed_molecules = [
        {"smiles": "CCO"},
        {"smiles": "invalid"}
    ]
    valid, invalid = processor.validate_and_standardize(mixed_molecules)
    assert len(valid) == 1
    assert len(invalid) == 1
    
    # Test the generate_summary method with processed molecules
    summary = processor.generate_summary(molecules)
    assert summary["total_molecules"] == 2
    assert summary["valid_molecules"] == 2
    
    # Test the clear_cache method
    processor.clear_cache()


def test_molecule_processor_caching():
    """Tests the caching functionality of MoleculeProcessor."""
    # Create a MoleculeProcessor instance
    processor = MoleculeProcessor()
    
    # Process the same molecule twice
    with patch('../../app.molecular.processor.process_molecule') as mock_process:
        # Mock the process_molecule function to return a valid result
        mock_process.return_value = {
            "smiles": "CCO",
            "validation": {"structure_valid": True}
        }
        
        # First call should use the real function
        processor.process("CCO")
        
        # Second call should use the cached result
        processor.process("CCO")
        
        # Verify that process_molecule was called only once
        assert mock_process.call_count == 1
        
        # Clear the cache and process again
        processor.clear_cache()
        processor.process("CCO")
        
        # Verify that process_molecule was called again after cache clearing
        assert mock_process.call_count == 2


def test_error_handling():
    """Tests error handling in the processor functions."""
    # Test with None SMILES and verify appropriate error handling
    result = process_molecule(None)
    assert "smiles" in result
    assert result["smiles"] is None
    assert "validation" in result
    assert result["validation"]["structure_valid"] is False
    
    # Test with empty SMILES and verify appropriate error handling
    result = process_molecule("")
    assert result["smiles"] == ""
    assert result["validation"]["structure_valid"] is False
    
    # Test with invalid property names
    result = process_molecule("CCO", ["invalid_property"])
    assert "properties" in result
    assert len(result["properties"]) == 0
    assert "property_error" in result
    
    # Test with invalid filter parameters
    molecules = [{"smiles": "CCO", "validation": {"structure_valid": True}}]
    result = filter_molecules(molecules, {"invalid_property": {"min": 0}})
    assert len(result) == 0  # No molecules match the filter