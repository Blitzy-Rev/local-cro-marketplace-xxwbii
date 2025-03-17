"""
Unit tests for the CSV utility functions in the Molecular Data Management and CRO Integration Platform.

This module tests the functionality of CSV file processing, validation, mapping, and data extraction 
utilities used for molecular data ingestion.
"""

import pytest
import os
import io
import tempfile
import pandas as pd
import numpy as np
from unittest.mock import patch, MagicMock

from src.backend.app.exceptions import ValidationException, CSVProcessingException
from src.backend.app.constants import MAX_CSV_SIZE_MB, CSV_CHUNK_SIZE, DEFAULT_MOLECULE_PROPERTIES
from src.backend.app.utils import csv_utils


# Helper functions to create test CSV files
def create_test_csv_file(content, filename=None):
    """Create a temporary CSV file with the given content."""
    if filename is None:
        fd, filename = tempfile.mkstemp(suffix='.csv')
        os.close(fd)
    
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
    
    return filename


def create_valid_csv_file():
    """Create a valid CSV file with test data."""
    content = "SMILES,MW,LogP,Activity\nCCO,46.07,-0.14,78.5\nCCCCO,74.12,0.88,65.2\nc1ccccc1,78.11,1.90,45.0"
    return create_test_csv_file(content)


def create_empty_csv_file():
    """Create an empty CSV file."""
    return create_test_csv_file("")


def create_headers_only_csv_file():
    """Create a CSV file with only headers."""
    content = "SMILES,MW,LogP,Activity"
    return create_test_csv_file(content)


def test_read_csv_file():
    """Test the read_csv_file function with various CSV files."""
    # Create a temporary CSV file with test data
    temp_csv = create_valid_csv_file()
    non_csv_file = None
    
    try:
        # Test with a valid CSV file
        df = csv_utils.read_csv_file(temp_csv)
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 3
        assert list(df.columns) == ['SMILES', 'MW', 'LogP', 'Activity']
        
        # Test with custom pandas options
        df = csv_utils.read_csv_file(temp_csv, {'usecols': ['SMILES', 'MW']})
        assert list(df.columns) == ['SMILES', 'MW']
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.read_csv_file('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.read_csv_file(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
        
        # Test with file exceeding size limit
        with patch('src.backend.app.utils.csv_utils.validate_file_size', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.read_csv_file(temp_csv)
            assert "File size" in str(excinfo.value)
            assert "exceeds maximum allowed size" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(temp_csv):
                os.remove(temp_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_get_csv_headers():
    """Test the get_csv_headers function with various CSV files."""
    # Create a temporary CSV file with test headers
    temp_csv = create_valid_csv_file()
    non_csv_file = None
    
    try:
        # Test with a valid CSV file
        headers = csv_utils.get_csv_headers(temp_csv)
        assert headers == ['SMILES', 'MW', 'LogP', 'Activity']
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.get_csv_headers('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.get_csv_headers(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(temp_csv):
                os.remove(temp_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_validate_csv_structure():
    """Test the validate_csv_structure function with various CSV files."""
    # Create temporary CSV files with different structures
    valid_csv = create_valid_csv_file()
    headers_only_csv = create_headers_only_csv_file()
    empty_csv = create_empty_csv_file()
    non_csv_file = None
    
    try:
        # Test with valid CSV file containing headers and data
        result = csv_utils.validate_csv_structure(valid_csv)
        assert result['is_valid'] is True
        assert result['headers'] == ['SMILES', 'MW', 'LogP', 'Activity']
        assert result['row_count'] == 3
        
        # Test with CSV file containing only headers (no data rows)
        result = csv_utils.validate_csv_structure(headers_only_csv)
        assert result['is_valid'] is False
        assert result['headers'] == ['SMILES', 'MW', 'LogP', 'Activity']
        assert result['row_count'] == 0
        
        # Test with empty CSV file
        with pytest.raises(CSVProcessingException):
            csv_utils.validate_csv_structure(empty_csv)
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.validate_csv_structure('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.validate_csv_structure(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
        
        # Test with file exceeding size limit
        with patch('src.backend.app.utils.csv_utils.validate_file_size', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.validate_csv_structure(valid_csv)
            assert "File size" in str(excinfo.value)
            assert "exceeds maximum allowed size" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(valid_csv):
                os.remove(valid_csv)
        except:
            pass
        
        try:
            if os.path.exists(headers_only_csv):
                os.remove(headers_only_csv)
        except:
            pass
        
        try:
            if os.path.exists(empty_csv):
                os.remove(empty_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_validate_mapping():
    """Test the validate_mapping function with various mapping configurations."""
    # Test with valid mapping containing SMILES column
    headers = ['Column1', 'Column2', 'Column3']
    valid_mapping = {'smiles': 'Column1', 'mw': 'Column2', 'logp': 'Column3'}
    is_valid, error_messages = csv_utils.validate_mapping(headers, valid_mapping)
    assert is_valid is True
    assert len(error_messages) == 0
    
    # Test with mapping missing SMILES column
    invalid_mapping = {'mw': 'Column2', 'logp': 'Column3'}
    is_valid, error_messages = csv_utils.validate_mapping(headers, invalid_mapping)
    assert is_valid is False
    assert "SMILES column must be mapped" in error_messages
    
    # Test with mapping referencing non-existent CSV headers
    invalid_mapping = {'smiles': 'Column1', 'mw': 'NonExistentColumn'}
    is_valid, error_messages = csv_utils.validate_mapping(headers, invalid_mapping)
    assert is_valid is False
    assert "Mapped column 'NonExistentColumn' does not exist in CSV headers" in error_messages
    
    # Test with empty mapping
    is_valid, error_messages = csv_utils.validate_mapping(headers, {})
    assert is_valid is False
    assert "Column mapping cannot be empty" in error_messages
    
    # Test with None mapping
    is_valid, error_messages = csv_utils.validate_mapping(headers, None)
    assert is_valid is False
    assert "Column mapping cannot be empty" in error_messages


def test_map_csv_columns():
    """Test the map_csv_columns function with various mapping configurations."""
    # Create test DataFrame with sample columns
    df = pd.DataFrame({
        'Column1': ['CCO', 'CCCCO', 'c1ccccc1'],
        'Column2': [46.07, 74.12, 78.11],
        'Column3': [-0.14, 0.88, 1.90]
    })
    
    # Test with valid mapping for all columns
    valid_mapping = {'smiles': 'Column1', 'mw': 'Column2', 'logp': 'Column3'}
    mapped_df = csv_utils.map_csv_columns(df, valid_mapping)
    assert list(mapped_df.columns) == ['smiles', 'mw', 'logp']
    assert mapped_df.shape[0] == 3
    
    # Test with partial mapping (some columns mapped, others unchanged)
    partial_mapping = {'smiles': 'Column1', 'mw': 'Column2'}
    mapped_df = csv_utils.map_csv_columns(df, partial_mapping)
    assert list(mapped_df.columns) == ['smiles', 'mw']
    assert mapped_df.shape[0] == 3
    
    # Test with empty mapping
    mapped_df = csv_utils.map_csv_columns(df, {})
    assert mapped_df.empty
    
    # Test with None mapping
    mapped_df = csv_utils.map_csv_columns(df, None)
    assert mapped_df.empty


def test_get_csv_preview():
    """Test the get_csv_preview function with various CSV files."""
    # Create a temporary CSV file with test data
    temp_csv = create_valid_csv_file()
    non_csv_file = None
    
    try:
        # Test with default preview rows
        preview = csv_utils.get_csv_preview(temp_csv)
        assert preview['headers'] == ['SMILES', 'MW', 'LogP', 'Activity']
        assert len(preview['rows']) == min(csv_utils.DEFAULT_PREVIEW_ROWS, 3)
        
        # Test with custom number of preview rows
        preview = csv_utils.get_csv_preview(temp_csv, num_rows=1)
        assert len(preview['rows']) == 1
        
        # Test with CSV file having fewer rows than requested
        preview = csv_utils.get_csv_preview(temp_csv, num_rows=10)
        assert len(preview['rows']) == 3  # Only 3 rows in test file
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.get_csv_preview('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.get_csv_preview(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(temp_csv):
                os.remove(temp_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_process_csv_in_batches():
    """Test the process_csv_in_batches function with batch processing."""
    # Create a temporary CSV file with multiple rows of test data
    content = "SMILES,MW,LogP\n" + "\n".join([f"C{'C'*i}O,{46.07+i*14},{-0.14+i*0.5}" for i in range(20)])
    temp_csv = create_test_csv_file(content)
    empty_csv = create_headers_only_csv_file()
    non_csv_file = None
    
    try:
        # Create mock batch processor function
        mock_processor = MagicMock()
        mock_processor.return_value = {'processed_rows': 5, 'successful_rows': 4, 'failed_rows': 1}
        
        # Test processing with default batch size
        result = csv_utils.process_csv_in_batches(temp_csv, mock_processor)
        assert result['file_path'] == temp_csv
        assert result['total_rows'] == 20
        assert mock_processor.call_count > 0
        
        # Reset mock
        mock_processor.reset_mock()
        
        # Test processing with custom batch size
        result = csv_utils.process_csv_in_batches(temp_csv, mock_processor, batch_size=5)
        assert result['file_path'] == temp_csv
        assert result['total_rows'] == 20
        assert mock_processor.call_count == 4  # 20 rows / 5 batch size = 4 calls
        
        # Test with empty CSV file
        mock_processor.reset_mock()
        result = csv_utils.process_csv_in_batches(empty_csv, mock_processor)
        assert result['total_rows'] == 0
        assert mock_processor.call_count == 0
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.process_csv_in_batches('nonexistent.csv', mock_processor)
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.process_csv_in_batches(non_csv_file, mock_processor)
            assert "File is not a CSV file" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(temp_csv):
                os.remove(temp_csv)
        except:
            pass
        
        try:
            if os.path.exists(empty_csv):
                os.remove(empty_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_extract_molecules_from_dataframe():
    """Test the extract_molecules_from_dataframe function with molecular data."""
    # Create test DataFrame with valid and invalid SMILES data
    df = pd.DataFrame({
        'smiles': ['CCO', 'CCCCO', 'INVALID', 'c1ccccc1'],
        'mw': [46.07, 74.12, None, 78.11],
        'logp': [-0.14, 0.88, None, 1.90]
    })
    
    # Mock validate_smiles and calculate_molecular_properties functions
    with patch('src.backend.app.utils.validation_utils.validate_smiles', side_effect=lambda s: s != 'INVALID'):
        with patch('src.backend.app.utils.molecular_utils.calculate_molecular_properties', return_value={'calculated_prop': 100}):
            
            # Test with DataFrame containing valid SMILES
            valid_df = df[df['smiles'] != 'INVALID'].reset_index(drop=True)
            result = csv_utils.extract_molecules_from_dataframe(valid_df)
            assert result['valid_count'] == 3
            assert result['error_count'] == 0
            assert len(result['valid_molecules']) == 3
            
            # Test with DataFrame containing mix of valid and invalid SMILES
            result = csv_utils.extract_molecules_from_dataframe(df)
            assert result['valid_count'] == 3
            assert result['error_count'] == 1
            assert len(result['valid_molecules']) == 3
            assert len(result['errors']) == 1
            assert result['errors'][0]['smiles'] == 'INVALID'
            
            # Test with DataFrame missing required columns
            invalid_df = pd.DataFrame({'not_smiles': ['CCO']})
            with pytest.raises(ValidationException) as excinfo:
                csv_utils.extract_molecules_from_dataframe(invalid_df)
            assert "Missing required columns" in str(excinfo.value)
            
            # Test with custom required columns
            result = csv_utils.extract_molecules_from_dataframe(df, required_columns=['smiles', 'mw'])
            assert result['valid_count'] == 3
            assert result['error_count'] == 1


def test_validate_csv_file_with_exception():
    """Test the validate_csv_file_with_exception function with various CSV files."""
    # Create temporary CSV files with different structures
    valid_csv = create_valid_csv_file()
    headers_only_csv = create_headers_only_csv_file()
    non_csv_file = None
    
    try:
        # Test with valid CSV file
        result = csv_utils.validate_csv_file_with_exception(valid_csv)
        assert result['is_valid'] is True
        
        # Test with invalid CSV file and verify CSVProcessingException is raised
        with patch('src.backend.app.utils.csv_utils.validate_csv_structure', 
                   return_value={'is_valid': False, 'headers': [], 'row_count': 0}):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.validate_csv_file_with_exception(headers_only_csv)
            assert "CSV validation failed" in str(excinfo.value)
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.validate_csv_file_with_exception('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.validate_csv_file_with_exception(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
        
        # Test with file exceeding size limit
        with patch('src.backend.app.utils.csv_utils.validate_file_size', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.validate_csv_file_with_exception(valid_csv)
            assert "File size" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(valid_csv):
                os.remove(valid_csv)
        except:
            pass
        
        try:
            if os.path.exists(headers_only_csv):
                os.remove(headers_only_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_count_csv_rows():
    """Test the count_csv_rows function with various CSV files."""
    # Create temporary CSV files with known number of rows
    content = "SMILES,MW,LogP\n" + "\n".join([f"C{'C'*i}O,{46.07+i*14},{-0.14+i*0.5}" for i in range(10)])
    multi_row_csv = create_test_csv_file(content)
    header_only_csv = create_headers_only_csv_file()
    empty_csv = create_empty_csv_file()
    non_csv_file = None
    
    try:
        # Test with CSV file containing multiple rows
        count = csv_utils.count_csv_rows(multi_row_csv)
        assert count == 10
        
        # Test with CSV file containing only header row
        count = csv_utils.count_csv_rows(header_only_csv)
        assert count == 0
        
        # Test with empty CSV file
        with pytest.raises(CSVProcessingException):
            csv_utils.count_csv_rows(empty_csv)
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.count_csv_rows('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.count_csv_rows(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(multi_row_csv):
                os.remove(multi_row_csv)
        except:
            pass
        
        try:
            if os.path.exists(header_only_csv):
                os.remove(header_only_csv)
        except:
            pass
        
        try:
            if os.path.exists(empty_csv):
                os.remove(empty_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_detect_csv_delimiter():
    """Test the detect_csv_delimiter function with various CSV files."""
    # Create temporary CSV files with different delimiters
    comma_csv = create_test_csv_file("a,b,c\n1,2,3")
    tab_csv = create_test_csv_file("a\tb\tc\n1\t2\t3")
    semicolon_csv = create_test_csv_file("a;b;c\n1;2;3")
    non_csv_file = None
    
    try:
        # Test with comma-delimited CSV file
        delimiter = csv_utils.detect_csv_delimiter(comma_csv)
        assert delimiter == ','
        
        # Test with tab-delimited CSV file
        delimiter = csv_utils.detect_csv_delimiter(tab_csv)
        assert delimiter == '\t'
        
        # Test with semicolon-delimited CSV file
        delimiter = csv_utils.detect_csv_delimiter(semicolon_csv)
        assert delimiter == ';'
        
        # Test with non-existent file
        with pytest.raises(CSVProcessingException) as excinfo:
            csv_utils.detect_csv_delimiter('nonexistent.csv')
        assert "File does not exist" in str(excinfo.value)
        
        # Test with non-CSV file
        non_csv_file = tempfile.mkstemp(suffix='.txt')[1]
        with open(non_csv_file, 'w') as f:
            f.write("Not a CSV file")
        
        with patch('src.backend.app.utils.csv_utils.is_csv_file', return_value=False):
            with pytest.raises(CSVProcessingException) as excinfo:
                csv_utils.detect_csv_delimiter(non_csv_file)
            assert "File is not a CSV file" in str(excinfo.value)
    
    finally:
        # Clean up
        try:
            if os.path.exists(comma_csv):
                os.remove(comma_csv)
        except:
            pass
        
        try:
            if os.path.exists(tab_csv):
                os.remove(tab_csv)
        except:
            pass
        
        try:
            if os.path.exists(semicolon_csv):
                os.remove(semicolon_csv)
        except:
            pass
        
        try:
            if non_csv_file and os.path.exists(non_csv_file):
                os.remove(non_csv_file)
        except:
            pass


def test_validate_smiles_column():
    """Test the validate_smiles_column function with various SMILES data."""
    # Create test DataFrame with valid and invalid SMILES data
    df = pd.DataFrame({
        'smiles': ['CCO', 'CCCCO', 'INVALID', 'c1ccccc1'],
        'mw': [46.07, 74.12, None, 78.11],
        'logp': [-0.14, 0.88, None, 1.90]
    })
    
    # Mock validate_smiles function
    with patch('src.backend.app.utils.validation_utils.validate_smiles', side_effect=lambda s: s != 'INVALID'):
        
        # Test with DataFrame containing all valid SMILES
        valid_df = df[df['smiles'] != 'INVALID'].reset_index(drop=True)
        result = csv_utils.validate_smiles_column(valid_df)
        assert result['is_valid'] is True
        assert result['valid_count'] == 3
        assert result['invalid_count'] == 0
        
        # Test with DataFrame containing mix of valid and invalid SMILES
        result = csv_utils.validate_smiles_column(df)
        assert result['is_valid'] is False
        assert result['valid_count'] == 3
        assert result['invalid_count'] == 1
        assert len(result['invalid_rows']) == 1
        assert result['invalid_rows'][0]['smiles'] == 'INVALID'
        
        # Test with DataFrame missing SMILES column
        invalid_df = pd.DataFrame({'not_smiles': ['CCO']})
        with pytest.raises(ValidationException) as excinfo:
            csv_utils.validate_smiles_column(invalid_df)
        assert "SMILES column 'smiles' not found" in str(excinfo.value)
        
        # Test with custom SMILES column name
        custom_df = pd.DataFrame({'custom_smiles': ['CCO', 'CCCCO', 'c1ccccc1']})
        result = csv_utils.validate_smiles_column(custom_df, smiles_column='custom_smiles')
        assert result['is_valid'] is True
        assert result['valid_count'] == 3


def test_generate_csv_import_summary():
    """Test the generate_csv_import_summary function with various processing results."""
    # Create sample processing results with different statistics
    successful_results = {
        'total_rows': 100,
        'processed_rows': 100,
        'successful_rows': 95,
        'failed_rows': 5,
        'batch_results': [
            {
                'batch_index': 0,
                'processed_rows': 50,
                'successful_rows': 48,
                'failed_rows': 2,
                'errors': [
                    {'row_index': 10, 'error': 'Invalid SMILES'},
                    {'row_index': 20, 'error': 'Invalid SMILES'}
                ]
            },
            {
                'batch_index': 1,
                'processed_rows': 50,
                'successful_rows': 47,
                'failed_rows': 3,
                'errors': [
                    {'row_index': 60, 'error': 'Invalid SMILES'},
                    {'row_index': 70, 'error': 'Missing property'},
                    {'row_index': 80, 'error': 'Invalid SMILES'}
                ]
            }
        ]
    }
    
    partial_results = {
        'total_rows': 100,
        'processed_rows': 100,
        'successful_rows': 70,
        'failed_rows': 30,
        'batch_results': [
            {
                'batch_index': 0,
                'processed_rows': 50,
                'successful_rows': 30,
                'failed_rows': 20,
                'errors': [{'row_index': i, 'error': 'Invalid SMILES'} for i in range(20)]
            },
            {
                'batch_index': 1,
                'processed_rows': 50,
                'successful_rows': 40,
                'failed_rows': 10,
                'errors': [{'row_index': 50+i, 'error': 'Missing property'} for i in range(10)]
            }
        ]
    }
    
    failed_results = {
        'total_rows': 100,
        'processed_rows': 100,
        'successful_rows': 10,
        'failed_rows': 90,
        'batch_results': [
            {
                'batch_index': 0,
                'processed_rows': 50,
                'successful_rows': 5,
                'failed_rows': 45,
                'errors': [{'row_index': i, 'error': 'Invalid SMILES'} for i in range(45)]
            },
            {
                'batch_index': 1,
                'processed_rows': 50,
                'successful_rows': 5,
                'failed_rows': 45,
                'errors': [{'row_index': 50+i, 'error': 'Missing property'} for i in range(45)]
            }
        ]
    }
    
    empty_results = {
        'total_rows': 0,
        'processed_rows': 0,
        'successful_rows': 0,
        'failed_rows': 0,
        'batch_results': []
    }
    
    # Test with successful processing results (high success rate)
    summary = csv_utils.generate_csv_import_summary(successful_results)
    assert summary['total_rows'] == 100
    assert summary['successful_rows'] == 95
    assert summary['failed_rows'] == 5
    assert summary['success_rate'] == 95.0
    assert 'Invalid SMILES' in summary['error_categories']
    assert 'Missing property' in summary['error_categories']
    
    # Test with partially successful processing results (some failures)
    summary = csv_utils.generate_csv_import_summary(partial_results)
    assert summary['total_rows'] == 100
    assert summary['successful_rows'] == 70
    assert summary['failed_rows'] == 30
    assert summary['success_rate'] == 70.0
    assert 'Invalid SMILES' in summary['error_categories']
    assert 'Missing property' in summary['error_categories']
    
    # Test with failed processing results (high failure rate)
    summary = csv_utils.generate_csv_import_summary(failed_results)
    assert summary['total_rows'] == 100
    assert summary['successful_rows'] == 10
    assert summary['failed_rows'] == 90
    assert summary['success_rate'] == 10.0
    assert 'Invalid SMILES' in summary['error_categories']
    assert 'Missing property' in summary['error_categories']
    
    # Test with empty processing results
    summary = csv_utils.generate_csv_import_summary(empty_results)
    assert summary['total_rows'] == 0
    assert summary['successful_rows'] == 0
    assert summary['failed_rows'] == 0
    assert summary['success_rate'] == 0.0
    assert len(summary['error_categories']) == 0