"""
Utility module for CSV file processing in the Molecular Data Management and CRO Integration Platform.

This module provides functions for reading, validating, mapping, and processing CSV files
containing molecular data, with support for batch processing of large datasets.
"""

import os  # standard library
import io  # standard library
import logging  # standard library
from typing import List, Dict, Optional, Any, Tuple, Union, Callable  # standard library

import pandas as pd  # pandas 2.0+
import numpy as np  # numpy 1.24+

from ..exceptions import ValidationException, CSVProcessingException
from ..constants import MAX_CSV_SIZE_MB, CSV_CHUNK_SIZE, DEFAULT_MOLECULE_PROPERTIES
from .validation_utils import validate_smiles
from .file_utils import validate_file_size, is_csv_file
from .molecular_utils import calculate_molecular_properties

# Set up logger
logger = logging.getLogger(__name__)

# Constants
MAX_CSV_SIZE_BYTES = MAX_CSV_SIZE_MB * 1024 * 1024
DEFAULT_PREVIEW_ROWS = 5
REQUIRED_MAPPING_FIELDS = ['smiles']


def read_csv_file(file_path: str, options: Dict[str, Any] = None) -> pd.DataFrame:
    """
    Read a CSV file into a pandas DataFrame.
    
    Args:
        file_path: Path to the CSV file
        options: Additional options for pandas.read_csv
        
    Returns:
        DataFrame containing CSV data
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    # Get file size and validate it's within limits
    file_size = os.path.getsize(file_path)
    if not validate_file_size(file_size, MAX_CSV_SIZE_BYTES):
        max_size_mb = MAX_CSV_SIZE_MB
        current_size_mb = file_size / (1024 * 1024)
        raise CSVProcessingException(
            f"File size ({current_size_mb:.2f} MB) exceeds maximum allowed size ({max_size_mb} MB)",
            {
                "file_path": file_path,
                "file_size_bytes": file_size,
                "file_size_mb": current_size_mb,
                "max_size_bytes": MAX_CSV_SIZE_BYTES,
                "max_size_mb": max_size_mb
            }
        )
    
    # Set default options if not provided
    if options is None:
        options = {}
    
    # Default options for CSV reading
    default_options = {
        'encoding': 'utf-8',
        'engine': 'python',
        'on_bad_lines': 'warn',  # For pandas 2.0+
        'low_memory': False,     # Allow better dtype detection
    }
    
    # Update with user-provided options
    for key, value in default_options.items():
        if key not in options:
            options[key] = value
    
    try:
        df = pd.read_csv(file_path, **options)
        return df
    except Exception as e:
        raise CSVProcessingException(
            f"Error reading CSV file: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def get_csv_headers(file_path: str) -> List[str]:
    """
    Extract headers from a CSV file.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        List of header column names
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    try:
        # Use pandas to read only the header row
        df = pd.read_csv(file_path, nrows=0)
        return list(df.columns)
    except Exception as e:
        raise CSVProcessingException(
            f"Error extracting CSV headers: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def validate_csv_structure(file_path: str) -> Dict[str, Any]:
    """
    Validate the structure and content of a CSV file.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        Validation results with headers, row_count, and is_valid
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    # Get file size and validate it's within limits
    file_size = os.path.getsize(file_path)
    if not validate_file_size(file_size, MAX_CSV_SIZE_BYTES):
        max_size_mb = MAX_CSV_SIZE_MB
        current_size_mb = file_size / (1024 * 1024)
        raise CSVProcessingException(
            f"File size ({current_size_mb:.2f} MB) exceeds maximum allowed size ({max_size_mb} MB)",
            {
                "file_path": file_path,
                "file_size_bytes": file_size,
                "file_size_mb": current_size_mb,
                "max_size_bytes": MAX_CSV_SIZE_BYTES,
                "max_size_mb": max_size_mb
            }
        )
    
    try:
        # Get headers
        headers = get_csv_headers(file_path)
        
        # Count the number of rows in the file
        row_count = 0
        with open(file_path, 'r', encoding='utf-8') as f:
            # Skip header row
            next(f)
            for _ in f:
                row_count += 1
        
        # Check if file has at least one data row
        is_valid = row_count > 0
        
        validation_results = {
            "headers": headers,
            "row_count": row_count,
            "is_valid": is_valid,
            "file_path": file_path,
            "file_size_bytes": file_size,
            "file_size_mb": file_size / (1024 * 1024)
        }
        
        if not is_valid:
            validation_results["error"] = "CSV file contains no data rows"
        
        return validation_results
    except Exception as e:
        raise CSVProcessingException(
            f"Error validating CSV structure: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def validate_mapping(headers: List[str], mapping: Dict[str, str]) -> Tuple[bool, List[str]]:
    """
    Validate column mapping configuration.
    
    Args:
        headers: List of header column names from the CSV
        mapping: Dictionary mapping CSV column names to system property names
        
    Returns:
        Tuple of (is_valid, error_messages)
    """
    error_messages = []
    
    # Check if mapping is empty
    if not mapping:
        error_messages.append("Column mapping cannot be empty")
        return False, error_messages
    
    # Verify that SMILES column is mapped
    smiles_mapped = False
    for system_property, csv_column in mapping.items():
        if system_property.lower() == 'smiles':
            smiles_mapped = True
            break
    
    if not smiles_mapped:
        error_messages.append("SMILES column must be mapped")
    
    # Check that all mapped CSV columns exist in headers
    for system_property, csv_column in mapping.items():
        if csv_column not in headers:
            error_messages.append(f"Mapped column '{csv_column}' does not exist in CSV headers")
    
    # Validate that system properties are valid identifiers
    for system_property in mapping.keys():
        if not system_property.isidentifier():
            error_messages.append(f"Invalid system property name: '{system_property}'")
    
    return len(error_messages) == 0, error_messages


def map_csv_columns(df: pd.DataFrame, mapping: Dict[str, str]) -> pd.DataFrame:
    """
    Map CSV columns to system properties based on mapping configuration.
    
    Args:
        df: DataFrame containing CSV data
        mapping: Dictionary mapping CSV column names to system property names
        
    Returns:
        DataFrame with renamed columns according to mapping
    """
    # Create a copy of the DataFrame to avoid modifying the original
    mapped_df = df.copy()
    
    # Create a column renaming dictionary
    rename_dict = {}
    for system_property, csv_column in mapping.items():
        if csv_column in df.columns:
            rename_dict[csv_column] = system_property
    
    # Rename columns
    mapped_df = mapped_df.rename(columns=rename_dict)
    
    # Select only the mapped columns
    mapped_df = mapped_df[list(rename_dict.values())]
    
    return mapped_df


def get_csv_preview(file_path: str, num_rows: int = DEFAULT_PREVIEW_ROWS) -> Dict[str, Any]:
    """
    Get a preview of CSV data with a limited number of rows.
    
    Args:
        file_path: Path to the CSV file
        num_rows: Number of rows to include in preview
        
    Returns:
        Preview data with headers and rows
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    try:
        # Read a limited number of rows
        df = pd.read_csv(file_path, nrows=num_rows)
        
        # Extract headers
        headers = list(df.columns)
        
        # Convert DataFrame rows to list of dictionaries
        rows = df.replace({np.nan: None}).to_dict('records')
        
        return {
            "headers": headers,
            "rows": rows,
            "file_path": file_path,
            "preview_rows": len(rows),
            "total_rows": count_csv_rows(file_path)
        }
    except Exception as e:
        raise CSVProcessingException(
            f"Error generating CSV preview: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def process_csv_in_batches(
    file_path: str, 
    batch_processor: Callable[[pd.DataFrame, int, int], Dict[str, Any]],
    batch_size: int = CSV_CHUNK_SIZE,
    options: Dict[str, Any] = None
) -> Dict[str, Any]:
    """
    Process a CSV file in batches to handle large datasets efficiently.
    
    Args:
        file_path: Path to the CSV file
        batch_processor: Function to process each batch, accepting DataFrame, batch index, and total batches
        batch_size: Number of rows in each batch
        options: Additional options for pandas.read_csv
        
    Returns:
        Processing results with statistics
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    # Set default options if not provided
    if options is None:
        options = {}
    
    # Default options for CSV reading
    default_options = {
        'encoding': 'utf-8',
        'engine': 'python',
        'on_bad_lines': 'warn',  # For pandas 2.0+
        'low_memory': False,     # Better dtype detection
        'chunksize': batch_size  # Read in chunks
    }
    
    # Update with user-provided options
    for key, value in default_options.items():
        if key not in options:
            options[key] = value
    
    # Ensure chunksize is set correctly
    options['chunksize'] = batch_size
    
    try:
        # Initialize counters
        total_rows = 0
        processed_rows = 0
        successful_rows = 0
        failed_rows = 0
        batch_results = []
        
        # Get approximate total row count for progress reporting
        total_row_count = count_csv_rows(file_path)
        
        # Process file in chunks
        reader = pd.read_csv(file_path, **options)
        for i, chunk in enumerate(reader):
            logger.info(f"Processing batch {i+1} with {len(chunk)} rows")
            
            # Process the batch
            try:
                batch_result = batch_processor(chunk, i, total_row_count)
                
                # Update counters
                total_rows += len(chunk)
                processed_rows += batch_result.get('processed_rows', len(chunk))
                successful_rows += batch_result.get('successful_rows', 0)
                failed_rows += batch_result.get('failed_rows', 0)
                
                # Append batch result
                batch_results.append(batch_result)
                
                logger.info(f"Batch {i+1} processed successfully: "
                           f"{batch_result.get('successful_rows', 0)} successful, "
                           f"{batch_result.get('failed_rows', 0)} failed")
            except Exception as e:
                logger.error(f"Error processing batch {i+1}: {str(e)}")
                failed_rows += len(chunk)
                batch_results.append({
                    'batch_index': i,
                    'processed_rows': len(chunk),
                    'successful_rows': 0,
                    'failed_rows': len(chunk),
                    'error': str(e)
                })
        
        # Compile final results
        success_rate = successful_rows / total_rows if total_rows > 0 else 0
        
        return {
            'file_path': file_path,
            'total_rows': total_rows,
            'processed_rows': processed_rows,
            'successful_rows': successful_rows,
            'failed_rows': failed_rows,
            'success_rate': success_rate,
            'batch_count': len(batch_results),
            'batch_results': batch_results
        }
    except Exception as e:
        raise CSVProcessingException(
            f"Error processing CSV in batches: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def extract_molecules_from_dataframe(
    df: pd.DataFrame, 
    required_columns: List[str] = ['smiles']
) -> Dict[str, Any]:
    """
    Extract molecule data from DataFrame with validation.
    
    Args:
        df: DataFrame containing molecular data
        required_columns: List of required column names
        
    Returns:
        Extraction results with molecules and errors
        
    Raises:
        ValidationException: If required columns are missing
    """
    # Validate that all required columns exist in the DataFrame
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValidationException(
            f"Missing required columns in DataFrame: {', '.join(missing_columns)}",
            {"missing_columns": missing_columns, "required_columns": required_columns}
        )
    
    # Initialize lists for valid molecules and errors
    valid_molecules = []
    errors = []
    
    # Process each row
    for index, row in df.iterrows():
        try:
            # Extract SMILES string and validate
            smiles = row['smiles']
            if validate_smiles(smiles):
                # Extract other properties
                properties = {}
                for col in df.columns:
                    if col != 'smiles':
                        properties[col] = row[col]
                
                # Calculate additional molecular properties
                try:
                    calculated_properties = calculate_molecular_properties(smiles)
                    # Merge calculated properties with extracted properties
                    for prop, value in calculated_properties.items():
                        if prop not in properties:
                            properties[prop] = value
                except Exception as e:
                    logger.warning(f"Error calculating properties for {smiles}: {str(e)}")
                
                # Add molecule to valid list
                molecule = {
                    'smiles': smiles,
                    'properties': properties,
                    'row_index': index
                }
                valid_molecules.append(molecule)
            else:
                errors.append({
                    'row_index': index,
                    'smiles': smiles,
                    'error': 'Invalid SMILES structure'
                })
        except Exception as e:
            errors.append({
                'row_index': index,
                'error': str(e)
            })
    
    return {
        'valid_molecules': valid_molecules,
        'errors': errors,
        'total_rows': len(df),
        'valid_count': len(valid_molecules),
        'error_count': len(errors)
    }


def validate_csv_file_with_exception(file_path: str) -> Dict[str, Any]:
    """
    Validate a CSV file and raise an exception if invalid.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        Validation results if valid
        
    Raises:
        CSVProcessingException: If file is invalid
    """
    try:
        validation_results = validate_csv_structure(file_path)
        if not validation_results['is_valid']:
            raise CSVProcessingException(
                "CSV validation failed",
                validation_results
            )
        return validation_results
    except Exception as e:
        if isinstance(e, CSVProcessingException):
            raise
        else:
            raise CSVProcessingException(
                f"Error validating CSV file: {str(e)}",
                {"file_path": file_path, "error": str(e)}
            )


def count_csv_rows(file_path: str) -> int:
    """
    Count the number of rows in a CSV file.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        Number of data rows (excluding header)
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    # Validate file is a CSV
    if not is_csv_file(file_path):
        raise CSVProcessingException(
            f"File is not a CSV file: {file_path}",
            {"file_path": file_path}
        )
    
    try:
        # Use a more efficient approach for large files
        row_count = 0
        with open(file_path, 'r', encoding='utf-8') as f:
            # Skip header row
            next(f)
            for _ in f:
                row_count += 1
        return row_count
    except Exception as e:
        raise CSVProcessingException(
            f"Error counting CSV rows: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def detect_csv_delimiter(file_path: str) -> str:
    """
    Detect the delimiter used in a CSV file.
    
    Args:
        file_path: Path to the CSV file
        
    Returns:
        Detected delimiter character
        
    Raises:
        CSVProcessingException: If file doesn't exist, is not readable, or is not a valid CSV
    """
    # Check if file exists and is readable
    if not os.path.exists(file_path):
        raise CSVProcessingException(
            f"File does not exist: {file_path}",
            {"file_path": file_path}
        )
    
    if not os.access(file_path, os.R_OK):
        raise CSVProcessingException(
            f"File is not readable: {file_path}",
            {"file_path": file_path}
        )
    
    try:
        # Read a sample of the file
        with open(file_path, 'r', encoding='utf-8') as f:
            sample = f.readline()
        
        # Common delimiters to check
        delimiters = [',', ';', '\t', '|']
        delimiter_counts = {}
        
        # Count occurrences of each delimiter
        for delimiter in delimiters:
            delimiter_counts[delimiter] = sample.count(delimiter)
        
        # Find delimiter with the highest count
        detected_delimiter = max(delimiter_counts, key=delimiter_counts.get)
        
        # If no delimiter found or count is zero, default to comma
        if delimiter_counts[detected_delimiter] == 0:
            return ','
        
        return detected_delimiter
    except Exception as e:
        raise CSVProcessingException(
            f"Error detecting CSV delimiter: {str(e)}",
            {"file_path": file_path, "error": str(e)}
        )


def validate_smiles_column(df: pd.DataFrame, smiles_column: str = 'smiles') -> Dict[str, Any]:
    """
    Validate that a DataFrame contains a valid SMILES column.
    
    Args:
        df: DataFrame to validate
        smiles_column: Name of the SMILES column
        
    Returns:
        Validation results with valid_count, invalid_count, and invalid_rows
        
    Raises:
        ValidationException: If SMILES column doesn't exist
    """
    # Check if the specified SMILES column exists
    if smiles_column not in df.columns:
        raise ValidationException(
            f"SMILES column '{smiles_column}' not found in DataFrame",
            {"smiles_column": smiles_column, "available_columns": list(df.columns)}
        )
    
    # Initialize counters
    valid_count = 0
    invalid_count = 0
    invalid_rows = []
    
    # Validate each SMILES string
    for index, row in df.iterrows():
        smiles = row[smiles_column]
        if validate_smiles(smiles):
            valid_count += 1
        else:
            invalid_count += 1
            invalid_rows.append({
                'row_index': index,
                'smiles': smiles
            })
    
    return {
        'valid_count': valid_count,
        'invalid_count': invalid_count,
        'total_count': len(df),
        'invalid_rows': invalid_rows,
        'is_valid': invalid_count == 0
    }


def generate_csv_import_summary(processing_results: Dict[str, Any]) -> Dict[str, Any]:
    """
    Generate a summary of CSV import results.
    
    Args:
        processing_results: Results from CSV processing
        
    Returns:
        Formatted summary with statistics and error categories
    """
    # Extract key metrics
    total_rows = processing_results.get('total_rows', 0)
    processed_rows = processing_results.get('processed_rows', 0)
    successful_rows = processing_results.get('successful_rows', 0)
    failed_rows = processing_results.get('failed_rows', 0)
    
    # Calculate percentages
    success_rate = (successful_rows / total_rows) * 100 if total_rows > 0 else 0
    failure_rate = (failed_rows / total_rows) * 100 if total_rows > 0 else 0
    
    # Categorize errors if available
    error_categories = {}
    if 'batch_results' in processing_results:
        for batch in processing_results['batch_results']:
            if 'errors' in batch:
                for error in batch['errors']:
                    error_type = error.get('error', 'Unknown error')
                    if error_type in error_categories:
                        error_categories[error_type] += 1
                    else:
                        error_categories[error_type] = 1
    
    # Create summary
    summary = {
        'total_rows': total_rows,
        'processed_rows': processed_rows,
        'successful_rows': successful_rows,
        'failed_rows': failed_rows,
        'success_rate': success_rate,
        'failure_rate': failure_rate,
        'error_categories': error_categories,
        'status': 'complete',
        'timestamp': pd.Timestamp.now().isoformat()
    }
    
    return summary