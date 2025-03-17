"""
Utility functions for file operations in the Molecular Data Management and CRO Integration Platform.

This module provides helper functions for file validation, manipulation, and management
used throughout the application, particularly for CSV uploads, result files, and molecular images.
"""

import os  # standard library
import pathlib  # standard library
import uuid  # standard library
import datetime  # standard library
import mimetypes  # standard library
from typing import Optional, Set, Tuple  # standard library

import magic  # python-magic 0.4.27

from ..logging_config import logger
from ..exceptions import ValidationException, FileStorageException
from ..constants import MAX_FILE_SIZE_MB, ALLOWED_EXTENSIONS

# Calculate maximum file size in bytes
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024

# Convert allowed extensions to lowercase set with dots
ALLOWED_FILE_EXTENSIONS = {f'.{ext.lower()}' for ext in ALLOWED_EXTENSIONS}

# Define file type extension sets
IMAGE_EXTENSIONS = {'.png', '.jpg', '.jpeg', '.svg'}
DOCUMENT_EXTENSIONS = {'.pdf', '.docx', '.xlsx'}
DATA_EXTENSIONS = {'.csv', '.xlsx'}


def get_file_extension(filename: str) -> str:
    """
    Extract the file extension from a filename.
    
    Args:
        filename: The filename to extract extension from
        
    Returns:
        File extension including the dot (e.g., '.csv')
    """
    filename = filename.lower()  # Convert to lowercase for case-insensitive comparison
    return os.path.splitext(filename)[1]


def is_allowed_file(filename: str, allowed_extensions: Optional[Set[str]] = None) -> bool:
    """
    Check if a file has an allowed extension.
    
    Args:
        filename: The filename to check
        allowed_extensions: Set of allowed extensions (if None, uses ALLOWED_FILE_EXTENSIONS)
        
    Returns:
        True if file extension is allowed, False otherwise
    """
    if allowed_extensions is None:
        allowed_extensions = ALLOWED_FILE_EXTENSIONS
        
    extension = get_file_extension(filename)
    return extension.lower() in allowed_extensions


def validate_file_size(file_size: int, max_size: Optional[int] = None) -> bool:
    """
    Validate that a file size is within the allowed limit.
    
    Args:
        file_size: Size of the file in bytes
        max_size: Maximum allowed size in bytes (if None, uses MAX_FILE_SIZE_BYTES)
        
    Returns:
        True if file size is within limit, False otherwise
    """
    if max_size is None:
        max_size = MAX_FILE_SIZE_BYTES
        
    return file_size <= max_size


def get_mime_type(file_path: str) -> str:
    """
    Determine the MIME type of a file.
    
    Args:
        file_path: Path to the file
        
    Returns:
        MIME type of the file
    """
    try:
        # Try to detect MIME type using python-magic
        mime_type = magic.from_file(file_path, mime=True)
        if mime_type:
            return mime_type
    except Exception as e:
        logger.warning(f"Error detecting MIME type with magic: {str(e)}")
        
    # Fall back to mimetypes module
    mime_type, _ = mimetypes.guess_type(file_path)
    if mime_type:
        return mime_type
    
    # Default to binary if detection fails
    return 'application/octet-stream'


def generate_unique_filename(original_filename: str, prefix: str = '') -> str:
    """
    Generate a unique filename with the original extension.
    
    Args:
        original_filename: Original filename with extension
        prefix: Optional prefix to add to the filename
        
    Returns:
        Unique filename with original extension
    """
    # Get the extension from the original filename
    extension = get_file_extension(original_filename)
    
    # Generate a UUID for uniqueness
    unique_id = str(uuid.uuid4())
    
    # Format current timestamp
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    
    # Combine components to create unique filename
    unique_filename = f"{prefix}{'_' if prefix else ''}{timestamp}_{unique_id}{extension}"
    
    return unique_filename


def sanitize_filename(filename: str) -> str:
    """
    Sanitize a filename to remove invalid characters.
    
    Args:
        filename: Filename to sanitize
        
    Returns:
        Sanitized filename
    """
    # Define invalid characters
    invalid_chars = ['\\', '/', ':', '*', '?', '"', '<', '>', '|']
    
    # Replace invalid characters with underscores
    sanitized = filename
    for char in invalid_chars:
        sanitized = sanitized.replace(char, '_')
    
    # Remove any leading/trailing whitespace
    sanitized = sanitized.strip()
    
    # Ensure filename is not empty after sanitization
    if not sanitized:
        sanitized = "file"
        
    return sanitized


def create_directory_if_not_exists(directory_path: str) -> bool:
    """
    Create a directory if it doesn't exist.
    
    Args:
        directory_path: Path to the directory
        
    Returns:
        True if directory exists or was created successfully, False otherwise
    """
    try:
        path = pathlib.Path(directory_path)
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        return True
    except Exception as e:
        logger.error(f"Error creating directory {directory_path}: {str(e)}")
        return False


def is_csv_file(filename: str) -> bool:
    """
    Check if a file is a CSV file based on extension.
    
    Args:
        filename: Filename to check
        
    Returns:
        True if file is a CSV file, False otherwise
    """
    extension = get_file_extension(filename)
    return extension.lower() == '.csv'


def is_image_file(filename: str) -> bool:
    """
    Check if a file is an image file based on extension.
    
    Args:
        filename: Filename to check
        
    Returns:
        True if file is an image file, False otherwise
    """
    extension = get_file_extension(filename)
    return extension.lower() in IMAGE_EXTENSIONS


def is_document_file(filename: str) -> bool:
    """
    Check if a file is a document file based on extension.
    
    Args:
        filename: Filename to check
        
    Returns:
        True if file is a document file, False otherwise
    """
    extension = get_file_extension(filename)
    return extension.lower() in DOCUMENT_EXTENSIONS


def get_file_size(file_path: str) -> int:
    """
    Get the size of a file in bytes.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Size of the file in bytes, or -1 if file doesn't exist or is inaccessible
    """
    try:
        return os.path.getsize(file_path)
    except (OSError, FileNotFoundError) as e:
        logger.error(f"Error getting file size for {file_path}: {str(e)}")
        return -1


def format_file_size(size_in_bytes: int) -> str:
    """
    Format a file size in bytes to a human-readable string.
    
    Args:
        size_in_bytes: Size of the file in bytes
        
    Returns:
        Human-readable file size (e.g., '2.5 MB')
    """
    # Define units
    units = ['B', 'KB', 'MB', 'GB', 'TB']
    
    # Handle negative or zero size
    if size_in_bytes <= 0:
        return "0 B"
    
    # Calculate appropriate unit
    i = 0
    while size_in_bytes >= 1024 and i < len(units) - 1:
        size_in_bytes /= 1024
        i += 1
    
    # Format to 2 decimal places if needed
    if i > 0:
        return f"{size_in_bytes:.2f} {units[i]}"
    else:
        return f"{size_in_bytes} {units[i]}"


def validate_file(filename: str, file_size: int, allowed_extensions: Optional[Set[str]] = None, 
                 max_size: Optional[int] = None) -> bool:
    """
    Validate a file's extension and size.
    
    Args:
        filename: Filename to validate
        file_size: Size of the file in bytes
        allowed_extensions: Set of allowed extensions (if None, uses ALLOWED_FILE_EXTENSIONS)
        max_size: Maximum allowed size in bytes (if None, uses MAX_FILE_SIZE_BYTES)
        
    Returns:
        True if file is valid, False otherwise
    """
    # Check file extension
    if not is_allowed_file(filename, allowed_extensions):
        return False
    
    # Check file size
    if not validate_file_size(file_size, max_size):
        return False
    
    return True


def validate_file_with_exception(filename: str, file_size: int, allowed_extensions: Optional[Set[str]] = None, 
                                max_size: Optional[int] = None) -> bool:
    """
    Validate a file's extension and size, raising an exception if invalid.
    
    Args:
        filename: Filename to validate
        file_size: Size of the file in bytes
        allowed_extensions: Set of allowed extensions (if None, uses ALLOWED_FILE_EXTENSIONS)
        max_size: Maximum allowed size in bytes (if None, uses MAX_FILE_SIZE_BYTES)
        
    Returns:
        True if file is valid
        
    Raises:
        ValidationException: If file extension or size is invalid
    """
    # Use default values if None provided
    if allowed_extensions is None:
        allowed_extensions = ALLOWED_FILE_EXTENSIONS
    
    if max_size is None:
        max_size = MAX_FILE_SIZE_BYTES
    
    # Check file extension
    if not is_allowed_file(filename, allowed_extensions):
        # Format extensions for display
        formatted_extensions = [ext if ext.startswith('.') else f'.{ext}' for ext in allowed_extensions]
        extensions_str = ", ".join(formatted_extensions)
        raise ValidationException(
            f"Invalid file extension for {filename}. Allowed extensions: {extensions_str}",
            {"filename": filename, "allowed_extensions": list(allowed_extensions)}
        )
    
    # Check file size
    if not validate_file_size(file_size, max_size):
        max_size_formatted = format_file_size(max_size)
        raise ValidationException(
            f"File size ({format_file_size(file_size)}) exceeds maximum allowed size ({max_size_formatted})",
            {"filename": filename, "file_size": file_size, "max_size": max_size}
        )
    
    return True


def split_bucket_and_object_name(file_path: str) -> Tuple[Optional[str], str]:
    """
    Split a file path into bucket name and object name.
    
    Args:
        file_path: File path in format 'bucket_name/object_name'
        
    Returns:
        Tuple of (bucket_name, object_name)
    """
    # Check if path contains a bucket separator
    if '/' in file_path:
        # Split by first '/' to get bucket name and object name
        parts = file_path.split('/', 1)
        return (parts[0], parts[1])
    else:
        # If no separator, assume it's just an object name (default bucket)
        return (None, file_path)