"""
Utility functions for validating various types of data in the Molecular Data Management
and CRO Integration Platform.

This module provides validation functions for SMILES strings, molecular properties,
CSV data, and other input types to ensure data integrity throughout the application.
"""

import os  # standard library
import re  # standard library
import logging  # standard library
from typing import List, Dict, Optional, Union, Any, Callable  # standard library
from datetime import datetime  # standard library
from email_validator import validate_email as validate_email_format, EmailNotValidError  # version 1.3+
from rdkit import Chem  # version 2023.03+

from ..exceptions import ValidationException, MolecularProcessingException
from ..constants import (
    DEFAULT_MOLECULE_PROPERTIES,
    ALLOWED_EXTENSIONS,
    MAX_CSV_SIZE_MB,
    MAX_FILE_SIZE_MB
)

# Set up logging
logger = logging.getLogger(__name__)

# Global constants
EMAIL_REGEX = re.compile(r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$')
PASSWORD_REGEX = re.compile(r'^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[@$!%*?&])[A-Za-z\d@$!%*?&]{10,}$')
MAX_CSV_SIZE_BYTES = MAX_CSV_SIZE_MB * 1024 * 1024
MAX_FILE_SIZE_BYTES = MAX_FILE_SIZE_MB * 1024 * 1024

# Molecular property validation ranges
PROPERTY_RANGES = {
    "MW": {"min": 0, "max": 2000},
    "LogP": {"min": -10, "max": 10},
    "TPSA": {"min": 0, "max": 500},
    "HBA": {"min": 0, "max": 20},
    "HBD": {"min": 0, "max": 20},
    "RotBonds": {"min": 0, "max": 20}
}


def validate_smiles_structure(smiles: str) -> bool:
    """
    Validates if a SMILES string represents a valid molecular structure using RDKit.
    
    Args:
        smiles: The SMILES string to validate.
        
    Returns:
        True if the SMILES string represents a valid molecule, False otherwise.
    """
    if smiles is None or smiles.strip() == "":
        return False
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.debug(f"Invalid SMILES structure: {smiles}")
            return False
        return True
    except Exception as e:
        logger.error(f"Error validating SMILES structure: {str(e)}")
        return False


def validate_property_value(property_name: str, property_value: Any) -> bool:
    """
    Validates if a molecular property value is within the acceptable range.
    
    Args:
        property_name: The name of the property to validate.
        property_value: The value to validate.
        
    Returns:
        True if the property value is within acceptable range, False otherwise.
    """
    if property_name not in PROPERTY_RANGES:
        # No validation for custom properties
        return True
    
    try:
        value = float(property_value)
        min_value = PROPERTY_RANGES[property_name]["min"]
        max_value = PROPERTY_RANGES[property_name]["max"]
        
        return min_value <= value <= max_value
    except (ValueError, TypeError):
        return False


def validate_smiles(smiles: str) -> bool:
    """
    Validates a SMILES string using the molecular validator.
    
    Args:
        smiles: The SMILES string to validate.
        
    Returns:
        True if the SMILES string is valid, False otherwise.
    """
    if smiles is None or smiles.strip() == "":
        return False
    
    return validate_smiles_structure(smiles)


def validate_smiles_with_exception(smiles: str) -> bool:
    """
    Validates a SMILES string and raises an exception if invalid.
    
    Args:
        smiles: The SMILES string to validate.
        
    Returns:
        True if the SMILES string is valid.
        
    Raises:
        ValidationException: If the SMILES string is invalid.
    """
    if smiles is None or smiles.strip() == "":
        raise ValidationException(
            "SMILES string cannot be empty",
            {"smiles": smiles}
        )
    
    if not validate_smiles_structure(smiles):
        raise ValidationException(
            "Invalid SMILES structure",
            {"smiles": smiles}
        )
    
    return True


def validate_molecular_property(property_name: str, property_value: Any) -> bool:
    """
    Validates a molecular property value against expected ranges.
    
    Args:
        property_name: The name of the property to validate.
        property_value: The value to validate.
        
    Returns:
        True if the property value is valid, False otherwise.
    """
    return validate_property_value(property_name, property_value)


def validate_molecular_property_with_exception(property_name: str, property_value: Any) -> bool:
    """
    Validates a molecular property value and raises an exception if invalid.
    
    Args:
        property_name: The name of the property to validate.
        property_value: The value to validate.
        
    Returns:
        True if the property value is valid.
        
    Raises:
        ValidationException: If the property value is invalid.
    """
    if property_name not in PROPERTY_RANGES:
        # No validation for custom properties
        return True
    
    try:
        value = float(property_value)
    except (ValueError, TypeError):
        raise ValidationException(
            f"Property {property_name} must be a numeric value",
            {"property_name": property_name, "property_value": property_value}
        )
    
    min_value = PROPERTY_RANGES[property_name]["min"]
    max_value = PROPERTY_RANGES[property_name]["max"]
    
    if not (min_value <= value <= max_value):
        raise ValidationException(
            f"Property {property_name} must be between {min_value} and {max_value}",
            {
                "property_name": property_name,
                "property_value": value,
                "min_value": min_value,
                "max_value": max_value
            }
        )
    
    return True


def validate_email(email: str) -> bool:
    """
    Validates an email address format.
    
    Args:
        email: The email address to validate.
        
    Returns:
        True if the email is valid, False otherwise.
    """
    if email is None or email.strip() == "":
        return False
    
    try:
        validate_email_format(email)
        return True
    except EmailNotValidError:
        return False


def validate_email_with_exception(email: str) -> bool:
    """
    Validates an email address and raises an exception if invalid.
    
    Args:
        email: The email address to validate.
        
    Returns:
        True if the email is valid.
        
    Raises:
        ValidationException: If the email is invalid.
    """
    if email is None or email.strip() == "":
        raise ValidationException(
            "Email address cannot be empty",
            {"email": email}
        )
    
    try:
        validate_email_format(email)
        return True
    except EmailNotValidError as e:
        raise ValidationException(
            "Invalid email address format",
            {"email": email, "details": str(e)}
        )


def validate_password(password: str) -> bool:
    """
    Validates a password against security requirements.
    
    Password must be at least 10 characters long and contain at least:
    - One lowercase letter
    - One uppercase letter
    - One digit
    - One special character (@$!%*?&)
    
    Args:
        password: The password to validate.
        
    Returns:
        True if the password meets requirements, False otherwise.
    """
    if password is None or password.strip() == "":
        return False
    
    return bool(PASSWORD_REGEX.match(password))


def validate_password_with_exception(password: str) -> bool:
    """
    Validates a password and raises an exception if requirements are not met.
    
    Args:
        password: The password to validate.
        
    Returns:
        True if the password is valid.
        
    Raises:
        ValidationException: If the password does not meet security requirements.
    """
    if password is None or password.strip() == "":
        raise ValidationException(
            "Password cannot be empty",
            {"password": "empty"}
        )
    
    if not PASSWORD_REGEX.match(password):
        raise ValidationException(
            "Password does not meet security requirements",
            {
                "requirements": "Password must be at least 10 characters long and contain at least one "
                               "lowercase letter, one uppercase letter, one digit, and one special character (@$!%*?&)"
            }
        )
    
    return True


def validate_file_extension(filename: str, allowed_extensions: List[str] = None) -> bool:
    """
    Validates if a file has an allowed extension.
    
    Args:
        filename: The filename to validate.
        allowed_extensions: List of allowed file extensions (without dot).
                           If None, uses ALLOWED_EXTENSIONS from constants.
        
    Returns:
        True if the file extension is allowed, False otherwise.
    """
    if filename is None or filename.strip() == "":
        return False
    
    if allowed_extensions is None:
        allowed_extensions = ALLOWED_EXTENSIONS
    
    extension = os.path.splitext(filename)[1][1:].lower()
    return extension.lower() in [ext.lower() for ext in allowed_extensions]


def validate_file_extension_with_exception(filename: str, allowed_extensions: List[str] = None) -> bool:
    """
    Validates file extension and raises an exception if not allowed.
    
    Args:
        filename: The filename to validate.
        allowed_extensions: List of allowed file extensions (without dot).
                           If None, uses ALLOWED_EXTENSIONS from constants.
        
    Returns:
        True if the file extension is allowed.
        
    Raises:
        ValidationException: If the file extension is not allowed.
    """
    if filename is None or filename.strip() == "":
        raise ValidationException(
            "Filename cannot be empty",
            {"filename": filename}
        )
    
    if allowed_extensions is None:
        allowed_extensions = ALLOWED_EXTENSIONS
    
    extension = os.path.splitext(filename)[1][1:].lower()
    if extension.lower() not in [ext.lower() for ext in allowed_extensions]:
        raise ValidationException(
            f"File extension '{extension}' is not allowed",
            {
                "filename": filename,
                "extension": extension,
                "allowed_extensions": allowed_extensions
            }
        )
    
    return True


def validate_file_size(file_size: int, max_size_bytes: int = None) -> bool:
    """
    Validates if a file size is within allowed limits.
    
    Args:
        file_size: The file size in bytes.
        max_size_bytes: Maximum allowed file size in bytes.
                        If None, uses MAX_FILE_SIZE_BYTES from constants.
        
    Returns:
        True if the file size is within limits, False otherwise.
    """
    if max_size_bytes is None:
        max_size_bytes = MAX_FILE_SIZE_BYTES
    
    return file_size <= max_size_bytes


def validate_file_size_with_exception(file_size: int, max_size_bytes: int = None) -> bool:
    """
    Validates file size and raises an exception if too large.
    
    Args:
        file_size: The file size in bytes.
        max_size_bytes: Maximum allowed file size in bytes.
                        If None, uses MAX_FILE_SIZE_BYTES from constants.
        
    Returns:
        True if the file size is within limits.
        
    Raises:
        ValidationException: If the file size exceeds the maximum allowed.
    """
    if max_size_bytes is None:
        max_size_bytes = MAX_FILE_SIZE_BYTES
    
    if file_size > max_size_bytes:
        max_size_mb = max_size_bytes / (1024 * 1024)
        file_size_mb = file_size / (1024 * 1024)
        raise ValidationException(
            f"File size ({file_size_mb:.2f} MB) exceeds maximum allowed size ({max_size_mb:.2f} MB)",
            {
                "file_size_bytes": file_size,
                "file_size_mb": file_size_mb,
                "max_size_bytes": max_size_bytes,
                "max_size_mb": max_size_mb
            }
        )
    
    return True


def validate_csv_file(filename: str, file_size: int) -> bool:
    """
    Validates a CSV file for extension and size.
    
    Args:
        filename: The filename to validate.
        file_size: The file size in bytes.
        
    Returns:
        True if the CSV file is valid, False otherwise.
    """
    return (
        validate_file_extension(filename, ["csv"]) and
        validate_file_size(file_size, MAX_CSV_SIZE_BYTES)
    )


def validate_csv_file_with_exception(filename: str, file_size: int) -> bool:
    """
    Validates a CSV file and raises an exception if invalid.
    
    Args:
        filename: The filename to validate.
        file_size: The file size in bytes.
        
    Returns:
        True if the CSV file is valid.
        
    Raises:
        ValidationException: If the CSV file is invalid.
    """
    validate_file_extension_with_exception(filename, ["csv"])
    validate_file_size_with_exception(file_size, MAX_CSV_SIZE_BYTES)
    return True


def validate_required_fields(data: Dict[str, Any], required_fields: List[str]) -> bool:
    """
    Validates that all required fields are present in a dictionary.
    
    Args:
        data: The dictionary to check.
        required_fields: List of field names that must be present and not None.
        
    Returns:
        True if all required fields are present, False otherwise.
    """
    if data is None:
        return False
    
    for field in required_fields:
        if field not in data or data[field] is None:
            return False
    
    return True


def validate_required_fields_with_exception(data: Dict[str, Any], required_fields: List[str]) -> bool:
    """
    Validates required fields and raises an exception if any are missing.
    
    Args:
        data: The dictionary to check.
        required_fields: List of field names that must be present and not None.
        
    Returns:
        True if all required fields are present.
        
    Raises:
        ValidationException: If any required fields are missing.
    """
    if data is None:
        raise ValidationException(
            "Data cannot be None",
            {"required_fields": required_fields}
        )
    
    missing_fields = []
    for field in required_fields:
        if field not in data or data[field] is None:
            missing_fields.append(field)
    
    if missing_fields:
        raise ValidationException(
            "Required fields are missing",
            {"missing_fields": missing_fields, "required_fields": required_fields}
        )
    
    return True


def validate_numeric_value(value: Any, min_value: float, max_value: float) -> bool:
    """
    Validates if a value is numeric and within a specified range.
    
    Args:
        value: The value to validate.
        min_value: The minimum allowed value (inclusive).
        max_value: The maximum allowed value (inclusive).
        
    Returns:
        True if the value is numeric and within range, False otherwise.
    """
    try:
        numeric_value = float(value)
        return min_value <= numeric_value <= max_value
    except (ValueError, TypeError):
        return False


def validate_numeric_value_with_exception(
    value: Any, min_value: float, max_value: float, field_name: str
) -> bool:
    """
    Validates numeric value and raises an exception if invalid.
    
    Args:
        value: The value to validate.
        min_value: The minimum allowed value (inclusive).
        max_value: The maximum allowed value (inclusive).
        field_name: The name of the field being validated (for error messages).
        
    Returns:
        True if the value is numeric and within range.
        
    Raises:
        ValidationException: If the value is not numeric or out of range.
    """
    try:
        numeric_value = float(value)
    except (ValueError, TypeError):
        raise ValidationException(
            f"{field_name} must be a numeric value",
            {"field_name": field_name, "value": value}
        )
    
    if not (min_value <= numeric_value <= max_value):
        raise ValidationException(
            f"{field_name} must be between {min_value} and {max_value}",
            {
                "field_name": field_name,
                "value": numeric_value,
                "min_value": min_value,
                "max_value": max_value
            }
        )
    
    return True


def validate_string_length(value: str, min_length: int, max_length: int) -> bool:
    """
    Validates if a string length is within specified limits.
    
    Args:
        value: The string to validate.
        min_length: The minimum allowed length (inclusive).
        max_length: The maximum allowed length (inclusive).
        
    Returns:
        True if the string length is within limits, False otherwise.
    """
    if value is None:
        return False
    
    length = len(value)
    return min_length <= length <= max_length


def validate_string_length_with_exception(
    value: str, min_length: int, max_length: int, field_name: str
) -> bool:
    """
    Validates string length and raises an exception if outside limits.
    
    Args:
        value: The string to validate.
        min_length: The minimum allowed length (inclusive).
        max_length: The maximum allowed length (inclusive).
        field_name: The name of the field being validated (for error messages).
        
    Returns:
        True if the string length is within limits.
        
    Raises:
        ValidationException: If the string length is outside specified limits.
    """
    if value is None:
        raise ValidationException(
            f"{field_name} cannot be None",
            {"field_name": field_name}
        )
    
    length = len(value)
    if not (min_length <= length <= max_length):
        raise ValidationException(
            f"{field_name} length must be between {min_length} and {max_length} characters",
            {
                "field_name": field_name,
                "current_length": length,
                "min_length": min_length,
                "max_length": max_length
            }
        )
    
    return True


def validate_enum_value(value: Any, enum_type: type) -> bool:
    """
    Validates if a value is a valid member of an enumeration.
    
    Args:
        value: The value to validate.
        enum_type: The enumeration type to check against.
        
    Returns:
        True if the value is a valid enum member, False otherwise.
    """
    if value is None:
        return False
    
    try:
        enum_type(value)
        return True
    except (ValueError, TypeError):
        return False


def validate_enum_value_with_exception(value: Any, enum_type: type, field_name: str) -> bool:
    """
    Validates enum value and raises an exception if invalid.
    
    Args:
        value: The value to validate.
        enum_type: The enumeration type to check against.
        field_name: The name of the field being validated (for error messages).
        
    Returns:
        True if the value is a valid enum member.
        
    Raises:
        ValidationException: If the value is not a valid enum member.
    """
    if value is None:
        raise ValidationException(
            f"{field_name} cannot be None",
            {"field_name": field_name}
        )
    
    try:
        enum_type(value)
        return True
    except (ValueError, TypeError):
        valid_values = [e.name for e in enum_type]
        raise ValidationException(
            f"Invalid value for {field_name}",
            {
                "field_name": field_name,
                "value": value,
                "valid_values": valid_values
            }
        )


def validate_date_format(date_string: str, format_string: str) -> bool:
    """
    Validates if a string is in the correct date format.
    
    Args:
        date_string: The date string to validate.
        format_string: The expected format (e.g., '%Y-%m-%d').
        
    Returns:
        True if the date string matches format, False otherwise.
    """
    if date_string is None or date_string.strip() == "":
        return False
    
    try:
        datetime.strptime(date_string, format_string)
        return True
    except ValueError:
        return False


def validate_date_format_with_exception(
    date_string: str, format_string: str, field_name: str
) -> bool:
    """
    Validates date format and raises an exception if invalid.
    
    Args:
        date_string: The date string to validate.
        format_string: The expected format (e.g., '%Y-%m-%d').
        field_name: The name of the field being validated (for error messages).
        
    Returns:
        True if the date string matches format.
        
    Raises:
        ValidationException: If the date string does not match the expected format.
    """
    if date_string is None or date_string.strip() == "":
        raise ValidationException(
            f"{field_name} cannot be empty",
            {"field_name": field_name}
        )
    
    try:
        datetime.strptime(date_string, format_string)
        return True
    except ValueError:
        raise ValidationException(
            f"{field_name} must be in format: {format_string}",
            {
                "field_name": field_name,
                "value": date_string,
                "expected_format": format_string
            }
        )