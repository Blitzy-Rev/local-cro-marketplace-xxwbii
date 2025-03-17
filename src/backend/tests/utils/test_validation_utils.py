"""
Unit tests for validation utility functions in the Molecular Data Management and CRO Integration Platform.

This module contains tests for all validation functions including SMILES validation, property validation,
email validation, password validation, file validation, and general validation utilities.
"""

import pytest
from unittest.mock import patch
import datetime
from enum import Enum

from app.exceptions import ValidationException, MolecularProcessingException
from app.constants import (
    DEFAULT_MOLECULE_PROPERTIES,
    ALLOWED_EXTENSIONS,
    MAX_CSV_SIZE_MB,
    MAX_FILE_SIZE_MB,
    UserRole
)
from app.utils import validation_utils
from app.utils.validation_utils import PROPERTY_RANGES


def test_validate_smiles():
    """Test the validate_smiles function with valid and invalid SMILES strings."""
    # Test valid SMILES
    assert validation_utils.validate_smiles("CCO") is True
    assert validation_utils.validate_smiles("c1ccccc1") is True
    assert validation_utils.validate_smiles("CC(=O)O") is True
    
    # Test invalid SMILES
    assert validation_utils.validate_smiles("X") is False
    assert validation_utils.validate_smiles("C(C") is False
    assert validation_utils.validate_smiles("invalid") is False
    
    # Test edge cases
    assert validation_utils.validate_smiles("") is False
    assert validation_utils.validate_smiles(None) is False


def test_validate_smiles_with_exception():
    """Test the validate_smiles_with_exception function with valid and invalid SMILES strings."""
    # Test valid SMILES
    assert validation_utils.validate_smiles_with_exception("CCO") is True
    assert validation_utils.validate_smiles_with_exception("c1ccccc1") is True
    assert validation_utils.validate_smiles_with_exception("CC(=O)O") is True
    
    # Test invalid SMILES - should raise ValidationException
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_smiles_with_exception("X")
    assert "Invalid SMILES structure" in str(exc_info.value)
    
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_smiles_with_exception("C(C")
    assert "Invalid SMILES structure" in str(exc_info.value)
    
    # Test empty string
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_smiles_with_exception("")
    assert "SMILES string cannot be empty" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_smiles_with_exception(None)
    assert "SMILES string cannot be empty" in str(exc_info.value)


def test_validate_molecular_property():
    """Test the validate_molecular_property function with valid and invalid property values."""
    # Test valid property values
    assert validation_utils.validate_molecular_property("MW", 500) is True
    assert validation_utils.validate_molecular_property("LogP", 2.5) is True
    assert validation_utils.validate_molecular_property("TPSA", 100) is True
    
    # Test out-of-range values
    assert validation_utils.validate_molecular_property("MW", 3000) is False
    assert validation_utils.validate_molecular_property("LogP", 15) is False
    
    # Test non-numeric values
    assert validation_utils.validate_molecular_property("MW", "not-a-number") is False
    
    # Test custom property (not in PROPERTY_RANGES)
    assert validation_utils.validate_molecular_property("CustomProperty", 100) is True
    assert validation_utils.validate_molecular_property("CustomProperty", "string-value") is True


def test_validate_molecular_property_with_exception():
    """Test the validate_molecular_property_with_exception function with valid and invalid property values."""
    # Test valid property values
    assert validation_utils.validate_molecular_property_with_exception("MW", 500) is True
    assert validation_utils.validate_molecular_property_with_exception("LogP", 2.5) is True
    
    # Test out-of-range values
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_molecular_property_with_exception("MW", 3000)
    assert "must be between" in str(exc_info.value)
    
    # Test non-numeric values
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_molecular_property_with_exception("MW", "not-a-number")
    assert "must be a numeric value" in str(exc_info.value)
    
    # Test custom property (not in PROPERTY_RANGES)
    assert validation_utils.validate_molecular_property_with_exception("CustomProperty", 100) is True
    assert validation_utils.validate_molecular_property_with_exception("CustomProperty", "string-value") is True


def test_validate_email():
    """Test the validate_email function with valid and invalid email addresses."""
    # Test valid emails
    assert validation_utils.validate_email("user@example.com") is True
    assert validation_utils.validate_email("name.surname@domain.co.uk") is True
    assert validation_utils.validate_email("user+tag@example.com") is True
    
    # Test invalid emails
    assert validation_utils.validate_email("user@") is False
    assert validation_utils.validate_email("@domain.com") is False
    assert validation_utils.validate_email("user@domain") is False
    assert validation_utils.validate_email("user.domain.com") is False
    
    # Test edge cases
    assert validation_utils.validate_email("") is False
    assert validation_utils.validate_email(None) is False


def test_validate_email_with_exception():
    """Test the validate_email_with_exception function with valid and invalid email addresses."""
    # Test valid emails
    assert validation_utils.validate_email_with_exception("user@example.com") is True
    assert validation_utils.validate_email_with_exception("name.surname@domain.co.uk") is True
    
    # Test invalid emails
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_email_with_exception("user@")
    assert "Invalid email address" in str(exc_info.value)
    
    # Test empty string
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_email_with_exception("")
    assert "Email address cannot be empty" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_email_with_exception(None)
    assert "Email address cannot be empty" in str(exc_info.value)


def test_validate_password():
    """Test the validate_password function with valid and invalid passwords."""
    # Test valid passwords (meets all requirements)
    assert validation_utils.validate_password("Password123!") is True
    assert validation_utils.validate_password("Secure@Password2023") is True
    
    # Test passwords missing lowercase letters
    assert validation_utils.validate_password("PASSWORD123!") is False
    
    # Test passwords missing uppercase letters
    assert validation_utils.validate_password("password123!") is False
    
    # Test passwords missing digits
    assert validation_utils.validate_password("Password!") is False
    
    # Test passwords missing special characters
    assert validation_utils.validate_password("Password123") is False
    
    # Test passwords shorter than minimum length
    assert validation_utils.validate_password("Pass1!") is False
    
    # Test edge cases
    assert validation_utils.validate_password("") is False
    assert validation_utils.validate_password(None) is False


def test_validate_password_with_exception():
    """Test the validate_password_with_exception function with valid and invalid passwords."""
    # Test valid passwords
    assert validation_utils.validate_password_with_exception("Password123!") is True
    
    # Test invalid passwords
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_password_with_exception("password123!")
    assert "does not meet security requirements" in str(exc_info.value)
    
    # Test empty string
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_password_with_exception("")
    assert "Password cannot be empty" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_password_with_exception(None)
    assert "Password cannot be empty" in str(exc_info.value)


def test_validate_file_extension():
    """Test the validate_file_extension function with various filenames."""
    # Test allowed extensions
    assert validation_utils.validate_file_extension("file.csv") is True
    assert validation_utils.validate_file_extension("document.pdf") is True
    assert validation_utils.validate_file_extension("data.xlsx") is True
    
    # Test disallowed extensions
    assert validation_utils.validate_file_extension("file.exe") is False
    assert validation_utils.validate_file_extension("script.js") is False
    
    # Test uppercase extensions
    assert validation_utils.validate_file_extension("FILE.CSV") is True
    
    # Test filename without extension
    assert validation_utils.validate_file_extension("file_without_extension") is False
    
    # Test custom allowed extensions
    assert validation_utils.validate_file_extension("file.txt", ["txt", "docx"]) is True
    assert validation_utils.validate_file_extension("file.csv", ["txt", "docx"]) is False
    
    # Test edge cases
    assert validation_utils.validate_file_extension("") is False
    assert validation_utils.validate_file_extension(None) is False


def test_validate_file_extension_with_exception():
    """Test the validate_file_extension_with_exception function with various filenames."""
    # Test allowed extensions
    assert validation_utils.validate_file_extension_with_exception("file.csv") is True
    
    # Test disallowed extensions
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_file_extension_with_exception("file.exe")
    assert "File extension 'exe' is not allowed" in str(exc_info.value)
    
    # Test custom allowed extensions
    assert validation_utils.validate_file_extension_with_exception("file.txt", ["txt", "docx"]) is True
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_file_extension_with_exception("file.csv", ["txt", "docx"])
    assert "File extension 'csv' is not allowed" in str(exc_info.value)
    
    # Test empty string
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_file_extension_with_exception("")
    assert "Filename cannot be empty" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_file_extension_with_exception(None)
    assert "Filename cannot be empty" in str(exc_info.value)


def test_validate_file_size():
    """Test the validate_file_size function with various file sizes."""
    max_size = 10 * 1024 * 1024  # 10 MB
    
    # Test file size under limit
    assert validation_utils.validate_file_size(5 * 1024 * 1024, max_size) is True
    
    # Test file size equal to limit
    assert validation_utils.validate_file_size(max_size, max_size) is True
    
    # Test file size over limit
    assert validation_utils.validate_file_size(15 * 1024 * 1024, max_size) is False
    
    # Test with default max size
    assert validation_utils.validate_file_size(1 * 1024 * 1024) is True


def test_validate_file_size_with_exception():
    """Test the validate_file_size_with_exception function with various file sizes."""
    max_size = 10 * 1024 * 1024  # 10 MB
    
    # Test file size under limit
    assert validation_utils.validate_file_size_with_exception(5 * 1024 * 1024, max_size) is True
    
    # Test file size equal to limit
    assert validation_utils.validate_file_size_with_exception(max_size, max_size) is True
    
    # Test file size over limit
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_file_size_with_exception(15 * 1024 * 1024, max_size)
    assert "File size" in str(exc_info.value)
    assert "exceeds maximum allowed size" in str(exc_info.value)


def test_validate_csv_file():
    """Test the validate_csv_file function with various CSV files."""
    # Test valid CSV file
    assert validation_utils.validate_csv_file("file.csv", 5 * 1024 * 1024) is True
    
    # Test invalid extension
    assert validation_utils.validate_csv_file("file.txt", 5 * 1024 * 1024) is False
    
    # Test file size over limit
    max_csv_size = MAX_CSV_SIZE_MB * 1024 * 1024
    assert validation_utils.validate_csv_file("file.csv", max_csv_size + 1) is False


def test_validate_csv_file_with_exception():
    """Test the validate_csv_file_with_exception function with various CSV files."""
    # Test valid CSV file
    assert validation_utils.validate_csv_file_with_exception("file.csv", 5 * 1024 * 1024) is True
    
    # Test invalid extension
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_csv_file_with_exception("file.txt", 5 * 1024 * 1024)
    assert "File extension 'txt' is not allowed" in str(exc_info.value)
    
    # Test file size over limit
    max_csv_size = MAX_CSV_SIZE_MB * 1024 * 1024
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_csv_file_with_exception("file.csv", max_csv_size + 1)
    assert "File size" in str(exc_info.value)
    assert "exceeds maximum allowed size" in str(exc_info.value)


def test_validate_required_fields():
    """Test the validate_required_fields function with various data dictionaries."""
    # Test data with all required fields
    data = {"field1": "value1", "field2": "value2", "field3": "value3"}
    required_fields = ["field1", "field2"]
    assert validation_utils.validate_required_fields(data, required_fields) is True
    
    # Test data missing some required fields
    data = {"field1": "value1", "field3": "value3"}
    required_fields = ["field1", "field2"]
    assert validation_utils.validate_required_fields(data, required_fields) is False
    
    # Test data with None values for required fields
    data = {"field1": "value1", "field2": None}
    required_fields = ["field1", "field2"]
    assert validation_utils.validate_required_fields(data, required_fields) is False
    
    # Test empty data dictionary
    assert validation_utils.validate_required_fields({}, ["field1"]) is False
    
    # Test None data
    assert validation_utils.validate_required_fields(None, ["field1"]) is False


def test_validate_required_fields_with_exception():
    """Test the validate_required_fields_with_exception function with various data dictionaries."""
    # Test data with all required fields
    data = {"field1": "value1", "field2": "value2", "field3": "value3"}
    required_fields = ["field1", "field2"]
    assert validation_utils.validate_required_fields_with_exception(data, required_fields) is True
    
    # Test data missing some required fields
    data = {"field1": "value1", "field3": "value3"}
    required_fields = ["field1", "field2"]
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_required_fields_with_exception(data, required_fields)
    assert "Required fields are missing" in str(exc_info.value)
    
    # Test data with None values for required fields
    data = {"field1": "value1", "field2": None}
    required_fields = ["field1", "field2"]
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_required_fields_with_exception(data, required_fields)
    assert "Required fields are missing" in str(exc_info.value)
    
    # Test empty data dictionary
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_required_fields_with_exception({}, ["field1"])
    assert "Required fields are missing" in str(exc_info.value)
    
    # Test None data
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_required_fields_with_exception(None, ["field1"])
    assert "Data cannot be None" in str(exc_info.value)


def test_validate_numeric_value():
    """Test the validate_numeric_value function with various values."""
    # Test numeric value within range
    assert validation_utils.validate_numeric_value(5, 0, 10) is True
    
    # Test numeric value at boundaries
    assert validation_utils.validate_numeric_value(0, 0, 10) is True
    assert validation_utils.validate_numeric_value(10, 0, 10) is True
    
    # Test numeric value outside range
    assert validation_utils.validate_numeric_value(-1, 0, 10) is False
    assert validation_utils.validate_numeric_value(11, 0, 10) is False
    
    # Test non-numeric value
    assert validation_utils.validate_numeric_value("not-a-number", 0, 10) is False
    
    # Test None value
    assert validation_utils.validate_numeric_value(None, 0, 10) is False


def test_validate_numeric_value_with_exception():
    """Test the validate_numeric_value_with_exception function with various values."""
    # Test numeric value within range
    assert validation_utils.validate_numeric_value_with_exception(5, 0, 10, "test_field") is True
    
    # Test numeric value at boundaries
    assert validation_utils.validate_numeric_value_with_exception(0, 0, 10, "test_field") is True
    assert validation_utils.validate_numeric_value_with_exception(10, 0, 10, "test_field") is True
    
    # Test numeric value outside range
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_numeric_value_with_exception(-1, 0, 10, "test_field")
    assert "test_field must be between 0 and 10" in str(exc_info.value)
    
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_numeric_value_with_exception(11, 0, 10, "test_field")
    assert "test_field must be between 0 and 10" in str(exc_info.value)
    
    # Test non-numeric value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_numeric_value_with_exception("not-a-number", 0, 10, "test_field")
    assert "test_field must be a numeric value" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_numeric_value_with_exception(None, 0, 10, "test_field")
    assert "test_field must be a numeric value" in str(exc_info.value)


def test_validate_string_length():
    """Test the validate_string_length function with various strings."""
    # Test string length within range
    assert validation_utils.validate_string_length("test", 1, 10) is True
    
    # Test string length at boundaries
    assert validation_utils.validate_string_length("a", 1, 10) is True
    assert validation_utils.validate_string_length("1234567890", 1, 10) is True
    
    # Test string length outside range
    assert validation_utils.validate_string_length("", 1, 10) is False
    assert validation_utils.validate_string_length("12345678901", 1, 10) is False
    
    # Test None value
    assert validation_utils.validate_string_length(None, 1, 10) is False


def test_validate_string_length_with_exception():
    """Test the validate_string_length_with_exception function with various strings."""
    # Test string length within range
    assert validation_utils.validate_string_length_with_exception("test", 1, 10, "test_field") is True
    
    # Test string length at boundaries
    assert validation_utils.validate_string_length_with_exception("a", 1, 10, "test_field") is True
    assert validation_utils.validate_string_length_with_exception("1234567890", 1, 10, "test_field") is True
    
    # Test string length outside range
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_string_length_with_exception("", 1, 10, "test_field")
    assert "test_field length must be between 1 and 10 characters" in str(exc_info.value)
    
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_string_length_with_exception("12345678901", 1, 10, "test_field")
    assert "test_field length must be between 1 and 10 characters" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_string_length_with_exception(None, 1, 10, "test_field")
    assert "test_field cannot be None" in str(exc_info.value)


def test_validate_enum_value():
    """Test the validate_enum_value function with various enum values."""
    # Test valid enum value (using UserRole enum)
    assert validation_utils.validate_enum_value(UserRole.PHARMA, UserRole) is True
    assert validation_utils.validate_enum_value("PHARMA", UserRole) is True
    
    # Test invalid enum value
    assert validation_utils.validate_enum_value("INVALID_ROLE", UserRole) is False
    
    # Test None value
    assert validation_utils.validate_enum_value(None, UserRole) is False


def test_validate_enum_value_with_exception():
    """Test the validate_enum_value_with_exception function with various enum values."""
    # Test valid enum value
    assert validation_utils.validate_enum_value_with_exception(UserRole.PHARMA, UserRole, "user_role") is True
    assert validation_utils.validate_enum_value_with_exception("PHARMA", UserRole, "user_role") is True
    
    # Test invalid enum value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_enum_value_with_exception("INVALID_ROLE", UserRole, "user_role")
    assert "Invalid value for user_role" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_enum_value_with_exception(None, UserRole, "user_role")
    assert "user_role cannot be None" in str(exc_info.value)


def test_validate_date_format():
    """Test the validate_date_format function with various date strings."""
    # Test valid date string
    assert validation_utils.validate_date_format("2023-06-01", "%Y-%m-%d") is True
    
    # Test valid date string with different format
    assert validation_utils.validate_date_format("01/06/2023", "%d/%m/%Y") is True
    
    # Test invalid date values
    assert validation_utils.validate_date_format("2023-13-01", "%Y-%m-%d") is False  # Invalid month
    assert validation_utils.validate_date_format("2023-06-32", "%Y-%m-%d") is False  # Invalid day
    
    # Test date string in wrong format
    assert validation_utils.validate_date_format("01/06/2023", "%Y-%m-%d") is False
    
    # Test empty string
    assert validation_utils.validate_date_format("", "%Y-%m-%d") is False
    
    # Test None value
    assert validation_utils.validate_date_format(None, "%Y-%m-%d") is False


def test_validate_date_format_with_exception():
    """Test the validate_date_format_with_exception function with various date strings."""
    # Test valid date string
    assert validation_utils.validate_date_format_with_exception("2023-06-01", "%Y-%m-%d", "date_field") is True
    
    # Test invalid date values
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_date_format_with_exception("2023-13-01", "%Y-%m-%d", "date_field")
    assert "date_field must be in format: %Y-%m-%d" in str(exc_info.value)
    
    # Test date string in wrong format
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_date_format_with_exception("01/06/2023", "%Y-%m-%d", "date_field")
    assert "date_field must be in format: %Y-%m-%d" in str(exc_info.value)
    
    # Test empty string
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_date_format_with_exception("", "%Y-%m-%d", "date_field")
    assert "date_field cannot be empty" in str(exc_info.value)
    
    # Test None value
    with pytest.raises(ValidationException) as exc_info:
        validation_utils.validate_date_format_with_exception(None, "%Y-%m-%d", "date_field")
    assert "date_field cannot be empty" in str(exc_info.value)


def test_validate_smiles_structure():
    """Test the validate_smiles_structure function directly."""
    # Test valid SMILES structures
    assert validation_utils.validate_smiles_structure("CCO") is True
    assert validation_utils.validate_smiles_structure("c1ccccc1") is True
    assert validation_utils.validate_smiles_structure("CC(=O)O") is True
    
    # Test invalid SMILES structures
    assert validation_utils.validate_smiles_structure("X") is False
    assert validation_utils.validate_smiles_structure("C(C") is False
    assert validation_utils.validate_smiles_structure("invalid") is False
    
    # Test edge cases
    assert validation_utils.validate_smiles_structure("") is False
    assert validation_utils.validate_smiles_structure(None) is False


def test_validate_property_value():
    """Test the validate_property_value function directly."""
    # Test property values within range
    assert validation_utils.validate_property_value("MW", 500) is True
    assert validation_utils.validate_property_value("LogP", 2.5) is True
    
    # Test property values outside range
    assert validation_utils.validate_property_value("MW", 3000) is False
    assert validation_utils.validate_property_value("LogP", 15) is False
    
    # Test non-numeric property values
    assert validation_utils.validate_property_value("MW", "not-a-number") is False
    
    # Test custom property not in PROPERTY_RANGES
    assert validation_utils.validate_property_value("CustomProperty", 100) is True
    assert validation_utils.validate_property_value("CustomProperty", "string-value") is True


def test_property_ranges_content():
    """Test that PROPERTY_RANGES contains expected properties and ranges."""
    # Verify PROPERTY_RANGES contains expected properties
    assert "MW" in PROPERTY_RANGES
    assert "LogP" in PROPERTY_RANGES
    assert "TPSA" in PROPERTY_RANGES
    assert "HBA" in PROPERTY_RANGES
    assert "HBD" in PROPERTY_RANGES
    assert "RotBonds" in PROPERTY_RANGES
    
    # Verify each property has min and max values
    for prop_name, prop_range in PROPERTY_RANGES.items():
        assert "min" in prop_range
        assert "max" in prop_range
        
        # Verify min is less than max
        assert prop_range["min"] < prop_range["max"]
        
    # Verify specific ranges for key properties
    assert PROPERTY_RANGES["MW"]["min"] == 0
    assert PROPERTY_RANGES["MW"]["max"] == 2000
    assert PROPERTY_RANGES["LogP"]["min"] == -10
    assert PROPERTY_RANGES["LogP"]["max"] == 10