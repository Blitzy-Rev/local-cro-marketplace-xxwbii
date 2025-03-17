"""
Initialization module for the utils package in the Molecular Data Management and CRO Integration Platform.

This module imports and re-exports utility functions and classes from various utility modules to provide
easy access to common functionality throughout the application. This centralized approach simplifies imports
in other modules and promotes code reusability.
"""

# Import all utilities from security_utils
from .security_utils import (
    validate_password_strength,
    get_password_validation_errors,
    generate_secure_token,
    is_common_password,
    calculate_password_strength,
    generate_random_password,
    secure_compare,
    hash_data,
    encode_base64,
    decode_base64,
    sanitize_input,
    PasswordValidator,
    TokenGenerator
)

# Import all utilities from file_utils
from .file_utils import (
    get_file_extension,
    is_allowed_file,
    validate_file_size,
    get_mime_type,
    generate_unique_filename,
    sanitize_filename,
    create_directory_if_not_exists,
    is_csv_file,
    is_image_file,
    is_document_file,
    get_file_size,
    format_file_size,
    validate_file,
    validate_file_with_exception,
    split_bucket_and_object_name
)

# Import all utilities from validation_utils
from .validation_utils import (
    validate_smiles,
    validate_smiles_with_exception,
    validate_molecular_property,
    validate_molecular_property_with_exception,
    validate_email,
    validate_email_with_exception,
    validate_password,
    validate_password_with_exception,
    validate_file_extension,
    validate_file_extension_with_exception,
    validate_csv_file,
    validate_csv_file_with_exception,
    validate_required_fields,
    validate_required_fields_with_exception,
    validate_numeric_value,
    validate_numeric_value_with_exception,
    validate_string_length,
    validate_string_length_with_exception,
    validate_enum_value,
    validate_enum_value_with_exception,
    validate_date_format,
    validate_date_format_with_exception
)

# Import all utilities from molecular_utils
from .molecular_utils import (
    is_valid_smiles,
    canonicalize_smiles,
    calculate_molecular_properties,
    is_property_in_range,
    filter_molecules_by_property,
    calculate_similarity,
    check_lipinski_rule_of_five,
    check_veber_rules,
    check_drug_likeness,
    has_substructure,
    get_molecular_formula
)

# Import all utilities from csv_utils
from .csv_utils import (
    read_csv_file,
    get_csv_headers,
    validate_csv_structure,
    validate_mapping,
    map_csv_columns,
    get_csv_preview,
    process_csv_in_batches,
    extract_molecules_from_dataframe,
    validate_csv_file_with_exception,
    count_csv_rows,
    detect_csv_delimiter,
    validate_smiles_column,
    generate_csv_import_summary
)

# Import all utilities from email_utils
from .email_utils import (
    send_email,
    send_verification_email,
    send_password_reset_email,
    send_notification_email,
    format_email_template
)

# Define __all__ to specify what is exported when using "from app.utils import *"
__all__ = [
    # Security utilities
    "validate_password_strength",
    "get_password_validation_errors",
    "generate_secure_token",
    "is_common_password",
    "calculate_password_strength",
    "generate_random_password",
    "secure_compare",
    "hash_data",
    "encode_base64",
    "decode_base64",
    "sanitize_input",
    "PasswordValidator",
    "TokenGenerator",
    
    # File utilities
    "get_file_extension",
    "is_allowed_file",
    "validate_file_size",
    "get_mime_type",
    "generate_unique_filename",
    "sanitize_filename",
    "create_directory_if_not_exists",
    "is_csv_file",
    "is_image_file",
    "is_document_file",
    "get_file_size",
    "format_file_size",
    "validate_file",
    "validate_file_with_exception",
    "split_bucket_and_object_name",
    
    # Validation utilities
    "validate_smiles",
    "validate_smiles_with_exception",
    "validate_molecular_property",
    "validate_molecular_property_with_exception",
    "validate_email",
    "validate_email_with_exception",
    "validate_password",
    "validate_password_with_exception",
    "validate_file_extension",
    "validate_file_extension_with_exception",
    "validate_csv_file",
    "validate_csv_file_with_exception",
    "validate_required_fields",
    "validate_required_fields_with_exception",
    "validate_numeric_value",
    "validate_numeric_value_with_exception",
    "validate_string_length",
    "validate_string_length_with_exception",
    "validate_enum_value",
    "validate_enum_value_with_exception",
    "validate_date_format",
    "validate_date_format_with_exception",
    
    # Molecular utilities
    "is_valid_smiles",
    "validate_smiles",
    "canonicalize_smiles",
    "calculate_molecular_properties",
    "is_property_in_range",
    "filter_molecules_by_property",
    "calculate_similarity",
    "check_lipinski_rule_of_five",
    "check_veber_rules",
    "check_drug_likeness",
    "has_substructure",
    "get_molecular_formula",
    
    # CSV utilities
    "read_csv_file",
    "get_csv_headers",
    "validate_csv_structure",
    "validate_mapping",
    "map_csv_columns",
    "get_csv_preview",
    "process_csv_in_batches",
    "extract_molecules_from_dataframe",
    "validate_csv_file_with_exception",
    "count_csv_rows",
    "detect_csv_delimiter",
    "validate_smiles_column",
    "generate_csv_import_summary",
    
    # Email utilities
    "send_email",
    "send_verification_email",
    "send_password_reset_email",
    "send_notification_email",
    "format_email_template"
]