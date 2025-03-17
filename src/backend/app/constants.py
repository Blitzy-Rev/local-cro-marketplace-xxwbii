"""
Constants module for the Molecular Data Management and CRO Integration Platform.

This module defines application-wide constants including API prefixes, 
user roles, status enumerations, security settings, and other configuration
values used across the application.
"""

import enum  # standard library

# Application information
VERSION = "1.0.0"
PROJECT_NAME = "Molecular Data Management and CRO Integration Platform"

# API endpoints
API_V1_PREFIX = "/api/v1"
HEALTH_PREFIX = "/health"

# Authentication settings
TOKEN_EXPIRY_MINUTES = 30
REFRESH_TOKEN_EXPIRY_DAYS = 7
ALGORITHM = "RS256"
PASSWORD_RESET_TOKEN_EXPIRY_HOURS = 1
EMAIL_VERIFICATION_TOKEN_EXPIRY_HOURS = 24

# Pagination settings
DEFAULT_PAGE_SIZE = 20
MAX_PAGE_SIZE = 100

# File upload settings
MAX_CSV_SIZE_MB = 50
MAX_FILE_SIZE_MB = 100
ALLOWED_EXTENSIONS = ["csv", "xlsx", "pdf", "zip"]
CSV_CHUNK_SIZE = 1000

# User roles and status enumerations
UserRole = enum.Enum('UserRole', ['PHARMA', 'CRO', 'ADMIN'])
UserStatus = enum.Enum('UserStatus', ['PENDING', 'ACTIVE', 'INACTIVE', 'LOCKED'])
ExperimentStatus = enum.Enum('ExperimentStatus', [
    'DRAFT', 'QUEUED', 'SUBMITTED', 'REJECTED', 
    'QUOTE_PENDING', 'QUOTE_REJECTED', 'IN_PROGRESS', 
    'RESULTS_PENDING', 'RESULTS_AVAILABLE', 'RESULTS_REJECTED', 
    'COMPLETED', 'CANCELLED'
])
SubmissionStatus = enum.Enum('SubmissionStatus', [
    'PENDING', 'REJECTED', 'QUOTE_PROVIDED', 'QUOTE_REJECTED', 
    'APPROVED', 'IN_PROGRESS', 'COMPLETED', 'CANCELLED'
])
ResultStatus = enum.Enum('ResultStatus', ['PENDING', 'UPLOADED', 'APPROVED', 'REJECTED'])
NotificationType = enum.Enum('NotificationType', [
    'EXPERIMENT_STATUS_CHANGE', 'SUBMISSION_CREATED', 
    'QUOTE_PROVIDED', 'RESULTS_UPLOADED', 
    'SYSTEM_ALERT', 'USER_MENTION'
])

# Role-based permissions mapping
ROLE_PERMISSIONS = {
    UserRole.PHARMA: [
        'view_molecules', 'create_molecules', 'edit_molecules', 'delete_molecules',
        'view_libraries', 'create_libraries', 'edit_libraries', 'delete_libraries',
        'view_experiments', 'create_experiments', 'edit_experiments', 'delete_experiments',
        'view_submissions', 'create_submissions',
        'view_results'
    ],
    UserRole.CRO: [
        'view_submissions', 'update_submissions', 'provide_quotes', 
        'upload_results', 'view_results'
    ],
    UserRole.ADMIN: [
        'view_molecules', 'create_molecules', 'edit_molecules', 'delete_molecules',
        'view_libraries', 'create_libraries', 'edit_libraries', 'delete_libraries',
        'view_experiments', 'create_experiments', 'edit_experiments', 'delete_experiments',
        'view_submissions', 'create_submissions',
        'view_results',
        'manage_users', 'view_system', 'configure_system'
    ]
}

# Storage bucket names
BUCKET_NAMES = {
    "CSV_UPLOADS": "csv-uploads",
    "MOLECULE_IMAGES": "molecule-images",
    "EXPERIMENT_FILES": "experiment-files",
    "RESULT_FILES": "result-files",
    "TEMP_FILES": "temp-files"
}

# Logging settings
LOG_LEVELS = {
    "DEBUG": "DEBUG",
    "INFO": "INFO",
    "WARNING": "WARNING",
    "ERROR": "ERROR",
    "CRITICAL": "CRITICAL"
}
DEFAULT_LOG_LEVEL = LOG_LEVELS["INFO"]
DEFAULT_LOG_FORMAT = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"