"""
Initializes the API v1 schemas package and exports all schema classes.

This module re-exports all Pydantic schema models from individual schema modules,
allowing consumers to import schemas directly from the schemas package rather than
from individual modules, simplifying imports throughout the application.

Example:
    from app.api.api_v1.schemas import LoginRequest, MoleculeRead

Version: 1.0.0
"""

# Authentication schemas
from .auth import (
    LoginRequest, LoginResponse, RefreshTokenRequest, RegistrationRequest, 
    RegistrationResponse, EmailVerificationRequest, PasswordResetRequest,
    PasswordResetConfirmRequest, ChangePasswordRequest, MessageResponse
)

# User management schemas
from .users import (
    UserProfileResponse, UserUpdateRequest, UserFilterParams, 
    UserListResponse, UserStatusUpdateRequest
)

# Molecule schemas
from .molecules import (
    MoleculeRead, MoleculeList, PropertyFilter, MoleculeFlagRequest, 
    MoleculeBulkOperation, BulkOperationResponse, SimilaritySearchRequest, 
    SubstructureSearchRequest, BULK_OPERATION_TYPE
)

# Library schemas
from .libraries import (
    LibraryRead, LibraryDetailRead, LibraryList, LibraryCreateRequest, 
    LibraryUpdateRequest, LibraryMoleculeOperationRequest, LibraryOperationResponse
)

# CSV processing schemas
from .csv import (
    CSVUploadResponse, PropertyMapping, CSVMappingRequest, CSVMappingResponse,
    CSVProcessingStatus, CSVImportSummary, CSVStatusResponse, SystemProperty,
    AvailablePropertiesResponse, CSV_PROCESSING_STATUS
)

# Experiment schemas
from .experiments import (
    ExperimentRead, ExperimentDetailRead, ExperimentList, ExperimentCreateRequest,
    ExperimentUpdateRequest, ExperimentMoleculeOperationRequest, ExperimentOperationResponse
)

# Submission schemas
from .submissions import (
    SubmissionRead, SubmissionDetailRead, SubmissionList, SubmissionCreateRequest,
    SubmissionUpdateRequest, SubmissionStatusUpdateRequest, QuoteRequest, QuoteResponse
)

# Result schemas
from .results import (
    ResultRead, ResultDetailRead, ResultList, ResultUploadRequest, ResultStatusUpdateRequest
)

# Admin schemas
from .admin import (
    AdminUserRead, AdminUserList, AdminUserCreateRequest, AdminUserUpdateRequest,
    SystemStatsResponse
)