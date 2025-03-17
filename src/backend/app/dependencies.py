"""
Provides centralized dependency injection functions for the Molecular Data Management and CRO Integration Platform.

This module defines reusable dependencies for database sessions, authentication, authorization, and service instances
that can be used throughout the application with FastAPI's dependency injection system.
"""

from typing import Dict, List, Any, Optional

from fastapi import Depends, HTTPException, status
from sqlalchemy.orm import Session

from .db.session import get_db
from .core.jwt import get_current_user
from .core.auth import check_permissions
from .constants import UserRole, UserStatus
from .exceptions import AuthenticationException, AuthorizationException
from .services.file_storage_service import FileStorageService
from .services.auth_service import AuthService

# Service singletons
_file_storage_service = None
_auth_service = None
_molecule_service = None
_library_service = None
_experiment_service = None
_submission_service = None
_result_service = None
_notification_service = None
_csv_service = None


def get_db_session():
    """
    Provides a database session dependency.
    
    Returns:
        Session: SQLAlchemy database session
    """
    yield from get_db()


def get_current_active_user(current_user: Dict[str, Any] = Depends(get_current_user)):
    """
    Ensures the current user is active.
    
    Args:
        current_user: Current user information from token
        
    Returns:
        Dict[str, Any]: Current active user information
        
    Raises:
        HTTPException: If user is not active
    """
    if current_user.get("status") != UserStatus.ACTIVE.name:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail={"message": "Inactive user", "details": None},
        )
    return current_user


def get_current_pharma_user(current_user: Dict[str, Any] = Depends(get_current_active_user)):
    """
    Ensures the current user has pharma role.
    
    Args:
        current_user: Current user information from token
        
    Returns:
        Dict[str, Any]: Current pharma user information
        
    Raises:
        HTTPException: If user does not have pharma role
    """
    if current_user.get("role") != UserRole.PHARMA.name:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail={"message": "Pharma role required", "details": None},
        )
    return current_user


def get_current_cro_user(current_user: Dict[str, Any] = Depends(get_current_active_user)):
    """
    Ensures the current user has CRO role.
    
    Args:
        current_user: Current user information from token
        
    Returns:
        Dict[str, Any]: Current CRO user information
        
    Raises:
        HTTPException: If user does not have CRO role
    """
    if current_user.get("role") != UserRole.CRO.name:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail={"message": "CRO role required", "details": None},
        )
    return current_user


def get_current_admin_user(current_user: Dict[str, Any] = Depends(get_current_active_user)):
    """
    Ensures the current user has admin role.
    
    Args:
        current_user: Current user information from token
        
    Returns:
        Dict[str, Any]: Current admin user information
        
    Raises:
        HTTPException: If user does not have admin role
    """
    if current_user.get("role") != UserRole.ADMIN.name:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail={"message": "Admin role required", "details": None},
        )
    return current_user


def require_permissions(required_permissions: List[str]):
    """
    Dependency factory that requires specific permissions for API endpoints.
    
    Args:
        required_permissions: List of required permission strings
        
    Returns:
        function: Dependency function that checks permissions
    """
    def dependency(current_user: Dict[str, Any] = Depends(get_current_active_user)):
        has_permissions = check_permissions(current_user, required_permissions)
        if not has_permissions:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail={
                    "message": "Not enough permissions",
                    "details": {"required": required_permissions, "user_role": current_user.get("role")}
                },
            )
        return current_user
    return dependency


def get_auth_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the authentication service.
    
    Args:
        db: Database session
        
    Returns:
        AuthService: Authentication service instance
    """
    global _auth_service
    if _auth_service is None:
        _auth_service = AuthService()
    return _auth_service


def get_file_storage_service():
    """
    Provides an instance of the file storage service.
    
    Returns:
        FileStorageService: File storage service instance
    """
    global _file_storage_service
    if _file_storage_service is None:
        _file_storage_service = FileStorageService()
    return _file_storage_service


def get_molecule_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the molecule service.
    
    Args:
        db: Database session
        
    Returns:
        MoleculeService: Molecule service instance
    """
    global _molecule_service
    if _molecule_service is None:
        from .services.molecule_service import MoleculeService
        _molecule_service = MoleculeService()
    return _molecule_service


def get_library_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the library service.
    
    Args:
        db: Database session
        
    Returns:
        LibraryService: Library service instance
    """
    global _library_service
    if _library_service is None:
        from .services.library_service import LibraryService
        _library_service = LibraryService()
    return _library_service


def get_experiment_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the experiment service.
    
    Args:
        db: Database session
        
    Returns:
        ExperimentService: Experiment service instance
    """
    global _experiment_service
    if _experiment_service is None:
        from .services.experiment_service import ExperimentService
        _experiment_service = ExperimentService()
    return _experiment_service


def get_submission_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the submission service.
    
    Args:
        db: Database session
        
    Returns:
        SubmissionService: Submission service instance
    """
    global _submission_service
    if _submission_service is None:
        from .services.submission_service import SubmissionService
        _submission_service = SubmissionService()
    return _submission_service


def get_result_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the result service.
    
    Args:
        db: Database session
        
    Returns:
        ResultService: Result service instance
    """
    global _result_service
    if _result_service is None:
        from .services.result_service import ResultService
        _result_service = ResultService()
    return _result_service


def get_notification_service(db: Session = Depends(get_db_session)):
    """
    Provides an instance of the notification service.
    
    Args:
        db: Database session
        
    Returns:
        NotificationService: Notification service instance
    """
    global _notification_service
    if _notification_service is None:
        from .services.notification_service import NotificationService
        _notification_service = NotificationService()
    return _notification_service


def get_csv_service(
    db: Session = Depends(get_db_session),
    molecule_service = Depends(get_molecule_service)
):
    """
    Provides an instance of the CSV service.
    
    Args:
        db: Database session
        molecule_service: Molecule service instance
        
    Returns:
        CSVService: CSV service instance
    """
    global _csv_service
    if _csv_service is None:
        from .services.csv_service import CSVService
        _csv_service = CSVService(molecule_service)
    return _csv_service