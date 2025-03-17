from typing import Optional, Dict, Any, Generator
from fastapi import Depends, HTTPException, status, SecurityScopes
from sqlalchemy.orm import Session

from ...core.auth import (
    get_current_user, 
    get_current_active_user,
    check_permissions,
    RoleChecker
)
from ...constants import USER_ROLES
from ...db.session import get_db
from ...crud.crud_user import user
from ...schemas.token import TokenData
from ...exceptions import AuthenticationException, AuthorizationException


def get_db_session() -> Generator[Session, None, None]:
    """
    Dependency that provides a database session
    
    Yields:
        Generator[Session, None, None]: Database session
    """
    # Call get_db() to get a database session
    # Yield the session to the endpoint function
    # Session is automatically closed after the endpoint function completes
    with get_db() as db:
        yield db


def get_current_user_from_token(db: Session = Depends(get_db_session)) -> Dict[str, Any]:
    """
    Dependency that extracts and validates the current user from a JWT token
    
    Args:
        db: Database session
        
    Returns:
        Dict[str, Any]: Current user information
    """
    try:
        # Call get_current_user() to extract user from token
        # If token is invalid, AuthenticationException is raised
        return get_current_user(SecurityScopes())
    except AuthenticationException as e:
        raise e


def get_current_active_user_from_token(
    current_user: Dict[str, Any] = Depends(get_current_user_from_token)
) -> Dict[str, Any]:
    """
    Dependency that ensures the current user is active
    
    Args:
        current_user: Current user information
        
    Returns:
        Dict[str, Any]: Current active user information
    """
    try:
        # Call get_current_active_user() to verify user is active
        # If user is not active, HTTPException is raised
        db = next(get_db())
        return get_current_active_user(current_user, db)
    except AuthenticationException as e:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=str(e)
        )


def get_current_user_from_db(
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token),
    db: Session = Depends(get_db_session)
) -> Any:
    """
    Dependency that retrieves the full user object from the database
    
    Args:
        current_user: Current user information
        db: Database session
        
    Returns:
        Any: User database object
    """
    # Extract user_id from current_user token data
    user_id = current_user.get("user_id")
    # Query the database for the user with the given ID
    db_user = user.get(db, id=user_id)
    # If user not found, raise HTTPException with 404 status
    if db_user is None:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    # Return the user database object
    return db_user


def get_current_pharma_user(
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
) -> Dict[str, Any]:
    """
    Dependency that ensures the current user has pharma role
    
    Args:
        current_user: Current user information
        
    Returns:
        Dict[str, Any]: Current pharma user information
    """
    # Check if user has pharma role using RoleChecker.is_pharma()
    # If not a pharma user, raise HTTPException with 403 status
    if not RoleChecker.is_pharma(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User does not have pharma role"
        )
    # Return current user information
    return current_user


def get_current_cro_user(
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
) -> Dict[str, Any]:
    """
    Dependency that ensures the current user has CRO role
    
    Args:
        current_user: Current user information
        
    Returns:
        Dict[str, Any]: Current CRO user information
    """
    # Check if user has CRO role using RoleChecker.is_cro()
    # If not a CRO user, raise HTTPException with 403 status
    if not RoleChecker.is_cro(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User does not have CRO role"
        )
    # Return current user information
    return current_user


def get_current_admin_user(
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
) -> Dict[str, Any]:
    """
    Dependency that ensures the current user has admin role
    
    Args:
        current_user: Current user information
        
    Returns:
        Dict[str, Any]: Current admin user information
    """
    # Check if user has admin role using RoleChecker.is_admin()
    # If not an admin user, raise HTTPException with 403 status
    if not RoleChecker.is_admin(current_user):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN,
            detail="User does not have admin role"
        )
    # Return current user information
    return current_user


def check_user_permissions(required_permissions: list):
    """
    Dependency factory that creates a dependency to check user permissions
    
    Args:
        required_permissions: List of required permissions
        
    Returns:
        function: Dependency function that checks permissions
    """
    # Define inner dependency function that takes current_user
    def permission_dependency(
        current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
    ) -> Dict[str, Any]:
        # Check if user has all required permissions using check_permissions()
        if not check_permissions(current_user, required_permissions):
            # If user lacks permissions, raise HTTPException with 403 status
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail=f"User does not have required permissions: {required_permissions}"
            )
        # Return current user if permissions are valid
        return current_user
    
    # Return the inner dependency function
    return permission_dependency