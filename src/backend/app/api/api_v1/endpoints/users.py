from typing import Dict, Any, List, Optional
from fastapi import APIRouter, Depends, HTTPException, status, Query
from sqlalchemy.orm import Session

from ..deps import get_db_session, get_current_user_from_db, get_current_admin_user
from ...crud.crud_user import user
from ..schemas.users import (
    UserProfileResponse, 
    UserUpdateRequest, 
    UserFilterParams,
    UserListResponse,
    UserStatusUpdateRequest,
    ChangePasswordRequest
)
from ..schemas.auth import MessageResponse
from ...services.auth_service import AuthService
from ...logging_config import logger

router = APIRouter(prefix='/users', tags=['users'])


@router.get('/me', response_model=UserProfileResponse, status_code=status.HTTP_200_OK)
def get_current_user_profile(current_user = Depends(get_current_user_from_db)):
    """
    Get the current authenticated user's profile information.
    
    Returns the profile information of the currently authenticated user.
    """
    return current_user


@router.put('/me', response_model=UserProfileResponse, status_code=status.HTTP_200_OK)
def update_current_user_profile(
    user_data: UserUpdateRequest,
    current_user = Depends(get_current_user_from_db),
    db: Session = Depends(get_db_session)
):
    """
    Update the current authenticated user's profile information.
    
    Allows a user to update their own profile information including email.
    If email is being changed, verifies that it's not already in use.
    """
    # Check if email is being updated and is different
    if user_data.email and user_data.email != current_user.email:
        # Check if email is already in use
        existing_user = user.get_by_email(db, email=user_data.email)
        if existing_user:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Email already registered"
            )
    
    # Update user with provided data
    updated_user = user.update(db, db_obj=current_user, obj_in=user_data.dict(exclude_unset=True))
    logger.info(f"User {updated_user.id} updated their profile", {"user_id": updated_user.id})
    
    return updated_user


@router.post('/me/change-password', response_model=MessageResponse, status_code=status.HTTP_200_OK)
def change_password(
    password_data: ChangePasswordRequest,
    current_user = Depends(get_current_user_from_db),
    db: Session = Depends(get_db_session)
):
    """
    Change the current authenticated user's password.
    
    Validates the current password before allowing the password to be changed.
    The new password must meet the system's password complexity requirements.
    """
    auth_service = AuthService()
    try:
        auth_service.change_password(
            db,
            current_user.id, 
            password_data.current_password, 
            password_data.new_password
        )
        logger.info(f"User {current_user.id} changed their password successfully", {"user_id": current_user.id})
        return {"success": True, "message": "Password changed successfully"}
    except Exception as e:
        logger.warning(f"Password change failed for user {current_user.id}: {str(e)}", 
                      {"user_id": current_user.id, "error": str(e)})
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Failed to change password. Please check your current password."
        )


@router.get('/{user_id}', response_model=UserProfileResponse, status_code=status.HTTP_200_OK)
def get_user_by_id(
    user_id: int,
    current_user = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """
    Get a user by ID (admin only).
    
    Retrieves detailed user information by ID. This endpoint is restricted to admin users only.
    """
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    
    return db_user


@router.put('/{user_id}', response_model=UserProfileResponse, status_code=status.HTTP_200_OK)
def update_user(
    user_id: int,
    user_data: UserUpdateRequest,
    current_user = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """
    Update a user's information (admin only).
    
    Allows an admin to update any user's information including email and role.
    If email is being changed, verifies that it's not already in use by another user.
    """
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    
    # Check if email is being updated and is different
    if user_data.email and user_data.email != db_user.email:
        # Check if email is already in use by another user
        existing_user = user.get_by_email(db, email=user_data.email)
        if existing_user and existing_user.id != user_id:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Email already registered"
            )
    
    # Update user with provided data
    updated_user = user.update(db, db_obj=db_user, obj_in=user_data.dict(exclude_unset=True))
    logger.info(f"Admin {current_user.id} updated user {user_id}", 
               {"admin_id": current_user.id, "user_id": user_id})
    
    return updated_user


@router.put('/{user_id}/status', response_model=UserProfileResponse, status_code=status.HTTP_200_OK)
def update_user_status(
    user_id: int,
    status_data: UserStatusUpdateRequest,
    current_user = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """
    Update a user's status (admin only).
    
    Allows an admin to update a user's status (active, inactive, pending, locked).
    The status change also affects the user's ability to login.
    """
    db_user = user.get(db, id=user_id)
    if not db_user:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail="User not found"
        )
    
    # Update user status
    update_data = {"status": status_data.status}
    
    # Set is_active based on status
    if status_data.status == 'active':
        update_data["is_active"] = True
    elif status_data.status in ['inactive', 'pending', 'locked']:
        update_data["is_active"] = False
    
    updated_user = user.update(db, db_obj=db_user, obj_in=update_data)
    logger.info(f"Admin {current_user.id} updated user {user_id} status to {status_data.status}", 
               {"admin_id": current_user.id, "user_id": user_id, "new_status": status_data.status})
    
    return updated_user


@router.get('/', response_model=UserListResponse, status_code=status.HTTP_200_OK)
def get_users(
    filters: UserFilterParams = Depends(),
    page: int = Query(1, ge=1),
    size: int = Query(50, ge=1, le=100),
    current_user = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """
    Get a paginated list of users with optional filtering (admin only).
    
    Retrieves a paginated list of users with optional filtering by role, status, etc.
    This endpoint is restricted to admin users only.
    """
    skip = (page - 1) * size
    
    # Apply filters based on role if provided
    if filters.role:
        users = user.get_users_by_role(db, filters.role, skip=skip, limit=size)
        total = user.count_by_role(db, role=filters.role)
    else:
        # Otherwise get all users with pagination
        users = user.get_multi(db, skip=skip, limit=size)
        total = user.count(db)
    
    return {
        "items": users,
        "total": total,
        "page": page,
        "size": size
    }