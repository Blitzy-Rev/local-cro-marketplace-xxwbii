from datetime import datetime
from typing import List, Dict, Optional, Any

from fastapi import APIRouter, Depends, HTTPException, status, Query, Path
from sqlalchemy.orm import Session

from ..deps import get_db_session, get_current_admin_user
from ...services.admin_service import AdminService
from ...services.file_storage_service import FileStorageService
from ..schemas.admin import (
    UserCreateRequest, UserAdminUpdateRequest, UserPasswordResetRequest,
    SystemStatsResponse, SystemResourcesResponse,
    ActivityLogResponse, ActivityLogFilterParams, ActivityLogListResponse,
    SystemAlertResponse, SystemAlertListResponse, SystemAlertUpdateRequest
)
from ...schemas.user import UserResponse
from ..schemas.auth import MessageResponse
from ...exceptions import AdminException

admin_router = APIRouter(prefix='/admin', tags=['admin'])


@admin_router.get('/stats', response_model=SystemStatsResponse)
def get_system_stats(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Get system statistics including user counts, data counts, and resource usage"""
    admin_service = AdminService(db, FileStorageService(db))
    stats = admin_service.get_system_stats()
    return stats


@admin_router.get('/resources', response_model=SystemResourcesResponse)
def get_system_resources(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Get detailed system resource usage information"""
    admin_service = AdminService(db, FileStorageService(db))
    resources = admin_service.get_system_resources()
    return resources


@admin_router.get('/health', response_model=Dict[str, Any])
def check_system_health(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Perform a comprehensive system health check"""
    admin_service = AdminService(db, FileStorageService(db))
    health = admin_service.check_system_health()
    return health


@admin_router.get('/users', response_model=List[UserResponse])
def get_users(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
    role: Optional[str] = Query(None, description='Filter by user role'),
    status: Optional[str] = Query(None, description='Filter by user status'),
    skip: int = Query(0, ge=0, description='Number of users to skip'),
    limit: int = Query(100, ge=1, le=100, description='Maximum number of users to return')
):
    """Get list of all users with optional filtering"""
    admin_service = AdminService(db, FileStorageService(db))
    users = admin_service.get_users(role=role, status=status, skip=skip, limit=limit)
    return users


@admin_router.get('/users/{user_id}', response_model=UserResponse)
def get_user(
    user_id: int = Path(..., gt=0, description='The ID of the user to get'),
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Get user by ID"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        user = admin_service.get_user(user_id)
        return user
    except AdminException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


@admin_router.post('/users', response_model=UserResponse, status_code=status.HTTP_201_CREATED)
def create_user(
    user_data: UserCreateRequest,
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Create a new user"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        new_user = admin_service.create_user(user_data.dict())
        return new_user
    except AdminException as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@admin_router.put('/users/{user_id}', response_model=UserResponse)
def update_user(
    user_id: int = Path(..., gt=0, description='The ID of the user to update'),
    user_data: UserAdminUpdateRequest,
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Update user information"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        updated_user = admin_service.update_user(user_id, user_data.dict(exclude_unset=True))
        return updated_user
    except AdminException as e:
        if "not found" in str(e).lower():
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
        elif "already exists" in str(e).lower():
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))


@admin_router.delete('/users/{user_id}', response_model=MessageResponse)
def delete_user(
    user_id: int = Path(..., gt=0, description='The ID of the user to delete'),
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Delete a user"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        admin_service.delete_user(user_id)
        return {"success": True, "message": f"User with ID {user_id} deleted successfully"}
    except AdminException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


@admin_router.post('/users/{user_id}/reset-password', response_model=MessageResponse)
def reset_user_password(
    user_id: int = Path(..., gt=0, description='The ID of the user to reset password'),
    password_data: UserPasswordResetRequest,
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Reset a user's password"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        admin_service.reset_user_password(user_id, password_data.new_password)
        return {"success": True, "message": f"Password for user with ID {user_id} reset successfully"}
    except AdminException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))


@admin_router.get('/activity-logs', response_model=ActivityLogListResponse)
def get_activity_logs(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
    filter_params: ActivityLogFilterParams = Depends(),
    page: int = Query(1, ge=1, description='Page number'),
    size: int = Query(20, ge=1, le=100, description='Page size')
):
    """Get activity logs with optional filtering"""
    admin_service = AdminService(db, FileStorageService(db))
    logs = admin_service.get_activity_logs(
        filter_params=filter_params.dict(exclude_unset=True),
        page=page,
        size=size
    )
    return logs


@admin_router.get('/alerts', response_model=SystemAlertListResponse)
def get_system_alerts(
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session),
    severity: Optional[str] = Query(None, description='Filter by alert severity'),
    alert_type: Optional[str] = Query(None, description='Filter by alert type'),
    resolved: Optional[bool] = Query(None, description='Filter by resolved status')
):
    """Get system alerts with optional filtering"""
    admin_service = AdminService(db, FileStorageService(db))
    alerts = admin_service.get_system_alerts(
        severity=severity,
        alert_type=alert_type,
        resolved=resolved
    )
    return alerts


@admin_router.put('/alerts/{alert_id}', response_model=SystemAlertResponse)
def update_system_alert(
    alert_id: int = Path(..., gt=0, description='The ID of the alert to update'),
    alert_data: SystemAlertUpdateRequest,
    current_user: Dict[str, Any] = Depends(get_current_admin_user),
    db: Session = Depends(get_db_session)
):
    """Update system alert status"""
    admin_service = AdminService(db, FileStorageService(db))
    try:
        updated_alert = admin_service.update_system_alert(alert_id, alert_data.dict())
        return updated_alert
    except AdminException as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))