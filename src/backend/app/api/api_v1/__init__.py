# fastapi==0.95+
from fastapi import APIRouter

# Import the main API router for version 1
from .api import api_router

# Import all dependency functions for API endpoints
from .deps import (
    get_db_session,
    get_current_user_from_token,
    get_current_active_user_from_token,
    get_current_user_from_db,
    get_current_pharma_user,
    get_current_cro_user,
    get_current_admin_user,
    check_user_permissions
)

# Export the main API router for use in the FastAPI application
__all__ = [
    "api_router",
    "get_db_session",
    "get_current_user_from_token",
    "get_current_active_user_from_token",
    "get_current_user_from_db",
    "get_current_pharma_user",
    "get_current_cro_user",
    "get_current_admin_user",
    "check_user_permissions"
]