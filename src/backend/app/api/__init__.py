# fastapi==0.95+
from fastapi import APIRouter

# Import the main API router from the v1 module
# src/backend/app/api/api_v1/api.py
from .api_v1.api import api_router

# Export the main API router for use in the FastAPI application
__all__ = ["api_router"]