# fastapi==0.95+
from fastapi import APIRouter

# Import health check endpoints router
# src/backend/app/api/api_v1/endpoints/health.py
from .endpoints.health import health_router

# Import authentication endpoints router
# src/backend/app/api/api_v1/endpoints/auth.py
from .endpoints.auth import router as auth_router

# Import user management endpoints router
# src/backend/app/api/api_v1/endpoints/users.py
from .endpoints.users import router as users_router

# Import molecule management endpoints router
# src/backend/app/api/api_v1/endpoints/molecules.py
from .endpoints.molecules import router as molecules_router

# Import library management endpoints router
# src/backend/app/api/api_v1/endpoints/libraries.py
from .endpoints.libraries import router as libraries_router

# Import CSV upload and processing endpoints router
# src/backend/app/api/api_v1/endpoints/csv.py
from .endpoints.csv import router as csv_router

# Import experiment management endpoints router
# src/backend/app/api/api_v1/endpoints/experiments.py
from .endpoints.experiments import router as experiments_router

# Import submission management endpoints router
# src/backend/app/api/api_v1/endpoints/submissions.py
from .endpoints.submissions import router as submissions_router

# Import result management endpoints router
# src/backend/app/api/api_v1/endpoints/results.py
from .endpoints.results import router as results_router

# Import admin endpoints router
# src/backend/app/api/api_v1/endpoints/admin.py
from .endpoints.admin import router as admin_router

# Create main API router with prefix
api_router = APIRouter(prefix='/api/v1')

# Include routers for different endpoints
api_router.include_router(health_router)
api_router.include_router(auth_router)
api_router.include_router(users_router)
api_router.include_router(molecules_router)
api_router.include_router(libraries_router)
api_router.include_router(csv_router)
api_router.include_router(experiments_router)
api_router.include_router(submissions_router)
api_router.include_router(results_router)
api_router.include_router(admin_router)