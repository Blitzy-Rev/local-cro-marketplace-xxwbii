# Import routers from their respective endpoint modules
from .health import router as health_router
from .auth import router as auth_router
from .users import router as users_router
from .molecules import router as molecules_router
from .libraries import router as libraries_router
from .csv import router as csv_router
from .experiments import router as experiments_router
from .submissions import router as submissions_router
from .results import router as results_router
from .admin import router as admin_router

# Export all routers for inclusion in the main API router
__all__ = [
    "health_router",
    "auth_router",
    "users_router",
    "molecules_router",
    "libraries_router",
    "csv_router",
    "experiments_router",
    "submissions_router",
    "results_router",
    "admin_router"
]