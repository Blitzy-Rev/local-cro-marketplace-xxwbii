"""
Main entry point for the FastAPI backend application of the Molecular Data Management and CRO Integration Platform. This file initializes the application, configures middleware, sets up database connections, and defines API routes. It implements a fully local deployment model without external dependencies as required by the technical specifications.
"""

import logging  # -- standard library
import os  # -- standard library
import sys  # -- standard library

from fastapi import FastAPI  # -- 0.95+
from starlette.middleware.cors import CORSMiddleware  # -- 0.27+
from starlette.middleware.trustedhost import TrustedHostMiddleware  # -- 0.27+
from starlette.middleware.gzip import GZipMiddleware  # -- 0.27+
import uvicorn  # -- 0.22+

from app import init_app  # src/backend/app/__init__.py
from app.core.config import get_settings  # src/backend/app/core/config.py
from app.api.api_v1.api import api_router  # src/backend/app/api/api_v1/api.py
from app.db.session import init_db  # src/backend/app/db/session.py
from app.db.init_db import init_db_data  # src/backend/app/db/init_db.py
from app.logging_config import setup_logging  # src/backend/app/logging_config.py

# Initialize logger
logger = logging.getLogger(__name__)

# Create FastAPI application instance
app = FastAPI(
    title=get_settings().PROJECT_NAME,
    openapi_url=f"{get_settings().API_V1_PREFIX}/openapi.json",
    docs_url=f"{get_settings().API_V1_PREFIX}/docs",
    redoc_url=f"{get_settings().API_V1_PREFIX}/redoc"
)

def configure_middleware():
    """
    Configures middleware for the FastAPI application
    """
    # Add CORS middleware with settings from configuration
    app.add_middleware(
        CORSMiddleware,
        allow_origins=get_settings().get_cors_origins,
        allow_credentials=get_settings().CORS_CREDENTIALS,
        allow_methods=get_settings().get_cors_methods,
        allow_headers=get_settings().get_cors_headers,
    )

    # Add TrustedHost middleware for security
    app.add_middleware(TrustedHostMiddleware, allowed_hosts=["*"])

    # Add GZip middleware for response compression
    app.add_middleware(GZipMiddleware, compresslevel=5)

    # Add custom rate limiting middleware
    # TODO: Implement rate limiting middleware

    # Add authentication middleware
    # TODO: Implement authentication middleware

def configure_routers():
    """
    Configures API routers for the FastAPI application
    """
    # Include the main API router with prefix from settings
    app.include_router(api_router)

    # Add root endpoint for health check and API information
    @app.get("/")
    def root_endpoint():
        """Root endpoint for API health check and information"""
        return {"message": "API is running", "version": "1.0.0"}

def configure_exception_handlers():
    """
    Configures global exception handlers for the FastAPI application
    """
    # Add handler for validation exceptions
    # TODO: Implement validation exception handler

    # Add handler for authentication exceptions
    # TODO: Implement authentication exception handler

    # Add handler for authorization exceptions
    # TODO: Implement authorization exception handler

    # Add handler for not found exceptions
    # TODO: Implement not found exception handler

    # Add handler for internal server errors
    # TODO: Implement internal server error handler

def configure_events():
    """
    Configures startup and shutdown event handlers for the FastAPI application
    """
    # Add startup event handler to initialize database
    @app.on_event("startup")
    def startup_db_init():
        """Initializes database connection and schema on application startup"""
        logger.info("Initializing database...")
        init_db()
        logger.info("Database initialized")

    # Add startup event handler to initialize database data
    @app.on_event("startup")
    def startup_db_data_init():
        """Initializes database with seed data on application startup"""
        logger.info("Initializing database data...")
        init_db_data()
        logger.info("Database data initialized")

    # Add shutdown event handler to close database connections
    @app.on_event("shutdown")
    def shutdown_db():
        """Closes database connections on application shutdown"""
        logger.info("Shutting down database...")
        # TODO: Close database connections
        logger.info("Database shutdown")

def initialize_application():
    """
    Initializes the FastAPI application with all configurations
    """
    # Set up application logging
    setup_logging()

    # Initialize application components
    init_app()

    # Configure middleware
    configure_middleware()

    # Configure routers
    configure_routers()

    # Configure exception handlers
    configure_exception_handlers()

    # Configure events
    configure_events()

    return app

if __name__ == "__main__":
    """
    Main entry point for running the application directly
    """
    # Check if script is being run directly
    if __name__ == "__main__":
        # Get application settings
        settings = get_settings()

        # Run the application with uvicorn
        uvicorn.run(
            "src.backend.main:app",
            host=settings.SERVER_HOST,
            port=settings.SERVER_PORT,
            reload=settings.DEBUG,
            workers=settings.SERVER_WORKERS
        )