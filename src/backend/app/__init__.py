"""
Molecular Data Management and CRO Integration Platform - Backend Package

This module serves as the entry point for the backend application package,
providing core functionality, constants, and utilities while setting up
version information and initialization mechanisms.

This package implements a modular architecture for the platform, enabling
fully local deployment without external dependencies by centralizing imports
and providing a clean interface to the application's components.
"""

import logging  # standard library

from .constants import VERSION, PROJECT_NAME
from .logging_config import setup_logging
from .core.config import get_settings, initialize_settings

# Expose version information
__version__ = VERSION
__app_name__ = PROJECT_NAME

# Module-level logger
logger = logging.getLogger(__name__)

def init_app():
    """
    Initializes the application by setting up logging and configuration.
    
    This function should be called once at application startup to ensure
    proper configuration of logging and application settings.
    
    Returns:
        None
    """
    # Set up logging
    setup_logging()
    
    # Initialize application settings
    initialize_settings()
    
    # Log application initialization
    logger.info(
        f"Initializing {__app_name__} v{__version__}",
        {"version": __version__, "app_name": __app_name__}
    )

# Define exported symbols
__all__ = ["__version__", "__app_name__", "init_app", "get_settings", "setup_logging"]