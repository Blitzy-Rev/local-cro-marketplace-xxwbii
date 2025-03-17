"""
Celery application module for the Molecular Data Management and CRO Integration Platform.

This module initializes and configures the Celery application for background task processing,
providing asynchronous execution for resource-intensive operations like CSV processing,
molecular calculations, and notification delivery.
"""

# Standard library imports
import os
import logging

# External package imports
from celery import Celery  # celery 5.2+

# Internal imports
from . import celery_config
from ..logging_config import setup_logging

# Configure module logger
logger = logging.getLogger(__name__)

# Initialize Celery application with configuration from celery_config
app = Celery('molecular-platform', config_source=celery_config)

def init_celery():
    """
    Initializes the Celery application with task modules and logging configuration.
    
    This function sets up logging for the Celery workers, ensures task modules
    are imported to register tasks with Celery, and returns the configured
    Celery application instance.
    
    Returns:
        celery.Celery: Configured Celery application instance
    """
    # Configure logging
    setup_logging()
    logger.info("Initializing Celery application for background task processing")
    
    # Note: Task modules should be imported at application startup to register tasks
    # with the Celery app. This typically happens in the main application entry point.
    # 
    # For example:
    #   - from app.worker.tasks import csv_tasks, molecule_tasks, notification_tasks
    #
    # These imports register the tasks with the Celery application without creating
    # circular import dependencies that would occur if imported here.
    
    return app

# If this module is run directly, initialize Celery
if __name__ == '__main__':
    init_celery()