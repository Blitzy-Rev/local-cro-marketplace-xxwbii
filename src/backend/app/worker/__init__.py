"""
Initializes the worker module for the Molecular Data Management and CRO Integration Platform. This module serves as the entry point for the Celery-based background task processing system, exposing the Celery application instance and task modules for asynchronous operations like CSV processing, molecular calculations, and notifications.
"""

# Standard library imports
import logging  # Configure logging for the worker module

# Internal imports
from .celery_app import app, init_celery  # Import the Celery application instance
from .tasks import notification_tasks  # Import notification task module
from .tasks import csv_tasks  # Import CSV processing task module
from .tasks import molecule_tasks  # Import molecule processing task module
from .tasks import result_tasks  # Import result processing task module
from .tasks import report_tasks  # Import report generation task module

# Configure module logger
logger = logging.getLogger(__name__)

# Initialize Celery application
celery_app = init_celery()

# Export the Celery application instance for task definition and execution
__all__ = ['app', 'notification_tasks', 'csv_tasks', 'molecule_tasks', 'result_tasks', 'report_tasks']