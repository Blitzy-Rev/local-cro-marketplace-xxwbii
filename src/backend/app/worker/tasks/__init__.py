"""
Initialization module for Celery tasks in the Molecular Data Management and CRO Integration Platform.

This module imports and exposes all task modules for background processing of resource-intensive
operations such as CSV processing, molecular calculations, notification delivery, result processing,
and report generation.
"""

# Standard library imports
import logging  # standard library

# Internal imports
from ..celery_app import app
from ...logging_config import logger

# Configure module logger
logger = logging.getLogger(__name__)

# Import task modules and expose their functions
from .csv_tasks import process_csv, cleanup_csv_files, cancel_csv_processing, retry_csv_processing, get_csv_processing_stats
from .molecule_tasks import calculate_molecule_properties, batch_process_molecules, stream_process_molecules, calculate_batch_properties, calculate_similarity, cluster_by_similarity, select_diverse_molecules, standardize_molecule_set, check_molecules_drug_likeness, bulk_update_molecule_properties
from .notification_tasks import deliver_notification, batch_deliver_notifications, mark_notification_read, cleanup_old_notifications, send_experiment_status_notification, send_submission_notification, send_quote_notification, send_results_notification, send_system_notification
from .result_tasks import process_result_file, update_result_status, batch_process_result_files, generate_result_report, cleanup_temporary_result_files
from .report_tasks import generate_experiment_report, generate_result_report, generate_molecule_comparison_report, generate_library_report, generate_activity_report