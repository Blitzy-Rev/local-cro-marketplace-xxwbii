"""
Celery configuration for the Molecular Data Management and CRO Integration Platform.

This module defines all settings for the Celery task queue system, including broker
connection settings, task routing rules, worker configurations, and other Celery-specific
settings to optimize background task processing for CSV imports, molecular calculations,
notifications, and report generation.
"""

import os  # standard library
from kombu import Queue  # kombu 5.2+

from ..core.config import get_settings, get_redis_settings

# Get application settings
settings = get_settings()
redis_settings = get_redis_settings()

# Redis broker and result backend settings
broker_url = f"redis://{redis_settings['host']}:{redis_settings['port']}/{redis_settings['db']}"
result_backend = f"redis://{redis_settings['host']}:{redis_settings['port']}/{redis_settings['db']}"

# Task serialization settings
task_serializer = "json"
accept_content = ["json"]
result_serializer = "json"
enable_utc = True

# Worker settings
worker_hijack_root_logger = False  # Don't hijack root logger
worker_prefetch_multiplier = 4  # Number of tasks to prefetch per worker process
task_acks_late = True  # Acknowledge tasks after execution
task_reject_on_worker_lost = True  # Reject tasks when worker connection is lost

# Queue definitions
task_queues = [
    Queue('default'),
    Queue('csv-processing'),  # Queue for CSV import and processing tasks
    Queue('molecular-tasks'),  # Queue for molecular calculations and processing
    Queue('notifications'),    # Queue for notification tasks
    Queue('results'),          # Queue for result processing tasks
    Queue('report-generation') # Queue for report generation tasks
]

# Default queue
task_default_queue = 'default'

# Task routing based on module names
task_routes = {
    'csv_tasks.*': {'queue': 'csv-processing'},
    'molecule_tasks.*': {'queue': 'molecular-tasks'},
    'notification_tasks.*': {'queue': 'notifications'},
    'result_tasks.*': {'queue': 'results'},
    'report_tasks.*': {'queue': 'report-generation'}
}

# Task execution time limits (in seconds)
task_time_limit = 3600  # Hard limit: 1 hour
task_soft_time_limit = 3300  # Soft limit: 55 minutes (allows for cleanup)

# Worker concurrency (number of parallel processes)
worker_concurrency = os.cpu_count() or 4

# Worker task limits
worker_max_tasks_per_child = 1000  # Restart worker process after processing this many tasks

# Broker connection settings
broker_connection_retry = True  # Retry connecting to broker if connection fails
broker_connection_retry_on_startup = True  # Retry connecting to broker on startup
broker_connection_max_retries = 10  # Maximum number of retry attempts
broker_pool_limit = 10  # Connection pool size

# Task creation settings
task_create_missing_queues = True  # Create queues if they don't exist
task_default_rate_limit = "100/m"  # Default rate limit for tasks

# Task tracking settings
task_track_started = True  # Track when tasks are started
task_send_sent_event = True  # Send event when task is sent
worker_send_task_events = True  # Send task-related events

# Result settings
task_ignore_result = False  # Store task results
result_expires = 86400  # Results expire after 1 day (in seconds)

# Scheduled tasks (Celery Beat)
beat_schedule = {
    # Cleanup old notifications (run daily, delete notifications older than 90 days)
    'cleanup_old_notifications': {
        'task': 'notification_tasks.cleanup_old_notifications',
        'schedule': 86400,  # Once per day (in seconds)
        'args': (90,)  # Keep notifications for 90 days
    },
    # Cleanup temporary CSV files (run hourly)
    'cleanup_temporary_files': {
        'task': 'csv_tasks.cleanup_csv_files',
        'schedule': 3600,  # Once per hour (in seconds)
        'args': ([],)  # No specific arguments
    }
}