"""
Gunicorn WSGI server configuration for the Molecular Data Management and CRO Integration Platform.

This module contains settings for the Gunicorn production server that hosts the FastAPI
backend application. It configures worker processes, timeouts, server binding, 
logging, and more to ensure optimal performance and reliability.
"""

import multiprocessing  # standard library
import os  # standard library
import logging  # standard library

from app.core.config import get_settings
from app.logging_config import setup_logging

# Set up logger for this configuration module
logger = logging.getLogger('gunicorn.conf')

# Get application settings
settings = get_settings()

# Server socket binding
bind = f"{settings.SERVER_HOST}:{settings.SERVER_PORT}"

# Worker processes
# If SERVER_WORKERS is defined in settings, use it; otherwise calculate based on CPU cores
# Formula: (2 x num_cores) + 1, but capped at 8 to prevent excessive resource usage
workers = settings.SERVER_WORKERS if hasattr(settings, 'SERVER_WORKERS') else min(multiprocessing.cpu_count() * 2 + 1, 8)

# Worker class - using Uvicorn worker for ASGI support (required for FastAPI)
worker_class = "uvicorn.workers.UvicornWorker"

# Timeout for worker processes (in seconds)
# Set to 120 seconds to handle long-running operations like CSV processing
timeout = 120

# Keep-alive timeout for inactive clients (in seconds)
keepalive = 5

# Logging configuration
errorlog = "-"  # stderr
accesslog = "-"  # stdout
loglevel = settings.LOG_LEVEL.lower()
access_log_format = "%({X-Forwarded-For}i)s %(l)s %(u)s %(t)s \"%(r)s\" %(s)s %(b)s \"%(f)s\" \"%(a)s\" %(L)s"

# Worker temporary directory - use /dev/shm for faster temporary storage
worker_tmp_dir = "/dev/shm"

# Graceful timeout - wait this many seconds for workers to finish their current requests before forcefully killing them
graceful_timeout = 30

# Worker recycling to prevent memory leaks - restart workers after handling this many requests
max_requests = 1000
max_requests_jitter = 50  # Add randomness to max_requests to prevent all workers from restarting at the same time

# Preload application code before forking worker processes to improve performance
preload_app = True


def on_starting(server):
    """
    Handler called when Gunicorn is starting up.
    
    Args:
        server: Gunicorn server instance
    """
    logger.info(f"Starting Gunicorn server for {settings.PROJECT_NAME}")
    
    # Configure application logging
    setup_logging()
    
    logger.info(f"Server configuration: bind={bind}, workers={workers}, worker_class={worker_class}")


def post_fork(server, worker):
    """
    Handler called after a worker has been forked.
    
    Args:
        server: Gunicorn server instance
        worker: Worker instance
    """
    logger.info(f"Worker {worker.pid} forked and initialized")
    
    # Reset the random seed for each worker to ensure different random sequences
    import random
    random.seed()


def worker_exit(server, worker):
    """
    Handler called when a worker is exiting.
    
    Args:
        server: Gunicorn server instance
        worker: Worker instance
    """
    logger.info(f"Worker {worker.pid} exiting")
    
    # Perform any worker-specific cleanup if needed


def on_exit(server):
    """
    Handler called when Gunicorn is shutting down.
    
    Args:
        server: Gunicorn server instance
    """
    logger.info(f"Shutting down Gunicorn server for {settings.PROJECT_NAME}")
    
    # Perform any final cleanup operations if needed