"""
Worker Tests Package for the Molecular Data Management and CRO Integration Platform

This package contains tests for the background worker tasks and asynchronous
processing components of the platform. It includes unit tests for individual
task functions, integration tests for worker queues, and tests for error
handling and recovery mechanisms.
"""

import logging  # standard library

from .. import logger

# Configure worker-specific test logger
worker_test_logger = logging.getLogger(__name__)