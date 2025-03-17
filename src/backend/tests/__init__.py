"""
Backend Test Package for the Molecular Data Management and CRO Integration Platform

This package contains all tests for the backend components, including unit tests,
integration tests, and functional tests. It provides the foundation for test discovery
and execution within the testing framework.
"""

import logging  # standard library

from ..app.constants import VERSION

# Test package metadata
__test_package__ = True
__version__ = VERSION

# Configure logger for test package
logger = logging.getLogger(__name__)