"""
Logging configuration module for the Molecular Data Management and CRO Integration Platform.

This module configures the application's logging system with structured JSON logging,
log rotation, and different log levels based on the application environment.
It provides centralized logging configuration that ensures consistent log formats
and appropriate handling across different parts of the application.
"""

import logging  # standard library
import logging.handlers  # standard library
import os  # standard library
import sys  # standard library
import pathlib  # standard library
import json  # standard library
from pathlib import Path  # standard library

from .core.config import get_settings
from .constants import LOG_LEVELS, DEFAULT_LOG_LEVEL, DEFAULT_LOG_FORMAT

# Module logger
logger = logging.getLogger(__name__)

# Logging configuration constants
LOG_DIR = pathlib.Path('logs')
MAX_LOG_SIZE_BYTES = 10 * 1024 * 1024  # 10MB
BACKUP_COUNT = 5


class JsonFormatter(logging.Formatter):
    """
    Custom log formatter that outputs log records as JSON strings.
    
    This formatter creates structured JSON logs with consistent fields:
    timestamp, service_name, level, message, and additional context
    from the log record.
    """
    
    def __init__(self):
        """Initializes the JsonFormatter with no format string (uses JSON instead)."""
        super().__init__()
        
    def format(self, record):
        """
        Formats the log record as a JSON string.
        
        Args:
            record: The log record to format
            
        Returns:
            str: JSON-formatted log message
        """
        log_object = {
            "timestamp": self.formatTime(record, "%Y-%m-%d %H:%M:%S,%f")[:-3],
            "service_name": record.name,
            "level": record.levelname,
            "message": record.getMessage(),
            "path": record.pathname,
            "line": record.lineno,
            "function": record.funcName
        }
        
        # Add exception info if present
        if record.exc_info:
            log_object["exception"] = self.formatException(record.exc_info)
            
        # Add correlation_id, user_id, and request_path if available
        if hasattr(record, 'correlation_id'):
            log_object['correlation_id'] = record.correlation_id
            
        if hasattr(record, 'user_id'):
            log_object['user_id'] = record.user_id
            
        if hasattr(record, 'request_path'):
            log_object['request_path'] = record.request_path
            
        # Add any extra fields from the record args if it's a dict
        if isinstance(record.args, dict):
            for key, value in record.args.items():
                if key not in log_object:
                    log_object[key] = value
                    
        return json.dumps(log_object)


def get_log_level():
    """
    Determines the appropriate log level based on application settings.
    
    Returns:
        int: Logging level constant from the logging module
    """
    settings = get_settings()
    log_level_str = settings.LOG_LEVEL
    
    # Convert string log level to logging module constant
    try:
        return getattr(logging, log_level_str)
    except (AttributeError, TypeError):
        logger.warning(
            f"Invalid log level: {log_level_str}. Defaulting to INFO.",
            {"invalid_level": log_level_str}
        )
        return logging.INFO


def create_console_handler(log_level, log_format):
    """
    Creates and configures a console logging handler.
    
    Args:
        log_level: The log level for the handler
        log_format: The format string for log messages
        
    Returns:
        logging.StreamHandler: Configured console logging handler
    """
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    
    if log_format.lower() == 'json':
        formatter = JsonFormatter()
    else:
        formatter = logging.Formatter(log_format)
        
    console_handler.setFormatter(formatter)
    return console_handler


def create_file_handler(log_level, log_format):
    """
    Creates and configures a rotating file logging handler.
    
    Args:
        log_level: The log level for the handler
        log_format: The format string for log messages
        
    Returns:
        logging.handlers.RotatingFileHandler: Configured file logging handler
    """
    # Ensure the logs directory exists
    if not LOG_DIR.exists():
        LOG_DIR.mkdir(parents=True, exist_ok=True)
        
    log_file = LOG_DIR / 'application.log'
    file_handler = logging.handlers.RotatingFileHandler(
        log_file,
        maxBytes=MAX_LOG_SIZE_BYTES,
        backupCount=BACKUP_COUNT
    )
    file_handler.setLevel(log_level)
    
    if log_format.lower() == 'json':
        formatter = JsonFormatter()
    else:
        formatter = logging.Formatter(log_format)
        
    file_handler.setFormatter(formatter)
    return file_handler


def setup_logging():
    """
    Configures the logging system for the application.
    
    This function sets up logging handlers for both console and file output,
    configures formatters, and sets the appropriate log level based on
    application settings.
    
    Returns:
        None
    """
    settings = get_settings()
    log_level = get_log_level()
    log_format = settings.LOG_FORMAT
    
    # Ensure logs directory exists
    if not LOG_DIR.exists():
        LOG_DIR.mkdir(parents=True, exist_ok=True)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)
    
    # Remove existing handlers if any
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Create and add console handler
    console_handler = create_console_handler(log_level, log_format)
    root_logger.addHandler(console_handler)
    
    # Create and add file handler
    file_handler = create_file_handler(log_level, log_format)
    root_logger.addHandler(file_handler)
    
    logger.info(
        f"Logging configured with level: {logging.getLevelName(log_level)}",
        {
            "log_level": logging.getLevelName(log_level),
            "log_format": log_format,
            "log_file": str(LOG_DIR / 'application.log')
        }
    )