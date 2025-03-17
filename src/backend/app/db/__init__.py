"""
Database package initialization for the Molecular Data Management and CRO Integration Platform.

This module serves as the entry point for the database package, providing centralized access
to database components and utilities. It supports the platform's fully local deployment
requirement by offering comprehensive database functionality without external dependencies.

The module exports essential database components including the SQLAlchemy Base class,
session management functions, and initialization utilities used throughout the application.

These components enable:
- Database model definition and relationship management
- Database connection pooling and session handling
- Database initialization and seed data creation
- Transaction management for data consistency
"""

# Import SQLAlchemy base class for model definitions
from .base_class import Base

# Import session management and database initialization utilities
from .session import get_engine, init_db, get_db, SessionLocal

# Import function to populate initial database data
from .init_db import init_db_data

# Export all necessary database components
__all__ = [
    "Base",
    "get_engine",
    "init_db",
    "get_db",
    "SessionLocal",
    "init_db_data"
]