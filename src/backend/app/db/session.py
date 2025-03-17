"""
Manages database session creation and connection pooling for the Molecular Data Management and CRO Integration Platform.

This module provides a SQLAlchemy engine and session factory with optimized connection pooling settings,
ensuring efficient database access while supporting the requirement for fully local deployment.
"""

from contextlib import contextmanager
from typing import Generator, Any

from sqlalchemy import create_engine  # sqlalchemy v2.0+
from sqlalchemy.engine import Engine  # sqlalchemy v2.0+
from sqlalchemy.orm import sessionmaker, Session  # sqlalchemy v2.0+

from ..core.config import get_database_url, Settings, get_settings

# Global variables for SQLAlchemy engine and session factory
engine = None
SessionLocal = None


def get_engine() -> Engine:
    """
    Creates or returns the SQLAlchemy engine instance with connection pooling.
    
    Returns:
        Engine: SQLAlchemy engine instance
    """
    global engine
    if engine is None:
        # Get database URL from configuration
        database_url = get_database_url()
        
        # Get database pool settings from configuration
        settings = get_settings()
        pool_size = settings.DATABASE_POOL_SIZE
        max_overflow = settings.DATABASE_MAX_OVERFLOW
        
        # Create engine with connection pooling settings
        engine = create_engine(
            database_url,
            pool_size=pool_size,
            max_overflow=max_overflow,
            pool_pre_ping=True,  # Check connections before usage to prevent stale connections
            pool_recycle=3600    # Recycle connections after 1 hour to avoid stale connections
        )
    
    return engine


def init_db() -> None:
    """
    Initializes the database engine and session factory.
    """
    global SessionLocal, engine
    
    # Ensure engine is initialized
    engine = get_engine()
    
    # Create session factory
    SessionLocal = sessionmaker(
        autocommit=False,  # Transactions must be explicitly committed
        autoflush=False,   # Changes are not automatically flushed to the database
        bind=engine
    )


@contextmanager
def get_db() -> Generator[Session, Any, None]:
    """
    Provides a database session context manager with automatic cleanup.
    
    Yields:
        Generator[Session, Any, None]: Database session generator
    """
    if SessionLocal is None:
        init_db()
        
    db = SessionLocal()
    try:
        yield db
    except Exception:
        # Ensure transaction is rolled back on exception
        db.rollback()
        raise
    finally:
        # Always close the session to return connection to pool
        db.close()