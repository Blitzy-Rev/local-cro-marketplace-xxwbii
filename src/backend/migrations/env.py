"""
Alembic migration environment configuration for the Molecular Data Management and CRO Integration Platform.

This script configures the Alembic migration environment by setting up the database connection,
loading all SQLAlchemy models, and providing context for running migrations. It supports both
'offline' mode (generating SQL scripts) and 'online' mode (directly applying migrations to the database).

This environment ensures that all database schema changes are properly tracked and can be applied
or rolled back as needed, supporting the migration strategy requirements of the platform.
"""

import logging
import os
import sys
from pathlib import Path

# Add the parent directory to Python path to allow imports from application
BASE_PATH = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(BASE_PATH))

from alembic import context
from sqlalchemy import engine_from_config, pool, create_engine

# Import app-specific modules
from app.core.config import get_database_url
from app.db.base_class import Base
# Import all models to ensure they're included in migrations
import app.db.base

# Configure Alembic logging
logger = logging.getLogger('alembic.env')

# Get Alembic configuration from alembic.ini
config = context.config

# Set metadata for Alembic to use for migrations
# This is the SQLAlchemy metadata object that contains all model definitions
target_metadata = Base.metadata

def run_migrations_offline():
    """
    Run migrations in 'offline' mode.
    
    This mode generates SQL scripts without connecting to the database, which can
    be useful for:
    - Reviewing changes before applying them
    - Generating scripts for DBA to run manually
    - Generating scripts for environments without direct database access
    
    Returns:
        None
    """
    # Get database URL from app config or from Alembic config
    url = get_database_url() or config.get_main_option("sqlalchemy.url")
    
    if not url:
        logger.error("Database URL not configured. Cannot run migrations.")
        raise ValueError("Database URL is required for migrations")
    
    logger.info("Running offline migrations")
    
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        compare_type=True,
        compare_server_default=True,
        include_schemas=True,
    )
    
    with context.begin_transaction():
        context.run_migrations()
    
    logger.info("Offline migrations completed successfully")

def run_migrations_online():
    """
    Run migrations in 'online' mode.
    
    This mode directly applies migrations to the connected database, which is
    the typical approach for development and automated deployment scenarios.
    
    Returns:
        None
    """
    # Get database URL from app config or from Alembic config
    url = get_database_url() or config.get_main_option("sqlalchemy.url")
    
    if not url:
        logger.error("Database URL not configured. Cannot run migrations.")
        raise ValueError("Database URL is required for migrations")
    
    # Set sqlalchemy.url in the Alembic config
    config.set_main_option("sqlalchemy.url", url)
    
    logger.info("Running online migrations")
    
    # Create engine with connection pooling settings
    # Using NullPool to ensure connections are closed after migration
    connectable = engine_from_config(
        config.get_section(config.config_ini_section),
        prefix="sqlalchemy.",
        poolclass=pool.NullPool,
    )
    
    try:
        with connectable.connect() as connection:
            context.configure(
                connection=connection, 
                target_metadata=target_metadata,
                compare_type=True,
                compare_server_default=True,
                include_schemas=True,
            )
            
            with context.begin_transaction():
                context.run_migrations()
        
        logger.info("Online migrations completed successfully")
    except Exception as e:
        logger.error(f"Error during migration: {str(e)}")
        raise

# Determine which function to run based on Alembic's directives
if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()