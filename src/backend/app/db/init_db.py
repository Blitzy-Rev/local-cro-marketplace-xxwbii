"""
Database initialization module for the Molecular Data Management and CRO Integration Platform.

This module handles initial database setup including user creation and populating
reference data for experiment types. It supports the requirement for completely
local deployment without external dependencies.
"""

import logging  # standard library
from datetime import datetime  # standard library
from sqlalchemy.orm import Session  # sqlalchemy v2.0+

from .base import Base  # Import all models to ensure they're registered with SQLAlchemy metadata
from .session import get_db, init_db
from ..core.config import get_settings
from ..core.security import get_password_hash
from ..schemas.user import UserCreate
from ..models.user import User
from ..models.experiment_type import ExperimentType
from ..constants import UserRole, UserStatus

# Set up logging
logger = logging.getLogger(__name__)

def create_initial_users(db_session: Session) -> None:
    """
    Creates initial admin, pharma, and CRO users if they don't exist.
    
    Args:
        db_session: SQLAlchemy database session
    """
    settings = get_settings()
    
    # Check if admin user exists
    admin_user = db_session.query(User).filter(User.email == "admin@example.com").first()
    if not admin_user:
        logger.info("Creating initial admin user")
        # Create admin user
        admin_password_hash = get_password_hash("Admin@123456")
        admin_user = User(
            email="admin@example.com",
            password_hash=admin_password_hash,
            role=UserRole.ADMIN,
            status=UserStatus.ACTIVE,
            is_active=True,
            email_verified=True,
            created_at=datetime.utcnow(),
            updated_at=datetime.utcnow(),
            password_history=[admin_password_hash]
        )
        db_session.add(admin_user)
    
    # Check if demo pharma user exists
    pharma_user = db_session.query(User).filter(User.email == "pharma@example.com").first()
    if not pharma_user:
        logger.info("Creating demo pharma user")
        # Create pharma user
        pharma_password_hash = get_password_hash("Pharma@123456")
        pharma_user = User(
            email="pharma@example.com",
            password_hash=pharma_password_hash,
            role=UserRole.PHARMA,
            status=UserStatus.ACTIVE,
            is_active=True,
            email_verified=True,
            created_at=datetime.utcnow(),
            updated_at=datetime.utcnow(),
            password_history=[pharma_password_hash]
        )
        db_session.add(pharma_user)
    
    # Check if demo CRO user exists
    cro_user = db_session.query(User).filter(User.email == "cro@example.com").first()
    if not cro_user:
        logger.info("Creating demo CRO user")
        # Create CRO user
        cro_password_hash = get_password_hash("CRO@123456")
        cro_user = User(
            email="cro@example.com",
            password_hash=cro_password_hash,
            role=UserRole.CRO,
            status=UserStatus.ACTIVE,
            is_active=True,
            email_verified=True,
            created_at=datetime.utcnow(),
            updated_at=datetime.utcnow(),
            password_history=[cro_password_hash]
        )
        db_session.add(cro_user)
    
    # Commit changes to the database
    db_session.commit()

def create_initial_experiment_types(db_session: Session) -> None:
    """
    Creates initial experiment types if they don't exist.
    
    Args:
        db_session: SQLAlchemy database session
    """
    # Define default experiment types
    default_experiment_types = [
        {
            "name": "Binding Assay",
            "description": "Measures the binding affinity of compounds to target proteins",
            "category": "Binding"
        },
        {
            "name": "IC50 Determination",
            "description": "Measures the concentration of compound required for 50% inhibition",
            "category": "Binding"
        },
        {
            "name": "ADME Panel",
            "description": "Comprehensive assessment of absorption, distribution, metabolism, and excretion properties",
            "category": "ADME"
        },
        {
            "name": "Solubility Assessment",
            "description": "Measures compound solubility in various media",
            "category": "ADME"
        },
        {
            "name": "Metabolic Stability",
            "description": "Evaluates stability of compounds in liver microsomes or hepatocytes",
            "category": "ADME"
        },
        {
            "name": "Permeability Assessment",
            "description": "Measures compound permeability across cell membranes",
            "category": "ADME"
        },
        {
            "name": "Plasma Protein Binding",
            "description": "Measures binding of compounds to plasma proteins",
            "category": "ADME"
        },
        {
            "name": "CYP Inhibition",
            "description": "Assesses inhibition of cytochrome P450 enzymes",
            "category": "ADME"
        },
        {
            "name": "Toxicity Assessment",
            "description": "Evaluates potential toxic effects of compounds",
            "category": "Toxicity"
        },
        {
            "name": "Cardiac Safety",
            "description": "Evaluates potential cardiac toxicity via hERG channel inhibition",
            "category": "Toxicity"
        },
        {
            "name": "Genotoxicity",
            "description": "Assesses potential for genetic damage",
            "category": "Toxicity"
        },
        {
            "name": "Cytotoxicity Panel",
            "description": "Measures toxicity across multiple cell lines",
            "category": "Toxicity"
        },
        {
            "name": "Formulation Development",
            "description": "Development of suitable formulations for compounds",
            "category": "Formulation"
        },
        {
            "name": "Stability Testing",
            "description": "Evaluates compound stability under various conditions",
            "category": "Formulation"
        },
        {
            "name": "Pharmacokinetics",
            "description": "Evaluates compound absorption, distribution, and elimination in animal models",
            "category": "In Vivo"
        },
        {
            "name": "Efficacy Model",
            "description": "Evaluates compound efficacy in disease models",
            "category": "In Vivo"
        },
    ]
    
    # Add experiment types if they don't exist
    for exp_type in default_experiment_types:
        # Check if experiment type exists
        existing_type = db_session.query(ExperimentType).filter(
            ExperimentType.name == exp_type["name"]
        ).first()
        
        if not existing_type:
            logger.info(f"Creating experiment type: {exp_type['name']}")
            new_type = ExperimentType(
                name=exp_type["name"],
                description=exp_type["description"],
                category=exp_type["category"]
            )
            db_session.add(new_type)
    
    # Commit changes to the database
    db_session.commit()

def init_db_data() -> None:
    """
    Initializes the database with required seed data.
    
    This function is called during application startup to ensure the database
    contains all required initial data.
    """
    try:
        # Initialize database engine and session factory
        init_db()
        
        # Create a database session
        with get_db() as db:
            # Create initial users
            create_initial_users(db)
            
            # Create initial experiment types
            create_initial_experiment_types(db)
            
        logger.info("Database initialization completed successfully")
    except Exception as e:
        logger.error(f"Error initializing database: {str(e)}")
        raise