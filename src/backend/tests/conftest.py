"""
Defines pytest fixtures and configuration for the Molecular Data Management and CRO Integration Platform backend test suite.

This file provides reusable test components such as database sessions, test users with different roles,
authentication tokens, and mock services to facilitate comprehensive testing of the application.
"""

import pytest
import os
from typing import Dict, Generator, Any
from datetime import datetime

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from fastapi import FastAPI
from fastapi.testclient import TestClient

from app.db.session import get_db
from app.db.base import Base
from app.models.user import User
from app.core.security import get_password_hash
from app.core.jwt import create_access_token, create_refresh_token
from app.constants import UserRole, UserStatus

# Test database URL (use SQLite in-memory for tests by default)
TEST_DATABASE_URL = os.environ.get('TEST_DATABASE_URL', 'sqlite:///./test.db')

# Create test database engine
engine = create_engine(
    TEST_DATABASE_URL,
    # Add connect_args for SQLite only
    connect_args={"check_same_thread": False} if TEST_DATABASE_URL.startswith("sqlite") else {}
)

# Create test session factory
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def pytest_configure(config):
    """
    Pytest hook to configure the test environment.
    
    Args:
        config: Pytest config object
    """
    # Register custom markers
    config.addinivalue_line("markers", "unit: mark a test as a unit test")
    config.addinivalue_line("markers", "integration: mark a test as an integration test")
    config.addinivalue_line("markers", "database: mark a test as requiring database access")
    config.addinivalue_line("markers", "api: mark a test as an API test")
    config.addinivalue_line("markers", "auth: mark a test as requiring authentication")


@pytest.fixture(scope="session")
def test_db():
    """
    Creates test database tables for the test session and drops them after tests.
    
    Returns:
        None: This fixture doesn't return a value, it just sets up and tears down the test database.
    """
    # Create all tables
    Base.metadata.create_all(bind=engine)
    
    yield
    
    # Drop all tables after tests complete
    Base.metadata.drop_all(bind=engine)


@pytest.fixture(scope="function")
def db_session(test_db):
    """
    Provides a database session for each test with automatic transaction rollback.
    
    Args:
        test_db: Session-scoped fixture that sets up the test database
        
    Returns:
        Generator[Session, Any, None]: Database session that will be rolled back after the test
    """
    # Create a connection and begin a transaction
    connection = engine.connect()
    transaction = connection.begin()
    
    # Create a session bound to the connection
    session = TestingSessionLocal(bind=connection)
    
    yield session
    
    # Roll back the transaction and close the session
    session.close()
    transaction.rollback()
    connection.close()


@pytest.fixture(scope="function")
def client(db_session):
    """
    Provides a FastAPI test client with an overridden database session.
    
    Args:
        db_session: Database session fixture
        
    Returns:
        TestClient: FastAPI test client
    """
    from app.main import app

    # Override the get_db dependency to use the test database session
    def override_get_db():
        try:
            yield db_session
        finally:
            pass

    app.dependency_overrides[get_db] = override_get_db
    
    with TestClient(app) as test_client:
        yield test_client
        
    # Clear dependency overrides after test
    app.dependency_overrides.clear()


@pytest.fixture(scope="function")
def test_password() -> str:
    """
    Provides a test password for user creation.
    
    Returns:
        str: A password string that meets complexity requirements
    """
    return "Password123!"


@pytest.fixture(scope="function")
def test_pharma_user(db_session, test_password):
    """
    Creates a test pharma user.
    
    Args:
        db_session: Database session fixture
        test_password: Test password fixture
        
    Returns:
        User: A User instance with PHARMA role
    """
    user = User(
        email="pharma@example.com",
        password_hash=get_password_hash(test_password),
        role=UserRole.PHARMA,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True,
        password_history=[]  # Empty history for new user
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    
    return user


@pytest.fixture(scope="function")
def test_cro_user(db_session, test_password):
    """
    Creates a test CRO user.
    
    Args:
        db_session: Database session fixture
        test_password: Test password fixture
        
    Returns:
        User: A User instance with CRO role
    """
    user = User(
        email="cro@example.com",
        password_hash=get_password_hash(test_password),
        role=UserRole.CRO,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True,
        password_history=[]
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    
    return user


@pytest.fixture(scope="function")
def test_admin_user(db_session, test_password):
    """
    Creates a test admin user.
    
    Args:
        db_session: Database session fixture
        test_password: Test password fixture
        
    Returns:
        User: A User instance with ADMIN role
    """
    user = User(
        email="admin@example.com",
        password_hash=get_password_hash(test_password),
        role=UserRole.ADMIN,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True,
        password_history=[]
    )
    db_session.add(user)
    db_session.commit()
    db_session.refresh(user)
    
    return user


@pytest.fixture(scope="function")
def pharma_token_headers(test_pharma_user):
    """
    Provides authentication headers for pharma user.
    
    Args:
        test_pharma_user: Test pharma user fixture
        
    Returns:
        Dict: Headers containing Bearer token for authentication
    """
    access_token = create_access_token(
        data={"user_id": test_pharma_user.id, "email": test_pharma_user.email, "role": test_pharma_user.role.name}
    )
    return {"Authorization": f"Bearer {access_token}"}


@pytest.fixture(scope="function")
def cro_token_headers(test_cro_user):
    """
    Provides authentication headers for CRO user.
    
    Args:
        test_cro_user: Test CRO user fixture
        
    Returns:
        Dict: Headers containing Bearer token for authentication
    """
    access_token = create_access_token(
        data={"user_id": test_cro_user.id, "email": test_cro_user.email, "role": test_cro_user.role.name}
    )
    return {"Authorization": f"Bearer {access_token}"}


@pytest.fixture(scope="function")
def admin_token_headers(test_admin_user):
    """
    Provides authentication headers for admin user.
    
    Args:
        test_admin_user: Test admin user fixture
        
    Returns:
        Dict: Headers containing Bearer token for authentication
    """
    access_token = create_access_token(
        data={"user_id": test_admin_user.id, "email": test_admin_user.email, "role": test_admin_user.role.name}
    )
    return {"Authorization": f"Bearer {access_token}"}


@pytest.fixture(scope="function")
def create_test_molecule(db_session):
    """
    Factory fixture for creating test molecules.
    
    Args:
        db_session: Database session fixture
        
    Returns:
        Callable: Function to create test molecules
    """
    def _create_test_molecule(smiles, created_by, properties=None, flag_status=None):
        from app.models.molecule import Molecule
        from app.models.molecule_property import MoleculeProperty
        
        # Create the molecule
        molecule = Molecule(
            smiles=smiles,
            created_by=created_by,
            flag_status=flag_status,
            created_at=datetime.utcnow()
        )
        db_session.add(molecule)
        db_session.flush()
        
        # Add properties if provided
        if properties:
            for name, value in properties.items():
                prop = MoleculeProperty(
                    molecule_id=molecule.id,
                    property_name=name,
                    property_value=value,
                    property_unit=None,
                    is_calculated=False
                )
                db_session.add(prop)
        
        db_session.commit()
        db_session.refresh(molecule)
        return molecule
    
    return _create_test_molecule


@pytest.fixture(scope="function")
def create_test_library(db_session):
    """
    Factory fixture for creating test libraries.
    
    Args:
        db_session: Database session fixture
        
    Returns:
        Callable: Function to create test libraries
    """
    def _create_test_library(name, created_by, description=None, molecules=None):
        from app.models.library import Library
        
        # Create the library
        library = Library(
            name=name,
            description=description,
            created_by=created_by,
            created_at=datetime.utcnow(),
            updated_at=datetime.utcnow()
        )
        db_session.add(library)
        db_session.flush()
        
        # Add molecules if provided
        if molecules:
            library.molecules = molecules
        
        db_session.commit()
        db_session.refresh(library)
        return library
    
    return _create_test_library


@pytest.fixture(scope="function")
def create_test_experiment(db_session):
    """
    Factory fixture for creating test experiments.
    
    Args:
        db_session: Database session fixture
        
    Returns:
        Callable: Function to create test experiments
    """
    def _create_test_experiment(name, created_by, experiment_type_id, status=None, description=None, molecules=None, parameters=None):
        from app.models.experiment import Experiment
        from app.models.experiment_parameter import ExperimentParameter
        from app.constants import ExperimentStatus
        
        # Create the experiment
        experiment = Experiment(
            name=name,
            type_id=experiment_type_id,
            status=status or ExperimentStatus.DRAFT,
            created_by=created_by,
            description=description,
            created_at=datetime.utcnow(),
            updated_at=datetime.utcnow()
        )
        db_session.add(experiment)
        db_session.flush()
        
        # Add molecules if provided
        if molecules:
            experiment.molecules = molecules
        
        # Add parameters if provided
        if parameters:
            for name, value in parameters.items():
                param = ExperimentParameter(
                    experiment_id=experiment.id,
                    parameter_name=name,
                    parameter_value=value
                )
                db_session.add(param)
        
        db_session.commit()
        db_session.refresh(experiment)
        return experiment
    
    return _create_test_experiment


@pytest.fixture(scope="function")
def mock_file_storage():
    """
    Provides a mock file storage service for testing.
    
    Returns:
        MockFileStorage: Mock implementation of file storage service
    """
    class MockFileStorage:
        def __init__(self):
            self.files = {}
            self.buckets = {}
            
        async def create_bucket(self, bucket_name):
            self.buckets[bucket_name] = True
            return True
            
        async def bucket_exists(self, bucket_name):
            return bucket_name in self.buckets
            
        async def upload_file(self, bucket_name, file_path, file_obj, content_type=None):
            # Create bucket if it doesn't exist
            if bucket_name not in self.buckets:
                self.buckets[bucket_name] = True
                
            # Store file content
            self.files[(bucket_name, file_path)] = {
                'content': file_obj,
                'content_type': content_type,
                'size': len(file_obj) if hasattr(file_obj, '__len__') else 0,
                'created_at': datetime.utcnow()
            }
            return file_path
            
        async def download_file(self, bucket_name, file_path):
            file_data = self.files.get((bucket_name, file_path))
            return file_data['content'] if file_data else None
            
        async def get_file_url(self, bucket_name, file_path, expires=3600):
            return f"http://mock-storage/{bucket_name}/{file_path}"
            
        async def delete_file(self, bucket_name, file_path):
            if (bucket_name, file_path) in self.files:
                del self.files[(bucket_name, file_path)]
                return True
            return False
            
        async def list_files(self, bucket_name, prefix=None):
            files = []
            for (b, p), data in self.files.items():
                if b == bucket_name and (prefix is None or p.startswith(prefix)):
                    files.append({
                        'path': p,
                        'size': data['size'],
                        'last_modified': data['created_at']
                    })
            return files
    
    return MockFileStorage()


@pytest.fixture(scope="function")
def mock_molecular_processor():
    """
    Provides a mock molecular processor service for testing.
    
    Returns:
        MockMolecularProcessor: Mock implementation of molecular processor
    """
    class MockMolecularProcessor:
        def validate_smiles(self, smiles):
            """Validate SMILES string"""
            if not smiles or not isinstance(smiles, str):
                return False
                
            # Check for basic SMILES format indicators
            has_atoms = any(atom in smiles for atom in ['C', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I'])
            has_valid_chars = all(c in 'CONSPFIBrClc[]()=#@+-1234567890' for c in smiles)
            
            return has_atoms and has_valid_chars
            
        def normalize_smiles(self, smiles):
            """Normalize SMILES string"""
            if not self.validate_smiles(smiles):
                return smiles
                
            # Simple "normalization" for testing
            if 'CC' in smiles:
                return smiles.replace('CC', 'C-C')
            return smiles
            
        def calculate_properties(self, smiles):
            """Calculate molecular properties"""
            length = len(smiles)
            return {
                "MW": 100.0 + length * 10,
                "LogP": 1.5 if 'O' in smiles else 2.5,
                "TPSA": 20.0 + smiles.count('O') * 10 + smiles.count('N') * 8,
                "HBA": smiles.count('O') + smiles.count('N'),
                "HBD": smiles.count('OH') + smiles.count('NH')
            }
            
        def generate_image(self, smiles, width=300, height=200, format="svg"):
            """Generate molecular structure image"""
            if format.lower() == "svg":
                return f'<svg width="{width}" height="{height}"><text x="10" y="20">{smiles}</text></svg>'.encode()
            else:
                return b"Mock molecule image data " + smiles.encode()
            
        def search_substructure(self, smiles, pattern):
            """Search for substructure pattern"""
            return pattern in smiles
            
        def calculate_similarity(self, smiles1, smiles2):
            """Calculate similarity between molecules"""
            len_shorter = min(len(smiles1), len(smiles2))
            len_longer = max(len(smiles1), len(smiles2))
            
            matches = sum(s1 == s2 for s1, s2 in zip(smiles1, smiles2))
            
            return matches / len_longer if len_longer > 0 else 0.0
    
    return MockMolecularProcessor()