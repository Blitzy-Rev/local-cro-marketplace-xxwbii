# Introduction

This document provides comprehensive guidelines for developers working on the Molecular Data Management and CRO Integration Platform. It covers development environment setup, codebase structure, coding standards, testing procedures, and contribution guidelines.

The platform is designed for local deployment with no external dependencies, enabling small to mid-cap pharmaceutical companies to streamline their molecular data management and CRO interactions. The development process follows this principle, with all components designed to run locally during development.

# Development Environment Setup

## Prerequisites

Before setting up the development environment, ensure you have the following prerequisites installed:

- **Git**: Version control system
- **Docker**: Version 23.0 or higher
- **Docker Compose**: Version 2.17 or higher
- **Python**: Version 3.10 or higher (for backend development without Docker)
- **Poetry**: Version 1.4 or higher (for Python dependency management)
- **Node.js**: Version 16.0 or higher (for frontend development without Docker)
- **npm**: Version 8.0 or higher (for JavaScript package management)
- **Visual Studio Code** (recommended) or another IDE with Python and TypeScript support

## Repository Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/your-organization/molecular-platform.git
   cd molecular-platform
   ```

2. Create development branches from the main branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```

3. Set up Git hooks for code quality checks (optional):
   ```bash
   # Backend
   cd src/backend
   poetry install
   poetry run pre-commit install
   
   # Frontend
   cd src/web
   npm install
   npx husky install
   ```

## Docker-based Development Environment

The recommended approach for development is to use Docker Compose, which sets up all required services in isolated containers:

1. Navigate to the infrastructure directory:
   ```bash
   cd infrastructure
   ```

2. Create a `.env` file based on the example:
   ```bash
   cp .env.example .env
   # Edit .env with your preferred configuration
   ```

3. Start the development environment:
   ```bash
   docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d
   ```

This will start all services in development mode with the following features:
- Hot reloading for both frontend and backend code
- Development-specific configurations
- Exposed debugging ports
- Mounted source code volumes

The development environment will be accessible at:
- Frontend: http://localhost:3000
- Backend API: http://localhost:8000
- API Documentation: http://localhost:8000/docs
- MinIO Console: http://localhost:9001
- Grafana: http://localhost:3001
- Kibana: http://localhost:5601

## Local Development Setup (Without Docker)

For developers who prefer to run services directly on their local machine, follow these steps:

**Backend Setup:**

1. Navigate to the backend directory:
   ```bash
   cd src/backend
   ```

2. Install dependencies using Poetry:
   ```bash
   poetry install
   ```

3. Create a `.env` file based on the example:
   ```bash
   cp .env.example .env
   # Edit .env with your local configuration
   ```

4. Start the backend service:
   ```bash
   poetry run uvicorn main:app --reload
   ```

The backend API will be accessible at http://localhost:8000

**Frontend Setup:**

1. Navigate to the frontend directory:
   ```bash
   cd src/web
   ```

2. Install dependencies using npm:
   ```bash
   npm install
   ```

3. Create a `.env` file based on the example:
   ```bash
   cp .env.development .env
   # Edit .env with your local configuration
   ```

4. Start the frontend development server:
   ```bash
   npm run dev
   ```

The frontend will be accessible at http://localhost:5173

**Note:** When running services locally, you'll need to set up and configure PostgreSQL, Redis, and MinIO separately, or connect to the containerized versions of these services.

## IDE Configuration

**Visual Studio Code (Recommended):**

Install the following extensions for an optimal development experience:

- Python extension (Microsoft)
- Pylance for Python language server
- ESLint for JavaScript/TypeScript linting
- Prettier for code formatting
- Docker extension for container management
- GitLens for enhanced Git integration

Workspace settings (`.vscode/settings.json`):
```json
{
  "python.linting.enabled": true,
  "python.linting.pylintEnabled": false,
  "python.linting.flake8Enabled": true,
  "python.formatting.provider": "black",
  "editor.formatOnSave": true,
  "editor.codeActionsOnSave": {
    "source.organizeImports": true,
    "source.fixAll.eslint": true
  },
  "[python]": {
    "editor.formatOnSave": true,
    "editor.defaultFormatter": "ms-python.python"
  },
  "[typescript]": {
    "editor.formatOnSave": true,
    "editor.defaultFormatter": "esbenp.prettier-vscode"
  },
  "[typescriptreact]": {
    "editor.formatOnSave": true,
    "editor.defaultFormatter": "esbenp.prettier-vscode"
  }
}
```

Launch configurations (`.vscode/launch.json`) for debugging:
```json
{
  "version": "0.2.0",
  "configurations": [
    {
      "name": "Python: FastAPI",
      "type": "python",
      "request": "launch",
      "module": "uvicorn",
      "args": ["main:app", "--reload"],
      "cwd": "${workspaceFolder}/src/backend",
      "env": {
        "PYTHONPATH": "${workspaceFolder}/src/backend"
      }
    },
    {
      "name": "Chrome: React",
      "type": "chrome",
      "request": "launch",
      "url": "http://localhost:5173",
      "webRoot": "${workspaceFolder}/src/web"
    }
  ]
}
```

# Project Structure

## Repository Organization

The repository is organized into the following main directories:

```
/
├── docs/                  # Documentation
├── infrastructure/        # Docker and deployment configuration
├── scripts/              # Utility scripts for deployment and maintenance
├── src/                  # Source code
│   ├── backend/          # Backend API and services
│   └── web/              # Frontend application
└── .github/              # GitHub workflows and templates
```

This organization separates the application code from infrastructure and deployment concerns, making it easier to understand and maintain the codebase.

## Backend Structure

The backend code is organized following a modular structure with clear separation of concerns:

```
src/backend/
├── app/                    # Application code
│   ├── api/                # API endpoints and schemas
│   │   └── api_v1/         # API version 1
│   │       ├── endpoints/  # API endpoint implementations
│   │       └── schemas/    # Request/response schemas
│   ├── core/               # Core functionality
│   │   ├── config.py       # Application configuration
│   │   ├── security.py     # Security utilities
│   │   └── jwt.py          # JWT handling
│   ├── crud/               # Database CRUD operations
│   ├── db/                 # Database setup and session
│   ├── models/             # SQLAlchemy models
│   ├── schemas/            # Pydantic schemas
│   ├── services/           # Business logic services
│   ├── utils/              # Utility functions
│   ├── worker/             # Celery worker tasks
│   └── molecular/          # Molecular processing
├── migrations/             # Alembic database migrations
├── tests/                  # Test suite
├── docker/                 # Docker configuration
├── scripts/                # Utility scripts
├── .env.example            # Example environment variables
├── alembic.ini             # Alembic configuration
├── Dockerfile              # Docker build configuration
├── entrypoint.sh           # Container entrypoint script
├── gunicorn_conf.py        # Gunicorn configuration
├── main.py                 # Application entry point
├── pyproject.toml          # Poetry configuration
└── requirements.txt        # Python dependencies
```

This structure follows a layered architecture pattern:

1. **API Layer**: Handles HTTP requests and responses (app/api)
2. **Service Layer**: Implements business logic (app/services)
3. **Data Access Layer**: Manages database operations (app/crud, app/models)
4. **Domain Layer**: Defines core domain models and logic (app/schemas, app/molecular)

The separation of concerns makes the codebase more maintainable and testable.

## Frontend Structure

The frontend code is organized following a feature-based structure with shared components:

```
src/web/
├── src/                  # Source code
│   ├── api/              # API client functions
│   ├── components/       # Shared UI components
│   │   ├── common/       # Generic UI components
│   │   ├── molecular/    # Molecule-specific components
│   │   └── data-visualization/ # Charts and visualizations
│   ├── features/         # Feature modules
│   │   ├── auth/         # Authentication
│   │   ├── dashboard/    # Dashboard
│   │   ├── molecules/    # Molecule management
│   │   ├── libraries/    # Library management
│   │   ├── experiments/  # Experiment management
│   │   ├── submissions/  # CRO submissions
│   │   ├── results/      # Experimental results
│   │   ├── cro-interface/# CRO user interface
│   │   ├── admin/        # Admin interface
│   │   └── communications/ # User communications
│   ├── hooks/            # Custom React hooks
│   ├── layouts/          # Page layouts
│   ├── store/            # Redux store configuration
│   ├── theme/            # UI theme configuration
│   ├── types/            # TypeScript type definitions
│   ├── utils/            # Utility functions
│   ├── __tests__/        # Test setup and mocks
│   ├── App.tsx           # Main application component
│   ├── index.tsx         # Application entry point
│   └── routes.tsx        # Routing configuration
├── public/               # Static assets
├── docker/               # Docker configuration
├── .env.development      # Development environment variables
├── .env.production       # Production environment variables
├── package.json          # Dependencies and scripts
├── tsconfig.json         # TypeScript configuration
├── vite.config.ts        # Vite configuration
└── README.md             # Documentation
```

This structure follows a feature-based organization pattern:

1. **Features**: Self-contained modules for specific functionality
2. **Shared Components**: Reusable UI components
3. **Infrastructure**: Core application setup and configuration

Each feature contains its own components, hooks, and pages, promoting separation of concerns and code organization.

## Infrastructure Structure

The infrastructure code is organized to support both development and production deployments:

```
infrastructure/
├── docker-compose.yml           # Production container configuration
├── docker-compose.dev.yml       # Development overrides
├── .env.example                 # Example environment variables
├── postgres/                    # PostgreSQL configuration
│   ├── init.sql                 # Database initialization script
│   ├── backup.sh                # Backup script
│   └── restore.sh               # Restore script
├── redis/                       # Redis configuration
│   └── redis.conf               # Redis configuration file
├── minio/                       # MinIO configuration
│   └── init.sh                  # Bucket initialization script
├── nginx/                       # Nginx configuration
│   ├── nginx.conf               # Main configuration
│   └── default.conf             # Default site configuration
├── prometheus/                  # Prometheus configuration
│   └── prometheus.yml           # Prometheus configuration file
├── grafana/                     # Grafana configuration
│   ├── datasources/             # Data source configurations
│   └── dashboards/              # Dashboard definitions
├── fluentd/                     # Fluentd configuration
│   └── fluent.conf              # Fluentd configuration file
├── elasticsearch/               # Elasticsearch configuration
│   └── elasticsearch.yml        # Elasticsearch configuration file
└── kibana/                      # Kibana configuration
    └── kibana.yml               # Kibana configuration file
```

This structure separates configuration for each service, making it easier to manage and customize the deployment environment.

# Backend Development

## Technology Stack

The backend uses the following key technologies:

- **Python 3.10+**: Core programming language
- **FastAPI 0.95+**: Web framework for building APIs
- **SQLAlchemy 2.0+**: ORM for database interactions
- **Alembic**: Database migration tool
- **Pydantic 2.0+**: Data validation and settings management
- **RDKit 2023.03+**: Cheminformatics toolkit for molecular processing
- **Pandas 2.0+**: Data analysis and manipulation library
- **Celery 5.2+**: Distributed task queue for background processing
- **Redis 4.5+**: Cache and message broker
- **MinIO 7.1+**: S3-compatible object storage
- **PyJWT 2.7+**: JWT token handling
- **Passlib 1.7.4+**: Password hashing

These technologies were selected to provide a robust, high-performance backend that can handle the specific requirements of molecular data management while maintaining the ability to deploy locally without external dependencies.

## Development Workflow

The recommended workflow for backend development is:

1. **Create a Feature Branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Implement Changes**:
   - Add/modify API endpoints in `app/api/api_v1/endpoints/`
   - Implement business logic in `app/services/`
   - Define data models in `app/models/` and schemas in `app/schemas/`
   - Add database operations in `app/crud/`

3. **Run Tests**:
   ```bash
   # Using Poetry
   poetry run pytest
   
   # Using Docker
   docker-compose exec backend pytest
   ```

4. **Format and Lint Code**:
   ```bash
   # Using Poetry
   poetry run black app tests
   poetry run isort app tests
   poetry run flake8 app tests
   
   # Using Docker
   docker-compose exec backend black app tests
   docker-compose exec backend isort app tests
   docker-compose exec backend flake8 app tests
   ```

5. **Create Database Migrations** (if needed):
   ```bash
   # Using Poetry
   poetry run alembic revision --autogenerate -m "description of changes"
   
   # Using Docker
   docker-compose exec backend alembic revision --autogenerate -m "description of changes"
   ```

6. **Commit Changes**:
   ```bash
   git add .
   git commit -m "feat: add your feature description"
   ```

7. **Push Changes and Create Pull Request**:
   ```bash
   git push origin feature/your-feature-name
   # Create PR through GitHub interface
   ```

## API Development

When developing new API endpoints, follow these guidelines:

1. **Create a new endpoint file** in `app/api/api_v1/endpoints/` if needed

2. **Define request and response schemas** in `app/api/api_v1/schemas/`

3. **Implement the endpoint** using FastAPI's dependency injection pattern:

```python
from fastapi import APIRouter, Depends, HTTPException, status
from sqlalchemy.orm import Session

from app.api.api_v1.schemas import ItemCreate, ItemResponse
from app.core.auth import get_current_user
from app.db.session import get_db
from app.models.user import User
from app.services.item_service import ItemService

router = APIRouter()

@router.post("/", response_model=ItemResponse, status_code=status.HTTP_201_CREATED)
def create_item(
    item_in: ItemCreate,
    db: Session = Depends(get_db),
    current_user: User = Depends(get_current_user),
):
    """Create a new item."""
    item_service = ItemService(db)
    item = item_service.create_item(item_in=item_in, user_id=current_user.id)
    return item
```

4. **Register the router** in `app/api/api_v1/api.py`:

```python
from fastapi import APIRouter

from app.api.api_v1.endpoints import items, users, auth

api_router = APIRouter()
api_router.include_router(auth.router, prefix="/auth", tags=["auth"])
api_router.include_router(users.router, prefix="/users", tags=["users"])
api_router.include_router(items.router, prefix="/items", tags=["items"])
```

5. **Document the endpoint** with comprehensive docstrings and parameter descriptions for automatic OpenAPI documentation

6. **Implement proper error handling** using appropriate HTTP status codes and error messages

7. **Add tests** for the new endpoint in `tests/api/api_v1/`

## Database Operations

When working with the database, follow these guidelines:

1. **Define models** in `app/models/` using SQLAlchemy's declarative syntax:

```python
from sqlalchemy import Column, ForeignKey, Integer, String, DateTime, func
from sqlalchemy.orm import relationship

from app.db.base_class import Base

class Item(Base):
    __tablename__ = "items"

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    description = Column(String)
    owner_id = Column(Integer, ForeignKey("users.id"))
    created_at = Column(DateTime, server_default=func.now())
    updated_at = Column(DateTime, onupdate=func.now())

    owner = relationship("User", back_populates="items")
```

2. **Create CRUD operations** in `app/crud/` using the base CRUD class:

```python
from typing import List, Optional

from sqlalchemy.orm import Session

from app.crud.base import CRUDBase
from app.models.item import Item
from app.schemas.item import ItemCreate, ItemUpdate

class CRUDItem(CRUDBase[Item, ItemCreate, ItemUpdate]):
    def get_by_name(self, db: Session, *, name: str) -> Optional[Item]:
        return db.query(self.model).filter(self.model.name == name).first()

    def get_multi_by_owner(
        self, db: Session, *, owner_id: int, skip: int = 0, limit: int = 100
    ) -> List[Item]:
        return (
            db.query(self.model)
            .filter(self.model.owner_id == owner_id)
            .offset(skip)
            .limit(limit)
            .all()
        )

item = CRUDItem(Item)
```

3. **Create database migrations** when changing models:

```bash
# Generate migration
alembic revision --autogenerate -m "add item table"

# Apply migration
alembic upgrade head
```

4. **Use transactions** for operations that require atomicity:

```python
from sqlalchemy.orm import Session

def transfer_item(db: Session, item_id: int, from_user_id: int, to_user_id: int) -> None:
    try:
        # Start transaction
        item = db.query(Item).filter(Item.id == item_id, Item.owner_id == from_user_id).first()
        if not item:
            raise ValueError("Item not found or not owned by from_user")
        
        item.owner_id = to_user_id
        db.add(item)
        
        # Create transfer record
        transfer = ItemTransfer(item_id=item_id, from_user_id=from_user_id, to_user_id=to_user_id)
        db.add(transfer)
        
        # Commit transaction
        db.commit()
    except Exception as e:
        # Rollback on error
        db.rollback()
        raise e
```

## Molecular Processing

When working with molecular data using RDKit, follow these guidelines:

1. **Use the MoleculeProcessor class** for high-level operations:

```python
from app.molecular.processor import MoleculeProcessor

processor = MoleculeProcessor()

# Validate SMILES
valid, error = processor.validate_smiles("CCO")

# Calculate properties
properties = processor.calculate_properties("CCO")

# Process multiple molecules in batch
results = processor.process_molecules(["CCO", "CCCCO", "c1ccccc1"])
```

2. **Implement efficient batch processing** for large datasets:

```python
from app.molecular.processor import MoleculeProcessor
import pandas as pd

def process_csv(csv_path: str, smiles_column: str):
    # Read CSV in chunks
    chunk_size = 1000
    results = []
    
    for chunk in pd.read_csv(csv_path, chunksize=chunk_size):
        # Extract SMILES
        smiles_list = chunk[smiles_column].tolist()
        
        # Process in batch
        processor = MoleculeProcessor()
        batch_results = processor.process_molecules(smiles_list)
        
        # Combine with original data
        chunk_with_results = pd.concat([chunk, pd.DataFrame(batch_results)], axis=1)
        results.append(chunk_with_results)
    
    # Combine all results
    return pd.concat(results)
```

3. **Use caching for repeated operations** on the same molecules:

```python
from app.molecular.processor import MoleculeProcessor
from functools import lru_cache

class CachedMoleculeProcessor(MoleculeProcessor):
    @lru_cache(maxsize=1000)
    def get_molecule(self, smiles: str):
        return self._create_molecule(smiles)
    
    def calculate_properties(self, smiles: str):
        mol = self.get_molecule(smiles)
        if mol is None:
            return None
        return self._calculate_properties(mol)
```

4. **Implement error handling** for invalid molecular structures:

```python
from app.molecular.processor import MoleculeProcessor
from app.molecular.exceptions import InvalidSmilesError

processor = MoleculeProcessor()

try:
    properties = processor.calculate_properties("invalid_smiles")
except InvalidSmilesError as e:
    print(f"Invalid SMILES: {e}")
```

## Background Tasks

For long-running operations, use Celery tasks to process them in the background:

1. **Define a task** in `app/worker/tasks/`:

```python
from app.worker.celery_app import celery_app
from app.molecular.processor import MoleculeProcessor

@celery_app.task(bind=True)
def process_molecules_task(self, smiles_list: list[str]):
    """Process a list of molecules in the background."""
    processor = MoleculeProcessor()
    
    results = []
    total = len(smiles_list)
    
    for i, smiles in enumerate(smiles_list):
        try:
            # Process molecule
            result = processor.process_molecule(smiles)
            results.append(result)
            
            # Update progress
            self.update_state(
                state="PROGRESS",
                meta={"current": i + 1, "total": total, "status": "Processing"}
            )
        except Exception as e:
            results.append({"smiles": smiles, "error": str(e)})
    
    return {"status": "Complete", "results": results}
```

2. **Submit the task** from a service or API endpoint:

```python
from app.worker.tasks.molecule_tasks import process_molecules_task

def process_molecules_async(smiles_list: list[str]):
    # Submit task to Celery
    task = process_molecules_task.delay(smiles_list)
    return {"task_id": task.id}
```

3. **Check task status** and retrieve results:

```python
from celery.result import AsyncResult

def get_task_status(task_id: str):
    task_result = AsyncResult(task_id)
    
    if task_result.state == "PENDING":
        response = {
            "state": task_result.state,
            "status": "Pending"
        }
    elif task_result.state == "FAILURE":
        response = {
            "state": task_result.state,
            "status": "Error",
            "error": str(task_result.info)
        }
    elif task_result.state == "PROGRESS":
        response = {
            "state": task_result.state,
            "status": task_result.info.get("status", ""),
            "current": task_result.info.get("current", 0),
            "total": task_result.info.get("total", 1),
            "percent": int(task_result.info.get("current", 0) / task_result.info.get("total", 1) * 100)
        }
    else:
        response = {
            "state": task_result.state,
            "status": "Complete",
            "results": task_result.info.get("results", [])
        }
    
    return response
```

## Testing

The backend uses pytest for testing. Follow these guidelines for effective testing:

1. **Organize tests** to mirror the application structure:

```
tests/
├── api/              # API endpoint tests
├── crud/             # Database operation tests
├── services/         # Service layer tests
├── molecular/        # Molecular processing tests
└── utils/            # Utility function tests
```

2. **Use fixtures** for common test setup:

```python
# tests/conftest.py
import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from fastapi.testclient import TestClient

from app.db.base import Base
from app.main import app
from app.db.session import get_db

@pytest.fixture
def db_engine():
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(bind=engine)
    return engine

@pytest.fixture
def db_session(db_engine):
    Session = sessionmaker(autocommit=False, autoflush=False, bind=db_engine)
    session = Session()
    try:
        yield session
    finally:
        session.close()

@pytest.fixture
def client(db_session):
    def override_get_db():
        try:
            yield db_session
        finally:
            pass
    
    app.dependency_overrides[get_db] = override_get_db
    with TestClient(app) as test_client:
        yield test_client
    app.dependency_overrides = {}
```

3. **Write comprehensive tests** for API endpoints:

```python
# tests/api/api_v1/test_items.py
def test_create_item(client, db_session, normal_user_token_headers):
    data = {"name": "Test Item", "description": "Test Description"}
    response = client.post(
        "/api/v1/items/",
        json=data,
        headers=normal_user_token_headers,
    )
    assert response.status_code == 201
    content = response.json()
    assert content["name"] == data["name"]
    assert content["description"] == data["description"]
    assert "id" in content
```

4. **Test database operations** with in-memory SQLite:

```python
# tests/crud/test_item.py
def test_create_item(db_session):
    item_in = ItemCreate(name="Test Item", description="Test Description")
    user_id = 1
    item = crud.item.create_with_owner(db=db_session, obj_in=item_in, owner_id=user_id)
    assert item.name == item_in.name
    assert item.description == item_in.description
    assert item.owner_id == user_id
```

5. **Test molecular processing** with known molecules:

```python
# tests/molecular/test_processor.py
def test_validate_smiles():
    processor = MoleculeProcessor()
    
    # Valid SMILES
    valid, error = processor.validate_smiles("CCO")
    assert valid is True
    assert error is None
    
    # Invalid SMILES
    valid, error = processor.validate_smiles("invalid_smiles")
    assert valid is False
    assert error is not None
```

6. **Run tests** with coverage reporting:

```bash
pytest --cov=app tests/
```

# Frontend Development

## Technology Stack

The frontend uses the following key technologies:

- **React 18.2+**: UI framework for building component-based interfaces
- **TypeScript 4.9+**: Type-safe JavaScript for improved developer experience
- **Redux Toolkit 1.9+**: State management for complex application state
- **React Query 4.28+**: Data fetching with caching for API requests
- **Material-UI 5.13+**: UI component library with comprehensive components
- **React Router 6.11+**: Client-side routing
- **Formik 2.2+**: Form handling with validation
- **React DnD 16.0+**: Drag and drop functionality for molecule organization
- **D3.js 7.8+**: Advanced data visualization
- **Chart.js 4.3+**: Charting library for data visualization
- **Vite 4.3+**: Build tool for fast development and optimized production builds

These technologies were selected to provide a responsive, interactive user interface that can handle complex molecular data visualization and organization while maintaining good performance.

## Development Workflow

The recommended workflow for frontend development is:

1. **Create a Feature Branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Implement Changes**:
   - Create/modify components in the appropriate feature directory
   - Update Redux state as needed
   - Add API client functions for new endpoints

3. **Run Tests**:
   ```bash
   # Using npm
   npm test
   
   # Using Docker
   docker-compose exec frontend npm test
   ```

4. **Format and Lint Code**:
   ```bash
   # Using npm
   npm run lint
   npm run format
   
   # Using Docker
   docker-compose exec frontend npm run lint
   docker-compose exec frontend npm run format
   ```

5. **Commit Changes**:
   ```bash
   git add .
   git commit -m "feat: add your feature description"
   ```

6. **Push Changes and Create Pull Request**:
   ```bash
   git push origin feature/your-feature-name
   # Create PR through GitHub interface
   ```

## Component Development

When developing new components, follow these guidelines:

1. **Create components in the appropriate feature directory**:

```tsx
// src/features/molecules/components/MoleculeCard.tsx
import React from 'react';
import { Card, CardContent, Typography, Box } from '@mui/material';
import { MoleculeViewer } from '../../../components/molecular/MoleculeViewer';
import { Molecule } from '../../../types/molecule';

interface MoleculeCardProps {
  molecule: Molecule;
  onClick?: (molecule: Molecule) => void;
}

export const MoleculeCard: React.FC<MoleculeCardProps> = ({ molecule, onClick }) => {
  const handleClick = () => {
    if (onClick) {
      onClick(molecule);
    }
  };

  return (
    <Card onClick={handleClick} sx={{ cursor: onClick ? 'pointer' : 'default' }}>
      <CardContent>
        <Box sx={{ height: 150, display: 'flex', justifyContent: 'center' }}>
          <MoleculeViewer smiles={molecule.smiles} width={120} height={120} />
        </Box>
        <Typography variant="subtitle1" gutterBottom>
          {molecule.smiles}
        </Typography>
        <Typography variant="body2" color="text.secondary">
          MW: {molecule.properties.molecular_weight?.toFixed(2)}
        </Typography>
        <Typography variant="body2" color="text.secondary">
          LogP: {molecule.properties.logp?.toFixed(2)}
        </Typography>
      </CardContent>
    </Card>
  );
};
```

2. **Create pages that use the components**:

```tsx
// src/features/molecules/pages/MoleculesListPage.tsx
import React, { useState } from 'react';
import { Container, Grid, Typography, Box } from '@mui/material';
import { MoleculeCard } from '../components/MoleculeCard';
import { MoleculeFilter } from '../components/MoleculeFilter';
import { useMolecules } from '../hooks/useMolecules';
import { Molecule } from '../../../types/molecule';
import { useNavigate } from 'react-router-dom';

export const MoleculesListPage: React.FC = () => {
  const [filters, setFilters] = useState({});
  const { molecules, isLoading, error } = useMolecules(filters);
  const navigate = useNavigate();

  const handleMoleculeClick = (molecule: Molecule) => {
    navigate(`/molecules/${molecule.id}`);
  };

  return (
    <Container maxWidth="lg">
      <Typography variant="h4" component="h1" gutterBottom>
        Molecules
      </Typography>
      
      <MoleculeFilter onFilterChange={setFilters} />
      
      {isLoading && <Typography>Loading...</Typography>}
      {error && <Typography color="error">Error: {error.message}</Typography>}
      
      <Box mt={3}>
        <Grid container spacing={3}>
          {molecules.map((molecule) => (
            <Grid item xs={12} sm={6} md={4} lg={3} key={molecule.id}>
              <MoleculeCard molecule={molecule} onClick={handleMoleculeClick} />
            </Grid>
          ))}\n        </Grid>
      </Box>
    </Container>
  );
};
```

3. **Create custom hooks for logic reuse**: