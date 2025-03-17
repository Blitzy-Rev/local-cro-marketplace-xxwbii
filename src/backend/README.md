# Molecular Data Management and CRO Integration Platform - Backend

Backend component of the Molecular Data Management and CRO Integration Platform, providing API endpoints, molecular data processing, and integration with CRO services. This system is designed for local deployment with no external dependencies, enabling small to mid-cap pharmaceutical companies to streamline their molecular data management and CRO interactions.

## Features

- User authentication and role-based access control (Pharma, CRO, Admin)
- CSV upload and molecular data ingestion with flexible mapping
- Molecular structure validation and property calculation using RDKit
- Interactive molecule sorting, filtering, and organization
- Experiment creation and submission to CROs
- Result management and visualization
- Background task processing for computationally intensive operations
- Comprehensive logging and monitoring

## Technology Stack

- **Language:** Python 3.10+
- **Web Framework:** FastAPI 0.95+
- **ORM:** SQLAlchemy 2.0+
- **Cheminformatics:** RDKit 2023.03+
- **Data Processing:** Pandas 2.0+
- **Validation:** Pydantic 2.0+
- **Authentication:** PyJWT 2.7+, Passlib 1.7.4+
- **Database:** PostgreSQL 15+
- **Caching & Queue:** Redis 7.0+
- **Object Storage:** MinIO (S3-compatible)
- **Task Queue:** Celery
- **Containerization:** Docker, Docker Compose

## Project Structure

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
├── requirements.txt        # Python dependencies
└── README.md               # This file
```

## Setup and Installation

### Prerequisites

- Docker and Docker Compose
- Git

### Local Development Setup

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. Create a `.env` file based on `.env.example`:
   ```bash
   cp .env.example .env
   # Edit .env with your configuration
   ```

3. Install dependencies with Poetry (optional, for development without Docker):
   ```bash
   poetry install
   ```

4. Run the application with Docker Compose:
   ```bash
   cd infrastructure
   docker-compose up -d
   ```

5. Access the API documentation at http://localhost/api/docs

### Running Tests

```bash
# Run with pytest
poetry run pytest

# Run with coverage report
poetry run pytest --cov=app tests/
```

## API Documentation

The API documentation is automatically generated using FastAPI's built-in Swagger UI and ReDoc:

- Swagger UI: http://localhost/api/docs
- ReDoc: http://localhost/api/redoc

### API Endpoints

The API is organized into the following endpoint groups:

- `/api/auth`: Authentication endpoints (login, register, token refresh)
- `/api/users`: User management endpoints
- `/api/admin`: Administrative endpoints
- `/api/molecules`: Molecule management endpoints
- `/api/libraries`: Library management endpoints
- `/api/experiments`: Experiment management endpoints
- `/api/submissions`: CRO submission endpoints
- `/api/results`: Result management endpoints
- `/api/csv`: CSV upload and processing endpoints
- `/api/health`: Health check endpoints

## Molecular Processing

The system uses RDKit for molecular processing, including:

- SMILES validation and normalization
- Molecular property calculation
- Substructure and similarity searching
- Molecule visualization

The `MoleculeProcessor` class in `app/molecular/processor.py` provides a high-level interface for these operations, with optimizations for handling large datasets efficiently.

## Database Models

The system uses SQLAlchemy ORM with the following key models:

- `User`: User accounts with role-based permissions
- `Molecule`: Molecular structures with SMILES representation
- `MoleculeProperty`: Properties associated with molecules
- `Library`: User-defined collections of molecules
- `Experiment`: Experimental configurations for testing molecules
- `Submission`: CRO submissions for experimental testing
- `Result`: Experimental results from CROs

Database migrations are managed using Alembic.

## Background Tasks

Computationally intensive operations are handled asynchronously using Celery tasks:

- CSV processing and molecule import
- Batch property calculation
- Report generation
- Notification delivery

The worker container runs these tasks in the background, communicating with the main application through Redis.

## Configuration

The application is configured using environment variables, with defaults defined in `app/core/config.py`. Key configuration options include:

- Database connection parameters
- Redis connection parameters
- MinIO connection parameters
- JWT authentication settings
- CORS settings
- Logging configuration

For local development, these can be set in the `.env` file. In Docker, they are passed through the Docker Compose configuration.

## Deployment

The application is designed for local deployment using Docker Compose. The `infrastructure/docker-compose.yml` file defines all required services:

- `nginx`: Reverse proxy for routing requests
- `frontend`: React frontend application
- `backend`: FastAPI backend application
- `worker`: Celery worker for background tasks
- `postgres`: PostgreSQL database
- `redis`: Redis for caching and message queue
- `minio`: MinIO object storage
- `prometheus`: Metrics collection
- `grafana`: Metrics visualization
- `fluentd`: Log collection
- `elasticsearch`: Log storage and indexing
- `kibana`: Log visualization

To deploy the application:

```bash
cd infrastructure
docker-compose up -d
```

## Development Guidelines

### Code Style

The project follows PEP 8 style guidelines with some modifications. Code formatting is enforced using Black and isort.

```bash
# Format code
poetry run black app tests
poetry run isort app tests

# Check code style
poetry run flake8 app tests
```

### Testing

All new features should include tests. The project uses pytest for testing.

### Documentation

- All functions, classes, and modules should include docstrings
- API endpoints should include comprehensive descriptions for Swagger documentation
- Complex logic should be explained with inline comments

### Git Workflow

- Create feature branches from `develop`
- Submit pull requests for review
- Ensure all tests pass before merging
- Follow conventional commit message format

## Monitoring and Observability

The system includes comprehensive monitoring and observability features:

- Health check endpoints at `/api/health`
- Prometheus metrics for performance monitoring
- Grafana dashboards for visualization
- Structured logging with correlation IDs
- ELK stack (Elasticsearch, Logstash, Kibana) for log analysis

Access the monitoring interfaces at:

- Grafana: http://localhost:3000
- Kibana: http://localhost:5601

## License

This project is licensed under the terms specified in the LICENSE file.