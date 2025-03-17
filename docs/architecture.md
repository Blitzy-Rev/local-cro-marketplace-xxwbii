# Introduction

This document provides a comprehensive overview of the architecture for the Molecular Data Management and CRO Integration Platform. The platform is designed to streamline the small molecule drug discovery process for small to mid-cap pharmaceutical companies by providing an intuitive interface for molecular data management, organization, and experimental workflow submission to Contract Research Organizations (CROs).

The architecture follows these key principles:

- **Modularity**: Independent services with well-defined responsibilities and interfaces
- **Local-first design**: All components operate without external dependencies
- **Stateful persistence**: Reliable data storage with transaction support
- **Layered security**: Authentication and authorization at multiple levels
- **Asynchronous processing**: Background tasks for computationally intensive operations

This document is intended for developers, system administrators, and technical stakeholders who need to understand the system's architecture, component interactions, and deployment model.

# System Overview

## Business Context

The Molecular Data Management and CRO Integration Platform addresses a critical gap in the small molecule drug discovery process. Small to mid-cap pharmaceutical companies often struggle with efficient molecular data management and seamless interactions with Contract Research Organizations (CROs) for experimental testing.

The platform serves as a bridge between computational chemistry and experimental validation by providing:

1. Efficient molecular data management through CSV import and organization
2. Interactive molecular library creation and management
3. Streamlined experiment definition and submission to CROs
4. Bidirectional communication with CROs for quotes, status updates, and results
5. Comprehensive result visualization and analysis

All of this functionality is delivered in a fully local deployment model, ensuring data security and eliminating external dependencies.

## User Roles

The system supports three primary user roles:

1. **Pharma User**: Researchers from pharmaceutical companies who manage molecular data, create libraries, define experiments, and review results

2. **CRO User**: Staff from Contract Research Organizations who receive experiment submissions, provide quotes, update experiment status, and upload results

3. **Administrator**: System administrators who manage users, monitor system health, and configure system settings

Each role has a dedicated user interface with role-specific functionality and access controls.

## Key Workflows

The platform supports these key workflows:

1. **Molecular Data Management**:
   - CSV upload with flexible column mapping
   - Molecular property calculation and validation
   - Interactive filtering, sorting, and organization
   - Custom library creation and management

2. **Experiment Management**:
   - Experiment definition with configurable parameters
   - Molecule selection for experiments
   - Experiment queuing and prioritization
   - Submission to selected CROs

3. **CRO Interaction**:
   - Experiment submission review by CROs
   - Quote provision and approval workflow
   - Experiment status tracking and updates
   - Result upload and notification

4. **Result Analysis**:
   - Result visualization and comparison
   - Data export for further analysis
   - Result organization and annotation
   - Historical result tracking

# High-Level Architecture

## Architecture Overview

The Molecular Data Management and CRO Integration Platform employs a containerized microservices architecture to enable full local deployment while maintaining separation of concerns. The system is composed of the following layers:

1. **Presentation Layer**: React-based frontend application with role-specific interfaces

2. **API Layer**: FastAPI backend providing RESTful endpoints for all functionality

3. **Service Layer**: Business logic implementation with domain-specific services

4. **Data Access Layer**: Database abstraction with repositories and data models

5. **Infrastructure Layer**: Supporting services including database, caching, file storage, and message queue

All components are containerized using Docker and orchestrated with Docker Compose for simplified deployment and management.

```
+---------------------+
|    User Interface   |
| (React, TypeScript) |
+---------------------+
           |
           v
+---------------------+
|     API Gateway     |
|      (Nginx)        |
+---------------------+
           |
           v
+---------------------+
|    Backend API      |
|     (FastAPI)       |
+---------------------+
           |
      +----+----+
      |         |
      v         v
+----------+ +----------+
| Database | | Storage  |
|(PostgreSQL)| (MinIO)   |
+----------+ +----------+
```

## Component Architecture

The system consists of the following core components:

1. **Frontend Application**: React-based single-page application (SPA) providing the user interface for both pharma and CRO users. Built with TypeScript, Material-UI, and Redux for state management.

2. **Backend API**: FastAPI application providing RESTful endpoints for all system functionality. Implements business logic, data validation, and service orchestration.

3. **Authentication Service**: Handles user authentication, authorization, and session handling using JWT tokens.

4. **Molecular Processing Engine**: Specialized component for molecular structure validation, property calculation, and analysis using RDKit.

5. **Database Service**: PostgreSQL database for persistent storage of all structured data including users, molecules, libraries, experiments, and results.

6. **File Storage Service**: MinIO object storage for managing CSV files, experimental results, and other unstructured data.

7. **Queue Service**: Redis-based message queue for handling asynchronous tasks such as CSV processing and property calculation.

8. **Worker Service**: Celery workers for processing background tasks from the queue.

9. **Notification Service**: Handles real-time notifications and alerts for users.

10. **Monitoring Services**: Prometheus, Grafana, Fluentd, Elasticsearch, and Kibana for system monitoring and observability.

## Containerization Strategy

The system is fully containerized using Docker with the following containers:

1. **Frontend Container**: Serves the React frontend application

2. **Backend API Container**: Runs the FastAPI application

3. **Worker Container**: Executes background tasks using Celery

4. **PostgreSQL Container**: Provides the relational database

5. **Redis Container**: Provides caching and message queue functionality

6. **MinIO Container**: Provides S3-compatible object storage

7. **Nginx Container**: Serves as a reverse proxy and API gateway

8. **Monitoring Containers**: Prometheus, Grafana, Fluentd, Elasticsearch, and Kibana

All containers are orchestrated using Docker Compose, with configuration defined in `infrastructure/docker-compose.yml`. This approach ensures consistent deployment across environments and eliminates external dependencies.

## Network Architecture

The containerized services are organized into the following networks:

1. **Frontend Network**: Contains the frontend container and Nginx reverse proxy

2. **Backend Network**: Contains the backend API container, worker container, and supporting services

3. **Database Network**: Contains the PostgreSQL container and services that need direct database access

4. **Monitoring Network**: Contains monitoring and observability services

5. **Logging Network**: Contains logging and log aggregation services

This network segmentation provides security isolation between components while allowing necessary communication paths.

# Frontend Architecture

## Component Structure

The frontend application follows a feature-based organization with shared components:

```
src/web/
├── src/
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
│   ├── App.tsx           # Main application component
│   ├── index.tsx         # Application entry point
│   └── routes.tsx        # Routing configuration
```

Each feature module contains its own components, hooks, and pages, promoting separation of concerns and code organization.

## State Management

The frontend uses Redux Toolkit for global state management with the following slices:

1. **Auth Slice**: Manages authentication state, user information, and permissions

2. **Molecules Slice**: Manages molecular data, filters, and selection state

3. **Libraries Slice**: Manages library data and organization

4. **Experiments Slice**: Manages experiment definitions and status

5. **Submissions Slice**: Manages CRO submission workflow state

6. **Results Slice**: Manages experimental result data

7. **UI Slice**: Manages UI state such as sidebar visibility, active view, and modal state

8. **Notifications Slice**: Manages user notifications and alerts

For server state management and data fetching, the application uses React Query, which provides caching, background updates, and optimistic UI updates.

## Routing and Navigation

The application uses React Router for client-side routing with the following route structure:

1. **Public Routes**:
   - `/login`: User login
   - `/register`: User registration
   - `/forgot-password`: Password recovery

2. **Pharma User Routes**:
   - `/dashboard`: User dashboard
   - `/molecules`: Molecule management
   - `/libraries`: Library management
   - `/experiments`: Experiment management
   - `/submissions`: CRO submissions
   - `/results`: Experimental results

3. **CRO User Routes**:
   - `/cro/dashboard`: CRO dashboard
   - `/cro/submissions`: Submission management
   - `/cro/experiments`: Experiment management
   - `/cro/results`: Result upload

4. **Admin Routes**:
   - `/admin/dashboard`: Admin dashboard
   - `/admin/users`: User management
   - `/admin/system`: System monitoring

All routes (except public routes) are protected by authentication and authorization middleware that verifies the user's role and permissions.

## API Integration

The frontend communicates with the backend API using a custom API client built on top of Axios. The API client handles:

1. **Authentication**: Automatically attaches JWT tokens to requests

2. **Error Handling**: Processes API errors and provides consistent error responses

3. **Request/Response Transformation**: Transforms data between API and application formats

4. **Caching**: Integrates with React Query for efficient caching

5. **Retry Logic**: Implements exponential backoff for failed requests

API endpoints are organized by domain (molecules, libraries, experiments, etc.) with consistent patterns for CRUD operations.

## UI Component Library

The frontend uses Material-UI as its component library, with a custom theme that provides:

1. **Consistent Styling**: Typography, colors, spacing, and component styling

2. **Responsive Design**: Adapts to different screen sizes and devices

3. **Accessibility**: WCAG 2.1 AA compliance with proper contrast and screen reader support

4. **Dark/Light Mode**: User-selectable theme preference

Custom components extend Material-UI to provide domain-specific functionality such as molecule visualization, property display, and experimental data visualization.

# Backend Architecture

## API Layer

The backend API is implemented using FastAPI with the following structure:

```
src/backend/app/
├── api/
│   └── api_v1/
│       ├── endpoints/    # API endpoint handlers
│       │   ├── auth.py
│       │   ├── users.py
│       │   ├── molecules.py
│       │   ├── libraries.py
│       │   ├── experiments.py
│       │   ├── submissions.py
│       │   └── results.py
│       ├── schemas/      # Request/response schemas
│       └── api.py        # API router configuration
```

The API follows RESTful principles with these characteristics:

1. **Versioned Endpoints**: All endpoints are prefixed with `/api/v1/`

2. **Resource-Based Routes**: Endpoints are organized around resources (molecules, libraries, etc.)

3. **Standard HTTP Methods**: Uses GET, POST, PUT, DELETE for CRUD operations

4. **JSON Payloads**: All requests and responses use JSON format

5. **Comprehensive Validation**: Input validation using Pydantic schemas

6. **Automatic Documentation**: OpenAPI documentation via Swagger UI and ReDoc

7. **Authentication**: JWT-based authentication with role-based access control

## Service Layer

The service layer implements the business logic of the application with the following services:

```
src/backend/app/services/
├── auth_service.py       # Authentication and authorization
├── admin_service.py      # Administrative functions
├── csv_service.py        # CSV processing
├── molecule_service.py   # Molecule management
├── library_service.py    # Library management
├── experiment_service.py # Experiment management
├── submission_service.py # CRO submission management
├── result_service.py     # Result management
├── notification_service.py # User notifications
└── file_storage_service.py # File storage operations
```

Each service encapsulates domain-specific business logic and coordinates interactions between different components of the system. Services are designed to be stateless and use dependency injection for testability.

## Data Access Layer

The data access layer is implemented using SQLAlchemy ORM with the following structure:

```
src/backend/app/
├── models/              # SQLAlchemy models
│   ├── user.py
│   ├── molecule.py
│   ├── molecule_property.py
│   ├── library.py
│   ├── library_molecule.py
│   ├── experiment.py
│   ├── experiment_molecule.py
│   ├── submission.py
│   └── result.py
├── crud/                # CRUD operations
│   ├── base.py
│   ├── crud_user.py
│   ├── crud_molecule.py
│   ├── crud_library.py
│   ├── crud_experiment.py
│   ├── crud_submission.py
│   └── crud_result.py
└── db/                  # Database configuration
    ├── base.py
    ├── session.py
    └── init_db.py
```

The data access layer follows the repository pattern, providing a clean abstraction over database operations. Each model has a corresponding CRUD module that implements standard operations (create, read, update, delete) as well as domain-specific queries.

## Molecular Processing Engine

The Molecular Processing Engine is a specialized component for handling molecular data processing:

```
src/backend/app/molecular/
├── molecule_converter.py  # SMILES conversion utilities
├── validator.py          # Molecular structure validation
├── processor.py          # Batch processing of molecules
├── property_calculator.py # Molecular property calculation
├── similarity_searcher.py # Molecular similarity search
└── substructure_searcher.py # Substructure search
```

This component leverages RDKit for cheminformatics operations and implements:

1. **SMILES Validation**: Verifies the correctness of molecular structures

2. **Property Calculation**: Computes molecular properties such as LogP, molecular weight, etc.

3. **Batch Processing**: Efficiently processes large sets of molecules in parallel

4. **Similarity Search**: Finds molecules similar to a query structure

5. **Substructure Search**: Identifies molecules containing specific substructures

The engine is designed for performance, with parallel processing capabilities and caching of computed properties.

## Background Processing

The system uses Celery with Redis as the message broker for background processing:

```
src/backend/app/worker/
├── celery_app.py        # Celery configuration
├── celery_config.py     # Celery settings
└── tasks/               # Task definitions
    ├── csv_tasks.py     # CSV processing tasks
    ├── molecule_tasks.py # Molecular processing tasks
    ├── notification_tasks.py # Notification delivery tasks
    └── report_tasks.py  # Report generation tasks
```

Background processing is used for computationally intensive or time-consuming operations such as:

1. **CSV Processing**: Parsing and importing large CSV files

2. **Property Calculation**: Computing molecular properties for large datasets

3. **Notification Delivery**: Sending email notifications

4. **Report Generation**: Creating complex reports and exports

Tasks are queued in Redis and processed by Celery worker containers, allowing the API to remain responsive during intensive operations.

## Security Implementation

The backend implements a comprehensive security model:

```
src/backend/app/core/
├── security.py          # Security utilities
├── auth.py              # Authentication logic
└── jwt.py               # JWT token handling
```

Security features include:

1. **Password Hashing**: Secure password storage using bcrypt

2. **JWT Authentication**: Token-based authentication with short-lived access tokens

3. **Role-Based Access Control**: Permission enforcement based on user roles

4. **Input Validation**: Comprehensive validation of all API inputs

5. **Rate Limiting**: Protection against brute force and DoS attacks

6. **CORS Protection**: Controlled cross-origin resource sharing

7. **Content Security Policy**: Protection against XSS attacks

All security mechanisms are implemented locally without external dependencies, ensuring the system can operate in isolated environments.

# Database Design

## Schema Overview

The database schema is designed to support the molecular data management and CRO integration workflows with the following key entities:

1. **Users**: User accounts and authentication data

2. **Molecules**: Molecular structures and basic properties

3. **Molecule Properties**: Flexible property storage for molecules

4. **Libraries**: User-defined collections of molecules

5. **Experiments**: Experimental definitions and configurations

6. **Submissions**: CRO submission records and status

7. **Results**: Experimental result data and metadata

The schema uses a relational model with PostgreSQL as the database management system, leveraging its support for JSON/JSONB data types for flexible property storage.

## Entity Relationships

The key entity relationships in the database are:

```
+-------+     +------------+     +---------+
| Users | 1:N | Molecules  | M:N | Library |
+-------+     +------------+     +---------+
    |              |  |             |
    |              |  |             |
    |              |  |             |
    | 1:N          |  |             | 1:N
    |              |  |             |
    v              |  |             v
+------------+     |  |      +------------+
| Experiment | M:N +--+      | LibraryMol |
+------------+     |         +------------+
    |              |
    | 1:N          | M:N
    |              |
    v              v
+------------+     +------------+
| Submission |     | ExperMol   |
+------------+     +------------+
    |
    | 1:N
    |
    v
+------------+
| Result     |
+------------+