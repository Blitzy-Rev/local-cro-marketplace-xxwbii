# Molecular Data Management and CRO Integration Platform

A comprehensive platform for streamlining the small molecule drug discovery process for small to mid-cap pharmaceutical companies. This system bridges the gap between computational data aggregation and Contract Research Organization (CRO) services by providing an intuitive interface for molecular data management, organization, and experimental workflow submission.

## Key Features
- **Fully Local Deployment**: Operates entirely within your organization's infrastructure with no external dependencies
- **CSV Upload & Molecular Data Ingestion**: Easily import molecular data with flexible mapping capabilities
- **Interactive Molecule Organization**: Sort, filter, and organize molecules into custom libraries
- **Experiment Definition & Queuing**: Create and configure experiments for submission to CROs
- **CRO Integration**: Streamlined submission process with bidirectional communication
- **Result Management**: Comprehensive tools for reviewing and analyzing experimental results
- **Role-Based Access Control**: Dedicated interfaces for pharma users, CRO users, and administrators
- **Containerized Architecture**: Docker-based deployment for simplified installation and maintenance

## System Architecture
The platform employs a containerized microservices architecture to enable full local deployment while maintaining separation of concerns:

- **Frontend**: React-based single-page application with TypeScript and Material-UI
- **Backend API**: FastAPI application providing RESTful endpoints
- **Database**: PostgreSQL for structured data storage
- **File Storage**: MinIO for S3-compatible object storage
- **Cache & Queue**: Redis for caching and asynchronous task processing
- **Monitoring**: Prometheus, Grafana, Fluentd, Elasticsearch, and Kibana

All components are containerized using Docker and orchestrated with Docker Compose, ensuring consistent deployment across environments.

For detailed architecture information, see [Architecture Documentation](docs/architecture.md).

## Prerequisites
- **Docker**: Version 23.0 or higher
- **Docker Compose**: Version 2.17 or higher
- **Hardware Requirements**:
  - Small Deployment (< 50,000 molecules): 4 cores, 8GB RAM, 100GB storage
  - Medium Deployment (50,000-250,000 molecules): 8 cores, 16GB RAM, 500GB storage
  - Large Deployment (250,000-1,000,000 molecules): 16 cores, 32GB RAM, 1TB storage
- **Network**: Local network with ports 80/443 available

## Quick Start
1. Clone this repository:
   ```bash
   git clone https://github.com/your-organization/molecular-platform.git
   cd molecular-platform
   ```

2. Run the setup script:
   ```bash
   ./scripts/setup.sh
   ```

3. Access the platform:
   - Open a web browser and navigate to `http://localhost`
   - Log in with the default admin credentials (see setup output)

For detailed deployment instructions, see [Deployment Guide](docs/deployment.md).

## User Roles
The platform supports three primary user roles:

1. **Pharma User**: Researchers from pharmaceutical companies who manage molecular data, create libraries, define experiments, and review results

2. **CRO User**: Staff from Contract Research Organizations who receive experiment submissions, provide quotes, update experiment status, and upload results

3. **Administrator**: System administrators who manage users, monitor system health, and configure system settings

Each role has a dedicated user interface with role-specific functionality and access controls.

## Documentation
- [User Guide](docs/user-guide.md): Comprehensive guide for all user roles
- [Architecture Documentation](docs/architecture.md): Detailed system architecture
- [Deployment Guide](docs/deployment.md): Installation and configuration instructions
- [API Documentation](docs/api-docs.md): API reference for programmatic access
- [Molecular Data Guidelines](docs/molecular-data-guidelines.md): Guidelines for molecular data formats
- [Monitoring Guide](docs/monitoring.md): Monitoring and observability documentation
- [Maintenance Guide](docs/maintenance.md): Routine maintenance procedures
- [Troubleshooting Guide](docs/troubleshooting.md): Solutions for common issues

## Development
For development setup and guidelines, see [Development Guide](docs/development.md).

The project is organized into the following main directories:

- `src/backend`: FastAPI backend application with Python
- `src/web`: React frontend application with TypeScript
- `infrastructure`: Docker Compose and service configurations
- `scripts`: Deployment and maintenance scripts
- `docs`: Documentation files

## Security
The platform implements comprehensive security measures:

- **Authentication**: JWT-based authentication with secure password storage
- **Authorization**: Role-based access control with fine-grained permissions
- **Data Protection**: Encryption of sensitive data at rest and in transit
- **Input Validation**: Comprehensive validation of all user inputs
- **Audit Logging**: Detailed logging of security-relevant events

For detailed security information, see [Security Documentation](docs/security.md).

## License
This project is licensed under the [MIT License](LICENSE).

## Contributing
Contributions are welcome! Please see [Contributing Guidelines](CONTRIBUTING.md) for details.