# Introduction

This document provides comprehensive instructions for deploying, configuring, and maintaining the Molecular Data Management and CRO Integration Platform. The platform is designed for fully local deployment without external dependencies, enabling small to mid-cap pharmaceutical companies to manage molecular data and interact with Contract Research Organizations (CROs) in a secure and efficient manner.

The deployment process is containerized using Docker and Docker Compose, ensuring consistent deployment across different environments while maintaining isolation between components. This document covers initial deployment, updates, scaling, backup and recovery, and routine maintenance procedures.

# System Requirements

Before deploying the platform, ensure that your environment meets the following requirements:

## Hardware Requirements

The platform can be deployed on various scales depending on the expected usage. The following table provides guidelines for resource allocation:

| Deployment Scale | Users | Molecules | CPU | Memory | Storage | Network |
|------------------|-------|-----------|-----|--------|---------|----------|
| Small | <10 | <50,000 | 4 cores | 8GB | 100GB | 1 Gbps |
| Medium | 10-25 | 50,000-250,000 | 8 cores | 16GB | 500GB | 1 Gbps |
| Large | 25-50 | 250,000-1,000,000 | 16 cores | 32GB | 1TB | 10 Gbps |

These requirements ensure that the system can handle the expected load while maintaining performance targets. For deployments exceeding these guidelines, additional resources should be allocated according to scaling recommendations.

## Software Requirements

The following software is required on the host system:

| Software | Version | Purpose |
|----------|---------|----------|
| Docker | 23.0+ | Container runtime for application deployment |
| Docker Compose | 2.17+ | Multi-container orchestration |
| Git | Latest | Source code management (for updates) |
| Bash | Latest | Shell for running deployment scripts |
| OpenSSL | Latest | Generation of security certificates and keys |
| Curl | Latest | Network requests for health checks |

All other dependencies are included in the Docker containers and do not need to be installed separately on the host system.

## Network Requirements

The platform requires the following network configuration:

| Requirement | Description |
|------------|-------------|
| Ports | 80/443 (HTTP/HTTPS) exposed to users |
| Internal Network | Communication between containers |
| DNS | Local hostname resolution |
| Firewall | Allow inbound connections to HTTP/HTTPS ports |

The platform is designed to operate in a local network environment without external dependencies. Internet access is not required for operation but may be needed for initial setup to pull Docker images.

## Security Requirements

The platform should be deployed in a secure environment with the following considerations:

| Requirement | Description |
|------------|-------------|
| Physical Security | Secure access to host servers |
| Network Security | Firewall protection, network segmentation |
| Access Control | Limited SSH access to host servers |
| Data Protection | Encrypted storage for sensitive data |
| Backup Security | Secure storage of backup files |

Implement the following security practices for the platform:

1. **Authentication and Authorization:**
   - Enforce strong password policies (minimum 10 characters, complexity requirements)
   - Implement role-based access control for all users
   - Configure JWT token expiration (15 minutes for access tokens)
   - Rotate refresh tokens regularly (7 days)

2. **Transport Security:**
   - Configure TLS 1.3 with strong cipher suites
   - Implement proper certificate management
   - Configure HTTP security headers (HSTS, CSP, X-Content-Type-Options)
   - Enable HTTP to HTTPS redirection

3. **Data Protection:**
   - Enable database encryption at rest
   - Secure sensitive files with appropriate permissions
   - Implement secure backup encryption
   - Configure secure deletion policies for temporary files

4. **Container Security:**
   - Use minimal base images to reduce attack surface
   - Run containers with non-root users where possible
   - Scan container images for vulnerabilities regularly
   - Implement proper resource limits to prevent DoS attacks

5. **Monitoring and Auditing:**
   - Configure security event logging and alerting
   - Implement failed login attempt monitoring
   - Enable audit logging for sensitive operations
   - Regularly review security logs for unusual activity

# System Architecture Overview

The Molecular Data Management and CRO Integration Platform employs a containerized microservices architecture to enable full local deployment while maintaining separation of concerns. This section provides an overview of the system architecture to help understand the deployment model.

## Core Components

The platform consists of the following core components:

| Component | Primary Responsibility | Key Dependencies |
|----------|------------------------|------------------|
| Frontend Application | User interface for molecule management and CRO interactions | Backend API, Authentication Service |
| Backend API | Core business logic, data validation, and service orchestration | Database, File Storage, Molecular Processing Engine |
| Authentication Service | User management, role-based access control | Database |
| Molecular Processing Engine | SMILES validation, property calculation, molecular analysis | RDKit |
| Database Service | Persistent storage for molecular data, user information, and experiment tracking | None |
| File Storage Service | Management of CSV files, experimental results, and documents | None |
| Queue Service | Asynchronous task processing for long-running operations | Backend API |
| Notification Service | User alerts for experiment status changes and communications | Backend API, Queue Service |

These components are deployed as Docker containers and orchestrated using Docker Compose, allowing for a fully local deployment without external dependencies.

## Container Architecture

The platform consists of the following containers:

| Container | Purpose | Base Image |
|-----------|---------|------------|
| Frontend | Serves the React frontend application | node:18-alpine |
| Backend API | Runs the FastAPI application | python:3.10-slim |
| Worker | Executes background tasks using Celery | python:3.10-slim |
| PostgreSQL | Provides the relational database | postgres:15-alpine |
| Redis | Provides caching and message queue | redis:7-alpine |
| MinIO | Provides S3-compatible object storage | minio/minio:latest |
| Nginx | Serves as a reverse proxy and API gateway | nginx:1.24-alpine |
| Prometheus | Collects and stores metrics | prom/prometheus:latest |
| Grafana | Provides visualization dashboards | grafana/grafana:latest |
| Fluentd | Collects and aggregates logs | fluent/fluentd:latest |
| Elasticsearch | Stores and indexes logs | elasticsearch:7.17.0 |
| Kibana | Provides log visualization | kibana:7.17.0 |

These containers are defined in the `infrastructure/docker-compose.yml` file and are orchestrated together to form the complete platform.

## Network Architecture

The containers are organized into the following networks:

| Network | Purpose | Connected Containers |
|---------|---------|----------------------|
| frontend_network | Frontend communication | nginx, frontend, grafana, kibana |
| backend_network | Backend services | nginx, backend, worker, redis, minio, prometheus |
| database_network | Database access | backend, worker, postgres, redis, minio |
| monitoring_network | Monitoring services | prometheus, grafana |
| logging_network | Logging services | fluentd, elasticsearch, kibana |

This network segmentation provides security isolation between components while allowing necessary communication paths.

## Data Flow

The system's data flow begins with CSV uploads containing molecular data. These files are processed by the Backend API, which validates the format and extracts SMILES strings and properties. The Molecular Processing Engine validates the SMILES structures and calculates any missing properties. Validated molecular data is then stored in the Database Service.

Users interact with molecules through the Frontend Application, which communicates with the Backend API to retrieve, filter, and organize molecular data. When molecules are selected for experimental testing, they are added to experiment queues in the Database Service. The Backend API manages the submission of these experiments to CRO users.

CRO users receive notifications about new submissions through the Notification Service. They interact with the Frontend Application to review submissions, provide pricing, and upload results. These results are stored in the File Storage Service, with metadata in the Database Service. The Notification Service alerts pharma users about new results.

Throughout these flows, the Authentication Service verifies user identity and permissions for each operation. The Queue Service handles resource-intensive tasks asynchronously, such as processing large CSV files or generating complex reports.

# Deployment Architecture

The platform is deployed as a set of Docker containers orchestrated by Docker Compose. This section provides an overview of the deployment architecture.

## Container Architecture

The platform consists of the following containers:

| Container | Purpose | Base Image |
|-----------|---------|------------|
| Frontend | Serves the React frontend application | node:18-alpine |
| Backend API | Runs the FastAPI application | python:3.10-slim |
| Worker | Executes background tasks using Celery | python:3.10-slim |
| PostgreSQL | Provides the relational database | postgres:15-alpine |
| Redis | Provides caching and message queue | redis:7-alpine |
| MinIO | Provides S3-compatible object storage | minio/minio:latest |
| Nginx | Serves as a reverse proxy and API gateway | nginx:1.24-alpine |
| Prometheus | Collects and stores metrics | prom/prometheus:latest |
| Grafana | Provides visualization dashboards | grafana/grafana:latest |
| Fluentd | Collects and aggregates logs | fluent/fluentd:latest |
| Elasticsearch | Stores and indexes logs | elasticsearch:7.17.0 |
| Kibana | Provides log visualization | kibana:7.17.0 |

These containers are defined in the `infrastructure/docker-compose.yml` file and are orchestrated together to form the complete platform.

## Network Architecture

The containers are organized into the following networks:

| Network | Purpose | Connected Containers |
|---------|---------|----------------------|
| frontend_network | Frontend communication | nginx, frontend, grafana, kibana |
| backend_network | Backend services | nginx, backend, worker, redis, minio, prometheus |
| database_network | Database access | backend, worker, postgres, redis, minio |
| monitoring_network | Monitoring services | prometheus, grafana |
| logging_network | Logging services | fluentd, elasticsearch, kibana |

This network segmentation provides security isolation between components while allowing necessary communication paths.

## Volume Architecture

The platform uses the following Docker volumes for persistent storage:

| Volume | Purpose | Connected Containers |
|--------|---------|----------------------|
| postgres_data | PostgreSQL database files | postgres |
| redis_data | Redis data | redis |
| minio_data | MinIO object storage | minio |
| backend_data | Backend application data | backend |
| worker_data | Worker application data | worker |
| nginx_logs | Nginx log files | nginx, fluentd |
| prometheus_data | Prometheus metrics | prometheus |
| grafana_data | Grafana dashboards and settings | grafana |
| elasticsearch_data | Elasticsearch indices | elasticsearch |
| fluentd_data | Fluentd logs | fluentd |

These volumes ensure that data persists across container restarts and system reboots.

## Port Mapping

The platform exposes the following ports to the host system:

| Container | Host Port | Container Port | Purpose |
|-----------|-----------|----------------|----------|
| Nginx | 80 | 80 | HTTP access |
| Nginx | 443 | 443 | HTTPS access |
| Prometheus | 9090 | 9090 | Metrics UI (optional) |
| Grafana | 3000 | 3000 | Monitoring dashboards (optional) |
| Kibana | 5601 | 5601 | Log analysis (optional) |

Only the HTTP/HTTPS ports (80/443) need to be exposed to users. The other ports are optional and can be restricted to administrators or not exposed at all.

# Initial Deployment

This section provides step-by-step instructions for the initial deployment of the platform.

## Prerequisites

Before starting the deployment, ensure that:

1. The host system meets all hardware and software requirements
2. You have administrative access to the host system
3. Docker and Docker Compose are properly installed and configured
4. Required ports are available and not used by other services
5. Sufficient disk space is available for the deployment

## Obtaining the Software

Clone the repository to the host system:

```bash
git clone https://github.com/your-organization/molecular-platform.git
cd molecular-platform
```

Alternatively, you can download and extract a release archive from the project's release page.

## Configuration

The platform uses environment variables for configuration. These are stored in `.env` files in different directories:

1. `infrastructure/.env`: Configuration for Docker Compose services
2. `src/backend/.env`: Configuration for the backend API
3. `src/web/.env`: Configuration for the frontend application

The setup script will create these files from templates if they don't exist. You can also create them manually by copying the corresponding `.env.example` files and adjusting the values as needed.

Key configuration parameters include:

| Parameter | Description | Default | Recommended |
|-----------|-------------|---------|-------------|
| POSTGRES_USER | PostgreSQL username | postgres | Custom secure username |
| POSTGRES_PASSWORD | PostgreSQL password | postgres | Strong random password |
| POSTGRES_DB | PostgreSQL database name | molecular_platform | No change needed |
| MINIO_ROOT_USER | MinIO root username | minioadmin | Custom secure username |
| MINIO_ROOT_PASSWORD | MinIO root password | minioadmin | Strong random password |
| SECRET_KEY | JWT signing key | Generated | Auto-generated secure key |
| APP_ENVIRONMENT | Environment name | production | No change needed |
| EXTERNAL_HTTP_PORT | HTTP port mapping | 80 | Adjust if port 80 is in use |
| EXTERNAL_HTTPS_PORT | HTTPS port mapping | 443 | Adjust if port 443 is in use |

Additional configuration options are documented in the example files.

## Running the Setup Script

The platform includes a setup script that automates the initial deployment process. Run the script with the following command:

```bash
./scripts/setup.sh
```

The setup script performs the following actions:

1. Checks for required dependencies
2. Creates necessary directories
3. Generates JWT keys for authentication
4. Configures environment variables
5. Builds Docker images
6. Initializes the database
7. Initializes MinIO object storage
8. Starts all services
9. Verifies the deployment

You can customize the setup process with the following options:

```bash
./scripts/setup.sh -h  # Show help
./scripts/setup.sh -e  # Include example data
./scripts/setup.sh -a admin@example.com -p StrongPassword123  # Set admin credentials
```

The setup process may take 10-15 minutes depending on your system's performance and internet connection speed.

## Verifying the Deployment

After the setup script completes, verify that the deployment was successful:

1. Check that all containers are running:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps
   ```

2. Verify that the web interface is accessible:
   - Open a web browser and navigate to `http://localhost` (or the configured HTTP port)
   - You should see the login page of the platform

3. Verify that the API is working:
   ```bash
   curl -s http://localhost/api/v1/health/live | grep UP
   ```

4. Log in with the admin credentials:
   - Use the email and password specified during setup
   - If not specified, check the setup log for the generated admin credentials

If any of these verification steps fail, check the troubleshooting section for guidance.

## Post-Installation Steps

After verifying the deployment, perform the following post-installation steps:

1. Change the default admin password:
   - Log in as admin
   - Navigate to the user profile
   - Change the password to a strong, unique password

2. Configure TLS/SSL for secure communication:
   - Generate or obtain SSL certificates
   - Place certificates in `infrastructure/nginx/ssl/`
   - Update Nginx configuration in `infrastructure/nginx/conf.d/default.conf`
   - Restart the Nginx container

3. Configure backup schedule:
   - Set up a cron job to run the backup script regularly
   - Example: `0 2 * * * /path/to/molecular-platform/scripts/backup.sh`

4. Review and adjust monitoring settings:
   - Access Grafana at `http://localhost:3000`
   - Review default dashboards
   - Configure alert notifications if needed

5. Document the deployment:
   - Record configuration choices
   - Document any customizations
   - Store credentials securely

# Updating the Platform

This section provides instructions for updating an existing deployment to a new version.

## Update Process Overview

The update process involves the following steps:

1. Backup the current deployment
2. Update the source code
3. Build new Docker images
4. Apply database migrations
5. Restart services
6. Verify the update

The platform includes a deployment script that automates this process.

## Preparing for the Update

Before updating, perform the following preparations:

1. Review the release notes for the new version
2. Check for any breaking changes or special update instructions
3. Ensure that you have sufficient disk space for the update
4. Schedule the update during a maintenance window
5. Notify users of the planned downtime

## Running the Update

Run the deployment script in update mode:

```bash
./scripts/deploy.sh -m update
```

The script performs the following actions:

1. Creates a backup of the current deployment
2. Pulls the latest code from the repository
3. Stops the running services
4. Builds new Docker images
5. Runs database migrations
6. Starts the services with the new images
7. Verifies the deployment

You can customize the update process with the following options:

```bash
./scripts/deploy.sh -h  # Show help
./scripts/deploy.sh -m update -s backup  # Skip backup
./scripts/deploy.sh -m update -s build  # Skip image building
./scripts/deploy.sh -m update -s migrations  # Skip database migrations
./scripts/deploy.sh -m update -f  # Force recreation of containers
```

The update process may take 5-10 minutes depending on your system's performance and the extent of the changes.

## Verifying the Update

After the update completes, verify that it was successful:

1. Check that all containers are running:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml ps
   ```

2. Verify that the web interface is accessible:
   - Open a web browser and navigate to `http://localhost` (or the configured HTTP port)
   - You should see the login page of the platform

3. Verify that the API is working:
   ```bash
   curl -s http://localhost/api/v1/health/live | grep UP
   ```

4. Log in and verify that your data is intact

If any of these verification steps fail, the deployment script will offer to roll back to the previous version using the backup created before the update.

## Rolling Back an Update

If you need to roll back to a previous version after an update, you can use the restore script:

```bash
./scripts/restore.sh -b backups/pre_deploy_YYYYMMDD_HHMMSS
```

Replace `YYYYMMDD_HHMMSS` with the timestamp of the backup you want to restore. The backup directory is created automatically by the deployment script before updating.

The restore process will:

1. Stop the current services
2. Restore the database from the backup
3. Restore the file storage from the backup
4. Restore the configuration from the backup
5. Start the services using the backed-up images

After the restore completes, verify that the rollback was successful using the same verification steps as for an update.

# Scaling the Platform

This section provides guidance on scaling the platform to accommodate growing usage.

## Vertical Scaling

The primary scaling approach for the platform is vertical scaling (increasing resources for existing containers). You can adjust resource allocation in the Docker Compose configuration:

```yaml
services:
  backend:
    deploy:
      resources:
        limits:
          cpus: '2'
          memory: 2G
        reservations:
          cpus: '1'
          memory: 1G
```

Recommended resource allocations for different components:

| Component | CPU Allocation | Memory Allocation | Disk I/O Priority |
|-----------|----------------|-------------------|-------------------|
| Frontend | Low (0.5-1 core) | Medium (1GB) | Low |
| Backend API | High (1-2 cores) | High (2GB) | Medium |
| Database | High (2-4 cores) | High (4-8GB) | High |
| Molecular Processing | Very High (2-4 cores) | High (2-4GB) | Medium |
| File Storage | Medium (1-2 cores) | Medium (1-2GB) | Very High |

After adjusting resource allocations, restart the affected services:

```bash
docker-compose -f infrastructure/docker-compose.yml up -d --no-deps <service_name>
```

## Horizontal Scaling

For stateless services (API, workers), you can implement horizontal scaling by increasing the number of container replicas. This requires additional configuration:

1. Add a load balancer (e.g., HAProxy) in front of the API containers
2. Configure the load balancer to distribute requests
3. Adjust the Docker Compose configuration to create multiple replicas

Example configuration for horizontal scaling:

```yaml
services:
  backend:
    deploy:
      replicas: 3
  worker:
    deploy:
      replicas: 2
```

Note that horizontal scaling requires Docker Swarm or Kubernetes for orchestration, which is beyond the scope of this basic deployment guide.

## Database Scaling

As your dataset grows, you may need to scale the database:

1. Increase database resources (CPU, memory)
2. Optimize database configuration for larger datasets
3. Implement database replication for read scaling

Database optimization parameters for PostgreSQL:

```
shared_buffers = 2GB
effective_cache_size = 6GB
maintenance_work_mem = 512MB
work_mem = 32MB
max_connections = 100
random_page_cost = 1.1
effective_io_concurrency = 200
```

These parameters should be adjusted based on your specific hardware and usage patterns. Add them to `infrastructure/postgres/postgresql.conf` and restart the PostgreSQL container.

## Storage Scaling

To scale storage capacity:

1. Increase volume sizes in Docker Compose configuration
2. Add additional storage volumes if needed
3. Configure MinIO for distributed storage (advanced)

Example of increasing volume size:

```bash
# Stop the services
docker-compose -f infrastructure/docker-compose.yml down

# Backup the data
./scripts/backup.sh

# Create a new, larger volume
docker volume create --name=minio_data_new

# Copy data from old volume to new volume (requires additional tools)
# ...

# Update docker-compose.yml to use the new volume
# ...

# Start the services
docker-compose -f infrastructure/docker-compose.yml up -d
```

For advanced storage scaling, consider implementing a distributed MinIO deployment with multiple nodes.

## Monitoring Scaling Requirements

Use the monitoring system to identify scaling requirements:

1. Monitor resource utilization (CPU, memory, disk, network)
2. Track performance metrics (response time, processing time)
3. Analyze usage patterns and growth trends

Key metrics to watch for scaling decisions:

| Metric | Warning Threshold | Action |
|--------|-------------------|--------|
| CPU Usage | >70% sustained | Increase CPU allocation |
| Memory Usage | >75% sustained | Increase memory allocation |
| Disk Usage | >70% of allocated | Increase storage allocation |
| API Response Time | >300ms average | Scale API resources |
| Database Connections | >70% of max | Increase connection limit |
| Queue Length | >100 tasks sustained | Add worker containers |

The Grafana dashboards provide visualizations of these metrics to help with scaling decisions.

# Backup and Recovery

This section provides instructions for backing up and recovering the platform.

## Backup Strategy

The platform implements a comprehensive backup strategy:

1. Database Backup: PostgreSQL data is backed up using pg_dump
2. File Storage Backup: MinIO data is backed up using mc mirror
3. Configuration Backup: Environment files and configurations are backed up
4. Container Image Backup: Docker images are saved for version consistency

Backups are stored in the `backups` directory with timestamped subdirectories. Each backup includes metadata about the version, configuration, and contents.

## Creating Backups

To create a backup, use the backup script:

```bash
./scripts/backup.sh
```

This script performs the following actions:

1. Creates a timestamped backup directory
2. Backs up the PostgreSQL database
3. Backs up MinIO object storage
4. Backs up configuration files
5. Saves Docker image information
6. Creates a backup manifest

You can customize the backup process with the following options:

```bash
./scripts/backup.sh -h  # Show help
./scripts/backup.sh -d /path/to/backup/dir  # Specify backup directory
./scripts/backup.sh -t full  # Full backup (default)
./scripts/backup.sh -t db  # Database-only backup
./scripts/backup.sh -t files  # File storage-only backup
```

It is recommended to schedule regular backups using cron or another scheduler.

## Backup Retention

Implement a backup retention policy to manage backup storage:

1. Daily backups: Keep for 7 days
2. Weekly backups: Keep for 4 weeks
3. Monthly backups: Keep for 3 months

You can implement this policy using a script that runs after the backup process:

```bash
# Example retention script
find /path/to/backups -name "backup_daily_*" -type d -mtime +7 -exec rm -rf {} \;
find /path/to/backups -name "backup_weekly_*" -type d -mtime +28 -exec rm -rf {} \;
find /path/to/backups -name "backup_monthly_*" -type d -mtime +90 -exec rm -rf {} \;
```

Adjust the retention periods based on your organization's requirements and available storage.

## Restoring from Backup

To restore the platform from a backup, use the restore script:

```bash
./scripts/restore.sh -b backups/backup_YYYYMMDD_HHMMSS
```

Replace `YYYYMMDD_HHMMSS` with the timestamp of the backup you want to restore.

The restore script performs the following actions:

1. Stops the current services
2. Restores the database from the backup
3. Restores the file storage from the backup
4. Restores the configuration from the backup
5. Starts the services

You can customize the restore process with the following options:

```bash
./scripts/restore.sh -h  # Show help
./scripts/restore.sh -b backups/backup_YYYYMMDD_HHMMSS -t full  # Full restore (default)
./scripts/restore.sh -b backups/backup_YYYYMMDD_HHMMSS -t db  # Database-only restore
./scripts/restore.sh -b backups/backup_YYYYMMDD_HHMMSS -t files  # File storage-only restore
```

After the restore completes, verify that the platform is functioning correctly.

## Disaster Recovery

For disaster recovery scenarios (complete system failure), follow these steps:

1. Install the operating system and prerequisites on the new host
2. Install Docker and Docker Compose
3. Clone the repository or extract the release archive
4. Copy the backup to the new host
5. Run the restore script to restore from the backup

The Recovery Time Objective (RTO) for this process is approximately 4 hours, depending on the size of the backup and the performance of the new host.

The Recovery Point Objective (RPO) depends on your backup frequency. With daily backups, the RPO is 24 hours.

# Routine Maintenance

This section provides guidance on routine maintenance tasks for the platform.

## Database Maintenance

Perform the following database maintenance tasks:

1. Regular Vacuuming: PostgreSQL's autovacuum should handle this automatically, but you can manually vacuum the database if needed:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "VACUUM ANALYZE;"
   ```

2. Index Maintenance: Rebuild indexes periodically to maintain performance:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "REINDEX DATABASE molecular_platform;"
   ```

3. Database Statistics: Update statistics to improve query planning:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml exec postgres psql -U postgres -d molecular_platform -c "ANALYZE;"
   ```

These tasks can be scheduled to run during off-peak hours.

## Log Management

Manage logs to prevent disk space issues:

1. Log Rotation: The platform uses Docker's built-in log rotation, but you can adjust the settings in `daemon.json`:
   ```json
   {
     "log-driver": "json-file",
     "log-opts": {
       "max-size": "10m",
       "max-file": "3"
     }
   }
   ```

2. Log Cleanup: Periodically clean up old logs:
   ```bash
   docker system prune -f --volumes
   ```
   Note: This command removes unused volumes, so use with caution.

3. Log Analysis: Use Kibana to analyze logs for patterns and issues:
   - Access Kibana at `http://localhost:5601`
   - Use the pre-configured dashboards for log analysis
   - Create custom queries for specific investigations

## Security Updates

Keep the platform secure with regular updates:

1. Container Base Images: Update base images regularly:
   ```bash
   docker-compose -f infrastructure/docker-compose.yml pull
   docker-compose -f infrastructure/docker-compose.yml up -d
   ```

2. Application Updates: Follow the update process described in the "Updating the Platform" section.

3. Dependency Scanning: Regularly scan for vulnerabilities in dependencies:
   ```bash
   # For Python dependencies
   docker-compose -f infrastructure/docker-compose.yml exec backend pip-audit
   
   # For JavaScript dependencies
   docker-compose -f infrastructure/docker-compose.yml exec frontend npm audit
   ```

4. Security Patches: Apply security patches promptly when they become available.

## Performance Tuning

Optimize performance based on monitoring data:

1. Database Query Optimization: Identify and optimize slow queries:
   - Review slow query logs
   - Add indexes for frequently queried columns
   - Optimize complex queries

2. API Performance Tuning: Optimize API endpoints:
   - Identify slow endpoints using monitoring data
   - Implement caching for frequently accessed data
   - Optimize database queries used by endpoints

3. Frontend Performance: Optimize frontend performance:
   - Analyze loading times using browser developer tools
   - Optimize asset sizes and loading
   - Implement client-side caching

Use the monitoring dashboards to identify performance bottlenecks and measure the impact of optimizations.

## Health Monitoring

Regularly monitor system health:

1. Automated Health Checks: Use the health monitoring script for regular checks:
   ```bash
   ./scripts/monitor-health.sh -c -i 300  # Run continuous monitoring at 5-minute intervals