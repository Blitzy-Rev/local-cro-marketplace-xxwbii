#!/bin/bash
# Molecular Data Management and CRO Integration Platform
# Deployment Script
# Version: 1.0.0

# Exit immediately if a command exits with a non-zero status
set -e

# Script directories and paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
INFRASTRUCTURE_DIR="$PROJECT_ROOT/infrastructure"
BACKUP_DIR="$PROJECT_ROOT/backups/pre_deploy_$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$SCRIPT_DIR/deploy_$(date +%Y%m%d_%H%M%S).log"

# Required dependencies version
DOCKER_REQUIRED_VERSION="23.0.0"
DOCKER_COMPOSE_REQUIRED_VERSION="2.17.0"

# Default configuration
DEPLOYMENT_MODE="update"  # 'initial' or 'update'
SKIP_BACKUP="false"
SKIP_BUILD="false"
SKIP_MIGRATIONS="false"
FORCE_RECREATE="false"
HEALTH_CHECK_TIMEOUT="300"  # in seconds

# Log messages to console and file
log() {
    local message="$1"
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    
    # Create log directory if needed
    mkdir -p "$(dirname "$LOG_FILE")"
    
    echo "[$timestamp] $message"
    echo "[$timestamp] $message" >> "$LOG_FILE"
}

# Print a deployment banner
print_banner() {
    echo "====================================================================="
    echo "      Molecular Data Management and CRO Integration Platform"
    echo "                     Deployment Script"
    echo ""
    echo "      Mode: $DEPLOYMENT_MODE"
    echo "      Date: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "====================================================================="
}

# Check if required dependencies are installed
check_dependencies() {
    log "Checking dependencies..."
    
    # Check if Docker is installed
    if ! command -v docker &> /dev/null; then
        log "ERROR: Docker is not installed. Please install Docker version $DOCKER_REQUIRED_VERSION or higher."
        return 1
    fi
    
    # Check Docker version
    DOCKER_VERSION=$(docker --version | awk '{print $3}' | sed 's/,//')
    if [ "$(printf '%s\n' "$DOCKER_REQUIRED_VERSION" "$DOCKER_VERSION" | sort -V | head -n1)" != "$DOCKER_REQUIRED_VERSION" ]; then
        log "ERROR: Docker version $DOCKER_VERSION is lower than required version $DOCKER_REQUIRED_VERSION"
        return 1
    else
        log "Docker version $DOCKER_VERSION is installed"
    fi
    
    # Check if Docker Compose is installed
    if ! command -v docker-compose &> /dev/null; then
        log "ERROR: Docker Compose is not installed. Please install Docker Compose version $DOCKER_COMPOSE_REQUIRED_VERSION or higher."
        return 1
    fi
    
    # Check Docker Compose version
    DOCKER_COMPOSE_VERSION=$(docker-compose --version | awk '{print $3}' | sed 's/,//')
    if [ "$(printf '%s\n' "$DOCKER_COMPOSE_REQUIRED_VERSION" "$DOCKER_COMPOSE_VERSION" | sort -V | head -n1)" != "$DOCKER_COMPOSE_REQUIRED_VERSION" ]; then
        log "ERROR: Docker Compose version $DOCKER_COMPOSE_VERSION is lower than required version $DOCKER_COMPOSE_REQUIRED_VERSION"
        return 1
    else
        log "Docker Compose version $DOCKER_COMPOSE_VERSION is installed"
    fi
    
    # Check if Git is installed
    if ! command -v git &> /dev/null; then
        log "WARNING: Git is not installed. Source code updates will not be available."
    else
        log "Git is installed"
    fi
    
    log "All dependencies satisfied."
    return 0
}

# Check if environment files exist and create them if needed
check_environment_files() {
    log "Checking environment files..."
    
    # Check infrastructure/.env file
    if [ ! -f "$INFRASTRUCTURE_DIR/.env" ]; then
        log "Environment file $INFRASTRUCTURE_DIR/.env not found. Creating from example..."
        if [ -f "$INFRASTRUCTURE_DIR/.env.example" ]; then
            cp "$INFRASTRUCTURE_DIR/.env.example" "$INFRASTRUCTURE_DIR/.env"
            log "Created $INFRASTRUCTURE_DIR/.env from example template"
            log "IMPORTANT: Please review and update the environment variables in $INFRASTRUCTURE_DIR/.env"
        else
            log "ERROR: Example file $INFRASTRUCTURE_DIR/.env.example not found"
            return 1
        fi
    else
        log "Environment file $INFRASTRUCTURE_DIR/.env exists"
    fi
    
    # Check backend .env file
    if [ ! -f "$PROJECT_ROOT/src/backend/.env" ]; then
        log "Backend environment file not found. Creating from example..."
        if [ -f "$PROJECT_ROOT/src/backend/.env.example" ]; then
            cp "$PROJECT_ROOT/src/backend/.env.example" "$PROJECT_ROOT/src/backend/.env"
            log "Created backend environment file from example template"
            log "IMPORTANT: Please review and update the backend environment variables"
        else
            log "WARNING: Backend example file not found, skipping"
        fi
    else
        log "Backend environment file exists"
    fi
    
    # Check frontend .env file
    if [ ! -f "$PROJECT_ROOT/src/web/.env" ]; then
        log "Frontend environment file not found. Creating from example..."
        if [ -f "$PROJECT_ROOT/src/web/.env.example" ]; then
            cp "$PROJECT_ROOT/src/web/.env.example" "$PROJECT_ROOT/src/web/.env"
            log "Created frontend environment file from example template"
            log "IMPORTANT: Please review and update the frontend environment variables"
        else
            log "WARNING: Frontend example file not found, skipping"
        fi
    else
        log "Frontend environment file exists"
    fi
    
    log "Environment files check completed"
    return 0
}

# Create backup before deployment
create_backup() {
    if [ "$SKIP_BACKUP" == "true" ]; then
        log "Skipping backup creation as requested"
        return 0
    fi
    
    log "Creating backup before deployment..."
    
    # Call backup.sh script
    if [ -f "$SCRIPT_DIR/backup.sh" ]; then
        bash "$SCRIPT_DIR/backup.sh" -o "$BACKUP_DIR"
        local result=$?
        
        if [ $result -eq 0 ]; then
            log "Backup created successfully at $BACKUP_DIR"
            return 0
        else
            log "ERROR: Backup creation failed with exit code $result"
            return $result
        fi
    else
        log "ERROR: Backup script not found at $SCRIPT_DIR/backup.sh"
        return 1
    fi
}

# Update source code if in update mode
update_source_code() {
    if [ "$DEPLOYMENT_MODE" != "update" ]; then
        log "Skipping source code update for $DEPLOYMENT_MODE mode"
        return 0
    fi
    
    log "Updating source code..."
    
    # Check if .git directory exists
    if [ ! -d "$PROJECT_ROOT/.git" ]; then
        log "WARNING: Not a git repository, skipping source code update"
        return 0
    fi
    
    # Update code from repository
    cd "$PROJECT_ROOT"
    
    # Stash any local changes
    log "Stashing any local changes..."
    git stash
    
    # Pull latest code
    log "Pulling latest code..."
    if git pull; then
        log "Source code updated successfully"
        return 0
    else
        log "ERROR: Failed to update source code"
        return 1
    fi
}

# Build Docker images
build_docker_images() {
    if [ "$SKIP_BUILD" == "true" ]; then
        log "Skipping Docker image build as requested"
        return 0
    fi
    
    log "Building Docker images..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    if docker-compose build; then
        log "Docker images built successfully"
        return 0
    else
        log "ERROR: Failed to build Docker images"
        return 1
    fi
}

# Stop running services
stop_services() {
    log "Stopping running services..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    if [ "$FORCE_RECREATE" == "true" ]; then
        log "Force recreation requested - removing volumes..."
        if docker-compose down -v; then
            log "Services stopped and volumes removed"
            return 0
        else
            log "ERROR: Failed to stop services and remove volumes"
            return 1
        fi
    else
        if docker-compose down; then
            log "Services stopped"
            return 0
        else
            log "ERROR: Failed to stop services"
            return 1
        fi
    fi
}

# Run database migrations
run_database_migrations() {
    if [ "$SKIP_MIGRATIONS" == "true" ]; then
        log "Skipping database migrations as requested"
        return 0
    fi
    
    log "Running database migrations..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Start only the database service
    log "Starting database service..."
    if ! docker-compose up -d postgres; then
        log "ERROR: Failed to start database service"
        return 1
    fi
    
    # Wait for database to be ready
    log "Waiting for database to be ready..."
    local retry_count=0
    local max_retries=30  # 30 * 2 seconds = 60 seconds timeout
    
    while ! docker-compose exec -T postgres pg_isready >/dev/null 2>&1; do
        retry_count=$((retry_count+1))
        if [ $retry_count -ge $max_retries ]; then
            log "ERROR: Database did not become ready within timeout period"
            return 1
        fi
        log "Database not ready yet, waiting... ($retry_count/$max_retries)"
        sleep 2
    done
    
    log "Database is ready"
    
    # Run migrations through the backend service
    log "Running migrations..."
    if docker-compose run --rm backend alembic upgrade head; then
        log "Database migrations completed successfully"
        return 0
    else
        log "ERROR: Database migrations failed"
        return 1
    fi
}

# Start all services
start_services() {
    log "Starting services..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    if [ "$FORCE_RECREATE" == "true" ]; then
        log "Starting services with force recreation..."
        if docker-compose up -d --force-recreate; then
            log "Services started successfully with force recreation"
            return 0
        else
            log "ERROR: Failed to start services with force recreation"
            return 1
        fi
    else
        if docker-compose up -d; then
            log "Services started successfully"
            return 0
        else
            log "ERROR: Failed to start services"
            return 1
        fi
    fi
}

# Verify deployment
verify_deployment() {
    log "Verifying deployment..."
    
    # Wait for services to be healthy
    log "Waiting for services to become healthy..."
    sleep 10
    
    # Check container status
    cd "$INFRASTRUCTURE_DIR"
    log "Checking container status..."
    
    # Run container status check
    container_status=$(docker-compose ps)
    log "Container status: \n$container_status"
    
    # Check required containers are running (nginx, backend, frontend, postgres, redis, minio)
    required_services=("nginx" "backend" "frontend" "postgres" "redis" "minio")
    failed_services=0
    
    for service in "${required_services[@]}"; do
        if ! echo "$container_status" | grep -q "$service.*Up"; then
            log "ERROR: $service container is not running"
            failed_services=$((failed_services+1))
        else
            log "$service container is running"
        fi
    done
    
    # Use health check script if available
    if [ -f "$SCRIPT_DIR/monitor-health.sh" ]; then
        log "Running health check script..."
        
        # Allow some startup time before health check
        log "Waiting $HEALTH_CHECK_TIMEOUT seconds before health check..."
        sleep 30
        
        bash "$SCRIPT_DIR/monitor-health.sh" -s -l deep
        local health_result=$?
        
        if [ $health_result -ne 0 ]; then
            log "WARNING: Health check failed with exit code $health_result"
            failed_services=$((failed_services+1))
        else
            log "Health check passed"
        fi
    else
        log "WARNING: Health check script not found at $SCRIPT_DIR/monitor-health.sh"
    fi
    
    # Try to access the web interface
    log "Checking web interface accessibility..."
    if curl -s -f -m 10 http://localhost:80 >/dev/null 2>&1; then
        log "Web interface is accessible"
    else
        log "WARNING: Web interface is not accessible"
        failed_services=$((failed_services+1))
    fi
    
    # Try to access the API
    log "Checking API health endpoint..."
    if curl -s -f -m 10 http://localhost:80/api/v1/health/live >/dev/null 2>&1; then
        log "API health endpoint is accessible"
    else
        log "WARNING: API health endpoint is not accessible"
        failed_services=$((failed_services+1))
    fi
    
    if [ $failed_services -eq 0 ]; then
        log "Deployment verification completed successfully"
        return 0
    else
        log "WARNING: Deployment verification completed with $failed_services issues"
        return 1
    fi
}

# Handle deployment failure
handle_deployment_failure() {
    log "ERROR: Deployment failed"
    
    # Ask if user wants to restore from backup
    echo ""
    read -p "Do you want to restore from the backup created before deployment? (y/n): " restore_choice
    
    if [[ "$restore_choice" =~ ^[Yy]$ ]]; then
        log "Attempting to restore from backup..."
        
        # Check if restore script exists
        if [ -f "$SCRIPT_DIR/restore.sh" ]; then
            bash "$SCRIPT_DIR/restore.sh" -i "$BACKUP_DIR"
            local restore_result=$?
            
            if [ $restore_result -eq 0 ]; then
                log "Restoration completed successfully"
                return 0
            else
                log "ERROR: Restoration failed with exit code $restore_result"
                return $restore_result
            fi
        else
            log "ERROR: Restore script not found at $SCRIPT_DIR/restore.sh"
            return 1
        fi
    else
        log "Skipping restoration as requested"
        return 0
    fi
}

# Print completion message
print_completion() {
    echo ""
    echo "====================================================================="
    echo "      Deployment completed successfully!"
    echo ""
    echo "      Access the application at: http://localhost:80"
    echo "      API documentation at: http://localhost:80/docs"
    echo ""
    echo "      Log file: $LOG_FILE"
    echo ""
    echo "      For monitoring use:"
    echo "      - Health check: bash $SCRIPT_DIR/monitor-health.sh"
    echo "      - Container status: cd $INFRASTRUCTURE_DIR && docker-compose ps"
    echo "      - Logs: cd $INFRASTRUCTURE_DIR && docker-compose logs -f"
    echo "====================================================================="
}

# Parse command line arguments
parse_arguments() {
    local OPTIND
    while getopts "m:bfhiusr" opt; do
        case ${opt} in
            m )
                DEPLOYMENT_MODE=$OPTARG
                if [[ ! "$DEPLOYMENT_MODE" =~ ^(initial|update)$ ]]; then
                    echo "ERROR: Invalid deployment mode. Must be 'initial' or 'update'."
                    exit 1
                fi
                ;;
            b )
                SKIP_BACKUP="true"
                ;;
            i )
                DEPLOYMENT_MODE="initial"
                ;;
            u )
                DEPLOYMENT_MODE="update"
                ;;
            s )
                SKIP_BUILD="true"
                ;;
            r )
                FORCE_RECREATE="true"
                ;;
            f )
                SKIP_MIGRATIONS="true"
                ;;
            h )
                echo "Usage: $(basename "$0") [OPTIONS]"
                echo ""
                echo "Deployment script for the Molecular Data Management and CRO Integration Platform"
                echo ""
                echo "Options:"
                echo "  -m MODE    Deployment mode: 'initial' or 'update' (default: update)"
                echo "  -i         Shorthand for '-m initial'"
                echo "  -u         Shorthand for '-m update'"
                echo "  -b         Skip backup creation"
                echo "  -s         Skip Docker image building"
                echo "  -f         Skip database migrations"
                echo "  -r         Force recreation of containers and volumes"
                echo "  -h         Display this help message and exit"
                echo ""
                echo "Examples:"
                echo "  $(basename "$0") -i              # Initial deployment"
                echo "  $(basename "$0") -u -b -s        # Update deployment, skip backup and build"
                echo "  $(basename "$0") -r              # Update with force recreation"
                exit 0
                ;;
            \? )
                echo "Invalid option: -$OPTARG" >&2
                exit 1
                ;;
        esac
    done
}

# Main function to orchestrate the deployment process
main() {
    # Parse command line arguments
    parse_arguments "$@"
    
    # Print banner
    print_banner
    
    # Check dependencies
    check_dependencies || exit 1
    
    # Check environment files
    check_environment_files || exit 1
    
    # Create backup
    create_backup || log "WARNING: Backup creation failed, continuing deployment..."
    
    # Update source code if in update mode
    update_source_code || log "WARNING: Source code update failed, continuing with current code..."
    
    # Stop running services
    stop_services || exit 1
    
    # Build Docker images
    build_docker_images || exit 1
    
    # Run database migrations
    run_database_migrations || exit 1
    
    # Start services
    start_services || exit 1
    
    # Verify deployment
    if ! verify_deployment; then
        handle_deployment_failure || exit 1
    fi
    
    # Print completion message
    print_completion
    
    log "Deployment completed successfully"
    return 0
}

# Execute main function with all arguments
main "$@"
exit $?