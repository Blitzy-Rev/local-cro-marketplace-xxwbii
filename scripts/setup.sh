#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Define script and project directories
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
INFRASTRUCTURE_DIR=$PROJECT_ROOT/infrastructure
BACKEND_DIR=$PROJECT_ROOT/src/backend
WEB_DIR=$PROJECT_ROOT/src/web
LOG_FILE=$SCRIPT_DIR/setup_$(date +%Y%m%d_%H%M%S).log
KEYS_DIR=$BACKEND_DIR/keys

# Required versions
DOCKER_REQUIRED_VERSION="23.0.0"
DOCKER_COMPOSE_REQUIRED_VERSION="2.17.0"

# Default values
DEFAULT_ADMIN_EMAIL="admin@example.com"
DEFAULT_ADMIN_PASSWORD=$(openssl rand -base64 12)
SECRET_KEY=$(openssl rand -base64 32)

# Logging function
log() {
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    echo -e "${timestamp} - $1"
    echo "${timestamp} - $1" | sed -r "s/\x1B\[([0-9]{1,3}(;[0-9]{1,3})*)?[mGK]//g" >> "$LOG_FILE"
}

# Print banner function
print_banner() {
    echo '╔═══════════════════════════════════════════════════════════════════════╗'
    echo '║                                                                       ║'
    echo '║   Molecular Data Management and CRO Integration Platform Setup        ║'
    echo '║                                                                       ║'
    echo '╚═══════════════════════════════════════════════════════════════════════╝'
    echo ""
    echo "Setup starting at: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Log file: $LOG_FILE"
    echo "----------------------------------------"
}

# Check if required dependencies are installed
check_dependencies() {
    log "Checking required dependencies..."
    
    # Check Docker
    if ! command -v docker &> /dev/null; then
        log "Error: Docker is not installed. Please install Docker (v$DOCKER_REQUIRED_VERSION+) first."
        return 1
    fi
    
    # Check if Docker daemon is running
    if ! docker info &> /dev/null; then
        log "Error: Docker daemon is not running. Please start the Docker service first."
        return 1
    fi
    
    # Check Docker version
    DOCKER_VERSION=$(docker --version | sed 's/Docker version \([0-9.]*\).*/\1/')
    if ! awk -v v1="$DOCKER_VERSION" -v v2="$DOCKER_REQUIRED_VERSION" 'BEGIN{if (v1 >= v2) exit 0; else exit 1}'; then
        log "Error: Docker version $DOCKER_VERSION is less than required version $DOCKER_REQUIRED_VERSION"
        return 1
    fi
    log "Docker version $DOCKER_VERSION is installed (✓)"
    
    # Check Docker Compose
    if ! command -v docker-compose &> /dev/null; then
        log "Error: Docker Compose is not installed. Please install Docker Compose (v$DOCKER_COMPOSE_REQUIRED_VERSION+) first."
        return 1
    fi
    
    # Check Docker Compose version
    DOCKER_COMPOSE_VERSION=$(docker-compose --version | sed 's/.*version \([0-9.]*\).*/\1/')
    if ! awk -v v1="$DOCKER_COMPOSE_VERSION" -v v2="$DOCKER_COMPOSE_REQUIRED_VERSION" 'BEGIN{if (v1 >= v2) exit 0; else exit 1}'; then
        log "Error: Docker Compose version $DOCKER_COMPOSE_VERSION is less than required version $DOCKER_COMPOSE_REQUIRED_VERSION"
        return 1
    fi
    log "Docker Compose version $DOCKER_COMPOSE_VERSION is installed (✓)"
    
    # Check OpenSSL
    if ! command -v openssl &> /dev/null; then
        log "Error: OpenSSL is not installed. Please install OpenSSL first."
        return 1
    fi
    log "OpenSSL is installed (✓)"
    
    # Check curl
    if ! command -v curl &> /dev/null; then
        log "Warning: curl is not installed. It may be needed for some operations."
    else
        log "curl is installed (✓)"
    fi
    
    log "All required dependencies are satisfied."
    return 0
}

# Create necessary directories
create_directories() {
    log "Creating necessary directories..."
    
    # Create logs directory
    mkdir -p $PROJECT_ROOT/logs
    log "Created logs directory (✓)"
    
    # Create backups directory
    mkdir -p $PROJECT_ROOT/backups
    log "Created backups directory (✓)"
    
    # Create keys directory for JWT tokens
    mkdir -p $KEYS_DIR
    log "Created keys directory (✓)"
    
    # Create data directories
    mkdir -p $INFRASTRUCTURE_DIR/data/postgres
    mkdir -p $INFRASTRUCTURE_DIR/data/redis
    mkdir -p $INFRASTRUCTURE_DIR/data/minio
    log "Created data directories (✓)"
    
    # Set appropriate permissions
    chmod 700 $KEYS_DIR
    log "Set appropriate permissions on directories (✓)"
    
    return 0
}

# Generate JWT keys
generate_jwt_keys() {
    log "Generating JWT keys..."
    
    # Check if keys directory exists
    if [ ! -d "$KEYS_DIR" ]; then
        mkdir -p "$KEYS_DIR"
        chmod 700 "$KEYS_DIR"
    fi
    
    # Generate private key
    openssl genrsa -out "$KEYS_DIR/jwt-private.pem" 2048
    chmod 600 "$KEYS_DIR/jwt-private.pem"
    
    # Generate public key
    openssl rsa -in "$KEYS_DIR/jwt-private.pem" -pubout -out "$KEYS_DIR/jwt-public.pem"
    chmod 644 "$KEYS_DIR/jwt-public.pem"
    
    log "JWT keys generated successfully (✓)"
    return 0
}

# Configure environment variables
configure_environment() {
    log "Configuring environment files..."
    
    # Infrastructure .env file
    if [ ! -f "$INFRASTRUCTURE_DIR/.env" ]; then
        cp "$INFRASTRUCTURE_DIR/.env.example" "$INFRASTRUCTURE_DIR/.env"
        
        # Update with secure random values
        sed -i.bak "s/POSTGRES_PASSWORD=changeme/POSTGRES_PASSWORD=$(openssl rand -base64 12)/g" "$INFRASTRUCTURE_DIR/.env"
        sed -i.bak "s/MINIO_ROOT_PASSWORD=changeme/MINIO_ROOT_PASSWORD=$(openssl rand -base64 12)/g" "$INFRASTRUCTURE_DIR/.env"
        sed -i.bak "s/SECRET_KEY=generate_a_secure_random_key_here/SECRET_KEY=$SECRET_KEY/g" "$INFRASTRUCTURE_DIR/.env"
        sed -i.bak "s/GRAFANA_ADMIN_PASSWORD=changeme/GRAFANA_ADMIN_PASSWORD=$(openssl rand -base64 12)/g" "$INFRASTRUCTURE_DIR/.env"
        
        log "Infrastructure .env file created and configured (✓)"
    else
        log "Infrastructure .env file already exists, skipping..."
    fi
    
    # Export variables from infrastructure .env for use in other configs
    export $(grep -v '^#' "$INFRASTRUCTURE_DIR/.env" | xargs)
    
    # Backend .env file
    if [ ! -f "$BACKEND_DIR/.env" ]; then
        cp "$BACKEND_DIR/.env.example" "$BACKEND_DIR/.env"
        
        # Update with values from infrastructure .env
        sed -i.bak "s|DATABASE_URL=postgresql://postgres:postgres@localhost:5432/molecular_platform|DATABASE_URL=postgresql://$POSTGRES_USER:$POSTGRES_PASSWORD@postgres:5432/$POSTGRES_DB|g" "$BACKEND_DIR/.env"
        sed -i.bak "s/REDIS_HOST=localhost/REDIS_HOST=redis/g" "$BACKEND_DIR/.env"
        sed -i.bak "s/MINIO_HOST=localhost/MINIO_HOST=minio/g" "$BACKEND_DIR/.env"
        sed -i.bak "s/MINIO_ACCESS_KEY=minioadmin/MINIO_ACCESS_KEY=$MINIO_ROOT_USER/g" "$BACKEND_DIR/.env"
        sed -i.bak "s/MINIO_SECRET_KEY=minioadmin/MINIO_SECRET_KEY=$MINIO_ROOT_PASSWORD/g" "$BACKEND_DIR/.env"
        
        log "Backend .env file created and configured (✓)"
    else
        log "Backend .env file already exists, skipping..."
    fi
    
    # Frontend .env file
    if [ ! -f "$WEB_DIR/.env" ]; then
        cp "$WEB_DIR/.env.production" "$WEB_DIR/.env"
        log "Frontend .env file created and configured (✓)"
    else
        log "Frontend .env file already exists, skipping..."
    fi
    
    # Remove backup files
    rm -f "$INFRASTRUCTURE_DIR/.env.bak"
    rm -f "$BACKEND_DIR/.env.bak"
    
    log "Environment configuration completed (✓)"
    return 0
}

# Build Docker images
build_docker_images() {
    log "Building Docker images..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Build images
    if docker-compose build; then
        log "Docker images built successfully (✓)"
        return 0
    else
        log "Error: Failed to build Docker images."
        return 1
    fi
}

# Initialize database with schema and initial data
initialize_database() {
    log "Initializing database..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Start only the database service
    docker-compose up -d postgres
    
    # Wait for database to be ready
    log "Waiting for database to be ready..."
    sleep 10
    
    # Check if database is ready
    retry_count=0
    max_retries=30
    
    until docker-compose exec -T postgres pg_isready -U $POSTGRES_USER -d $POSTGRES_DB > /dev/null 2>&1; do
        retry_count=$((retry_count+1))
        if [ $retry_count -ge $max_retries ]; then
            log "Error: Timed out waiting for database to be ready."
            return 1
        fi
        log "Database not ready yet, waiting... (Attempt $retry_count/$max_retries)"
        sleep 5
    done
    
    log "Database is ready (✓)"
    
    # Apply database migrations
    log "Running database migrations..."
    if docker-compose exec -T backend alembic upgrade head; then
        log "Database migrations completed successfully (✓)"
    else
        log "Error: Failed to run database migrations."
        return 1
    fi
    
    # Create initial admin user
    log "Creating initial admin user..."
    docker-compose exec -T backend python -c "
from app.db.session import SessionLocal
from app.models.user import User
from app.core.security import get_password_hash
from app.schemas.user_role import UserRole, UserStatus

db = SessionLocal()
admin = db.query(User).filter(User.email == '$DEFAULT_ADMIN_EMAIL').first()
if not admin:
    admin = User(
        email='$DEFAULT_ADMIN_EMAIL',
        password_hash=get_password_hash('$DEFAULT_ADMIN_PASSWORD'),
        role=UserRole.ADMIN,
        status=UserStatus.ACTIVE,
        is_active=True,
        email_verified=True
    )
    db.add(admin)
    db.commit()
    print('Admin user created successfully')
else:
    print('Admin user already exists')
db.close()
"
    
    log "Initial admin user created (✓)"
    log "Database initialization completed (✓)"
    return 0
}

# Initialize MinIO object storage with required buckets
initialize_minio() {
    log "Initializing MinIO..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Start only the MinIO service
    docker-compose up -d minio
    
    # Wait for MinIO to be ready
    log "Waiting for MinIO to be ready..."
    sleep 10
    
    # Check if MinIO is ready
    retry_count=0
    max_retries=30
    
    until docker-compose exec -T minio curl -s http://localhost:9000/minio/health/live > /dev/null; do
        retry_count=$((retry_count+1))
        if [ $retry_count -ge $max_retries ]; then
            log "Error: Timed out waiting for MinIO to be ready."
            return 1
        fi
        log "MinIO not ready yet, waiting... (Attempt $retry_count/$max_retries)"
        sleep 5
    done
    
    log "MinIO is ready (✓)"
    
    # Create required buckets
    log "Creating required buckets..."
    docker-compose exec -T minio mkdir -p /data/csv-uploads
    docker-compose exec -T minio mkdir -p /data/molecule-images
    docker-compose exec -T minio mkdir -p /data/experiment-files
    docker-compose exec -T minio mkdir -p /data/result-files
    docker-compose exec -T minio mkdir -p /data/temp-files
    docker-compose exec -T minio mkdir -p /data/system-backups
    
    log "MinIO buckets created (✓)"
    log "MinIO initialization completed (✓)"
    return 0
}

# Start all services
start_services() {
    log "Starting all services..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Start all services
    if docker-compose up -d; then
        log "All services started successfully (✓)"
        return 0
    else
        log "Error: Failed to start services."
        return 1
    fi
}

# Verify that the setup was successful
verify_setup() {
    log "Verifying setup..."
    
    cd "$INFRASTRUCTURE_DIR"
    
    # Wait for services to initialize
    log "Waiting for services to initialize..."
    sleep 20
    
    # Check if all containers are running
    if docker-compose ps | grep -q "Exit"; then
        log "Error: Some containers are not running. Please check docker-compose logs."
        return 1
    fi
    
    # Check if web interface is accessible
    log "Checking web interface..."
    if curl -s http://localhost > /dev/null; then
        log "Web interface is accessible (✓)"
    else
        log "Warning: Web interface is not accessible. It may still be initializing."
    fi
    
    # Check if API is accessible
    log "Checking API health..."
    if curl -s http://localhost/api/v1/health/live | grep -q "UP"; then
        log "API is accessible and healthy (✓)"
    else
        log "Warning: API health check failed. It may still be initializing."
    fi
    
    log "Setup verification completed. The system may need a few minutes to fully initialize."
    return 0
}

# Setup example data for demonstration
setup_example_data() {
    if [ "$SETUP_EXAMPLE_DATA" = true ]; then
        log "Setting up example data..."
        
        cd "$INFRASTRUCTURE_DIR"
        
        # Run example data script if it exists
        if [ -f "$SCRIPT_DIR/example-data/import-examples.sh" ]; then
            bash "$SCRIPT_DIR/example-data/import-examples.sh"
            log "Example data setup completed (✓)"
        else
            log "Warning: Example data script not found. Skipping example data setup."
        fi
    else
        log "Skipping example data setup (not requested)"
    fi
    
    return 0
}

# Print completion message
print_completion() {
    echo "----------------------------------------"
    echo "Setup Completed Successfully!"
    echo "----------------------------------------"
    echo "The Molecular Data Management and CRO Integration Platform has been set up successfully."
    echo ""
    echo "Access the platform at: http://localhost"
    echo "API Documentation: http://localhost/docs"
    echo ""
    echo "Admin credentials:"
    echo "  Email: $DEFAULT_ADMIN_EMAIL"
    echo "  Password: $DEFAULT_ADMIN_PASSWORD"
    echo ""
    echo "Please save these credentials in a secure location. You can change the password after logging in."
    echo ""
    echo "Setup log file: $LOG_FILE"
    echo ""
    echo "Note: It may take a few minutes for all services to fully initialize."
    echo "----------------------------------------"
}

# Parse command line arguments
parse_arguments() {
    while getopts "he:p:sd" opt; do
        case $opt in
            h)
                echo "Usage: $0 [options]"
                echo ""
                echo "Options:"
                echo "  -h                  Show this help message"
                echo "  -e EMAIL            Set admin email (default: admin@example.com)"
                echo "  -p PASSWORD         Set admin password (default: random)"
                echo "  -s                  Setup example data"
                echo "  -d                  Skip building Docker images (use existing)"
                echo ""
                exit 0
                ;;
            e)
                DEFAULT_ADMIN_EMAIL=$OPTARG
                ;;
            p)
                DEFAULT_ADMIN_PASSWORD=$OPTARG
                ;;
            s)
                SETUP_EXAMPLE_DATA=true
                ;;
            d)
                SKIP_DOCKER_BUILD=true
                ;;
            \?)
                echo "Invalid option: -$OPTARG" >&2
                exit 1
                ;;
        esac
    done
}

# Main function
main() {
    # Parse command-line arguments
    parse_arguments "$@"
    
    # Print banner
    print_banner
    
    # Check dependencies
    if ! check_dependencies; then
        log "Setup failed: Dependency check failed."
        return 1
    fi
    
    # Create directories
    if ! create_directories; then
        log "Setup failed: Directory creation failed."
        return 1
    fi
    
    # Generate JWT keys
    if ! generate_jwt_keys; then
        log "Setup failed: JWT key generation failed."
        return 1
    fi
    
    # Configure environment
    if ! configure_environment; then
        log "Setup failed: Environment configuration failed."
        return 1
    fi
    
    # Build Docker images
    if [ "$SKIP_DOCKER_BUILD" != true ]; then
        if ! build_docker_images; then
            log "Setup failed: Docker image build failed."
            return 1
        fi
    else
        log "Skipping Docker image build (as requested)"
    fi
    
    # Initialize database
    if ! initialize_database; then
        log "Setup failed: Database initialization failed."
        return 1
    fi
    
    # Initialize MinIO
    if ! initialize_minio; then
        log "Setup failed: MinIO initialization failed."
        return 1
    fi
    
    # Start services
    if ! start_services; then
        log "Setup failed: Service startup failed."
        return 1
    fi
    
    # Verify setup
    verify_setup
    
    # Setup example data if requested
    setup_example_data
    
    # Print completion message
    print_completion
    
    return 0
}

# Run main function with all arguments
main "$@"
exit $?