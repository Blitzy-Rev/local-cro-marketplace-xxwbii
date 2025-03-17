#!/bin/bash
# 
# restore.sh - System restoration script for the Molecular Data Management and CRO Integration Platform
#
# This script automates the restoration process from backups, restoring the database,
# file storage, and configuration to enable disaster recovery and system restoration.
#
# Usage: ./restore.sh [OPTIONS]
#

# Exit immediately if a command exits with a non-zero status
set -e

# Script directory and project paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
INFRASTRUCTURE_DIR=$PROJECT_ROOT/infrastructure
LOG_FILE=$SCRIPT_DIR/restore_$(date +%Y%m%d_%H%M%S).log

# Default configuration
BACKUP_DIR=${BACKUP_DIR:-$PROJECT_ROOT/backups}
INCLUDE_DB=true
INCLUDE_FILES=true
INCLUDE_CONFIG=true
FORCE_STOP=false

# Log messages to both console and log file
log() {
    local log_dir=$(dirname "$LOG_FILE")
    
    # Create log directory if it doesn't exist
    if [ ! -d "$log_dir" ]; then
        mkdir -p "$log_dir"
    fi
    
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    local message="[$timestamp] $1"
    
    echo "$message"
    echo "$message" >> "$LOG_FILE"
}

# Display usage information
print_usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo "Restore the Molecular Data Management and CRO Integration Platform from backup."
    echo ""
    echo "Options:"
    echo "  -b DIR    Specify backup directory (default: $BACKUP_DIR)"
    echo "  -s TIME   Specify backup timestamp (format: YYYYMMDD_HHMMSS)"
    echo "  -d        Include/exclude database restoration (default: $INCLUDE_DB)"
    echo "  -f        Include/exclude file storage restoration (default: $INCLUDE_FILES)"
    echo "  -c        Include/exclude configuration restoration (default: $INCLUDE_CONFIG)"
    echo "  -x        Force stop all services before restoration (default: $FORCE_STOP)"
    echo "  -h        Display this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") -b /path/to/backups -s 20230701_120000"
    echo "  $(basename "$0") -d false -f false -c true"
    echo ""
    echo "Note: If no timestamp is specified, the latest backup will be used."
}

# Parse command line arguments
parse_arguments() {
    local OPTIND opt
    local BACKUP_TIMESTAMP=""
    
    while getopts ":b:s:d:f:c:xh" opt; do
        case $opt in
            b)
                BACKUP_DIR="$OPTARG"
                ;;
            s)
                BACKUP_TIMESTAMP="$OPTARG"
                ;;
            d)
                if [[ "$OPTARG" == "false" ]]; then
                    INCLUDE_DB=false
                else
                    INCLUDE_DB=true
                fi
                ;;
            f)
                if [[ "$OPTARG" == "false" ]]; then
                    INCLUDE_FILES=false
                else
                    INCLUDE_FILES=true
                fi
                ;;
            c)
                if [[ "$OPTARG" == "false" ]]; then
                    INCLUDE_CONFIG=false
                else
                    INCLUDE_CONFIG=true
                fi
                ;;
            x)
                FORCE_STOP=true
                ;;
            h)
                print_usage
                return 0
                ;;
            \?)
                log "Error: Invalid option -$OPTARG"
                print_usage
                return 1
                ;;
            :)
                log "Error: Option -$OPTARG requires an argument"
                print_usage
                return 1
                ;;
        esac
    done
    
    # If specific backup timestamp was provided, construct the backup directory path
    if [ -n "$BACKUP_TIMESTAMP" ]; then
        BACKUP_DIR="$BACKUP_DIR/backup_$BACKUP_TIMESTAMP"
    fi
    
    # Validate backup directory
    if [ ! -d "$BACKUP_DIR" ]; then
        log "Error: Backup directory does not exist: $BACKUP_DIR"
        return 1
    fi
    
    return 0
}

# Find the latest backup directory
find_latest_backup() {
    if [ ! -d "$BACKUP_DIR" ]; then
        log "Error: Backup directory does not exist: $BACKUP_DIR"
        return 1
    fi
    
    # Find the most recent backup directory
    local latest_backup=$(find "$BACKUP_DIR" -maxdepth 1 -type d -name "backup_*" | sort -r | head -n 1)
    
    if [ -z "$latest_backup" ]; then
        log "Error: No backup directories found in $BACKUP_DIR"
        return 1
    fi
    
    echo "$latest_backup"
    return 0
}

# Validate backup directory contains required files
validate_backup() {
    local backup_dir="$1"
    local valid=true
    
    log "Validating backup directory: $backup_dir"
    
    if [ ! -d "$backup_dir" ]; then
        log "Error: Backup directory does not exist: $backup_dir"
        return 1
    fi
    
    # Check database backup
    if [ "$INCLUDE_DB" = true ]; then
        if [ ! -f "$backup_dir/database/backup.sql.gz" ] && [ ! -f "$backup_dir/database/backup.sql" ]; then
            log "Warning: Database backup not found in $backup_dir/database"
            valid=false
        else
            log "Database backup found"
        fi
    fi
    
    # Check file storage backup
    if [ "$INCLUDE_FILES" = true ]; then
        if [ ! -d "$backup_dir/file_storage" ]; then
            log "Warning: File storage backup not found in $backup_dir/file_storage"
            valid=false
        else
            log "File storage backup found"
        fi
    fi
    
    # Check configuration backup
    if [ "$INCLUDE_CONFIG" = true ]; then
        if [ ! -d "$backup_dir/config" ]; then
            log "Warning: Configuration backup not found in $backup_dir/config"
            valid=false
        else
            log "Configuration backup found"
        fi
    fi
    
    if [ "$valid" = false ]; then
        log "Backup validation failed: Some required components are missing"
        return 1
    fi
    
    log "Backup validation successful"
    return 0
}

# Stop running services before restoration
stop_services() {
    log "Stopping services..."
    
    # Navigate to the infrastructure directory
    cd "$INFRASTRUCTURE_DIR"
    
    # Stop all services using docker-compose
    if docker-compose down; then
        log "Services stopped successfully"
        return 0
    else
        log "Failed to stop services"
        return 1
    fi
}

# Restore PostgreSQL database from backup
restore_database() {
    local backup_dir="$1"
    
    if [ "$INCLUDE_DB" != true ]; then
        log "Skipping database restoration as requested"
        return 0
    fi
    
    log "Starting database restoration..."
    
    # Navigate to the infrastructure directory
    cd "$INFRASTRUCTURE_DIR"
    
    # Start only the PostgreSQL container if it's not already running
    if ! docker-compose ps postgres | grep -q "Up"; then
        log "Starting PostgreSQL container..."
        docker-compose up -d postgres
        
        # Wait for PostgreSQL to be ready
        log "Waiting for PostgreSQL to be ready..."
        local retry_count=0
        local max_retries=30
        
        until docker-compose exec postgres pg_isready -U postgres -d postgres &> /dev/null; do
            retry_count=$((retry_count+1))
            if [ $retry_count -ge $max_retries ]; then
                log "Error: PostgreSQL did not become ready in time"
                return 1
            fi
            log "PostgreSQL is not ready yet, waiting... ($retry_count/$max_retries)"
            sleep 5
        done
    fi
    
    # Find the database backup file
    local db_backup=""
    if [ -f "$backup_dir/database/backup.sql.gz" ]; then
        db_backup="$backup_dir/database/backup.sql.gz"
    elif [ -f "$backup_dir/database/backup.sql" ]; then
        db_backup="$backup_dir/database/backup.sql"
    else
        log "Error: Could not find database backup file"
        return 1
    fi
    
    log "Using database backup: $db_backup"
    
    # Call the PostgreSQL restore script
    if "$INFRASTRUCTURE_DIR/postgres/restore.sh" -f "$db_backup"; then
        log "Database restored successfully"
        return 0
    else
        log "Error: Failed to restore database"
        return 1
    fi
}

# Restore MinIO object storage from backup
restore_file_storage() {
    local backup_dir="$1"
    
    if [ "$INCLUDE_FILES" != true ]; then
        log "Skipping file storage restoration as requested"
        return 0
    fi
    
    log "Starting file storage restoration..."
    
    # Check if file storage backup exists
    if [ ! -d "$backup_dir/file_storage" ]; then
        log "Error: File storage backup directory not found: $backup_dir/file_storage"
        return 1
    fi
    
    # Navigate to the infrastructure directory
    cd "$INFRASTRUCTURE_DIR"
    
    # Extract MinIO credentials from .env file if available
    local minio_user="minioadmin"
    local minio_pass="minioadmin"
    
    if [ -f "$INFRASTRUCTURE_DIR/.env" ]; then
        # Source the .env file to get credentials
        source "$INFRASTRUCTURE_DIR/.env"
        minio_user="${MINIO_ROOT_USER:-minioadmin}"
        minio_pass="${MINIO_ROOT_PASSWORD:-minioadmin}"
    fi
    
    # Start only the MinIO container if it's not already running
    if ! docker-compose ps minio | grep -q "Up"; then
        log "Starting MinIO container..."
        docker-compose up -d minio
        
        # Wait for MinIO to be ready
        log "Waiting for MinIO to be ready..."
        local retry_count=0
        local max_retries=30
        
        until docker-compose exec minio curl -s http://localhost:9000/minio/health/live &> /dev/null; do
            retry_count=$((retry_count+1))
            if [ $retry_count -ge $max_retries ]; then
                log "Error: MinIO did not become ready in time"
                return 1
            fi
            log "MinIO is not ready yet, waiting... ($retry_count/$max_retries)"
            sleep 5
        done
    fi
    
    # Configure mc client inside the container
    log "Configuring MinIO client..."
    docker-compose exec minio mc alias set myminio http://localhost:9000 "$minio_user" "$minio_pass"
    
    # Check if buckets exist, create them if they don't
    log "Checking and creating required buckets..."
    local buckets=("csv-uploads" "molecule-images" "experiment-files" "result-files" "temp-files" "system-backups")
    
    for bucket in "${buckets[@]}"; do
        if ! docker-compose exec minio mc ls myminio | grep -q "$bucket"; then
            log "Creating bucket: $bucket"
            docker-compose exec minio mc mb myminio/$bucket
        fi
    done
    
    # Restore data from backup to each bucket
    log "Restoring file data to MinIO buckets..."
    for bucket in "${buckets[@]}"; do
        if [ -d "$backup_dir/file_storage/$bucket" ]; then
            log "Restoring data to bucket: $bucket"
            # Copy the backup data to a temporary location in the container
            docker-compose exec -T minio mkdir -p /tmp/restore/$bucket
            docker cp "$backup_dir/file_storage/$bucket/." $(docker-compose ps -q minio):/tmp/restore/$bucket/
            
            # Mirror data from temporary location to the bucket
            docker-compose exec minio mc mirror /tmp/restore/$bucket myminio/$bucket
            
            # Clean up temporary data
            docker-compose exec minio rm -rf /tmp/restore/$bucket
        else
            log "Warning: No backup data found for bucket: $bucket"
        fi
    done
    
    log "File storage restored successfully"
    return 0
}

# Restore configuration files from backup
restore_configuration() {
    local backup_dir="$1"
    
    if [ "$INCLUDE_CONFIG" != true ]; then
        log "Skipping configuration restoration as requested"
        return 0
    fi
    
    log "Starting configuration restoration..."
    
    # Check if configuration backup exists
    if [ ! -d "$backup_dir/config" ]; then
        log "Error: Configuration backup directory not found: $backup_dir/config"
        return 1
    fi
    
    # Restore environment files
    if [ -f "$backup_dir/config/.env" ]; then
        log "Restoring environment configuration..."
        cp "$backup_dir/config/.env" "$INFRASTRUCTURE_DIR/.env"
    fi
    
    # Restore docker-compose.yml if it exists in the backup
    if [ -f "$backup_dir/config/docker-compose.yml" ]; then
        log "Restoring docker-compose.yml..."
        cp "$backup_dir/config/docker-compose.yml" "$INFRASTRUCTURE_DIR/docker-compose.yml"
    fi
    
    # Restore nginx configuration
    if [ -d "$backup_dir/config/nginx" ]; then
        log "Restoring nginx configuration..."
        mkdir -p "$INFRASTRUCTURE_DIR/nginx"
        cp -r "$backup_dir/config/nginx/." "$INFRASTRUCTURE_DIR/nginx/"
    fi
    
    # Restore postgres configuration
    if [ -d "$backup_dir/config/postgres" ]; then
        log "Restoring postgres configuration..."
        mkdir -p "$INFRASTRUCTURE_DIR/postgres"
        cp -r "$backup_dir/config/postgres/." "$INFRASTRUCTURE_DIR/postgres/"
    fi
    
    # Restore redis configuration
    if [ -d "$backup_dir/config/redis" ]; then
        log "Restoring redis configuration..."
        mkdir -p "$INFRASTRUCTURE_DIR/redis"
        cp -r "$backup_dir/config/redis/." "$INFRASTRUCTURE_DIR/redis/"
    fi
    
    # Restore minio configuration
    if [ -d "$backup_dir/config/minio" ]; then
        log "Restoring minio configuration..."
        mkdir -p "$INFRASTRUCTURE_DIR/minio"
        cp -r "$backup_dir/config/minio/." "$INFRASTRUCTURE_DIR/minio/"
    fi
    
    log "Configuration restored successfully"
    return 0
}

# Start all services after restoration
start_services() {
    log "Starting services..."
    
    # Navigate to the infrastructure directory
    cd "$INFRASTRUCTURE_DIR"
    
    # Start all services using docker-compose
    if docker-compose up -d; then
        log "Services started successfully"
        return 0
    else
        log "Failed to start services"
        return 1
    fi
}

# Verify that the restoration was successful
verify_restoration() {
    log "Verifying restoration..."
    
    # Wait for services to be healthy
    log "Waiting for services to be healthy..."
    sleep 30
    
    # Navigate to the infrastructure directory
    cd "$INFRASTRUCTURE_DIR"
    
    # Check if containers are running
    local containers=("nginx" "frontend" "backend" "worker" "postgres" "redis" "minio")
    local all_running=true
    
    for container in "${containers[@]}"; do
        if ! docker-compose ps "$container" | grep -q "Up"; then
            log "Warning: Container $container is not running"
            all_running=false
        else
            log "Container $container is running"
        fi
    done
    
    # Verify database connectivity if database was restored
    if [ "$INCLUDE_DB" = true ]; then
        log "Verifying database connectivity..."
        if docker-compose exec postgres pg_isready -U postgres -d molecular_platform &> /dev/null; then
            log "Database connectivity verified"
            
            # Check if database has expected tables
            local table_count=$(docker-compose exec -T postgres psql -U postgres -d molecular_platform -t -c "SELECT COUNT(*) FROM information_schema.tables WHERE table_schema = 'public';" | tr -d ' ')
            log "Database contains $table_count tables"
            
            if [ "$table_count" -lt 5 ]; then
                log "Warning: Database may not be properly restored (only $table_count tables found)"
                all_running=false
            fi
        else
            log "Warning: Could not verify database connectivity"
            all_running=false
        fi
    fi
    
    # Verify MinIO buckets if file storage was restored
    if [ "$INCLUDE_FILES" = true ]; then
        log "Verifying MinIO buckets..."
        local expected_buckets=("csv-uploads" "molecule-images" "experiment-files" "result-files" "temp-files" "system-backups")
        local buckets_verified=true
        
        for bucket in "${expected_buckets[@]}"; do
            if ! docker-compose exec minio mc ls myminio | grep -q "$bucket"; then
                log "Warning: MinIO bucket '$bucket' not found"
                buckets_verified=false
            else
                log "MinIO bucket '$bucket' verified"
            fi
        done
        
        if [ "$buckets_verified" = false ]; then
            log "Warning: Some MinIO buckets are missing"
            all_running=false
        fi
    fi
    
    # Wait a bit longer for the API to be ready
    log "Waiting for API to be ready..."
    sleep 30
    
    # Extract HTTP port from .env file if available
    local http_port=80
    if [ -f "$INFRASTRUCTURE_DIR/.env" ]; then
        source "$INFRASTRUCTURE_DIR/.env"
        http_port="${EXTERNAL_HTTP_PORT:-80}"
    fi
    
    # Verify API health endpoint
    log "Verifying API health at http://localhost:$http_port/health..."
    if curl -s "http://localhost:$http_port/health" | grep -q "UP"; then
        log "API health verified"
    else
        log "Warning: Could not verify API health"
        all_running=false
    fi
    
    if [ "$all_running" = true ]; then
        log "All services appear to be running correctly"
        return 0
    else
        log "Warning: Some services may not be running correctly"
        return 1
    fi
}

# Create a summary of the restoration operation
create_restoration_summary() {
    local backup_dir="$1"
    local status="$2"
    
    log "Creating restoration summary..."
    
    # Create summary
    local summary="
==============================================
RESTORATION SUMMARY
==============================================
Timestamp: $(date '+%Y-%m-%d %H:%M:%S')
Backup directory: $backup_dir
Status: $([ $status -eq 0 ] && echo "SUCCESS" || echo "WARNING/ERROR")

Components Restored:
- Database: $INCLUDE_DB
- File Storage: $INCLUDE_FILES
- Configuration: $INCLUDE_CONFIG

For detailed information, see log file:
$LOG_FILE
==============================================
"
    
    # Log the summary
    log "$summary"
    
    # Also save to a separate summary file
    echo "$summary" > "$SCRIPT_DIR/restore_summary_$(date +%Y%m%d_%H%M%S).txt"
}

# Main function
main() {
    local start_time=$(date '+%Y-%m-%d %H:%M:%S')
    log "Starting system restoration at $start_time"
    
    # Parse command line arguments
    if ! parse_arguments "$@"; then
        return 1
    fi
    
    # If no specific backup path was provided, use the latest backup
    if [[ "$BACKUP_DIR" == *"/backup_"* ]]; then
        log "Using specified backup directory: $BACKUP_DIR"
    else
        log "No specific backup timestamp provided, finding latest backup..."
        backup_dir=$(find_latest_backup)
        if [ $? -ne 0 ]; then
            return 1
        fi
        BACKUP_DIR="$backup_dir"
        log "Using latest backup directory: $BACKUP_DIR"
    fi
    
    # Validate backup
    if ! validate_backup "$BACKUP_DIR"; then
        log "Error: Backup validation failed for $BACKUP_DIR"
        return 1
    fi
    
    # Stop services if requested or required
    if [ "$FORCE_STOP" = true ] || [ "$INCLUDE_DB" = true ] || [ "$INCLUDE_FILES" = true ]; then
        if ! stop_services; then
            log "Warning: Failed to stop services, proceeding with restoration anyway"
        fi
    fi
    
    # Perform restoration
    local db_status=0
    local files_status=0
    local config_status=0
    
    # Restore database
    if [ "$INCLUDE_DB" = true ]; then
        restore_database "$BACKUP_DIR"
        db_status=$?
    fi
    
    # Restore file storage
    if [ "$INCLUDE_FILES" = true ]; then
        restore_file_storage "$BACKUP_DIR"
        files_status=$?
    fi
    
    # Restore configuration
    if [ "$INCLUDE_CONFIG" = true ]; then
        restore_configuration "$BACKUP_DIR"
        config_status=$?
    fi
    
    # Start services
    if ! start_services; then
        log "Error: Failed to start services after restoration"
        return 1
    fi
    
    # Verify restoration
    local verify_status=0
    verify_restoration
    verify_status=$?
    
    # Calculate overall status
    local overall_status=0
    if [ $db_status -ne 0 ] || [ $files_status -ne 0 ] || [ $config_status -ne 0 ] || [ $verify_status -ne 0 ]; then
        overall_status=1
    fi
    
    # Create restoration summary
    create_restoration_summary "$BACKUP_DIR" $overall_status
    
    local end_time=$(date '+%Y-%m-%d %H:%M:%S')
    log "System restoration completed at $end_time"
    
    return $overall_status
}

# Call the main function with all script arguments
main "$@"