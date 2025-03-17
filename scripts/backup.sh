#!/bin/bash
# Molecular Data Management and CRO Integration Platform
# Comprehensive Backup Script
# Version 1.0.0

# Exit immediately if a command exits with a non-zero status
set -e

# Script directory and paths
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(dirname "$SCRIPT_DIR")
INFRASTRUCTURE_DIR=$PROJECT_ROOT/infrastructure
BACKUP_DIR=$PROJECT_ROOT/backups/$(date +%Y%m%d_%H%M%S)
LOG_FILE=$SCRIPT_DIR/backup_$(date +%Y%m%d_%H%M%S).log

# Default configuration
RETENTION_DAYS=${RETENTION_DAYS:-7}
BACKUP_TYPE=full
INCLUDE_DB=true
INCLUDE_FILES=true
INCLUDE_CONFIG=true

# Log messages to console and file
log() {
    local message=$1
    local timestamp=$(date "+%Y-%m-%d %H:%M:%S")
    
    # Create log directory if it doesn't exist
    mkdir -p "$(dirname "$LOG_FILE")"
    
    echo "[$timestamp] $message"
    echo "[$timestamp] $message" >> "$LOG_FILE"
}

# Display usage information
print_usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo
    echo "Comprehensive backup script for the Molecular Data Management and CRO Integration Platform"
    echo
    echo "Options:"
    echo "  -o DIR    Output directory for backup files (default: auto-generated)"
    echo "  -t TYPE   Backup type: full or incremental (default: full)"
    echo "  -r DAYS   Retention period in days (default: ${RETENTION_DAYS})"
    echo "  -d        Disable database backup"
    echo "  -f        Disable file storage backup"
    echo "  -c        Disable configuration backup"
    echo "  -h        Display this help message and exit"
    echo
    echo "Examples:"
    echo "  $(basename "$0") -o /data/backups -r 14"
    echo "  $(basename "$0") -t incremental"
}

# Parse command-line arguments
parse_arguments() {
    local OPTIND
    while getopts "o:t:r:dfch" opt; do
        case $opt in
            o)
                BACKUP_DIR="$OPTARG"
                ;;
            t)
                BACKUP_TYPE="$OPTARG"
                if [[ ! "$BACKUP_TYPE" =~ ^(full|incremental)$ ]]; then
                    log "Error: Backup type must be 'full' or 'incremental'"
                    return 1
                fi
                ;;
            r)
                RETENTION_DAYS="$OPTARG"
                if ! [[ "$RETENTION_DAYS" =~ ^[0-9]+$ ]]; then
                    log "Error: Retention days must be a positive integer"
                    return 1
                fi
                ;;
            d)
                INCLUDE_DB=false
                ;;
            f)
                INCLUDE_FILES=false
                ;;
            c)
                INCLUDE_CONFIG=false
                ;;
            h)
                print_usage
                exit 0
                ;;
            \?)
                log "Error: Invalid option -$OPTARG"
                print_usage
                return 1
                ;;
        esac
    done

    # Validate that at least one backup type is enabled
    if [[ "$INCLUDE_DB" == "false" && "$INCLUDE_FILES" == "false" && "$INCLUDE_CONFIG" == "false" ]]; then
        log "Error: At least one backup type must be enabled"
        return 1
    fi

    return 0
}

# Create backup directory structure
create_backup_directory() {
    log "Creating backup directory: $BACKUP_DIR"
    
    mkdir -p "$BACKUP_DIR"
    if [ $? -ne 0 ]; then
        log "Error: Failed to create backup directory '$BACKUP_DIR'"
        return 1
    fi
    
    # Create subdirectories
    if [[ "$INCLUDE_DB" == "true" ]]; then
        mkdir -p "$BACKUP_DIR/database"
    fi
    
    if [[ "$INCLUDE_FILES" == "true" ]]; then
        mkdir -p "$BACKUP_DIR/files"
    fi
    
    if [[ "$INCLUDE_CONFIG" == "true" ]]; then
        mkdir -p "$BACKUP_DIR/config"
    fi
    
    # Set appropriate permissions
    chmod -R 750 "$BACKUP_DIR"
    
    log "Backup directory created successfully"
    return 0
}

# Backup PostgreSQL database
backup_database() {
    if [[ "$INCLUDE_DB" != "true" ]]; then
        log "Database backup is disabled, skipping..."
        return 0
    fi
    
    log "Starting database backup..."
    
    # Ensure PostgreSQL container is running
    if ! docker ps | grep -q "postgres"; then
        log "Error: PostgreSQL container is not running"
        return 1
    fi
    
    # Use the existing postgres backup script if available
    if [ -f "$INFRASTRUCTURE_DIR/postgres/backup.sh" ]; then
        log "Using PostgreSQL backup script from infrastructure directory"
        
        # Set environment variables for the backup script
        export BACKUP_DIR="$BACKUP_DIR/database"
        export TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
        
        # Execute the PostgreSQL backup script
        bash "$INFRASTRUCTURE_DIR/postgres/backup.sh" -o "$BACKUP_DIR/database"
        local result=$?
        
        if [ $result -ne 0 ]; then
            log "Error: PostgreSQL backup script failed with exit code $result"
            return $result
        fi
        
        log "Database backup completed successfully using infrastructure script"
        return 0
    fi
    
    # If postgres backup script is not available, perform backup directly
    log "PostgreSQL backup script not found, performing direct backup..."
    
    # Get PostgreSQL container name
    local postgres_container=$(docker ps | grep "postgres" | awk '{print $NF}' | head -1)
    
    # Get database parameters
    local db_user=$(docker exec "$postgres_container" printenv POSTGRES_USER || echo "postgres")
    local db_name=$(docker exec "$postgres_container" printenv POSTGRES_DB || echo "molecular_platform")
    local backup_file="$BACKUP_DIR/database/db_backup_$(date +%Y%m%d_%H%M%S).sql.gz"
    
    # Ensure database backup directory exists
    mkdir -p "$BACKUP_DIR/database"
    
    # Execute backup using pg_dump inside the container
    log "Executing pg_dump inside PostgreSQL container..."
    docker exec "$postgres_container" pg_dump \
        -U "$db_user" \
        -d "$db_name" \
        --format=custom \
        --compress=9 \
        | gzip > "$backup_file"
    
    local result=$?
    
    if [ $result -ne 0 ]; then
        log "Error: Database backup failed"
        return $result
    fi
    
    # Verify backup file
    if [ ! -s "$backup_file" ]; then
        log "Error: Database backup file is empty"
        return 1
    fi
    
    local file_size=$(du -h "$backup_file" | cut -f1)
    log "Database backup completed successfully (Size: $file_size)"
    return 0
}

# Backup MinIO file storage
backup_file_storage() {
    if [[ "$INCLUDE_FILES" != "true" ]]; then
        log "File storage backup is disabled, skipping..."
        return 0
    fi
    
    log "Starting file storage backup..."
    
    # Ensure MinIO container is running
    if ! docker ps | grep -q "minio"; then
        log "Error: MinIO container is not running"
        return 1
    fi
    
    # Get MinIO container name
    local minio_container=$(docker ps | grep "minio" | awk '{print $NF}' | head -1)
    
    # Create files backup directory
    mkdir -p "$BACKUP_DIR/files"
    
    # If mc client is not configured in the container, we'll need to configure it
    log "Checking MinIO client configuration..."
    docker exec "$minio_container" sh -c "mc config host ls | grep -q local" || {
        log "Configuring MinIO client..."
        
        # Get MinIO credentials from environment variables
        local minio_user=$(docker exec "$minio_container" printenv MINIO_ROOT_USER || echo "minioadmin")
        local minio_pass=$(docker exec "$minio_container" printenv MINIO_ROOT_PASSWORD || echo "minioadmin")
        
        # Configure MinIO client
        docker exec "$minio_container" sh -c "mc config host add local http://localhost:9000 ${minio_user} ${minio_pass}"
        
        if [ $? -ne 0 ]; then
            log "Error: Failed to configure MinIO client"
            return 1
        fi
    }
    
    if [[ "$BACKUP_TYPE" == "full" ]]; then
        log "Performing full backup of MinIO storage..."
        
        # Get list of buckets
        local buckets=$(docker exec "$minio_container" sh -c "mc ls local/ | awk '{print \$5}'")
        
        if [ -z "$buckets" ]; then
            log "Warning: No buckets found in MinIO"
            return 0
        fi
        
        # Backup each bucket
        for bucket in $buckets; do
            log "Backing up bucket: $bucket"
            
            # Create temporary directory in container
            docker exec "$minio_container" sh -c "mkdir -p /tmp/backup"
            
            # Use MinIO client to backup bucket data
            docker exec "$minio_container" sh -c "mc cp --recursive local/$bucket/ /tmp/backup/"
            
            if [ $? -ne 0 ]; then
                log "Error: Failed to copy bucket $bucket to temporary location"
                return 1
            fi
            
            # Create bucket directory in backup location
            mkdir -p "$BACKUP_DIR/files/$bucket"
            
            # Copy files from container to host
            docker cp "$minio_container:/tmp/backup/." "$BACKUP_DIR/files/$bucket/"
            
            if [ $? -ne 0 ]; then
                log "Error: Failed to copy bucket data from container to backup directory"
                return 1
            fi
            
            # Clean up temporary directory in container
            docker exec "$minio_container" sh -c "rm -rf /tmp/backup"
        done
    else
        # Incremental backup (simplified implementation)
        log "Performing incremental backup of MinIO storage..."
        
        # Find the latest backup directory
        local backup_base_dir=$(dirname "$BACKUP_DIR")
        local last_backup=$(find "$backup_base_dir" -maxdepth 1 -type d -name "[0-9]*_[0-9]*" -mtime -$RETENTION_DAYS | sort | tail -1)
        
        if [ -z "$last_backup" ] || [ ! -d "$last_backup/files" ]; then
            log "No previous backup found, performing full backup instead"
            BACKUP_TYPE=full
            backup_file_storage
            return $?
        fi
        
        # Get list of buckets
        local buckets=$(docker exec "$minio_container" sh -c "mc ls local/ | awk '{print \$5}'")
        
        # Backup each bucket incrementally
        for bucket in $buckets; do
            log "Backing up bucket incrementally: $bucket"
            
            # Create temporary directory in container
            docker exec "$minio_container" sh -c "mkdir -p /tmp/backup"
            
            # Use MinIO client to get the last modified time
            local last_backup_time=$(stat -c %Y "$last_backup/files/$bucket" 2>/dev/null || echo "0")
            local current_time=$(date +%s)
            
            # Use find to identify files modified since last backup
            docker exec "$minio_container" sh -c "mc find local/$bucket --newer-than ${last_backup_time}s --exec 'mc cp \"{}\" /tmp/backup/'"
            
            # Create bucket directory in backup location
            mkdir -p "$BACKUP_DIR/files/$bucket"
            
            # Copy only changed files from container to host
            docker cp "$minio_container:/tmp/backup/." "$BACKUP_DIR/files/$bucket/"
            
            # Copy files that exist in the last backup but not in the current one
            if [ -d "$last_backup/files/$bucket" ]; then
                rsync -a --ignore-existing "$last_backup/files/$bucket/" "$BACKUP_DIR/files/$bucket/"
            fi
            
            # Clean up temporary directory in container
            docker exec "$minio_container" sh -c "rm -rf /tmp/backup"
        done
    fi
    
    log "File storage backup completed successfully"
    return 0
}

# Backup configuration files
backup_configuration() {
    if [[ "$INCLUDE_CONFIG" != "true" ]]; then
        log "Configuration backup is disabled, skipping..."
        return 0
    fi
    
    log "Starting configuration backup..."
    
    # Backup environment files
    local env_files=(.env .env.* docker-compose.env)
    for env_file in "${env_files[@]}"; do
        if [ -f "$PROJECT_ROOT/$env_file" ]; then
            cp "$PROJECT_ROOT/$env_file" "$BACKUP_DIR/config/"
            log "Backed up $env_file"
        fi
    done
    
    # Backup docker-compose files
    local compose_files=(docker-compose.yml docker-compose.*.yml)
    for compose_file in "${compose_files[@]}"; do
        if [ -f "$PROJECT_ROOT/$compose_file" ]; then
            cp "$PROJECT_ROOT/$compose_file" "$BACKUP_DIR/config/"
            log "Backed up $compose_file"
        fi
    done
    
    # Backup nginx configuration
    if [ -d "$PROJECT_ROOT/nginx" ]; then
        mkdir -p "$BACKUP_DIR/config/nginx"
        cp -r "$PROJECT_ROOT/nginx" "$BACKUP_DIR/config/"
        log "Backed up nginx configuration"
    fi
    
    # Backup other important configuration directories
    local config_dirs=(config conf configuration settings)
    for config_dir in "${config_dirs[@]}"; do
        if [ -d "$PROJECT_ROOT/$config_dir" ]; then
            mkdir -p "$BACKUP_DIR/config/$config_dir"
            cp -r "$PROJECT_ROOT/$config_dir" "$BACKUP_DIR/config/"
            log "Backed up $config_dir directory"
        fi
    done
    
    log "Configuration backup completed successfully"
    return 0
}

# Clean up old backups
cleanup_old_backups() {
    log "Cleaning up backups older than $RETENTION_DAYS days..."
    
    local backup_base_dir=$(dirname "$BACKUP_DIR")
    
    if [ ! -d "$backup_base_dir" ]; then
        log "Backup base directory does not exist, skipping cleanup"
        return 0
    fi
    
    # Find backup directories older than RETENTION_DAYS
    local old_backups=$(find "$backup_base_dir" -maxdepth 1 -type d -name "[0-9]*_[0-9]*" -mtime +$RETENTION_DAYS)
    
    if [ -z "$old_backups" ]; then
        log "No old backups to clean up"
        return 0
    fi
    
    # Remove old backup directories
    local count=0
    for old_backup in $old_backups; do
        log "Removing old backup: $old_backup"
        rm -rf "$old_backup"
        if [ $? -eq 0 ]; then
            count=$((count+1))
        else
            log "Warning: Failed to remove backup directory: $old_backup"
        fi
    done
    
    log "Cleaned up $count old backup(s)"
    return 0
}

# Create backup summary
create_backup_summary() {
    local summary_file="$BACKUP_DIR/backup_summary.txt"
    
    log "Creating backup summary..."
    
    echo "Backup Summary" > "$summary_file"
    echo "==============" >> "$summary_file"
    echo "" >> "$summary_file"
    echo "Timestamp: $(date)" >> "$summary_file"
    echo "Backup Type: $BACKUP_TYPE" >> "$summary_file"
    echo "Backup Location: $BACKUP_DIR" >> "$summary_file"
    echo "" >> "$summary_file"
    
    echo "Components:" >> "$summary_file"
    
    if [[ "$INCLUDE_DB" == "true" ]]; then
        local db_size=$(du -sh "$BACKUP_DIR/database" 2>/dev/null | cut -f1 || echo "N/A")
        echo "- Database: Backed up (Size: $db_size)" >> "$summary_file"
    else
        echo "- Database: Skipped" >> "$summary_file"
    fi
    
    if [[ "$INCLUDE_FILES" == "true" ]]; then
        local files_size=$(du -sh "$BACKUP_DIR/files" 2>/dev/null | cut -f1 || echo "N/A")
        echo "- File Storage: Backed up (Size: $files_size)" >> "$summary_file"
    else
        echo "- File Storage: Skipped" >> "$summary_file"
    fi
    
    if [[ "$INCLUDE_CONFIG" == "true" ]]; then
        local config_size=$(du -sh "$BACKUP_DIR/config" 2>/dev/null | cut -f1 || echo "N/A")
        echo "- Configuration: Backed up (Size: $config_size)" >> "$summary_file"
    else
        echo "- Configuration: Skipped" >> "$summary_file"
    fi
    
    echo "" >> "$summary_file"
    echo "Total Backup Size: $(du -sh "$BACKUP_DIR" | cut -f1)" >> "$summary_file"
    echo "Retention Period: $RETENTION_DAYS days" >> "$summary_file"
    
    log "Backup summary created: $summary_file"
}

# Verify backup integrity
verify_backup() {
    log "Verifying backup integrity..."
    
    local errors=0
    
    # Verify database backup
    if [[ "$INCLUDE_DB" == "true" ]]; then
        log "Verifying database backup..."
        
        # Check if database backup directory exists and contains files
        if [ ! -d "$BACKUP_DIR/database" ] || [ -z "$(ls -A "$BACKUP_DIR/database")" ]; then
            log "Error: Database backup directory is empty or does not exist"
            errors=$((errors+1))
        else
            # Check if at least one backup file is not empty
            if ! find "$BACKUP_DIR/database" -type f -name "*.sql.gz" -size +0 | grep -q .; then
                log "Error: No valid database backup files found"
                errors=$((errors+1))
            else
                log "Database backup verification passed"
            fi
        fi
    fi
    
    # Verify file storage backup
    if [[ "$INCLUDE_FILES" == "true" ]]; then
        log "Verifying file storage backup..."
        
        # Check if file storage backup directory exists
        if [ ! -d "$BACKUP_DIR/files" ]; then
            log "Error: File storage backup directory does not exist"
            errors=$((errors+1))
        else
            # Check if MinIO container is running
            if docker ps | grep -q "minio"; then
                local minio_container=$(docker ps | grep "minio" | awk '{print $NF}' | head -1)
                
                # Get bucket list
                local buckets=$(docker exec "$minio_container" sh -c "mc ls local/ | awk '{print \$5}'")
                
                # Verify each bucket has been backed up
                for bucket in $buckets; do
                    if [ ! -d "$BACKUP_DIR/files/$bucket" ]; then
                        log "Error: Bucket $bucket was not backed up"
                        errors=$((errors+1))
                    fi
                done
                
                if [ $errors -eq 0 ]; then
                    log "File storage backup verification passed"
                fi
            else
                log "Warning: Cannot verify file storage backup as MinIO container is not running"
            fi
        fi
    fi
    
    # Verify configuration backup
    if [[ "$INCLUDE_CONFIG" == "true" ]]; then
        log "Verifying configuration backup..."
        
        # Check if configuration backup directory exists
        if [ ! -d "$BACKUP_DIR/config" ]; then
            log "Error: Configuration backup directory does not exist"
            errors=$((errors+1))
        else
            # Check for essential configuration files
            local essential_files=("docker-compose.yml" ".env")
            local missing=0
            
            for file in "${essential_files[@]}"; do
                if [ ! -f "$BACKUP_DIR/config/$file" ] && [ -f "$PROJECT_ROOT/$file" ]; then
                    log "Error: Essential configuration file $file is missing from backup"
                    missing=$((missing+1))
                fi
            done
            
            if [ $missing -gt 0 ]; then
                log "Error: $missing essential configuration file(s) missing from backup"
                errors=$((errors+1))
            else
                log "Configuration backup verification passed"
            fi
        fi
    fi
    
    if [ $errors -eq 0 ]; then
        log "Backup verification completed successfully"
        return 0
    else
        log "Backup verification failed with $errors error(s)"
        return 1
    fi
}

# Main function
main() {
    log "===================================================="
    log "Molecular Data Platform Backup Process Starting"
    log "===================================================="
    
    # Parse command-line arguments
    parse_arguments "$@"
    local parse_result=$?
    
    if [ $parse_result -ne 0 ]; then
        log "Failed to parse command-line arguments"
        return $parse_result
    fi
    
    # Create backup directory structure
    create_backup_directory
    local dir_result=$?
    
    if [ $dir_result -ne 0 ]; then
        log "Failed to create backup directory"
        return $dir_result
    fi
    
    # Perform component backups
    local backup_result=0
    
    # Database backup
    if [[ "$INCLUDE_DB" == "true" ]]; then
        backup_database
        local db_result=$?
        
        if [ $db_result -ne 0 ]; then
            log "Database backup failed"
            backup_result=1
        fi
    fi
    
    # File storage backup
    if [[ "$INCLUDE_FILES" == "true" ]]; then
        backup_file_storage
        local files_result=$?
        
        if [ $files_result -ne 0 ]; then
            log "File storage backup failed"
            backup_result=1
        fi
    fi
    
    # Configuration backup
    if [[ "$INCLUDE_CONFIG" == "true" ]]; then
        backup_configuration
        local config_result=$?
        
        if [ $config_result -ne 0 ]; then
            log "Configuration backup failed"
            backup_result=1
        fi
    fi
    
    # Verify backup integrity
    verify_backup
    local verify_result=$?
    
    if [ $verify_result -ne 0 ]; then
        log "Backup verification failed"
        backup_result=1
    fi
    
    # Create backup summary
    create_backup_summary
    
    # Clean up old backups
    cleanup_old_backups
    
    if [ $backup_result -eq 0 ]; then
        log "Backup completed successfully"
    else
        log "Backup completed with errors"
    fi
    
    log "===================================================="
    log "Molecular Data Platform Backup Process Completed"
    log "===================================================="
    
    return $backup_result
}

# Execute main function with all arguments
main "$@"
exit $?