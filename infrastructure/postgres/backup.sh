#!/bin/bash
# PostgreSQL Database Backup Script
# For Molecular Data Management and CRO Integration Platform
# Version: 1.0.0

# Exit immediately if a command exits with a non-zero status
set -e

# Default configuration variables
BACKUP_DIR=${BACKUP_DIR:-/backups}
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
POSTGRES_USER=${POSTGRES_USER:-postgres}
POSTGRES_PASSWORD=${POSTGRES_PASSWORD:-postgres}
POSTGRES_DB=${POSTGRES_DB:-molecular_platform}
POSTGRES_HOST=${POSTGRES_HOST:-postgres}
POSTGRES_PORT=${POSTGRES_PORT:-5432}
RETENTION_DAYS=${RETENTION_DAYS:-7}
BACKUP_FILE=${BACKUP_DIR}/db_backup_${TIMESTAMP}.sql.gz

# Print usage information
print_usage() {
    echo "Usage: $(basename $0) [OPTIONS]"
    echo
    echo "Automated PostgreSQL database backup script for the Molecular Data Management Platform"
    echo
    echo "Options:"
    echo "  -o DIR    Output directory for backup files (default: ${BACKUP_DIR})"
    echo "  -r DAYS   Retention period in days (default: ${RETENTION_DAYS})"
    echo "  -h        Display this help message and exit"
    echo
    echo "Environment Variables:"
    echo "  POSTGRES_USER     Database user (default: postgres)"
    echo "  POSTGRES_PASSWORD Database password"
    echo "  POSTGRES_DB       Database name (default: molecular_platform)"
    echo "  POSTGRES_HOST     Database host (default: postgres)"
    echo "  POSTGRES_PORT     Database port (default: 5432)"
    echo "  BACKUP_DIR        Backup directory (default: /backups)"
    echo "  RETENTION_DAYS    Number of days to keep backups (default: 7)"
    echo
    echo "Examples:"
    echo "  $(basename $0) -o /data/backups -r 14"
    echo "  POSTGRES_DB=custom_db $(basename $0)"
}

# Parse command-line arguments
parse_arguments() {
    local OPTIND
    while getopts "o:r:h" opt; do
        case $opt in
            o)
                BACKUP_DIR="$OPTARG"
                BACKUP_FILE="${BACKUP_DIR}/db_backup_${TIMESTAMP}.sql.gz"
                ;;
            r)
                RETENTION_DAYS="$OPTARG"
                if ! [[ "$RETENTION_DAYS" =~ ^[0-9]+$ ]]; then
                    echo "Error: Retention days must be a positive integer"
                    return 1
                fi
                ;;
            h)
                print_usage
                exit 0
                ;;
            \?)
                echo "Error: Invalid option -$OPTARG"
                print_usage
                return 1
                ;;
        esac
    done

    # Validate parameters
    if [ -z "$POSTGRES_DB" ]; then
        echo "Error: Database name is required"
        return 1
    fi

    if [ -z "$POSTGRES_HOST" ]; then
        echo "Error: Database host is required"
        return 1
    fi

    return 0
}

# Create backup directory if it doesn't exist
create_backup_directory() {
    if [ ! -d "$BACKUP_DIR" ]; then
        echo "Creating backup directory: $BACKUP_DIR"
        mkdir -p "$BACKUP_DIR"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create backup directory '$BACKUP_DIR'"
            return 1
        fi
    fi

    # Check directory is writable
    if [ ! -w "$BACKUP_DIR" ]; then
        echo "Error: Backup directory '$BACKUP_DIR' is not writable"
        return 1
    fi

    return 0
}

# Perform database backup
perform_backup() {
    echo "Starting database backup at $(date)"
    echo "Backup file: $BACKUP_FILE"

    # Check available disk space
    local available_space
    available_space=$(df -k "$BACKUP_DIR" | awk 'NR==2 {print $4}')
    local min_required_space=102400  # 100MB in KB as a safe minimum

    if [ "$available_space" -lt "$min_required_space" ]; then
        echo "Error: Insufficient disk space. Available: ${available_space}KB, Required: ${min_required_space}KB"
        return 1
    fi

    # Export database password as environment variable
    export PGPASSWORD="$POSTGRES_PASSWORD"

    # Execute pg_dump with compression
    pg_dump \
        --host="$POSTGRES_HOST" \
        --port="$POSTGRES_PORT" \
        --username="$POSTGRES_USER" \
        --dbname="$POSTGRES_DB" \
        --format=custom \
        --compress=9 \
        --verbose \
        --file="$BACKUP_FILE"

    local result=$?
    unset PGPASSWORD  # Clear password from environment

    return $result
}

# Clean up old backup files
cleanup_old_backups() {
    echo "Cleaning up backup files older than $RETENTION_DAYS days"
    
    local old_files
    old_files=$(find "$BACKUP_DIR" -name "db_backup_*.sql.gz" -type f -mtime +$RETENTION_DAYS)
    
    if [ -n "$old_files" ]; then
        echo "$old_files" | xargs rm -f
        echo "Deleted $(echo "$old_files" | wc -l) old backup files"
    else
        echo "No old backup files to delete"
    fi
}

# Log backup status
log_backup_status() {
    local status=$1
    
    if [ $status -eq 0 ]; then
        local file_size
        file_size=$(du -h "$BACKUP_FILE" | cut -f1)
        echo "Backup completed successfully at $(date)"
        echo "Backup file: $BACKUP_FILE (Size: $file_size)"
        # Add entry to backup log
        echo "$(date +"%Y-%m-%d %H:%M:%S") - Backup successful - File: $BACKUP_FILE - Size: $file_size" \
            >> "$BACKUP_DIR/backup.log"
    else
        echo "Backup failed with status $status at $(date)"
        echo "$(date +"%Y-%m-%d %H:%M:%S") - Backup failed with status $status" \
            >> "$BACKUP_DIR/backup.log"
    fi
}

# Main function
main() {
    echo "=========================================================="
    echo "Database Backup Started: $(date)"
    echo "=========================================================="

    parse_arguments "$@"
    local parse_result=$?
    
    if [ $parse_result -ne 0 ]; then
        return $parse_result
    fi

    create_backup_directory
    local dir_result=$?
    
    if [ $dir_result -ne 0 ]; then
        return $dir_result
    fi

    perform_backup
    local backup_result=$?
    
    log_backup_status $backup_result
    
    if [ $backup_result -eq 0 ]; then
        cleanup_old_backups
    fi

    echo "=========================================================="
    echo "Database Backup Completed: $(date)"
    echo "=========================================================="
    
    return $backup_result
}

# Execute main function with all arguments
main "$@"
exit $?