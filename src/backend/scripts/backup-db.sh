#!/bin/bash

# backup-db.sh - PostgreSQL database backup script for Molecular Data Management and CRO Integration Platform
# This script creates compressed backups of the database with proper timestamps and implements retention policies

# Exit immediately if a command exits with a non-zero status
set -e

# Set default values
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
BACKUP_DIR=${SCRIPT_DIR}/../../backups/database
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
BACKUP_FILE=${BACKUP_DIR}/db_backup_${TIMESTAMP}.sql.gz
DB_HOST=localhost
DB_PORT=5432
DB_NAME=molecular_platform
DB_USER=postgres
DB_PASSWORD=postgres
RETENTION_DAYS=7

# Print usage information
print_usage() {
    echo "Usage: $(basename $0) [options]"
    echo ""
    echo "Create a compressed backup of the PostgreSQL database for the Molecular Data Management Platform."
    echo ""
    echo "Options:"
    echo "  -o DIR    Specify output directory for backups (default: ../../backups/database)"
    echo "  -r DAYS   Specify number of days to retain backups (default: 7)"
    echo "  -h        Display this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $(basename $0)"
    echo "  $(basename $0) -o /path/to/backups -r 14"
    echo ""
    echo "Environment variables (can be set in .env file):"
    echo "  POSTGRES_HOST      Database host (default: localhost)"
    echo "  POSTGRES_PORT      Database port (default: 5432)"
    echo "  POSTGRES_DB        Database name (default: molecular_platform)"
    echo "  POSTGRES_USER      Database user (default: postgres)"
    echo "  POSTGRES_PASSWORD  Database password (default: postgres)"
}

# Parse command line arguments
parse_arguments() {
    while getopts ":o:r:h" opt; do
        case ${opt} in
            o )
                BACKUP_DIR=$OPTARG
                ;;
            r )
                RETENTION_DAYS=$OPTARG
                ;;
            h )
                print_usage
                exit 0
                ;;
            \? )
                echo "Invalid option: $OPTARG" 1>&2
                print_usage
                exit 1
                ;;
            : )
                echo "Invalid option: $OPTARG requires an argument" 1>&2
                print_usage
                exit 1
                ;;
        esac
    done
    
    # Validate parameters
    if ! [[ "$RETENTION_DAYS" =~ ^[0-9]+$ ]]; then
        echo "Error: Retention days must be a positive number" 1>&2
        return 1
    fi
    
    # Update backup file path based on potentially changed BACKUP_DIR
    BACKUP_FILE=${BACKUP_DIR}/db_backup_${TIMESTAMP}.sql.gz
    
    return 0
}

# Load environment variables from .env file if it exists
load_environment_variables() {
    ENV_FILE="${SCRIPT_DIR}/../../.env"
    if [ -f "$ENV_FILE" ]; then
        echo "Loading environment variables from $ENV_FILE"
        source "$ENV_FILE"
    fi
    
    # Override defaults with environment variables if they exist
    DB_HOST=${POSTGRES_HOST:-$DB_HOST}
    DB_PORT=${POSTGRES_PORT:-$DB_PORT}
    DB_NAME=${POSTGRES_DB:-$DB_NAME}
    DB_USER=${POSTGRES_USER:-$DB_USER}
    DB_PASSWORD=${POSTGRES_PASSWORD:-$DB_PASSWORD}
}

# Create backup directory if it doesn't exist
create_backup_directory() {
    if [ ! -d "$BACKUP_DIR" ]; then
        echo "Creating backup directory: $BACKUP_DIR"
        mkdir -p "$BACKUP_DIR"
        
        # Set secure permissions on backup directory
        chmod 700 "$BACKUP_DIR"
    fi
}

# Perform the database backup
perform_backup() {
    echo "Starting database backup at $(date)"
    echo "Backup file: $BACKUP_FILE"
    
    # Export password for pg_dump
    export PGPASSWORD="$DB_PASSWORD"
    
    # Perform the backup using pg_dump and compress with gzip
    pg_dump -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" \
        -F p -v | gzip > "$BACKUP_FILE"
    
    # Get the exit code of the backup command
    local result=$?
    
    # Clear the password from environment
    unset PGPASSWORD
    
    return $result
}

# Clean up old backups beyond retention period
cleanup_old_backups() {
    echo "Cleaning up backups older than $RETENTION_DAYS days"
    
    # Find and delete files older than RETENTION_DAYS
    find "$BACKUP_DIR" -name "db_backup_*.sql.gz" -type f -mtime +$RETENTION_DAYS -delete -print | \
        while read file; do
            echo "Deleted: $file"
        done
}

# Main function
main() {
    local exit_code=0
    
    # Parse command line arguments
    parse_arguments "$@"
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    # Load environment variables
    load_environment_variables
    
    # Create backup directory
    create_backup_directory
    
    # Perform backup
    perform_backup
    exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        # Get file size in a human-readable format
        local file_size=$(du -h "$BACKUP_FILE" | cut -f1)
        
        # Print success message
        echo "Backup completed successfully at $(date)"
        echo "Backup file size: $file_size"
        
        # Clean up old backups
        cleanup_old_backups
    else
        echo "Error: Backup failed with exit code $exit_code" >&2
    fi
    
    return $exit_code
}

# Call main function with all script arguments
main "$@"
exit $?