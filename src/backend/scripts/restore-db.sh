#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Script directory path
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")

# Default values
BACKUP_DIR="${SCRIPT_DIR}/../../backups/database"
DB_HOST="localhost"
DB_PORT="5432"
DB_NAME="molecular_platform"
DB_USER="postgres"
DB_PASSWORD="postgres"
BACKUP_FILE=""

# Display usage information
print_usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo ""
    echo "Restore the PostgreSQL database from a backup file for the Molecular Data Management and CRO Integration Platform."
    echo ""
    echo "Options:"
    echo "  -f, --file FILENAME   Specify the backup file to restore (if not provided, latest backup will be used)"
    echo "  -h, --help            Display this help message and exit"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") -f /path/to/backup.sql         # Restore from specific backup file"
    echo "  $(basename "$0") -f /path/to/backup.sql.gz      # Restore from compressed backup file"
    echo "  $(basename "$0")                                # Restore from latest backup in default location"
    echo ""
    echo "Note: Environment variables can be used to override database connection parameters:"
    echo "  POSTGRES_HOST, POSTGRES_PORT, POSTGRES_DB, POSTGRES_USER, POSTGRES_PASSWORD, DB_BACKUP_DIR"
}

# Parse command line arguments
parse_arguments() {
    # Parse command line options
    while getopts ":f:h" opt; do
        case ${opt} in
            f)
                BACKUP_FILE=$OPTARG
                ;;
            h)
                print_usage
                return 0
                ;;
            \?)
                echo "Error: Invalid option -$OPTARG" >&2
                print_usage
                return 1
                ;;
            :)
                echo "Error: Option -$OPTARG requires an argument." >&2
                print_usage
                return 1
                ;;
        esac
    done

    return 0
}

# Load environment variables
load_environment_variables() {
    # Try to load environment variables from .env file if it exists
    ENV_FILE="${SCRIPT_DIR}/../../.env"
    if [[ -f "$ENV_FILE" ]]; then
        echo "Loading environment variables from $ENV_FILE"
        source "$ENV_FILE"
    fi

    # Override defaults with environment variables if they exist
    DB_HOST=${POSTGRES_HOST:-$DB_HOST}
    DB_PORT=${POSTGRES_PORT:-$DB_PORT}
    DB_NAME=${POSTGRES_DB:-$DB_NAME}
    DB_USER=${POSTGRES_USER:-$DB_USER}
    DB_PASSWORD=${POSTGRES_PASSWORD:-$DB_PASSWORD}
    
    # Allow specific backup directory override
    BACKUP_DIR=${DB_BACKUP_DIR:-$BACKUP_DIR}
}

# Find the latest backup file
find_latest_backup() {
    if [[ ! -d "$BACKUP_DIR" ]]; then
        echo "Error: Backup directory $BACKUP_DIR does not exist." >&2
        return 1
    fi

    # Find the most recent backup file (including .gz files)
    latest_backup=$(find "$BACKUP_DIR" -type f \( -name "*.sql" -o -name "*.sql.gz" \) -printf "%T@ %p\n" | sort -nr | head -1 | cut -d' ' -f2-)
    
    if [[ -z "$latest_backup" ]]; then
        echo "Error: No backup files found in $BACKUP_DIR" >&2
        return 1
    fi
    
    # Show the backup file modification time
    backup_time=$(stat -c "%y" "$latest_backup")
    echo "Latest backup found: $latest_backup (modified: $backup_time)"
    
    echo "$latest_backup"
    return 0
}

# Verify backup file exists and is readable
verify_backup_file() {
    if [[ ! -f "$BACKUP_FILE" ]]; then
        echo "Error: Backup file $BACKUP_FILE does not exist." >&2
        return 1
    fi

    if [[ ! -r "$BACKUP_FILE" ]]; then
        echo "Error: Backup file $BACKUP_FILE is not readable." >&2
        return 1
    fi

    if [[ ! -s "$BACKUP_FILE" ]]; then
        echo "Error: Backup file $BACKUP_FILE is empty." >&2
        return 1
    fi

    return 0
}

# Perform database restore
perform_restore() {
    echo "Starting database restore from $BACKUP_FILE at $(date)"
    
    # Check if the file is gzipped
    if [[ "$BACKUP_FILE" == *.gz ]]; then
        echo "Detected compressed backup file, decompressing on-the-fly..."
        PGPASSWORD="$DB_PASSWORD" gunzip -c "$BACKUP_FILE" | psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME"
        return $?
    else
        echo "Detected uncompressed backup file..."
        PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" < "$BACKUP_FILE"
        return $?
    fi
}

# Main function
main() {
    # Parse command line arguments
    parse_arguments "$@"
    if [ $? -eq 0 ] && [ "$1" = "-h" ]; then
        exit 0
    elif [ $? -ne 0 ]; then
        exit 1
    fi

    # Load environment variables
    load_environment_variables

    # If no backup file specified, find the latest one
    if [[ -z "$BACKUP_FILE" ]]; then
        echo "No backup file specified, searching for latest backup..."
        BACKUP_FILE=$(find_latest_backup)
        if [ $? -ne 0 ]; then
            exit 1
        fi
        echo "Using latest backup file: $BACKUP_FILE"
    fi

    # Verify backup file
    verify_backup_file
    if [ $? -ne 0 ]; then
        exit 1
    fi

    # Perform database restore
    echo "Restoring database $DB_NAME from backup file $BACKUP_FILE..."
    perform_restore
    RESTORE_RESULT=$?

    # Check restore result
    if [ $RESTORE_RESULT -eq 0 ]; then
        echo "Database restore completed successfully at $(date)"
        echo "Database $DB_NAME has been restored from $BACKUP_FILE"
        
        # Simple verification of database connectivity
        echo "Verifying database connectivity..."
        if PGPASSWORD="$DB_PASSWORD" psql -h "$DB_HOST" -p "$DB_PORT" -U "$DB_USER" -d "$DB_NAME" -c "SELECT 1 as test_result;" > /dev/null 2>&1; then
            echo "Verification successful: Database $DB_NAME is accessible."
        else
            echo "Warning: Database restore completed but verification failed. The database may not be fully operational." >&2
            exit 1
        fi
    else
        echo "Error: Database restore failed with exit code $RESTORE_RESULT" >&2
        exit $RESTORE_RESULT
    fi

    return 0
}

# Execute main function with all arguments
main "$@"