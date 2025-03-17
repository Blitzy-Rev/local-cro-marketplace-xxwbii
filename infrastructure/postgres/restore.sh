#!/bin/bash
#
# restore.sh - PostgreSQL database restoration script
#
# This script restores a PostgreSQL database from a backup file (compressed or uncompressed)
# for the Molecular Data Management and CRO Integration Platform. It handles connection
# verification, database creation if needed, and provides detailed logging throughout
# the restoration process.
#
# Usage: ./restore.sh -f backup_file.sql.gz [-d database_name]
#
# Environment variables can be used to configure database connection parameters:
# - POSTGRES_USER: Database username (default: postgres)
# - POSTGRES_PASSWORD: Database password (default: postgres)
# - POSTGRES_DB: Database name (default: molecular_platform)
# - POSTGRES_HOST: Database host (default: postgres)
# - POSTGRES_PORT: Database port (default: 5432)
#

# Exit immediately if a command exits with a non-zero status
set -e

# Default configuration
BACKUP_DIR="/backups"
POSTGRES_USER="${POSTGRES_USER:-postgres}"
POSTGRES_PASSWORD="${POSTGRES_PASSWORD:-postgres}"
POSTGRES_DB="${POSTGRES_DB:-molecular_platform}"
POSTGRES_HOST="${POSTGRES_HOST:-postgres}"
POSTGRES_PORT="${POSTGRES_PORT:-5432}"
TIMESTAMP_FORMAT="%Y-%m-%d %H:%M:%S"

# Display usage information
print_usage() {
    echo "Usage: $(basename "$0") [OPTIONS]"
    echo "Restore PostgreSQL database from a compressed backup file."
    echo ""
    echo "Options:"
    echo "  -f FILE    Specify the backup file to restore (required)"
    echo "  -d DB      Specify the target database name (default: $POSTGRES_DB)"
    echo "  -h         Display this help message and exit"
    echo ""
    echo "Environment variables:"
    echo "  POSTGRES_USER      Database user (default: postgres)"
    echo "  POSTGRES_PASSWORD  Database password (default: postgres)"
    echo "  POSTGRES_DB        Database name (default: molecular_platform)"
    echo "  POSTGRES_HOST      Database host (default: postgres)"
    echo "  POSTGRES_PORT      Database port (default: 5432)"
    echo ""
    echo "Examples:"
    echo "  $(basename "$0") -f /backups/backup_20230701.sql.gz"
    echo "  $(basename "$0") -f backup_20230701.sql.gz -d custom_db_name"
    echo ""
    echo "Note: This script requires the PostgreSQL client tools (psql) to be installed."
}

# Parse command line arguments
parse_arguments() {
    local OPTIND opt
    BACKUP_FILE=""
    HELP_REQUESTED=false
    
    while getopts ":f:d:h" opt; do
        case $opt in
            f)
                BACKUP_FILE="$OPTARG"
                ;;
            d)
                POSTGRES_DB="$OPTARG"
                ;;
            h)
                HELP_REQUESTED=true
                ;;
            \?)
                echo "Error: Invalid option -$OPTARG" >&2
                print_usage
                return 1
                ;;
            :)
                echo "Error: Option -$OPTARG requires an argument" >&2
                print_usage
                return 1
                ;;
        esac
    done
    
    # Check if help was requested
    if [ "$HELP_REQUESTED" = true ]; then
        return 0
    fi
    
    # Check if backup file was specified
    if [ -z "$BACKUP_FILE" ]; then
        echo "Error: Backup file (-f) is required" >&2
        print_usage
        return 1
    fi
    
    # If the backup file doesn't include a path, assume it's in the BACKUP_DIR
    if [[ "$BACKUP_FILE" != /* ]]; then
        BACKUP_FILE="$BACKUP_DIR/$BACKUP_FILE"
    fi
    
    return 0
}

# Check for required PostgreSQL tools
check_required_tools() {
    local tools=("psql" "pg_isready")
    local missing_tools=()
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        echo "Error: Required PostgreSQL tools are missing: ${missing_tools[*]}" >&2
        echo "Please install PostgreSQL client tools package." >&2
        return 1
    fi
    
    return 0
}

# Check if the backup file exists and is readable
check_backup_file() {
    local backup_file="$1"
    
    # Check if the file exists
    if [ ! -f "$backup_file" ]; then
        echo "Error: Backup file does not exist: $backup_file" >&2
        return 1
    fi
    
    # Check if the file is readable
    if [ ! -r "$backup_file" ]; then
        echo "Error: Cannot read backup file: $backup_file" >&2
        return 1
    fi
    
    # Get file size in bytes and convert to MB for display
    local size_bytes=$(stat -c%s "$backup_file")
    local size_mb=$(echo "scale=2; $size_bytes / 1048576" | bc)
    
    # Basic check to see if it's a gzipped file
    if [[ "$backup_file" == *.gz ]]; then
        if ! gzip -t "$backup_file" &>/dev/null; then
            echo "Error: Invalid gzip file: $backup_file" >&2
            return 1
        fi
        echo "Backup file validated: $backup_file (Compressed, $size_mb MB)"
    else
        echo "Backup file validated: $backup_file (Uncompressed, $size_mb MB)"
        # For uncompressed files, we could do a basic check if it looks like SQL
        if ! grep -q "PostgreSQL database dump" "$backup_file" &>/dev/null; then
            echo "Warning: File does not appear to be a PostgreSQL dump. Restore may fail." >&2
        fi
    fi
    
    return 0
}

# Check database connection
check_database_connection() {
    echo "Verifying database connection to $POSTGRES_HOST:$POSTGRES_PORT..."
    
    # Use pg_isready for initial connection check
    if ! pg_isready -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" &>/dev/null; then
        echo "Error: PostgreSQL server at $POSTGRES_HOST:$POSTGRES_PORT is not ready" >&2
        return 1
    fi
    
    # Try to execute a simple query
    export PGPASSWORD="$POSTGRES_PASSWORD"
    if ! psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "postgres" -c "SELECT 1;" &>/dev/null; then
        echo "Error: Could not connect to PostgreSQL server using provided credentials" >&2
        unset PGPASSWORD
        return 1
    fi
    unset PGPASSWORD
    
    echo "Database connection successful."
    return 0
}

# Check if database exists, create if it doesn't
ensure_database_exists() {
    local db_name="$1"
    local exists=false
    
    echo "Checking if database '$db_name' exists..."
    export PGPASSWORD="$POSTGRES_PASSWORD"
    
    # Check if database exists
    if psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "postgres" -t -c "SELECT 1 FROM pg_database WHERE datname = '$db_name';" | grep -q 1; then
        echo "Database '$db_name' exists."
        exists=true
    else
        echo "Database '$db_name' does not exist. Creating..."
        if ! psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "postgres" -c "CREATE DATABASE \"$db_name\";" &>/dev/null; then
            echo "Error: Failed to create database '$db_name'" >&2
            unset PGPASSWORD
            return 1
        fi
        echo "Database '$db_name' created successfully."
    fi
    
    unset PGPASSWORD
    return 0
}

# Terminate existing connections to the database
terminate_connections() {
    local db_name="$1"
    
    echo "Terminating existing connections to database '$db_name'..."
    export PGPASSWORD="$POSTGRES_PASSWORD"
    
    psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "postgres" -c "
        SELECT pg_terminate_backend(pg_stat_activity.pid)
        FROM pg_stat_activity
        WHERE pg_stat_activity.datname = '$db_name'
        AND pid <> pg_backend_pid();" &>/dev/null || true
    
    unset PGPASSWORD
    echo "Existing connections terminated."
    return 0
}

# Perform the database restore
perform_restore() {
    local backup_file="$1"
    local start_time=$(date +"$TIMESTAMP_FORMAT")
    
    echo "[$start_time] Starting database restore from: $backup_file"
    echo "Target database: $POSTGRES_DB on $POSTGRES_HOST:$POSTGRES_PORT"
    
    # Export password for PostgreSQL tools
    export PGPASSWORD="$POSTGRES_PASSWORD"
    
    # Terminate existing connections
    terminate_connections "$POSTGRES_DB"
    
    echo "Restoring database from backup..."
    local restore_start=$(date +%s)
    
    # Check if it's a gzipped file
    if [[ "$backup_file" == *.gz ]]; then
        echo "Restoring from compressed backup file..."
        gunzip -c "$backup_file" | psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "$POSTGRES_DB"
        restore_status=$?
    else
        echo "Restoring from uncompressed backup file..."
        psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "$POSTGRES_DB" -f "$backup_file"
        restore_status=$?
    fi
    
    local restore_end=$(date +%s)
    local restore_duration=$((restore_end - restore_start))
    
    # Clear the password environment variable
    unset PGPASSWORD
    
    local end_time=$(date +"$TIMESTAMP_FORMAT")
    
    if [ $restore_status -eq 0 ]; then
        echo "[$end_time] Database restore completed successfully in $restore_duration seconds."
    else
        echo "[$end_time] Database restore failed with status: $restore_status" >&2
    fi
    
    return $restore_status
}

# Log the restore status
log_restore_status() {
    local status="$1"
    local timestamp=$(date +"$TIMESTAMP_FORMAT")
    
    if [ $status -eq 0 ]; then
        echo "[$timestamp] SUCCESS: Database restore completed successfully from: $BACKUP_FILE"
        
        # Get the size of the restored database
        export PGPASSWORD="$POSTGRES_PASSWORD"
        local db_size=$(psql -h "$POSTGRES_HOST" -p "$POSTGRES_PORT" -U "$POSTGRES_USER" -d "postgres" -t -c "
            SELECT pg_size_pretty(pg_database_size('$POSTGRES_DB'));" | tr -d ' ')
        unset PGPASSWORD
        
        echo "[$timestamp] Restored database size: $db_size"
    else
        echo "[$timestamp] ERROR: Database restore failed with status: $status" >&2
        echo "[$timestamp] Please check PostgreSQL logs for more details." >&2
    fi
}

# Clean up resources and handle interrupts
cleanup() {
    echo "Cleaning up resources..."
    unset PGPASSWORD
    echo "Done."
}

# Main function
main() {
    echo "========================================================"
    echo "Database Restore Script - $(date +"$TIMESTAMP_FORMAT")"
    echo "========================================================"
    
    # Set up trap for cleanup on exit or interrupt
    trap cleanup EXIT INT TERM
    
    # Check for required tools
    check_required_tools
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    # Parse command line arguments
    parse_arguments "$@"
    local parse_status=$?
    
    # Display help and exit if requested
    if [ "$HELP_REQUESTED" = true ]; then
        print_usage
        return 0
    fi
    
    # Check if argument parsing was successful
    if [ $parse_status -ne 0 ]; then
        return $parse_status
    fi
    
    # Check the backup file
    check_backup_file "$BACKUP_FILE"
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    # Check database connection
    check_database_connection
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    # Ensure database exists
    ensure_database_exists "$POSTGRES_DB"
    if [ $? -ne 0 ]; then
        return 1
    fi
    
    # Perform the restore
    perform_restore "$BACKUP_FILE"
    local restore_status=$?
    
    # Log the restore status
    log_restore_status $restore_status
    
    # Return the status of the restore operation
    return $restore_status
}

# Call the main function with all script arguments
main "$@"