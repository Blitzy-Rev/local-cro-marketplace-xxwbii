#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

# Configuration
MAX_RETRIES=5
RETRY_INTERVAL=5  # seconds

# Function to check database connection
check_database_connection() {
    echo "Checking database connection..."
    export PGPASSWORD=$POSTGRES_PASSWORD
    
    # Try to connect to the database
    psql -h "$POSTGRES_SERVER" -U "$POSTGRES_USER" -d "$POSTGRES_DB" -p "$POSTGRES_PORT" -c "SELECT 1" > /dev/null 2>&1
    return $?
}

# Function to wait for the database with retries
wait_for_database() {
    echo "Waiting for database to be available..."
    local retry_count=0
    
    while [ $retry_count -lt $MAX_RETRIES ]; do
        if check_database_connection; then
            echo "Database is available."
            return 0
        fi
        
        retry_count=$((retry_count + 1))
        echo "Database not available yet. Retrying in $RETRY_INTERVAL seconds... ($retry_count/$MAX_RETRIES)"
        sleep $RETRY_INTERVAL
    done
    
    echo "Error: Could not connect to the database after $MAX_RETRIES attempts."
    return 1
}

# Function to run database migrations
run_migrations() {
    echo "Running database migrations..."
    alembic upgrade head
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "Database migrations completed successfully."
    else
        echo "Error: Database migrations failed with exit code $exit_code."
    fi
    
    return $exit_code
}

# Function to initialize database with seed data
initialize_data() {
    echo "Initializing database with seed data..."
    python -c "from app.db.init_db import init_db_data; init_db_data()"
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "Database initialization completed successfully."
    else
        echo "Error: Database initialization failed with exit code $exit_code."
    fi
    
    return $exit_code
}

# Main function
main() {
    echo "Starting prestart script for Molecular Data Management Platform..."
    
    # Check if required environment variables are set
    if [ -z "$POSTGRES_SERVER" ] || [ -z "$POSTGRES_USER" ] || [ -z "$POSTGRES_PASSWORD" ] || [ -z "$POSTGRES_DB" ] || [ -z "$POSTGRES_PORT" ]; then
        echo "Error: Required database environment variables are not set."
        echo "Required: POSTGRES_SERVER, POSTGRES_USER, POSTGRES_PASSWORD, POSTGRES_DB, POSTGRES_PORT"
        exit 1
    fi
    
    # Wait for database to be available
    if ! wait_for_database; then
        echo "Fatal error: Could not connect to the database."
        exit 1
    fi
    
    # Run database migrations
    if ! run_migrations; then
        echo "Fatal error: Failed to run database migrations."
        exit 1
    fi
    
    # Initialize database with seed data
    if ! initialize_data; then
        echo "Fatal error: Failed to initialize database."
        exit 1
    fi
    
    echo "Prestart tasks completed successfully."
    exit 0
}

# Execute main function
main