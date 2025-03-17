#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

# Function to check and set environment variables with defaults
check_environment_variables() {
    # APP_MODULE - Python module path to FastAPI application
    export APP_MODULE=${APP_MODULE:-main:app}
    
    # HOST - Host address to bind the server to
    export HOST=${HOST:-0.0.0.0}
    
    # PORT - Port to bind the server to
    export PORT=${PORT:-8000}
    
    # LOG_LEVEL - Logging level
    export LOG_LEVEL=${LOG_LEVEL:-info}
    
    # WORKERS_PER_CORE - Number of worker processes per CPU core
    export WORKERS_PER_CORE=${WORKERS_PER_CORE:-1}
    
    # MAX_WORKERS - Maximum number of worker processes (0 means no limit)
    export MAX_WORKERS=${MAX_WORKERS:-0}
    
    echo "Environment variables set:"
    echo "APP_MODULE=$APP_MODULE"
    echo "HOST=$HOST"
    echo "PORT=$PORT"
    echo "LOG_LEVEL=$LOG_LEVEL"
    echo "WORKERS_PER_CORE=$WORKERS_PER_CORE"
    echo "MAX_WORKERS=$MAX_WORKERS"
}

# Function to wait for the database to be available with retries
wait_for_database() {
    echo "Checking database connectivity..."
    
    # Check if psql command is available
    if ! command -v psql &> /dev/null; then
        echo "Warning: psql command not found, cannot check database connectivity directly."
        echo "Database connectivity will be checked by prestart.sh script if it exists."
        return 0
    fi
    
    # Maximum number of retries
    MAX_RETRIES=5
    # Interval between retries in seconds
    RETRY_INTERVAL=5
    
    # Initialize retry counter
    retry_count=0
    
    # Check if required database environment variables are set
    if [ -z "$POSTGRES_SERVER" ] || [ -z "$POSTGRES_USER" ] || [ -z "$POSTGRES_PASSWORD" ] || [ -z "$POSTGRES_DB" ] || [ -z "$POSTGRES_PORT" ]; then
        echo "Warning: Required database environment variables are not set."
        echo "Required: POSTGRES_SERVER, POSTGRES_USER, POSTGRES_PASSWORD, POSTGRES_DB, POSTGRES_PORT"
        echo "Database connectivity check will be skipped."
        return 0
    fi
    
    # Try to connect to the database with retries
    while [ $retry_count -lt $MAX_RETRIES ]; do
        echo "Attempting database connection ($((retry_count+1))/$MAX_RETRIES)..."
        
        # Try to connect to PostgreSQL
        export PGPASSWORD=$POSTGRES_PASSWORD
        if psql -h "$POSTGRES_SERVER" -U "$POSTGRES_USER" -d "$POSTGRES_DB" -p "$POSTGRES_PORT" -c "SELECT 1" > /dev/null 2>&1; then
            echo "Database is available."
            return 0
        else
            echo "Database not available yet. Retrying in $RETRY_INTERVAL seconds..."
            retry_count=$((retry_count+1))
            sleep $RETRY_INTERVAL
        fi
    done
    
    echo "Warning: Could not connect to the database after $MAX_RETRIES attempts."
    echo "Will continue execution, but the application may not work properly if database is required."
    return 0  # Return success to continue execution
}

# Function to run prestart scripts
run_prestart_scripts() {
    echo "Running prestart scripts..."
    
    # Execute prestart.sh if it exists
    if [ -f ./scripts/prestart.sh ]; then
        echo "Running ./scripts/prestart.sh"
        bash ./scripts/prestart.sh
    else
        echo "No prestart.sh script found"
    fi
    
    # Execute scripts in the prestart directory if it exists
    if [ -d ./scripts/prestart ]; then
        echo "Running scripts in ./scripts/prestart/"
        for script in ./scripts/prestart/*.sh; do
            if [ -f "$script" ]; then
                echo "Running $script"
                bash "$script"
            fi
        done
    fi
}

# Function to start Gunicorn
start_gunicorn() {
    echo "Starting Gunicorn with APP_MODULE=$APP_MODULE"
    
    # Start Gunicorn with preloaded configuration from gunicorn_conf.py
    exec gunicorn "$APP_MODULE" \
        --worker-class uvicorn.workers.UvicornWorker \
        --bind "$HOST:$PORT" \
        --log-level "$LOG_LEVEL" \
        --config ./gunicorn_conf.py
}

# Main function
main() {
    echo "Starting entrypoint script for Molecular Data Management Platform..."
    
    # Check and set environment variables
    check_environment_variables
    
    # Wait for database to be available
    wait_for_database
    
    # Run prestart scripts
    run_prestart_scripts
    
    # Start Gunicorn
    start_gunicorn
}

# Execute main function
main