#!/bin/bash
# Docker container entrypoint script for Celery workers in the Molecular Data Management and CRO Integration Platform.
# This script initializes the environment, performs health checks, and starts Celery worker processes.

# Exit immediately if a command exits with a non-zero status
set -e

# Check and set default values for required environment variables
check_environment_variables() {
    echo "Checking environment variables..."
    
    # Set default values for variables if not provided
    export APP_MODULE=${APP_MODULE:-"app.worker.celery_app:app"}
    export QUEUE=${QUEUE:-"default"}
    export CONCURRENCY=${CONCURRENCY:-$(nproc)}
    export LOG_LEVEL=${LOG_LEVEL:-"info"}
    export MAX_TASKS_PER_CHILD=${MAX_TASKS_PER_CHILD:-1000}
    export TASK_TIME_LIMIT=${TASK_TIME_LIMIT:-3600}
    export TASK_SOFT_TIME_LIMIT=${TASK_SOFT_TIME_LIMIT:-3300}
    
    # Make sure CELERY_BROKER_URL is set
    if [ -z "$CELERY_BROKER_URL" ]; then
        # Construct default Redis URL if not provided
        REDIS_HOST=${REDIS_HOST:-"redis"}
        REDIS_PORT=${REDIS_PORT:-6379}
        REDIS_DB=${REDIS_DB:-0}
        export CELERY_BROKER_URL="redis://${REDIS_HOST}:${REDIS_PORT}/${REDIS_DB}"
        echo "CELERY_BROKER_URL not set, using default: redis://${REDIS_HOST}:${REDIS_PORT}/${REDIS_DB}"
    else
        echo "Using CELERY_BROKER_URL from environment variables"
    fi
    
    echo "Environment variables configured."
}

# Verify Redis connection before starting workers
check_redis_connection() {
    echo "Checking Redis connection..."
    
    # Extract Redis host and port from CELERY_BROKER_URL
    if [[ $CELERY_BROKER_URL =~ redis://([^@]+@)?([^:/]+)(:([0-9]+))? ]]; then
        REDIS_HOST="${BASH_REMATCH[2]}"
        REDIS_PORT="${BASH_REMATCH[4]:-6379}"
    else
        # Default fallback if URL parsing fails
        REDIS_HOST=${REDIS_HOST:-"redis"}
        REDIS_PORT=${REDIS_PORT:-6379}
    fi
    
    echo "Attempting to connect to Redis at $REDIS_HOST:$REDIS_PORT..."
    
    MAX_RETRIES=${REDIS_MAX_RETRIES:-5}
    RETRY_DELAY=${REDIS_RETRY_DELAY:-5}
    
    # Try to connect to Redis with retries
    for i in $(seq 1 $MAX_RETRIES); do
        echo "Attempt $i of $MAX_RETRIES to connect to Redis..."
        
        # Try redis-cli ping if available
        if command -v redis-cli &> /dev/null; then
            if redis-cli -h $REDIS_HOST -p $REDIS_PORT ping | grep -q "PONG"; then
                echo "Successfully connected to Redis."
                return 0
            fi
        # Fallback to nc if available
        elif command -v nc &> /dev/null; then
            if nc -z $REDIS_HOST $REDIS_PORT; then
                echo "Successfully connected to Redis."
                return 0
            fi
        # Last resort, try a basic connection with bash
        else
            if timeout 5 bash -c "</dev/tcp/$REDIS_HOST/$REDIS_PORT" 2>/dev/null; then
                echo "Successfully connected to Redis."
                return 0
            fi
        fi
        
        echo "Failed to connect to Redis. Retrying in $RETRY_DELAY seconds..."
        sleep $RETRY_DELAY
        
        # Exponential backoff
        RETRY_DELAY=$((RETRY_DELAY * 2))
    done
    
    echo "ERROR: Failed to connect to Redis after $MAX_RETRIES attempts."
    return 1
}

# Setup Python path to ensure all modules are accessible
setup_python_path() {
    echo "Setting up Python path..."
    
    # Add app directory to Python path if not already included
    if [[ ":$PYTHONPATH:" != *":/app:"* ]]; then
        export PYTHONPATH="${PYTHONPATH:+${PYTHONPATH}:}/app"
    fi
    
    echo "Python path configured: $PYTHONPATH"
}

# Main function
main() {
    echo "Starting worker entrypoint script..."
    
    # Check and set environment variables
    check_environment_variables
    
    # Set up Python path
    setup_python_path
    
    # Check Redis connection
    if ! check_redis_connection; then
        echo "Redis connection check failed. Exiting."
        exit 1
    fi
    
    echo "All checks passed. Executing start-worker.sh..."
    
    # Execute the start-worker.sh script, passing all arguments
    exec ./start-worker.sh "$@"
}

# Execute main function
main "$@"