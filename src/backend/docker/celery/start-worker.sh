#!/bin/bash
# start-worker.sh - Script for starting Celery worker processes
# Celery version: 5.2+
#
# This script configures and starts Celery worker processes based on
# environment variables, optimizing settings for specific queue types
# and providing sensible defaults for a containerized environment.

# Exit immediately if a command exits with a non-zero status
set -e

# Function to configure Celery worker options based on environment variables
configure_worker_options() {
    # Set default values for options if not provided
    APP_MODULE=${APP_MODULE:-"app.worker.celery_app:app"}
    # Default queue selection based on QUEUE env var, or "default" if not specified
    QUEUE=${QUEUE:-"default"}
    
    # Determine optimal concurrency based on CPU cores if not specified
    if [ -z "$CONCURRENCY" ]; then
        # Use 2 workers per CPU core by default for I/O-bound tasks
        CONCURRENCY=$(($(nproc) * 2))
        # Cap at 8 workers maximum to prevent resource exhaustion
        if [ "$CONCURRENCY" -gt 8 ]; then
            CONCURRENCY=8
        fi
    fi
    
    # Other configuration options with sensible defaults
    LOG_LEVEL=${LOG_LEVEL:-"info"}
    MAX_TASKS_PER_CHILD=${MAX_TASKS_PER_CHILD:-1000}
    TASK_TIME_LIMIT=${TASK_TIME_LIMIT:-3600}
    TASK_SOFT_TIME_LIMIT=${TASK_SOFT_TIME_LIMIT:-3300}
    
    # Generate a unique hostname for this worker based on queue and PID
    HOSTNAME="worker-${QUEUE//,/_}-$$"
    
    # Build worker options
    WORKER_OPTS="-Q $QUEUE"
    WORKER_OPTS="$WORKER_OPTS -c $CONCURRENCY"
    WORKER_OPTS="$WORKER_OPTS -l $LOG_LEVEL"
    WORKER_OPTS="$WORKER_OPTS -n $HOSTNAME"
    WORKER_OPTS="$WORKER_OPTS --max-tasks-per-child=$MAX_TASKS_PER_CHILD"
    WORKER_OPTS="$WORKER_OPTS --time-limit=$TASK_TIME_LIMIT"
    WORKER_OPTS="$WORKER_OPTS --soft-time-limit=$TASK_SOFT_TIME_LIMIT"
    
    # Add options for specific queues
    if [[ "$QUEUE" == *"csv-processing"* ]]; then
        # CSV processing is memory-intensive, optimize accordingly
        WORKER_OPTS="$WORKER_OPTS --prefetch-multiplier=1"
    elif [[ "$QUEUE" == *"molecular-tasks"* ]]; then
        # Molecular tasks are CPU-intensive, optimize accordingly
        WORKER_OPTS="$WORKER_OPTS --prefetch-multiplier=1"
    elif [[ "$QUEUE" == *"notifications"* ]]; then
        # Notifications are I/O-bound, can process more at once
        WORKER_OPTS="$WORKER_OPTS --prefetch-multiplier=4"
    elif [[ "$QUEUE" == *"reports"* ]]; then
        # Report generation can be memory-intensive
        WORKER_OPTS="$WORKER_OPTS --prefetch-multiplier=1"
    fi
    
    # Add any additional options if provided
    if [ -n "$CELERY_WORKER_OPTS" ]; then
        WORKER_OPTS="$WORKER_OPTS $CELERY_WORKER_OPTS"
    fi
    
    echo "$WORKER_OPTS"
}

# Function to start the Celery worker
start_worker() {
    # Ensure environment variables are set
    APP_MODULE=${APP_MODULE:-"app.worker.celery_app:app"}
    QUEUE=${QUEUE:-"default"}
    CONCURRENCY=${CONCURRENCY:-$(nproc)}
    LOG_LEVEL=${LOG_LEVEL:-"info"}
    MAX_TASKS_PER_CHILD=${MAX_TASKS_PER_CHILD:-1000}
    TASK_TIME_LIMIT=${TASK_TIME_LIMIT:-3600}
    TASK_SOFT_TIME_LIMIT=${TASK_SOFT_TIME_LIMIT:-3300}
    
    echo "======================================================"
    echo "Starting Celery worker with the following configuration:"
    echo "  - Application: $APP_MODULE"
    echo "  - Queues: $QUEUE"
    echo "  - Concurrency: $CONCURRENCY"
    echo "  - Log level: $LOG_LEVEL"
    echo "  - Max tasks per child: $MAX_TASKS_PER_CHILD"
    echo "  - Task time limit: $TASK_TIME_LIMIT"
    echo "  - Task soft time limit: $TASK_SOFT_TIME_LIMIT"
    echo "======================================================"
    
    # Get worker options
    WORKER_OPTS=$(configure_worker_options)
    
    # Build the Celery worker command
    CMD="celery -A $APP_MODULE worker $WORKER_OPTS"
    
    echo "Executing: $CMD"
    
    # Execute the command
    # Using exec to replace the current process, ensuring signals are properly handled
    exec $CMD
}

# Main execution
echo "Initializing Celery worker..."

# Handle errors and provide useful error messages
if ! command -v celery &> /dev/null; then
    echo "ERROR: celery command not found. Please ensure Celery is installed."
    exit 1
fi

# Check if Redis is available (assuming Redis is the broker)
if [ -n "$CELERY_BROKER_URL" ]; then
    echo "Using broker URL from CELERY_BROKER_URL environment variable"
else
    echo "CELERY_BROKER_URL not set, using default Redis configuration"
fi

# Start the worker
start_worker