#!/bin/bash
# Shell script that starts the FastAPI backend application for the Molecular Data Management and CRO Integration Platform.
# It configures and launches the Gunicorn WSGI server with Uvicorn workers to serve the FastAPI application with optimal performance settings.

# Exit immediately if a command exits with a non-zero status
set -e

# Check and set default values for environment variables
# Python module path to the FastAPI application, defaults to app.main:app
APP_MODULE="${APP_MODULE:-app.main:app}"
# Host address to bind the server to, defaults to 0.0.0.0
HOST="${HOST:-0.0.0.0}"
# Port to bind the server to, defaults to 8000
PORT="${PORT:-8000}"
# Logging level, defaults to 'info'
LOG_LEVEL="${LOG_LEVEL:-info}"
# Number of worker processes per CPU core, defaults to 1
WORKERS_PER_CORE="${WORKERS_PER_CORE:-1}"
# Maximum number of worker processes, 0 means no limit
MAX_WORKERS="${MAX_WORKERS:-0}"
# Worker timeout in seconds, defaults to 120
TIMEOUT="${TIMEOUT:-120}"
# Graceful worker timeout in seconds, defaults to 120
GRACEFUL_TIMEOUT="${GRACEFUL_TIMEOUT:-120}"
# Connection keep-alive timeout in seconds, defaults to 5
KEEP_ALIVE="${KEEP_ALIVE:-5}"
# Enable auto-reload for development, defaults to false
RELOAD="${RELOAD:-false}"

check_environment_variables() {
  # Check if APP_MODULE is set, use default if not
  : ${APP_MODULE:="app.main:app"}
  # Check if HOST is set, use default if not
  : ${HOST:="0.0.0.0"}
  # Check if PORT is set, use default if not
  : ${PORT:="8000"}
  # Check if LOG_LEVEL is set, use default if not
  : ${LOG_LEVEL:="info"}
  # Check if WORKERS_PER_CORE is set, use default if not
  : ${WORKERS_PER_CORE:="1"}
  # Check if MAX_WORKERS is set, use default if not
  : ${MAX_WORKERS:="0"}
  # Check if TIMEOUT is set, use default if not
  : ${TIMEOUT:="120"}
  # Check if GRACEFUL_TIMEOUT is set, use default if not
  : ${GRACEFUL_TIMEOUT:="120"}
  # Check if KEEP_ALIVE is set, use default if not
  : ${KEEP_ALIVE:="5"}
  # Check if RELOAD is set, use default if not
  : ${RELOAD:="false"}
}

log_configuration() {
  # Log the APP_MODULE that will be used
  echo "Starting server with APP_MODULE: $APP_MODULE"
  # Log the HOST and PORT that will be used
  echo "Binding to host: $HOST and port: $PORT"
  # Log the LOG_LEVEL that will be used
  echo "Setting log level to: $LOG_LEVEL"
  # Log the worker configuration (WORKERS_PER_CORE, MAX_WORKERS)
  echo "Worker configuration: WORKERS_PER_CORE=$WORKERS_PER_CORE, MAX_WORKERS=$MAX_WORKERS"
  # Log the timeout settings (TIMEOUT, GRACEFUL_TIMEOUT, KEEP_ALIVE)
  echo "Timeout settings: TIMEOUT=$TIMEOUT, GRACEFUL_TIMEOUT=$GRACEFUL_TIMEOUT, KEEP_ALIVE=$KEEP_ALIVE"
  # Log whether auto-reload is enabled (RELOAD)
  echo "Auto-reload is enabled: $RELOAD"
}

start_gunicorn() {
  # Log that the server is starting
  echo "Starting Gunicorn..."

  # Construct Gunicorn command with appropriate arguments
  # Set worker class to UvicornWorker for ASGI support
  # Configure number of workers based on CPU cores and environment variables
  # Set bind address using HOST and PORT
  # Configure timeouts and keep-alive settings
  # Enable or disable auto-reload based on RELOAD environment variable
  
  CMD="gunicorn -k uvicorn.workers.UvicornWorker -b $HOST:$PORT --log-level $LOG_LEVEL"

  # Calculate number of workers
  NUM_CORES=$(nproc)
  WORKERS=$((NUM_CORES * WORKERS_PER_CORE))

  # Apply maximum worker limit if set
  if [[ "$MAX_WORKERS" -gt 0 ]] && [[ "$WORKERS" -gt "$MAX_WORKERS" ]]; then
    WORKERS=$MAX_WORKERS
  fi

  CMD="$CMD --workers $WORKERS"
  CMD="$CMD --timeout $TIMEOUT --graceful-timeout $GRACEFUL_TIMEOUT --keep-alive $KEEP_ALIVE"

  if [[ "$RELOAD" == "true" ]]; then
    CMD="$CMD --reload"
  fi

  CMD="$CMD $APP_MODULE"

  # Execute Gunicorn command to start the server
  echo "Executing: $CMD"
  eval $CMD

  # Capture and return the exit code
  EXIT_CODE=$?
  return $EXIT_CODE
}

main() {
  # Log that the start script is running
  echo "Starting Molecular Data Platform..."

  # Check and set environment variables
  check_environment_variables

  # Log the configuration
  log_configuration

  # Start Gunicorn with the configured settings
  start_gunicorn
  EXIT_CODE=$?

  # Handle any errors that occur during startup
  if [[ "$EXIT_CODE" -ne 0 ]]; then
    echo "Gunicorn failed to start. Exit code: $EXIT_CODE"
    exit $EXIT_CODE
  fi

  # Return the exit code from Gunicorn
  exit 0
}

# Execute the script
main