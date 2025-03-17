#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Default environment variables if not provided
APP_MODULE=${APP_MODULE:-app.main:app}
HOST=${HOST:-0.0.0.0}
PORT=${PORT:-8000}
LOG_LEVEL=${LOG_LEVEL:-info}
WORKERS_PER_CORE=${WORKERS_PER_CORE:-1}
MAX_WORKERS=${MAX_WORKERS:-0}
TIMEOUT=${TIMEOUT:-120}
GRACEFUL_TIMEOUT=${GRACEFUL_TIMEOUT:-120}
KEEP_ALIVE=${KEEP_ALIVE:-5}
RELOAD=${RELOAD:-false}

# Log startup information
echo "Starting Molecular Data Management and CRO Integration Platform backend"
echo "============================================================================"
echo "App module: $APP_MODULE"
echo "Host: $HOST"
echo "Port: $PORT"
echo "Log level: $LOG_LEVEL"
echo "Workers per CPU core: $WORKERS_PER_CORE"
echo "Max workers: $MAX_WORKERS"
echo "Timeout: $TIMEOUT"
echo "Graceful timeout: $GRACEFUL_TIMEOUT"
echo "Keep alive: $KEEP_ALIVE"
echo "Reload: $RELOAD"
echo "============================================================================"

# Build Gunicorn command
GUNICORN_CMD="gunicorn"

# Set the app module
GUNICORN_CMD="$GUNICORN_CMD $APP_MODULE"

# Set the Uvicorn worker class for ASGI support (required for FastAPI)
GUNICORN_CMD="$GUNICORN_CMD --worker-class uvicorn.workers.UvicornWorker"

# Calculate the number of workers based on CPU cores
if [ $MAX_WORKERS -gt 0 ]; then
    # If MAX_WORKERS is set, use it
    GUNICORN_CMD="$GUNICORN_CMD --workers $MAX_WORKERS"
else
    # Otherwise calculate based on CPU cores
    CORES=$(nproc)
    WORKERS=$(($CORES * $WORKERS_PER_CORE))
    # Ensure we have at least 1 worker
    if [ $WORKERS -lt 1 ]; then
        WORKERS=1
    fi
    echo "Auto-detected $CORES CPU cores. Using $WORKERS workers."
    GUNICORN_CMD="$GUNICORN_CMD --workers $WORKERS"
fi

# Set bind address
GUNICORN_CMD="$GUNICORN_CMD --bind $HOST:$PORT"

# Set timeouts
GUNICORN_CMD="$GUNICORN_CMD --timeout $TIMEOUT"
GUNICORN_CMD="$GUNICORN_CMD --graceful-timeout $GRACEFUL_TIMEOUT"
GUNICORN_CMD="$GUNICORN_CMD --keep-alive $KEEP_ALIVE"

# Set log level
GUNICORN_CMD="$GUNICORN_CMD --log-level $LOG_LEVEL"

# Set reload option for development
if [ "$RELOAD" = "true" ]; then
    GUNICORN_CMD="$GUNICORN_CMD --reload"
    echo "Auto-reload enabled (for development use only)"
fi

# Set access log format for structured logging
GUNICORN_CMD="$GUNICORN_CMD --access-logformat '%({X-Forwarded-For}i)s %(l)s %(u)s %(t)s \"%(r)s\" %(s)s %(b)s \"%(f)s\" \"%(a)s\" %(L)s'"

# Using --config to load additional configuration from gunicorn_conf.py
GUNICORN_CMD="$GUNICORN_CMD --config ../gunicorn_conf.py"

# Log the final command
echo "Executing: $GUNICORN_CMD"
echo "============================================================================"

# Execute Gunicorn
exec $GUNICORN_CMD