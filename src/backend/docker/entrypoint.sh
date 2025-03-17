#!/bin/bash
# Exit immediately if a command exits with a non-zero status
set -e

# Function to check and set default values for environment variables
check_environment_variables() {
    echo "Checking environment variables..."
    
    # Check if APP_MODULE is set, use default if not
    if [ -z "$APP_MODULE" ]; then
        echo "APP_MODULE not set, using default: app.main:app"
        export APP_MODULE="app.main:app"
    fi
    
    # Check if PORT is set, use default if not
    if [ -z "$PORT" ]; then
        echo "PORT not set, using default: 8000"
        export PORT="8000"
    fi
    
    # Check if LOG_LEVEL is set, use default if not
    if [ -z "$LOG_LEVEL" ]; then
        echo "LOG_LEVEL not set, using default: info"
        export LOG_LEVEL="info"
    fi
    
    # Check if WORKERS_PER_CORE is set, use default if not
    if [ -z "$WORKERS_PER_CORE" ]; then
        echo "WORKERS_PER_CORE not set, using default: 1"
        export WORKERS_PER_CORE="1"
    fi
    
    # Check if MAX_WORKERS is set, use default if not
    if [ -z "$MAX_WORKERS" ]; then
        echo "MAX_WORKERS not set, using default: 4"
        export MAX_WORKERS="4"
    fi
    
    echo "Environment variables checked."
}

# Function to run prestart scripts
run_prestart_scripts() {
    echo "Running prestart scripts..."
    
    # Run prestart.sh
    if [ -f /app/scripts/prestart.sh ]; then
        echo "Running prestart.sh..."
        bash /app/scripts/prestart.sh
        prestart_exit_code=$?
        if [ $prestart_exit_code -ne 0 ]; then
            echo "Error: prestart.sh failed with exit code $prestart_exit_code"
            return $prestart_exit_code
        fi
    else
        echo "Warning: /app/scripts/prestart.sh not found."
    fi
    
    # Check if prestart directory exists and run scripts in it
    if [ -d /app/scripts/prestart.d ]; then
        echo "Running scripts in /app/scripts/prestart.d..."
        for script in /app/scripts/prestart.d/*.sh; do
            if [ -f "$script" ]; then
                echo "Running $script..."
                bash "$script"
                script_exit_code=$?
                if [ $script_exit_code -ne 0 ]; then
                    echo "Error: $script failed with exit code $script_exit_code"
                    return $script_exit_code
                fi
            fi
        done
    fi
    
    echo "Prestart scripts executed successfully."
    return 0
}

# Main function to execute the entrypoint script
main() {
    echo "Starting entrypoint script for Molecular Data Management Platform..."
    
    # Check environment variables
    check_environment_variables
    
    # Run prestart scripts
    run_prestart_scripts
    prestart_exit_code=$?
    if [ $prestart_exit_code -ne 0 ]; then
        echo "Error: Prestart scripts failed with exit code $prestart_exit_code"
        exit $prestart_exit_code
    fi
    
    # Execute the command passed to the entrypoint
    echo "Executing: $@"
    exec "$@"
}

# Run the main function with all arguments passed to the script
main "$@"