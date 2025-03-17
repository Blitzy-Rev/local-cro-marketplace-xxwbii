#!/bin/bash

# Set script to exit immediately if a command exits with a non-zero status
set -e
# Set script to treat unset variables as an error
set -u

# Main function that runs tests on the Python codebase
main() {
    echo "Running tests..."
    python -m pytest src/backend/tests -v --cov=src/backend/app --cov-report=term-missing --cov-report=html --cov-report=xml
    echo "Tests completed successfully!"
    return 0
}

# Execute the main function
main
exit $?