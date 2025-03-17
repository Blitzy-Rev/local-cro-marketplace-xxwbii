#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e
# Treat unset variables as an error
set -u

# Main function to run linting tools
main() {
    echo "Linting code..."
    
    # Run flake8 for style checking
    flake8 src/backend/app src/backend/tests
    
    # Run mypy for type checking
    mypy src/backend/app
    
    # Add other linting tools as needed
    # black --check src/backend/app src/backend/tests
    # isort --check-only src/backend/app src/backend/tests
    # pylint src/backend/app
    
    echo "Linting complete! No issues found."
    return 0
}

# Execute main function
main
exit $?