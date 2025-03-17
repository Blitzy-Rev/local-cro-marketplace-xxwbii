#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Treat unset variables as an error
set -u

# Print formatting start message
echo "Formatting code..."

# Run Black code formatter
black src/backend/app src/backend/tests

# Run isort to sort imports
isort src/backend/app src/backend/tests

# Print success message
echo "Code formatting complete!"

# Exit with success code
exit 0