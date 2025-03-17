#!/bin/bash
# MinIO initialization script for Molecular Data Management and CRO Integration Platform
# Creates required buckets and sets appropriate policies
# Version: 1.0.0

# Exit immediately if a command exits with a non-zero status
set -e

# Function to configure MinIO client
setup_minio_client() {
  echo "Setting up MinIO client..."
  
  # Configure mc with MinIO endpoint and credentials
  echo "Configuring MinIO client with alias 'minio'..."
  mc alias set minio $MINIO_ENDPOINT $MINIO_ROOT_USER $MINIO_ROOT_PASSWORD
  
  # Wait for MinIO server to be ready with timeout
  echo "Waiting for MinIO server to be ready..."
  local retry_count=0
  local max_retries=30  # 5 minutes timeout (30 * 10 seconds)
  
  until mc ping minio --quiet; do
    retry_count=$((retry_count+1))
    if [ $retry_count -ge $max_retries ]; then
      echo "Timed out waiting for MinIO server to be ready after 5 minutes"
      return 1
    fi
    echo "MinIO server is not ready yet, waiting... (Attempt $retry_count/$max_retries)"
    sleep 10
  done
  
  # Verify connection
  if mc ping minio --quiet; then
    echo "Successfully connected to MinIO server"
    return 0
  else
    echo "Failed to connect to MinIO server"
    return 1
  fi
}

# Function to create a bucket if it doesn't exist
create_bucket_if_not_exists() {
  local bucket_name=$1
  local bucket_description=$2
  
  echo "Checking if bucket '$bucket_name' exists ($bucket_description)..."
  if ! mc ls minio | grep -q "$bucket_name"; then
    echo "Creating bucket '$bucket_name'..."
    if mc mb minio/$bucket_name; then
      # Set bucket policy to allow read/write access
      echo "Setting policy for bucket '$bucket_name'..."
      if mc policy set download minio/$bucket_name; then
        echo "Bucket '$bucket_name' created successfully with download policy"
      else
        echo "Warning: Failed to set policy for bucket '$bucket_name'"
      fi
    else
      echo "Error: Failed to create bucket '$bucket_name'"
      return 1
    fi
  else
    echo "Bucket '$bucket_name' already exists"
  fi
  
  return 0
}

# Function to create all required buckets
create_required_buckets() {
  echo "Creating required buckets..."
  local failed=0
  
  # Create buckets according to storage architecture
  create_bucket_if_not_exists "csv-uploads" "For CSV files" || failed=1
  create_bucket_if_not_exists "molecule-images" "For molecular structure images" || failed=1
  create_bucket_if_not_exists "experiment-files" "For experiment specifications" || failed=1
  create_bucket_if_not_exists "result-files" "For experimental results" || failed=1
  create_bucket_if_not_exists "temp-files" "For temporary processing files" || failed=1
  create_bucket_if_not_exists "system-backups" "For database and configuration backups" || failed=1
  
  if [ $failed -eq 0 ]; then
    echo "All required buckets created successfully"
    return 0
  else
    echo "Warning: Some buckets could not be created"
    return 1
  fi
}

# Main function
main() {
  echo "----------------------------------------"
  echo "Starting MinIO initialization..."
  echo "----------------------------------------"
  
  # Check if environment variables are set
  if [ -z "$MINIO_ROOT_USER" ] || [ -z "$MINIO_ROOT_PASSWORD" ] || [ -z "$MINIO_ENDPOINT" ]; then
    echo "Error: Required environment variables are not set"
    echo "Please ensure MINIO_ROOT_USER, MINIO_ROOT_PASSWORD, and MINIO_ENDPOINT are defined"
    return 1
  fi
  
  # Set up MinIO client
  if setup_minio_client; then
    # Create required buckets
    if create_required_buckets; then
      echo "----------------------------------------"
      echo "MinIO initialization completed successfully"
      echo "----------------------------------------"
      return 0
    else
      echo "----------------------------------------"
      echo "MinIO initialization partially completed with warnings"
      echo "----------------------------------------"
      return 1
    fi
  else
    echo "----------------------------------------"
    echo "MinIO initialization failed: couldn't set up MinIO client"
    echo "----------------------------------------"
    return 1
  fi
}

# Run the main function
main
exit $?