"""
Service for managing file storage operations in the Molecular Data Management and CRO Integration Platform.

This module provides a high-level interface for storing, retrieving, and managing files
using MinIO as the underlying object storage system, supporting various file types including 
CSV uploads, molecule images, experiment files, and result files.
"""

import os  # standard library
import io  # standard library
import typing  # standard library
from typing import Dict, List, Optional, Union, BinaryIO, Tuple, Any  # standard library
from datetime import datetime  # standard library

from minio import Minio  # v7.1.0+
from minio.error import S3Error  # v7.1.0+

from ..logging_config import logger
from ..exceptions import FileStorageException, ValidationException
from ..core.config import get_minio_settings, get_bucket_names
from ..utils.file_utils import (
    get_file_extension,
    is_allowed_file,
    validate_file_size,
    generate_unique_filename,
    get_mime_type
)
from ..constants import BUCKET_NAMES, MAX_FILE_SIZE_MB

# MinIO client singleton
_minio_client = None


def get_minio_client() -> Minio:
    """
    Gets or initializes the MinIO client singleton.
    
    Returns:
        Minio: Initialized MinIO client instance
    """
    global _minio_client
    
    if _minio_client is None:
        # Get MinIO settings from configuration
        settings = get_minio_settings()
        
        # Create MinIO client
        _minio_client = Minio(
            f"{settings['host']}:{settings['port']}",
            access_key=settings['access_key'],
            secret_key=settings['secret_key'],
            secure=settings['secure'],
            region=settings['region']
        )
        
        logger.info("MinIO client initialized")
    
    return _minio_client


def initialize_buckets() -> bool:
    """
    Ensures all required storage buckets exist.
    
    Returns:
        bool: True if all buckets exist or were created successfully
    """
    client = get_minio_client()
    bucket_names = get_bucket_names()
    
    try:
        for bucket_type, bucket_name in bucket_names.items():
            # Check if bucket exists
            if not client.bucket_exists(bucket_name):
                # Create bucket if it doesn't exist
                client.make_bucket(bucket_name)
                logger.info(f"Created bucket: {bucket_name}")
            else:
                logger.info(f"Bucket already exists: {bucket_name}")
        
        return True
    except S3Error as e:
        logger.error(f"Error initializing buckets: {str(e)}")
        raise FileStorageException(f"Failed to initialize storage buckets: {str(e)}", 
                                  {"error": str(e)})


def upload_file(file_object: BinaryIO, filename: str, bucket_name: str, 
                content_type: str = None, metadata: Dict[str, str] = None) -> str:
    """
    Uploads a file to the specified bucket.
    
    Args:
        file_object: File-like object to upload
        filename: Original filename
        bucket_name: Target bucket name
        content_type: MIME type of the file (optional)
        metadata: Additional metadata for the file (optional)
    
    Returns:
        str: Object name (path) of the uploaded file
        
    Raises:
        ValidationException: If file validation fails
        FileStorageException: If upload fails
    """
    try:
        # Validate file extension
        if not is_allowed_file(filename):
            raise ValidationException(
                f"Invalid file extension for {filename}",
                {"filename": filename}
            )
        
        # Generate unique filename to prevent collisions
        object_name = generate_unique_filename(filename)
        
        # Determine content type if not provided
        if content_type is None:
            content_type = get_mime_type(filename)
        
        # Set metadata if provided or initialize empty dict
        file_metadata = metadata or {}
        
        # Add original filename to metadata
        file_metadata["original-filename"] = filename
        
        # Get file size from file object
        file_object.seek(0, os.SEEK_END)
        file_size = file_object.tell()
        file_object.seek(0)  # Reset file position
        
        # Validate file size
        max_size_bytes = MAX_FILE_SIZE_MB * 1024 * 1024
        if file_size > max_size_bytes:
            raise ValidationException(
                f"File size exceeds maximum allowed size of {MAX_FILE_SIZE_MB}MB",
                {"filename": filename, "size": file_size, "max_size": max_size_bytes}
            )
        
        # Get MinIO client
        client = get_minio_client()
        
        # Upload file to MinIO
        client.put_object(
            bucket_name=bucket_name,
            object_name=object_name,
            data=file_object,
            length=file_size,
            content_type=content_type,
            metadata=file_metadata
        )
        
        logger.info(f"File uploaded: {filename} to {bucket_name}/{object_name}")
        return object_name
    
    except ValidationException:
        # Re-raise validation exceptions
        raise
    except S3Error as e:
        logger.error(f"S3 error uploading file {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload file {filename}: {str(e)}",
            {"filename": filename, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error uploading file {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload file {filename}: {str(e)}",
            {"filename": filename, "bucket": bucket_name, "error": str(e)}
        )


def download_file(object_name: str, bucket_name: str) -> Tuple[io.BytesIO, str, str]:
    """
    Downloads a file from the specified bucket.
    
    Args:
        object_name: Object name (path) in the bucket
        bucket_name: Source bucket name
    
    Returns:
        Tuple[io.BytesIO, str, str]: Tuple containing file data, filename, and content type
        
    Raises:
        FileStorageException: If download fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # Get file stats to retrieve metadata and content type
        stat = client.stat_object(bucket_name, object_name)
        
        # Create BytesIO object to store file data
        file_data = io.BytesIO()
        
        # Download file from MinIO
        response = client.get_object(bucket_name, object_name)
        
        # Read all data from response
        for data in response.stream(32*1024):
            file_data.write(data)
        
        # Close response
        response.close()
        response.release_conn()
        
        # Reset file position to beginning
        file_data.seek(0)
        
        # Extract original filename from metadata or use object name
        filename = stat.metadata.get("original-filename", object_name)
        
        # Get content type
        content_type = stat.content_type
        
        logger.info(f"File downloaded: {bucket_name}/{object_name}")
        return file_data, filename, content_type
    
    except S3Error as e:
        logger.error(f"S3 error downloading file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download file {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error downloading file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download file {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )


def get_file_url(object_name: str, bucket_name: str, expires: int = 3600) -> str:
    """
    Generates a presigned URL for temporary access to a file.
    
    Args:
        object_name: Object name (path) in the bucket
        bucket_name: Source bucket name
        expires: Expiration time in seconds (default: 1 hour)
    
    Returns:
        str: Presigned URL for accessing the file
        
    Raises:
        FileStorageException: If URL generation fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # Generate presigned URL
        url = client.presigned_get_object(
            bucket_name=bucket_name,
            object_name=object_name,
            expires=expires
        )
        
        logger.info(f"Generated presigned URL for {bucket_name}/{object_name}")
        return url
    
    except S3Error as e:
        logger.error(f"S3 error generating URL for {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error generating URL for {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )


def delete_file(object_name: str, bucket_name: str) -> bool:
    """
    Deletes a file from the specified bucket.
    
    Args:
        object_name: Object name (path) in the bucket
        bucket_name: Source bucket name
    
    Returns:
        bool: True if file was deleted successfully
        
    Raises:
        FileStorageException: If deletion fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # Delete object from bucket
        client.remove_object(bucket_name, object_name)
        
        logger.info(f"File deleted: {bucket_name}/{object_name}")
        return True
    
    except S3Error as e:
        logger.error(f"S3 error deleting file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to delete file {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error deleting file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to delete file {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )


def list_files(bucket_name: str, prefix: str = "", recursive: bool = True) -> List[Dict[str, Any]]:
    """
    Lists files in the specified bucket with optional prefix.
    
    Args:
        bucket_name: Bucket name to list files from
        prefix: Prefix to filter objects (optional)
        recursive: Whether to recursively list objects (default: True)
    
    Returns:
        List[Dict[str, Any]]: List of file information dictionaries
        
    Raises:
        FileStorageException: If listing fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # List objects in bucket
        objects = client.list_objects(bucket_name, prefix=prefix, recursive=recursive)
        
        # Convert objects to dictionaries
        result = []
        for obj in objects:
            # Get object stats to retrieve metadata
            try:
                stat = client.stat_object(bucket_name, obj.object_name)
                metadata = stat.metadata
                content_type = stat.content_type
            except:
                metadata = {}
                content_type = ""
            
            # Create dictionary with object information
            file_info = {
                "name": obj.object_name,
                "size": obj.size,
                "last_modified": obj.last_modified,
                "content_type": content_type,
                "metadata": metadata
            }
            
            result.append(file_info)
        
        logger.info(f"Listed {len(result)} files in {bucket_name} with prefix '{prefix}'")
        return result
    
    except S3Error as e:
        logger.error(f"S3 error listing files in {bucket_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to list files in {bucket_name}: {str(e)}",
            {"bucket": bucket_name, "prefix": prefix, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error listing files in {bucket_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to list files in {bucket_name}: {str(e)}",
            {"bucket": bucket_name, "prefix": prefix, "error": str(e)}
        )


def file_exists(object_name: str, bucket_name: str) -> bool:
    """
    Checks if a file exists in the specified bucket.
    
    Args:
        object_name: Object name (path) in the bucket
        bucket_name: Bucket name to check
    
    Returns:
        bool: True if file exists, False otherwise
        
    Raises:
        FileStorageException: If check fails for reasons other than file not found
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # Try to get object stats
        client.stat_object(bucket_name, object_name)
        
        # If successful, file exists
        return True
    
    except S3Error as e:
        # Check if error is due to file not found
        if e.code == 'NoSuchKey':
            return False
        
        # Other S3 error
        logger.error(f"S3 error checking if file exists {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to check if file exists {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error checking if file exists {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to check if file exists {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )


def get_file_metadata(object_name: str, bucket_name: str) -> Dict[str, Any]:
    """
    Gets metadata for a file in the specified bucket.
    
    Args:
        object_name: Object name (path) in the bucket
        bucket_name: Bucket name to check
    
    Returns:
        Dict[str, Any]: Dictionary containing file metadata
        
    Raises:
        FileStorageException: If metadata retrieval fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # Get object stats
        stat = client.stat_object(bucket_name, object_name)
        
        # Create metadata dictionary
        metadata = {
            "name": object_name,
            "size": stat.size,
            "last_modified": stat.last_modified,
            "content_type": stat.content_type,
            "metadata": stat.metadata
        }
        
        logger.info(f"Retrieved metadata for {bucket_name}/{object_name}")
        return metadata
    
    except S3Error as e:
        logger.error(f"S3 error getting metadata for {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to get metadata for {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )
    except Exception as e:
        logger.error(f"Error getting metadata for {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to get metadata for {object_name}: {str(e)}",
            {"object_name": object_name, "bucket": bucket_name, "error": str(e)}
        )


def copy_file(source_bucket: str, source_object: str, 
             destination_bucket: str, destination_object: str = None) -> str:
    """
    Copies a file from one location to another within MinIO.
    
    Args:
        source_bucket: Source bucket name
        source_object: Source object name
        destination_bucket: Destination bucket name
        destination_object: Destination object name (if None, uses source_object)
    
    Returns:
        str: Destination object name
        
    Raises:
        FileStorageException: If copy fails
    """
    try:
        # Get MinIO client
        client = get_minio_client()
        
        # If destination_object is not provided, use source_object
        if destination_object is None:
            destination_object = source_object
        
        # Copy object
        client.copy_object(
            bucket_name=destination_bucket,
            object_name=destination_object,
            source_bucket_name=source_bucket,
            source_object_name=source_object
        )
        
        logger.info(f"Copied file from {source_bucket}/{source_object} to {destination_bucket}/{destination_object}")
        return destination_object
    
    except S3Error as e:
        logger.error(f"S3 error copying file {source_object}: {str(e)}")
        raise FileStorageException(
            f"Failed to copy file {source_object}: {str(e)}",
            {
                "source_bucket": source_bucket,
                "source_object": source_object,
                "destination_bucket": destination_bucket,
                "destination_object": destination_object,
                "error": str(e)
            }
        )
    except Exception as e:
        logger.error(f"Error copying file {source_object}: {str(e)}")
        raise FileStorageException(
            f"Failed to copy file {source_object}: {str(e)}",
            {
                "source_bucket": source_bucket,
                "source_object": source_object,
                "destination_bucket": destination_bucket,
                "destination_object": destination_object,
                "error": str(e)
            }
        )


# Specialized file type functions

def upload_csv_file(file_object: BinaryIO, filename: str, metadata: Dict[str, str] = None) -> str:
    """
    Uploads a CSV file to the CSV uploads bucket.
    
    Args:
        file_object: File-like object to upload
        filename: Original filename
        metadata: Additional metadata for the file (optional)
    
    Returns:
        str: Object name (path) of the uploaded CSV file
        
    Raises:
        ValidationException: If file validation fails
        FileStorageException: If upload fails
    """
    try:
        # Get bucket name for CSV uploads
        bucket_name = BUCKET_NAMES["CSV_UPLOADS"]
        
        # Upload file with CSV content type
        return upload_file(
            file_object=file_object,
            filename=filename,
            bucket_name=bucket_name,
            content_type="text/csv",
            metadata=metadata
        )
    
    except (ValidationException, FileStorageException):
        # Re-raise exceptions
        raise
    except Exception as e:
        logger.error(f"Error uploading CSV file {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload CSV file {filename}: {str(e)}",
            {"filename": filename, "error": str(e)}
        )


def upload_result_file(file_object: BinaryIO, filename: str, metadata: Dict[str, str] = None) -> str:
    """
    Uploads a result file to the result files bucket.
    
    Args:
        file_object: File-like object to upload
        filename: Original filename
        metadata: Additional metadata for the file (optional)
    
    Returns:
        str: Object name (path) of the uploaded result file
        
    Raises:
        ValidationException: If file validation fails
        FileStorageException: If upload fails
    """
    try:
        # Get bucket name for result files
        bucket_name = BUCKET_NAMES["RESULT_FILES"]
        
        # Upload file
        return upload_file(
            file_object=file_object,
            filename=filename,
            bucket_name=bucket_name,
            content_type=None,  # Auto-detect content type
            metadata=metadata
        )
    
    except (ValidationException, FileStorageException):
        # Re-raise exceptions
        raise
    except Exception as e:
        logger.error(f"Error uploading result file {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload result file {filename}: {str(e)}",
            {"filename": filename, "error": str(e)}
        )


def upload_experiment_file(file_object: BinaryIO, filename: str, metadata: Dict[str, str] = None) -> str:
    """
    Uploads an experiment file to the experiment files bucket.
    
    Args:
        file_object: File-like object to upload
        filename: Original filename
        metadata: Additional metadata for the file (optional)
    
    Returns:
        str: Object name (path) of the uploaded experiment file
        
    Raises:
        ValidationException: If file validation fails
        FileStorageException: If upload fails
    """
    try:
        # Get bucket name for experiment files
        bucket_name = BUCKET_NAMES["EXPERIMENT_FILES"]
        
        # Upload file
        return upload_file(
            file_object=file_object,
            filename=filename,
            bucket_name=bucket_name,
            content_type=None,  # Auto-detect content type
            metadata=metadata
        )
    
    except (ValidationException, FileStorageException):
        # Re-raise exceptions
        raise
    except Exception as e:
        logger.error(f"Error uploading experiment file {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload experiment file {filename}: {str(e)}",
            {"filename": filename, "error": str(e)}
        )


def upload_molecule_image(file_object: BinaryIO, filename: str, metadata: Dict[str, str] = None) -> str:
    """
    Uploads a molecule image to the molecule images bucket.
    
    Args:
        file_object: File-like object to upload
        filename: Original filename
        metadata: Additional metadata for the file (optional)
    
    Returns:
        str: Object name (path) of the uploaded molecule image
        
    Raises:
        ValidationException: If file validation fails
        FileStorageException: If upload fails
    """
    try:
        # Get bucket name for molecule images
        bucket_name = BUCKET_NAMES["MOLECULE_IMAGES"]
        
        # Determine image content type based on extension
        extension = get_file_extension(filename).lower()
        content_type = None
        
        if extension == '.png':
            content_type = 'image/png'
        elif extension in ['.jpg', '.jpeg']:
            content_type = 'image/jpeg'
        elif extension == '.svg':
            content_type = 'image/svg+xml'
        
        # Upload file
        return upload_file(
            file_object=file_object,
            filename=filename,
            bucket_name=bucket_name,
            content_type=content_type,
            metadata=metadata
        )
    
    except (ValidationException, FileStorageException):
        # Re-raise exceptions
        raise
    except Exception as e:
        logger.error(f"Error uploading molecule image {filename}: {str(e)}")
        raise FileStorageException(
            f"Failed to upload molecule image {filename}: {str(e)}",
            {"filename": filename, "error": str(e)}
        )


def get_csv_file(object_name: str) -> Tuple[io.BytesIO, str, str]:
    """
    Downloads a CSV file from the CSV uploads bucket.
    
    Args:
        object_name: Object name (path) in the bucket
    
    Returns:
        Tuple[io.BytesIO, str, str]: Tuple containing file data, filename, and content type
        
    Raises:
        FileStorageException: If download fails
    """
    try:
        # Get bucket name for CSV uploads
        bucket_name = BUCKET_NAMES["CSV_UPLOADS"]
        
        # Download file
        return download_file(object_name, bucket_name)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error downloading CSV file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download CSV file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_result_file(object_name: str) -> Tuple[io.BytesIO, str, str]:
    """
    Downloads a result file from the result files bucket.
    
    Args:
        object_name: Object name (path) in the bucket
    
    Returns:
        Tuple[io.BytesIO, str, str]: Tuple containing file data, filename, and content type
        
    Raises:
        FileStorageException: If download fails
    """
    try:
        # Get bucket name for result files
        bucket_name = BUCKET_NAMES["RESULT_FILES"]
        
        # Download file
        return download_file(object_name, bucket_name)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error downloading result file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download result file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_experiment_file(object_name: str) -> Tuple[io.BytesIO, str, str]:
    """
    Downloads an experiment file from the experiment files bucket.
    
    Args:
        object_name: Object name (path) in the bucket
    
    Returns:
        Tuple[io.BytesIO, str, str]: Tuple containing file data, filename, and content type
        
    Raises:
        FileStorageException: If download fails
    """
    try:
        # Get bucket name for experiment files
        bucket_name = BUCKET_NAMES["EXPERIMENT_FILES"]
        
        # Download file
        return download_file(object_name, bucket_name)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error downloading experiment file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download experiment file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_molecule_image(object_name: str) -> Tuple[io.BytesIO, str, str]:
    """
    Downloads a molecule image from the molecule images bucket.
    
    Args:
        object_name: Object name (path) in the bucket
    
    Returns:
        Tuple[io.BytesIO, str, str]: Tuple containing file data, filename, and content type
        
    Raises:
        FileStorageException: If download fails
    """
    try:
        # Get bucket name for molecule images
        bucket_name = BUCKET_NAMES["MOLECULE_IMAGES"]
        
        # Download file
        return download_file(object_name, bucket_name)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error downloading molecule image {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to download molecule image {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_csv_file_url(object_name: str, expires: int = 3600) -> str:
    """
    Generates a presigned URL for a CSV file.
    
    Args:
        object_name: Object name (path) in the bucket
        expires: Expiration time in seconds (default: 1 hour)
    
    Returns:
        str: Presigned URL for accessing the CSV file
        
    Raises:
        FileStorageException: If URL generation fails
    """
    try:
        # Get bucket name for CSV uploads
        bucket_name = BUCKET_NAMES["CSV_UPLOADS"]
        
        # Generate URL
        return get_file_url(object_name, bucket_name, expires)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error generating URL for CSV file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for CSV file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_result_file_url(object_name: str, expires: int = 3600) -> str:
    """
    Generates a presigned URL for a result file.
    
    Args:
        object_name: Object name (path) in the bucket
        expires: Expiration time in seconds (default: 1 hour)
    
    Returns:
        str: Presigned URL for accessing the result file
        
    Raises:
        FileStorageException: If URL generation fails
    """
    try:
        # Get bucket name for result files
        bucket_name = BUCKET_NAMES["RESULT_FILES"]
        
        # Generate URL
        return get_file_url(object_name, bucket_name, expires)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error generating URL for result file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for result file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_experiment_file_url(object_name: str, expires: int = 3600) -> str:
    """
    Generates a presigned URL for an experiment file.
    
    Args:
        object_name: Object name (path) in the bucket
        expires: Expiration time in seconds (default: 1 hour)
    
    Returns:
        str: Presigned URL for accessing the experiment file
        
    Raises:
        FileStorageException: If URL generation fails
    """
    try:
        # Get bucket name for experiment files
        bucket_name = BUCKET_NAMES["EXPERIMENT_FILES"]
        
        # Generate URL
        return get_file_url(object_name, bucket_name, expires)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error generating URL for experiment file {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for experiment file {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )


def get_molecule_image_url(object_name: str, expires: int = 3600) -> str:
    """
    Generates a presigned URL for a molecule image.
    
    Args:
        object_name: Object name (path) in the bucket
        expires: Expiration time in seconds (default: 1 hour)
    
    Returns:
        str: Presigned URL for accessing the molecule image
        
    Raises:
        FileStorageException: If URL generation fails
    """
    try:
        # Get bucket name for molecule images
        bucket_name = BUCKET_NAMES["MOLECULE_IMAGES"]
        
        # Generate URL
        return get_file_url(object_name, bucket_name, expires)
    
    except FileStorageException:
        # Re-raise exception
        raise
    except Exception as e:
        logger.error(f"Error generating URL for molecule image {object_name}: {str(e)}")
        raise FileStorageException(
            f"Failed to generate URL for molecule image {object_name}: {str(e)}",
            {"object_name": object_name, "error": str(e)}
        )