import pytest  # version 0.7.0+
import os
import tempfile
import uuid
from pathlib import Path
from unittest.mock import patch, MagicMock

from app.exceptions import ValidationException, FileStorageException
from app.constants import MAX_FILE_SIZE_MB, ALLOWED_EXTENSIONS
from app.utils.file_utils import (
    get_file_extension,
    is_allowed_file,
    validate_file_size,
    get_mime_type,
    generate_unique_filename,
    sanitize_filename,
    create_directory_if_not_exists,
    is_csv_file,
    is_image_file,
    is_document_file,
    get_file_size,
    format_file_size,
    validate_file,
    validate_file_with_exception,
    split_bucket_and_object_name
)


def test_get_file_extension():
    """Test the get_file_extension function with various filenames."""
    # Test with standard filename
    assert get_file_extension("file.txt") == ".txt"
    
    # Test with filename containing multiple dots
    assert get_file_extension("file.name.with.dots.csv") == ".csv"
    
    # Test with filename without extension
    assert get_file_extension("filename") == ""
    
    # Test with uppercase extension
    assert get_file_extension("file.TXT") == ".txt"


def test_is_allowed_file():
    """Test the is_allowed_file function with allowed and disallowed extensions."""
    # Test with allowed extensions from constants
    for ext in ALLOWED_EXTENSIONS:
        assert is_allowed_file(f"file.{ext}") is True
    
    # Test with custom allowed extensions
    custom_extensions = {".txt", ".doc"}
    assert is_allowed_file("file.txt", custom_extensions) is True
    assert is_allowed_file("file.jpg", custom_extensions) is False
    
    # Test with disallowed extensions
    assert is_allowed_file("file.exe") is False
    
    # Test with filename without extension
    assert is_allowed_file("filename") is False


def test_validate_file_size():
    """Test the validate_file_size function with various file sizes."""
    max_size = 10 * 1024 * 1024  # 10MB in bytes
    
    # Test with file size under the limit
    assert validate_file_size(5 * 1024 * 1024, max_size) is True
    
    # Test with file size equal to the limit
    assert validate_file_size(max_size, max_size) is True
    
    # Test with file size over the limit
    assert validate_file_size(max_size + 1, max_size) is False
    
    # Test with custom size limit
    custom_limit = 1 * 1024  # 1KB
    assert validate_file_size(512, custom_limit) is True
    assert validate_file_size(2 * 1024, custom_limit) is False


def test_get_mime_type():
    """Test the get_mime_type function with different file types."""
    # Test with magic library success cases
    with patch('magic.from_file') as mock_magic:
        # Test with CSV file
        mock_magic.return_value = "text/csv"
        assert get_mime_type("test.csv") == "text/csv"
        
        # Test with image file
        mock_magic.return_value = "image/png"
        assert get_mime_type("test.png") == "image/png"
        
        # Test with document file
        mock_magic.return_value = "application/pdf"
        assert get_mime_type("test.pdf") == "application/pdf"
    
    # Test with magic library failure, fallback to mimetypes
    with patch('magic.from_file', side_effect=Exception("Magic library error")):
        with patch('mimetypes.guess_type') as mock_mimetypes:
            # Mimetypes successful fallback
            mock_mimetypes.return_value = ("text/plain", None)
            assert get_mime_type("test.txt") == "text/plain"
            
            # Both magic and mimetypes fail
            mock_mimetypes.return_value = (None, None)
            assert get_mime_type("unknown.xyz") == "application/octet-stream"


def test_generate_unique_filename():
    """Test the generate_unique_filename function."""
    # Mock uuid and datetime for predictable output
    test_uuid = "12345678-1234-5678-1234-567812345678"
    test_timestamp = "20230601_123456"
    
    with patch('uuid.uuid4') as mock_uuid, \
         patch('datetime.datetime') as mock_datetime:
        
        mock_uuid.return_value = uuid.UUID(test_uuid)
        mock_datetime.now.return_value.strftime.return_value = test_timestamp
        
        # Test with standard filename
        unique_name = generate_unique_filename("test.txt")
        assert unique_name.endswith(".txt")
        assert test_uuid in unique_name
        assert test_timestamp in unique_name
        
        # Test with prefix
        unique_name = generate_unique_filename("test.csv", prefix="upload")
        assert unique_name.startswith("upload_")
        assert unique_name.endswith(".csv")
        assert test_uuid in unique_name
        assert test_timestamp in unique_name
        
        # Test with filename without extension
        unique_name = generate_unique_filename("test")
        assert test_uuid in unique_name
        assert "." not in unique_name


def test_sanitize_filename():
    """Test the sanitize_filename function with various filenames."""
    # Test with filename containing invalid characters
    assert sanitize_filename("file:name*?.txt") == "file_name__.txt"
    
    # Test with filename containing spaces
    assert sanitize_filename(" filename.txt ") == "filename.txt"
    
    # Test with filename containing path separators
    assert sanitize_filename("path/to/file.txt") == "path_to_file.txt"
    assert sanitize_filename("path\\to\\file.txt") == "path_to_file.txt"
    
    # Test with empty filename
    assert sanitize_filename("") == "file"
    
    # Test with only invalid characters
    assert sanitize_filename(":::**??") == "______"


def test_create_directory_if_not_exists():
    """Test the create_directory_if_not_exists function."""
    # Test with temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        test_dir = os.path.join(temp_dir, "test_dir")
        
        # Test creating directory that doesn't exist
        assert create_directory_if_not_exists(test_dir) is True
        assert os.path.exists(test_dir)
        
        # Test with directory that already exists
        assert create_directory_if_not_exists(test_dir) is True
        
        # Test with nested directory path
        nested_dir = os.path.join(test_dir, "nested", "dir")
        assert create_directory_if_not_exists(nested_dir) is True
        assert os.path.exists(nested_dir)
    
    # Test with invalid directory path (should return False)
    with patch('pathlib.Path.mkdir', side_effect=PermissionError("Permission denied")):
        assert create_directory_if_not_exists("/nonexistent/path/that/should/fail") is False


def test_is_csv_file():
    """Test the is_csv_file function with various filenames."""
    # Test with CSV filename
    assert is_csv_file("file.csv") is True
    
    # Test with non-CSV filename
    assert is_csv_file("file.txt") is False
    
    # Test with uppercase CSV extension
    assert is_csv_file("file.CSV") is True
    
    # Test with filename without extension
    assert is_csv_file("filename") is False


def test_is_image_file():
    """Test the is_image_file function with various filenames."""
    # Test with PNG filename
    assert is_image_file("file.png") is True
    
    # Test with JPG filename
    assert is_image_file("file.jpg") is True
    
    # Test with JPEG filename
    assert is_image_file("file.jpeg") is True
    
    # Test with SVG filename
    assert is_image_file("file.svg") is True
    
    # Test with non-image filename
    assert is_image_file("file.txt") is False
    
    # Test with uppercase image extension
    assert is_image_file("file.PNG") is True


def test_is_document_file():
    """Test the is_document_file function with various filenames."""
    # Test with PDF filename
    assert is_document_file("file.pdf") is True
    
    # Test with DOCX filename
    assert is_document_file("file.docx") is True
    
    # Test with XLSX filename
    assert is_document_file("file.xlsx") is True
    
    # Test with non-document filename
    assert is_document_file("file.txt") is False
    
    # Test with uppercase document extension
    assert is_document_file("file.PDF") is True


def test_get_file_size():
    """Test the get_file_size function with existing and non-existing files."""
    # Create a temporary file with known content
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        content = b"Test content for file size testing"
        temp_file.write(content)
        temp_file.flush()
        
        # Test with existing file
        file_size = get_file_size(temp_file.name)
        assert file_size == len(content)
    
    # Clean up the temporary file
    os.unlink(temp_file.name)
    
    # Test with non-existent file
    assert get_file_size("nonexistent_file.txt") == -1


def test_format_file_size():
    """Test the format_file_size function with various file sizes."""
    # Test with bytes (< 1KB)
    assert format_file_size(512) == "512 B"
    
    # Test with kilobytes (< 1MB)
    assert format_file_size(1024) == "1.00 KB"
    assert format_file_size(1536) == "1.50 KB"
    
    # Test with megabytes (< 1GB)
    assert format_file_size(1024 * 1024) == "1.00 MB"
    assert format_file_size(int(1.5 * 1024 * 1024)) == "1.50 MB"
    
    # Test with gigabytes (< 1TB)
    assert format_file_size(1024 * 1024 * 1024) == "1.00 GB"
    
    # Test with terabytes
    assert format_file_size(1024 * 1024 * 1024 * 1024) == "1.00 TB"
    
    # Test with zero size
    assert format_file_size(0) == "0 B"
    
    # Test with negative size (should be treated as zero)
    assert format_file_size(-100) == "0 B"


def test_validate_file():
    """Test the validate_file function with various file scenarios."""
    # Test with valid file (allowed extension and size)
    assert validate_file("file.csv", 1024 * 1024) is True
    
    # Test with invalid extension
    assert validate_file("file.exe", 1024 * 1024) is False
    
    # Test with file size over limit
    max_size_bytes = MAX_FILE_SIZE_MB * 1024 * 1024
    assert validate_file("file.csv", max_size_bytes + 1) is False
    
    # Test with custom allowed extensions
    custom_extensions = {".txt"}
    assert validate_file("file.txt", 1024, custom_extensions) is True
    assert validate_file("file.jpg", 1024, custom_extensions) is False
    
    # Test with custom size limit
    custom_size = 1 * 1024  # 1KB
    assert validate_file("file.csv", 512, None, custom_size) is True
    assert validate_file("file.csv", 2 * 1024, None, custom_size) is False


def test_validate_file_with_exception():
    """Test the validate_file_with_exception function with various file scenarios."""
    # Test with valid file (allowed extension and size)
    assert validate_file_with_exception("file.csv", 1024 * 1024) is True
    
    # Test with invalid extension
    with pytest.raises(ValidationException) as exc_info:
        validate_file_with_exception("file.exe", 1024 * 1024)
    assert "Invalid file extension" in str(exc_info.value)
    
    # Test with file size over limit
    max_size_bytes = MAX_FILE_SIZE_MB * 1024 * 1024
    with pytest.raises(ValidationException) as exc_info:
        validate_file_with_exception("file.csv", max_size_bytes + 1)
    assert "exceeds maximum allowed size" in str(exc_info.value)
    
    # Test with custom allowed extensions
    custom_extensions = {".txt"}
    assert validate_file_with_exception("file.txt", 1024, custom_extensions) is True
    with pytest.raises(ValidationException) as exc_info:
        validate_file_with_exception("file.jpg", 1024, custom_extensions)
    assert "Invalid file extension" in str(exc_info.value)
    
    # Test with custom size limit
    custom_size = 1 * 1024  # 1KB
    assert validate_file_with_exception("file.csv", 512, None, custom_size) is True
    with pytest.raises(ValidationException) as exc_info:
        validate_file_with_exception("file.csv", 2 * 1024, None, custom_size)
    assert "exceeds maximum allowed size" in str(exc_info.value)


def test_split_bucket_and_object_name():
    """Test the split_bucket_and_object_name function with various file paths."""
    # Test with standard bucket/object path
    bucket, object_name = split_bucket_and_object_name("bucket/object.txt")
    assert bucket == "bucket"
    assert object_name == "object.txt"
    
    # Test with path containing multiple slashes
    bucket, object_name = split_bucket_and_object_name("bucket/path/to/object.txt")
    assert bucket == "bucket"
    assert object_name == "path/to/object.txt"
    
    # Test with path without bucket separator
    bucket, object_name = split_bucket_and_object_name("object.txt")
    assert bucket is None
    assert object_name == "object.txt"
    
    # Test with empty path
    bucket, object_name = split_bucket_and_object_name("")
    assert bucket is None
    assert object_name == ""