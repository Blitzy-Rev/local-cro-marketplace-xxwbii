"""
Unit tests for the security utility functions in the Molecular Data Management and CRO Integration Platform.
This module tests password validation, token generation, secure comparison, hashing,
and other security-related utilities to ensure they meet the security requirements of the application.
"""

import pytest
import re
from unittest.mock import MagicMock, patch

from ../../app.utils.security_utils import (
    validate_password_strength,
    get_password_validation_errors,
    generate_secure_token,
    is_common_password,
    calculate_password_strength,
    generate_random_password,
    secure_compare,
    hash_data,
    encode_base64,
    decode_base64,
    sanitize_input,
    PasswordValidator,
    TokenGenerator,
    PASSWORD_MIN_LENGTH
)
from ../../app.exceptions import AuthenticationException, ValidationException


def test_validate_password_strength_valid():
    """Tests that valid passwords pass the strength validation."""
    # Define valid passwords with different combinations of characters
    valid_passwords = [
        "Abcdef1!ghij",  # Meets all requirements
        "P@ssw0rd123",   # Meets all requirements
        "C0mpl3x!P@ssw0rd",  # More complex password
        "a".ljust(PASSWORD_MIN_LENGTH-3, "a") + "A1!",  # Minimum length with required chars
        "ThisIs@V3ryStrongP@ssword",  # Long password with all requirements
    ]
    
    for password in valid_passwords:
        assert validate_password_strength(password) is True, f"Password {password} should be valid"


def test_validate_password_strength_invalid():
    """Tests that invalid passwords fail the strength validation."""
    # Define invalid passwords missing different requirements
    invalid_passwords = [
        "short",  # Too short
        "abcdefghij",  # Missing uppercase, digit, and special character
        "ABCDEFGHIJ",  # Missing lowercase, digit, and special character
        "abcdEFGHIJ",  # Missing digit and special character
        "abcd1234567",  # Missing uppercase and special character
        "ABCD1234567",  # Missing lowercase and special character
        "abcdEFGH12",   # Missing special character
        "abcdEFGH!@",   # Missing digit
        "1234!@#$%^",   # Missing letters
    ]
    
    for password in invalid_passwords:
        assert validate_password_strength(password) is False, f"Password {password} should be invalid"


def test_get_password_validation_errors():
    """Tests that appropriate error messages are returned for invalid passwords."""
    # Define invalid passwords missing different requirements
    test_cases = [
        {
            "password": "short",
            "expected_errors": [f"Password must be at least {PASSWORD_MIN_LENGTH} characters long"]
        },
        {
            "password": "abcdefghijklm",
            "expected_errors": [
                "Password must contain at least one uppercase letter",
                "Password must contain at least one digit",
                "Password must contain at least one special character"
            ]
        },
        {
            "password": "ABCDEFGHIJKLM",
            "expected_errors": [
                "Password must contain at least one lowercase letter",
                "Password must contain at least one digit",
                "Password must contain at least one special character"
            ]
        },
        {
            "password": "abcdEFGHIJKLM",
            "expected_errors": [
                "Password must contain at least one digit",
                "Password must contain at least one special character"
            ]
        }
    ]
    
    for case in test_cases:
        errors = get_password_validation_errors(case["password"])
        for expected_error in case["expected_errors"]:
            assert expected_error in errors, f"Expected error '{expected_error}' not found in errors for password '{case['password']}'"
        assert len(errors) == len(case["expected_errors"]), f"Error count mismatch for password '{case['password']}'"


def test_generate_secure_token():
    """Tests that generated tokens have the correct length and format."""
    # Generate tokens of different lengths
    token_lengths = [8, 16, 32, 64]
    
    for length in token_lengths:
        token = generate_secure_token(length)
        assert len(token) == length, f"Token length should be {length}, got {len(token)}"
        
        # Verify that token contains only URL-safe characters
        url_safe_pattern = r'^[A-Za-z0-9_-]*$'
        assert re.match(url_safe_pattern, token), f"Token {token} contains non-URL-safe characters"
    
    # Ensure that multiple calls generate different tokens (randomness check)
    tokens = [generate_secure_token(16) for _ in range(10)]
    assert len(set(tokens)) == 10, "Multiple token generations should produce unique tokens"


@patch('../../app.utils.security_utils.get_settings')
@patch('../../app.utils.security_utils.os.path.exists')
@patch('builtins.open')
def test_is_common_password(mock_open, mock_path_exists, mock_get_settings):
    """Tests detection of common passwords."""
    # Mock the settings and file operations
    mock_settings = MagicMock()
    mock_settings.JWT_PRIVATE_KEY_PATH = "/path/to/private_key.pem"
    mock_get_settings.return_value = mock_settings
    
    # Setup the mock for path existence
    mock_path_exists.return_value = True
    
    # Setup mock for file reading
    mock_file = MagicMock()
    mock_file.__enter__.return_value = mock_file
    mock_file.readlines.return_value = ["password\n", "123456\n", "qwerty\n"]
    mock_open.return_value = mock_file
    
    # Test common passwords
    common_passwords = ["password", "123456", "qwerty"]
    for password in common_passwords:
        assert is_common_password(password) is True, f"Password '{password}' should be identified as common"
    
    # Test uncommon passwords
    uncommon_passwords = ["ComplexP@ssw0rd123", "N0tInTheL1st!"]
    for password in uncommon_passwords:
        assert is_common_password(password) is False, f"Password '{password}' should not be identified as common"
    
    # Test edge cases
    # File doesn't exist
    mock_path_exists.return_value = False
    assert is_common_password("password") is False, "Should not identify common passwords when file doesn't exist"
    
    # Exception during file reading
    mock_path_exists.return_value = True
    mock_open.side_effect = Exception("File error")
    assert is_common_password("password") is False, "Should handle exceptions gracefully"


def test_calculate_password_strength():
    """Tests password strength calculation algorithm."""
    # Define passwords with different strength characteristics
    test_cases = [
        {"password": "", "expected_min": 0, "expected_max": 10},  # Empty password
        {"password": "a", "expected_min": 0, "expected_max": 15},  # Very short password
        {"password": "password", "expected_min": 10, "expected_max": 40},  # Common password
        {"password": "Password1", "expected_min": 30, "expected_max": 60},  # Some complexity
        {"password": "Password1!", "expected_min": 50, "expected_max": 80},  # Good complexity
        {"password": "P@ssw0rd123!Abcd", "expected_min": 70, "expected_max": 100},  # Very strong
    ]
    
    for case in test_cases:
        strength = calculate_password_strength(case["password"])
        assert case["expected_min"] <= strength <= case["expected_max"], \
            f"Password '{case['password']}' strength {strength} should be between {case['expected_min']} and {case['expected_max']}"
    
    # Test relative strengths
    passwords = ["a", "password", "Password1", "Password1!", "P@ssw0rd123!Abcd"]
    strengths = [calculate_password_strength(p) for p in passwords]
    
    # Verify that strengths are non-decreasing
    for i in range(1, len(strengths)):
        assert strengths[i] >= strengths[i-1], \
            f"Password '{passwords[i]}' ({strengths[i]}) should be stronger than '{passwords[i-1]}' ({strengths[i-1]})"


def test_generate_random_password():
    """Tests generation of random passwords that meet requirements."""
    # Generate random passwords of different lengths
    password = generate_random_password()
    assert len(password) >= PASSWORD_MIN_LENGTH, f"Generated password should be at least {PASSWORD_MIN_LENGTH} characters"
    assert validate_password_strength(password), "Generated password should meet strength requirements"
    
    # Test with specific lengths
    lengths = [12, 16, 20, 32]
    for length in lengths:
        password = generate_random_password(length)
        assert len(password) == length, f"Generated password should be {length} characters long"
        assert validate_password_strength(password), "Generated password should meet strength requirements"
    
    # Test randomness
    passwords = [generate_random_password() for _ in range(10)]
    assert len(set(passwords)) == 10, "Multiple password generations should produce unique passwords"
    
    # Test minimum length enforcement
    short_length = PASSWORD_MIN_LENGTH - 2
    password = generate_random_password(short_length)
    assert len(password) >= PASSWORD_MIN_LENGTH, f"Password length should be adjusted to minimum {PASSWORD_MIN_LENGTH}"


def test_secure_compare():
    """Tests secure string comparison function."""
    # Compare identical strings and assert they are equal
    assert secure_compare("test", "test") is True, "Identical strings should compare as equal"
    assert secure_compare("", "") is True, "Empty strings should compare as equal"
    assert secure_compare("complex!password123", "complex!password123") is True, "Complex identical strings should compare as equal"
    
    # Compare different strings and assert they are not equal
    assert secure_compare("test", "different") is False, "Different strings should not compare as equal"
    assert secure_compare("test", "Test") is False, "Case differences should not compare as equal"
    assert secure_compare("", "notempty") is False, "Empty and non-empty strings should not compare as equal"
    
    # Test edge cases like strings of different lengths
    assert secure_compare("short", "shortlong") is False, "Strings of different lengths should not compare as equal"


def test_hash_data():
    """Tests data hashing functionality."""
    # Hash sample data strings
    data = "test data"
    hash_result = hash_data(data)
    
    # Assert that hashes are in the expected format (hex string)
    assert isinstance(hash_result, str), "Hash result should be a string"
    assert len(hash_result) == 64, "SHA-256 hash should be 64 characters long"
    assert re.match(r'^[0-9a-f]+$', hash_result), "Hash should be a hexadecimal string"
    
    # Assert that identical inputs produce identical hashes
    hash_result2 = hash_data(data)
    assert hash_result == hash_result2, "Identical inputs should produce identical hashes"
    
    # Assert that different inputs produce different hashes
    data2 = "different data"
    hash_result3 = hash_data(data2)
    assert hash_result != hash_result3, "Different inputs should produce different hashes"
    
    # Test edge cases like empty strings
    empty_hash = hash_data("")
    assert len(empty_hash) == 64, "Empty string hash should still be 64 characters"
    
    # Test non-string inputs
    binary_data = b"binary data"
    binary_hash = hash_data(binary_data)
    assert len(binary_hash) == 64, "Binary data hash should be 64 characters"


def test_encode_decode_base64():
    """Tests base64 encoding and decoding functions."""
    # Encode sample data strings to base64
    test_strings = [
        "hello world",
        "special characters: !@#$%^&*()_+",
        "",  # Empty string
        "unicode characters: ñáéíóú",
        "1234567890"
    ]
    
    for test_str in test_strings:
        encoded = encode_base64(test_str)
        assert isinstance(encoded, str), "Encoded result should be a string"
        
        # Assert that encoded strings are in valid base64 format
        base64_pattern = r'^[A-Za-z0-9+/]*={0,2}$'
        assert re.match(base64_pattern, encoded), f"Encoded string {encoded} should be valid base64"
        
        # Decode the encoded strings
        decoded = decode_base64(encoded)
        assert decoded == test_str, f"Decoded string {decoded} should match original {test_str}"
    
    # Test edge cases like binary data
    binary_data = b"\x00\x01\x02\x03\x04"
    encoded = encode_base64(binary_data)
    decoded = decode_base64(encoded)
    assert decoded.encode() == binary_data, "Binary data should encode and decode correctly"


def test_sanitize_input():
    """Tests input sanitization for preventing injection attacks."""
    # Define input strings with potentially dangerous content
    test_cases = [
        {
            "input": "<script>alert('XSS')</script>",
            "expected": "alert('XSS')"  # Removed script tags
        },
        {
            "input": "Normal text",
            "expected": "Normal text"  # No change for safe text
        },
        {
            "input": "Text with <iframe src='evil.com'>",
            "expected": "Text with  src='evil.com'"  # Removed iframe tags
        },
        {
            "input": 'User input with quotes "double" and \'single\'',
            "expected": 'User input with quotes &quot;double&quot; and &#x27;single&#x27;'  # Escaped quotes
        },
        {
            "input": "Input with & ampersand",
            "expected": "Input with &amp; ampersand"  # Escaped ampersand
        },
        {
            "input": "",
            "expected": ""  # Empty string remains empty
        }
    ]
    
    for case in test_cases:
        sanitized = sanitize_input(case["input"])
        assert sanitized == case["expected"], \
            f"Input '{case['input']}' should sanitize to '{case['expected']}', got '{sanitized}'"


def test_password_validator_class():
    """Tests the PasswordValidator class functionality."""
    # Create PasswordValidator instances with different configurations
    validator = PasswordValidator()
    
    # Test validation of passwords against different requirements
    assert validator.validate("ValidP@ss123") is True, "Valid password should pass validation"
    
    # Test invalid passwords
    assert validator.validate("short") is False, "Short password should fail validation"
    assert validator.validate("nouppercase123!") is False, "Password without uppercase should fail validation"
    assert validator.validate("NOLOWERCASE123!") is False, "Password without lowercase should fail validation"
    assert validator.validate("NoDigits!") is False, "Password without digits should fail validation"
    assert validator.validate("NoSpecialChars123") is False, "Password without special chars should fail validation"
    
    # Test getting validation errors for invalid passwords
    errors = validator.get_validation_errors("short")
    assert any("length" in error.lower() for error in errors), "Should report length error"
    
    errors = validator.get_validation_errors("nouppercase123!")
    assert any("uppercase" in error.lower() for error in errors), "Should report uppercase error"
    
    # Test with custom configuration
    custom_validator = PasswordValidator(
        min_length=6,
        require_lowercase=True,
        require_uppercase=False,
        require_digit=True,
        require_special=False,
        check_common_passwords=True
    )
    
    # Should pass with the custom configuration
    assert custom_validator.validate("simple1") is True, "Should pass with custom validator settings"
    
    # Should fail even with custom configuration
    assert custom_validator.validate("short") is False, "Too short for custom validator"
    assert custom_validator.validate("nodigit") is False, "No digit for custom validator"
    
    # Test calculating password strength
    strong_password = "StrongP@ssw0rd123"
    weak_password = "weak123"
    
    strong_score = validator.calculate_strength(strong_password)
    weak_score = validator.calculate_strength(weak_password)
    
    assert strong_score > weak_score, "Strong password should have higher strength score"


def test_token_generator_class():
    """Tests the TokenGenerator class functionality."""
    # Create a TokenGenerator instance
    token_gen = TokenGenerator()
    
    # Test generating tokens of different lengths
    token_lengths = [8, 16, 32]
    for length in token_lengths:
        token = token_gen.generate_token(length)
        assert len(token) == length, f"Token should be {length} characters"
        
        # Verify that token contains only URL-safe characters
        url_safe_pattern = r'^[A-Za-z0-9_-]*$'
        assert re.match(url_safe_pattern, token), f"Token {token} should be URL-safe"
    
    # Test generating verification codes
    code_lengths = [4, 6, 8]
    for length in code_lengths:
        code = token_gen.generate_verification_code(length)
        assert len(code) == length, f"Verification code should be {length} digits"
        assert code.isdigit(), "Verification code should contain only digits"
    
    # Test generating API keys
    api_key = token_gen.generate_api_key()
    assert api_key.startswith("api_"), "API key should have 'api_' prefix"
    assert len(api_key) > 5, "API key should have sufficient length"