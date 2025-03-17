"""
Utility module providing security-related functions for password validation, token generation,
and other security operations for the Molecular Data Management and CRO Integration Platform.
This module supports the application's security requirements without external dependencies.
"""

import re  # standard library
import typing
from typing import List, Dict, Optional, Any, Tuple  # standard library
import secrets  # standard library
import string  # standard library
import base64  # standard library
import hashlib  # standard library

from ..exceptions import AuthenticationException, ValidationException
from ..core.config import get_settings

# Constants for password validation
PASSWORD_MIN_LENGTH = 10
PASSWORD_REGEX_LOWERCASE = r'[a-z]'
PASSWORD_REGEX_UPPERCASE = r'[A-Z]'
PASSWORD_REGEX_DIGIT = r'[0-9]'
PASSWORD_REGEX_SPECIAL = r'[!@#$%^&*(),.?":{}|<>]'
COMMON_PASSWORDS_FILE = "common_passwords.txt"

def validate_password_strength(password: str) -> bool:
    """
    Validates a password against security requirements.
    
    Args:
        password: The password to validate
        
    Returns:
        True if password meets requirements, False otherwise
    """
    # Check minimum length
    if len(password) < PASSWORD_MIN_LENGTH:
        return False
        
    # Check for presence of character types
    if not re.search(PASSWORD_REGEX_LOWERCASE, password):
        return False
        
    if not re.search(PASSWORD_REGEX_UPPERCASE, password):
        return False
        
    if not re.search(PASSWORD_REGEX_DIGIT, password):
        return False
        
    if not re.search(PASSWORD_REGEX_SPECIAL, password):
        return False
    
    # Check if it's a common password
    if is_common_password(password):
        return False
        
    return True

def get_password_validation_errors(password: str) -> List[str]:
    """
    Returns a list of validation errors for a password.
    
    Args:
        password: The password to validate
        
    Returns:
        List of validation error messages
    """
    errors = []
    
    # Check minimum length
    if len(password) < PASSWORD_MIN_LENGTH:
        errors.append(f"Password must be at least {PASSWORD_MIN_LENGTH} characters long")
        
    # Check for presence of character types
    if not re.search(PASSWORD_REGEX_LOWERCASE, password):
        errors.append("Password must contain at least one lowercase letter")
        
    if not re.search(PASSWORD_REGEX_UPPERCASE, password):
        errors.append("Password must contain at least one uppercase letter")
        
    if not re.search(PASSWORD_REGEX_DIGIT, password):
        errors.append("Password must contain at least one digit")
        
    if not re.search(PASSWORD_REGEX_SPECIAL, password):
        errors.append("Password must contain at least one special character")
        
    # Check if it's a common password
    if is_common_password(password):
        errors.append("Password is too common and easily guessable")
        
    return errors

def generate_secure_token(length: int) -> str:
    """
    Generates a cryptographically secure random token.
    
    Args:
        length: The desired length of the token
        
    Returns:
        Secure random token
    """
    # Generate a token with sufficient entropy
    token = secrets.token_urlsafe(length)
    
    # Ensure exact length by truncating or padding
    if len(token) > length:
        return token[:length]
    
    return token

def is_common_password(password: str) -> bool:
    """
    Checks if a password is in a list of common passwords.
    
    Args:
        password: The password to check
        
    Returns:
        True if password is common, False otherwise
    """
    try:
        settings = get_settings()
        common_passwords_path = settings.JWT_PRIVATE_KEY_PATH.rsplit("/", 1)[0] + "/" + COMMON_PASSWORDS_FILE
        
        # If common passwords file doesn't exist, consider no passwords common
        if not os.path.exists(common_passwords_path):
            return False
            
        with open(common_passwords_path, "r") as f:
            common_passwords = [line.strip() for line in f]
            
        return password.lower() in common_passwords
    except Exception:
        # If there's any error reading the file, fail safe (assume not common)
        return False

def calculate_password_strength(password: str) -> int:
    """
    Calculates a password strength score based on complexity.
    
    Args:
        password: The password to evaluate
        
    Returns:
        Password strength score (0-100)
    """
    score = 0
    
    # Length contribution (up to 25 points)
    length_score = min(len(password) * 2, 25)
    score += length_score
    
    # Character type diversity (up to 35 points)
    if re.search(PASSWORD_REGEX_LOWERCASE, password):
        score += 10
    if re.search(PASSWORD_REGEX_UPPERCASE, password):
        score += 10
    if re.search(PASSWORD_REGEX_DIGIT, password):
        score += 5
    if re.search(PASSWORD_REGEX_SPECIAL, password):
        score += 10
        
    # Non-repeating character bonus (up to 15 points)
    unique_chars = len(set(password))
    unique_ratio = unique_chars / len(password)
    score += int(unique_ratio * 15)
    
    # Penalize common patterns (up to -25 points)
    # Check for sequential characters
    has_sequential = False
    for i in range(len(password) - 2):
        if (ord(password[i]) + 1 == ord(password[i+1]) and
            ord(password[i+1]) + 1 == ord(password[i+2])):
            has_sequential = True
            break
    
    if has_sequential:
        score -= 10
        
    # Check for repeated patterns
    has_repeated = False
    for pattern_length in range(2, min(5, len(password) // 2 + 1)):
        for i in range(len(password) - pattern_length * 2 + 1):
            pattern = password[i:i+pattern_length]
            if pattern in password[i+pattern_length:]:
                has_repeated = True
                break
        if has_repeated:
            break
            
    if has_repeated:
        score -= 15
        
    # Check if it's a common password
    if is_common_password(password):
        score = min(score, 25)  # Cap at 25 if it's common
        
    # Ensure score is within 0-100 range
    score = max(0, min(score, 100))
    
    return score

def generate_random_password(length: int = 12) -> str:
    """
    Generates a secure random password that meets complexity requirements.
    
    Args:
        length: The desired length of the password (default: 12)
        
    Returns:
        Secure random password
    """
    # Ensure minimum length
    if length < PASSWORD_MIN_LENGTH:
        length = PASSWORD_MIN_LENGTH
    
    # Define character sets
    lowercase = string.ascii_lowercase
    uppercase = string.ascii_uppercase
    digits = string.digits
    special = "!@#$%^&*(),.?\":{}|<>"
    
    # Ensure at least one character from each set
    password = [
        secrets.choice(lowercase),
        secrets.choice(uppercase),
        secrets.choice(digits),
        secrets.choice(special)
    ]
    
    # Fill the rest with random characters from all sets
    all_chars = lowercase + uppercase + digits + special
    password.extend(secrets.choice(all_chars) for _ in range(length - 4))
    
    # Shuffle the characters for randomness
    secrets.SystemRandom().shuffle(password)
    
    return ''.join(password)

def secure_compare(a: str, b: str) -> bool:
    """
    Performs a constant-time comparison of two strings to prevent timing attacks.
    
    Args:
        a: First string
        b: Second string
        
    Returns:
        True if strings are equal, False otherwise
    """
    # Check if lengths are equal (this is not constant time, but necessary)
    if len(a) != len(b):
        return False
        
    # Use secrets.compare_digest for constant-time comparison
    return secrets.compare_digest(a.encode(), b.encode())

def hash_data(data: str) -> str:
    """
    Creates a secure hash of data using SHA-256.
    
    Args:
        data: Data to hash
        
    Returns:
        Hexadecimal hash string
    """
    # Convert to bytes if it's a string
    if isinstance(data, str):
        data = data.encode()
        
    # Create SHA-256 hash
    hash_obj = hashlib.sha256(data)
    
    # Return hexadecimal digest
    return hash_obj.hexdigest()

def encode_base64(data: str) -> str:
    """
    Encodes data as base64.
    
    Args:
        data: Data to encode
        
    Returns:
        Base64 encoded string
    """
    # Convert to bytes if it's a string
    if isinstance(data, str):
        data = data.encode()
        
    # Encode as base64
    encoded = base64.b64encode(data)
    
    # Return as string
    return encoded.decode()

def decode_base64(data: str) -> str:
    """
    Decodes base64 data.
    
    Args:
        data: Base64 data to decode
        
    Returns:
        Decoded string
    """
    # Decode base64
    decoded = base64.b64decode(data)
    
    # Return as string
    return decoded.decode()

def sanitize_input(input_string: str) -> str:
    """
    Sanitizes user input to prevent injection attacks.
    
    Args:
        input_string: The input string to sanitize
        
    Returns:
        Sanitized string
    """
    # Remove potentially dangerous patterns
    sanitized = re.sub(r'[<>]', '', input_string)
    
    # Escape special characters
    sanitized = sanitized.replace('&', '&amp;')
    sanitized = sanitized.replace('"', '&quot;')
    sanitized = sanitized.replace("'", '&#x27;')
    
    return sanitized

class PasswordValidator:
    """
    Class that validates passwords against security requirements.
    """
    
    def __init__(
        self,
        min_length: int = PASSWORD_MIN_LENGTH,
        require_lowercase: bool = True,
        require_uppercase: bool = True,
        require_digit: bool = True,
        require_special: bool = True,
        check_common_passwords: bool = True
    ):
        """
        Initializes the password validator with configurable settings.
        
        Args:
            min_length: Minimum password length
            require_lowercase: Require at least one lowercase letter
            require_uppercase: Require at least one uppercase letter
            require_digit: Require at least one digit
            require_special: Require at least one special character
            check_common_passwords: Check if password is in common passwords list
        """
        self.min_length = min_length
        self.require_lowercase = require_lowercase
        self.require_uppercase = require_uppercase
        self.require_digit = require_digit
        self.require_special = require_special
        self.check_common_passwords = check_common_passwords
    
    def validate(self, password: str) -> bool:
        """
        Validates a password against the configured requirements.
        
        Args:
            password: The password to validate
            
        Returns:
            True if password meets requirements, False otherwise
        """
        # Check minimum length
        if len(password) < self.min_length:
            return False
            
        # Check for character types if required
        if self.require_lowercase and not re.search(PASSWORD_REGEX_LOWERCASE, password):
            return False
            
        if self.require_uppercase and not re.search(PASSWORD_REGEX_UPPERCASE, password):
            return False
            
        if self.require_digit and not re.search(PASSWORD_REGEX_DIGIT, password):
            return False
            
        if self.require_special and not re.search(PASSWORD_REGEX_SPECIAL, password):
            return False
            
        # Check if it's a common password
        if self.check_common_passwords and is_common_password(password):
            return False
            
        return True
    
    def get_validation_errors(self, password: str) -> List[str]:
        """
        Returns a list of validation errors for a password.
        
        Args:
            password: The password to validate
            
        Returns:
            List of validation error messages
        """
        errors = []
        
        # Check minimum length
        if len(password) < self.min_length:
            errors.append(f"Password must be at least {self.min_length} characters long")
            
        # Check for character types if required
        if self.require_lowercase and not re.search(PASSWORD_REGEX_LOWERCASE, password):
            errors.append("Password must contain at least one lowercase letter")
            
        if self.require_uppercase and not re.search(PASSWORD_REGEX_UPPERCASE, password):
            errors.append("Password must contain at least one uppercase letter")
            
        if self.require_digit and not re.search(PASSWORD_REGEX_DIGIT, password):
            errors.append("Password must contain at least one digit")
            
        if self.require_special and not re.search(PASSWORD_REGEX_SPECIAL, password):
            errors.append("Password must contain at least one special character")
            
        # Check if it's a common password
        if self.check_common_passwords and is_common_password(password):
            errors.append("Password is too common and easily guessable")
            
        return errors
    
    def calculate_strength(self, password: str) -> int:
        """
        Calculates the strength of a password.
        
        Args:
            password: The password to evaluate
            
        Returns:
            Password strength score (0-100)
        """
        return calculate_password_strength(password)

class TokenGenerator:
    """
    Class that generates secure tokens for various purposes.
    """
    
    def generate_token(self, length: int) -> str:
        """
        Generates a secure random token.
        
        Args:
            length: The desired length of the token
            
        Returns:
            Secure random token
        """
        return generate_secure_token(length)
    
    def generate_verification_code(self, length: int = 6) -> str:
        """
        Generates a numeric verification code.
        
        Args:
            length: The desired length of the code (default: 6)
            
        Returns:
            Numeric verification code
        """
        # Generate a random digits-only code
        return ''.join(secrets.choice(string.digits) for _ in range(length))
    
    def generate_api_key(self) -> str:
        """
        Generates a secure API key.
        
        Args:
            None
            
        Returns:
            Secure API key
        """
        # Generate a random token with prefix
        token = secrets.token_hex(16)
        return f"api_{token}"