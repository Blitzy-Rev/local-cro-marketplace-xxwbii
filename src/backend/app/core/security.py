"""
Core security module that implements password hashing, verification, and other security utilities
for the Molecular Data Management and CRO Integration Platform. This module provides essential
security functions without external dependencies, supporting the platform's requirement for fully
local deployment.
"""

import re  # standard library
import secrets  # standard library
import string  # standard library
import os  # standard library
from typing import Optional, List, Dict, Any  # standard library

from passlib.context import CryptContext  # version 1.7.4+

from .config import get_settings
from ..exceptions import AuthenticationException
from ..utils.security_utils import validate_password_strength, get_password_validation_errors
from ..constants import PASSWORD_MIN_LENGTH

# Configure password hashing context with bcrypt and appropriate security settings
pwd_context = CryptContext(schemes=["bcrypt"], deprecated="auto", bcrypt__rounds=12)

# Maximum number of previous passwords to check for history
PASSWORD_HISTORY_SIZE = a = 5

def verify_password(plain_password: str, hashed_password: str) -> bool:
    """
    Verifies a plain password against a hashed password.
    
    Args:
        plain_password: The plain password to verify
        hashed_password: The hashed password to verify against
        
    Returns:
        True if password matches hash, False otherwise
    """
    return pwd_context.verify(plain_password, hashed_password)

def get_password_hash(password: str) -> str:
    """
    Generates a secure hash of a password using bcrypt.
    
    Args:
        password: The password to hash
        
    Returns:
        Hashed password string
    """
    return pwd_context.hash(password)

def validate_password(password: str) -> bool:
    """
    Validates a password against security requirements.
    
    Args:
        password: The password to validate
        
    Returns:
        True if password meets requirements, False otherwise
    """
    return validate_password_strength(password)

def get_validation_errors(password: str) -> List[str]:
    """
    Gets a list of validation errors for a password.
    
    Args:
        password: The password to validate
        
    Returns:
        List of validation error messages
    """
    return get_password_validation_errors(password)

def check_password_history(password: str, password_history: List[str]) -> bool:
    """
    Checks if a password has been used before in the password history.
    
    Args:
        password: The password to check
        password_history: List of previous password hashes
        
    Returns:
        True if password is in history, False otherwise
    """
    for old_password_hash in password_history:
        if verify_password(password, old_password_hash):
            return True
    return False

def update_password_history(new_password_hash: str, password_history: List[str]) -> List[str]:
    """
    Updates password history with a new password hash.
    
    Args:
        new_password_hash: The new password hash to add
        password_history: Current password history list
        
    Returns:
        Updated password history list
    """
    # Add new password hash to the beginning
    updated_history = [new_password_hash] + password_history
    
    # Trim to maximum size
    if len(updated_history) > PASSWORD_HISTORY_SIZE:
        updated_history = updated_history[:PASSWORD_HISTORY_SIZE]
    
    return updated_history

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
    
    # Fill remaining length with random characters
    all_chars = lowercase + uppercase + digits + special
    for _ in range(length - 4):
        password.append(secrets.choice(all_chars))
    
    # Shuffle the password characters
    secrets.SystemRandom().shuffle(password)
    
    # Convert to string
    return ''.join(password)

def is_password_compromised(password: str) -> bool:
    """
    Checks if a password is in a list of known compromised passwords.
    
    Args:
        password: The password to check
        
    Returns:
        True if password is compromised, False otherwise
    """
    try:
        settings = get_settings()
        # Get the common passwords file path
        common_passwords_path = settings.JWT_PRIVATE_KEY_PATH.rsplit("/", 1)[0] + "/common_passwords.txt"
        
        # If common passwords file doesn't exist, consider no passwords compromised
        if not os.path.exists(common_passwords_path):
            return False
            
        with open(common_passwords_path, "r") as f:
            common_passwords = [line.strip() for line in f]
            
        return password.lower() in common_passwords
    except Exception:
        # If there's any error reading the file, fail safe (assume not compromised)
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
    if re.search(r'[a-z]', password):
        score += 10
    if re.search(r'[A-Z]', password):
        score += 10
    if re.search(r'\d', password):
        score += 5
    if re.search(r'[^a-zA-Z0-9]', password):
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
    if is_password_compromised(password):
        score = min(score, 25)  # Cap at 25 if it's common
        
    # Ensure score is within 0-100 range
    score = max(0, min(score, 100))
    
    return score

class PasswordPolicy:
    """
    Class that enforces password policies and validates passwords.
    """
    
    def __init__(
        self,
        min_length: int = PASSWORD_MIN_LENGTH,
        require_lowercase: bool = True,
        require_uppercase: bool = True,
        require_digit: bool = True,
        require_special: bool = True,
        history_size: int = PASSWORD_HISTORY_SIZE
    ):
        """
        Initializes the password policy with configurable settings.
        
        Args:
            min_length: Minimum password length
            require_lowercase: Require at least one lowercase letter
            require_uppercase: Require at least one uppercase letter
            require_digit: Require at least one digit
            require_special: Require at least one special character
            history_size: Number of previous passwords to check
        """
        self.min_length = min_length
        self.require_lowercase = require_lowercase
        self.require_uppercase = require_uppercase
        self.require_digit = require_digit
        self.require_special = require_special
        self.history_size = history_size
    
    def validate(self, password: str) -> bool:
        """
        Validates a password against the policy.
        
        Args:
            password: The password to validate
            
        Returns:
            True if password meets policy requirements, False otherwise
        """
        # Check length
        if len(password) < self.min_length:
            return False
        
        # Check character requirements
        if self.require_lowercase and not re.search(r'[a-z]', password):
            return False
        
        if self.require_uppercase and not re.search(r'[A-Z]', password):
            return False
        
        if self.require_digit and not re.search(r'\d', password):
            return False
        
        if self.require_special and not re.search(r'[^a-zA-Z0-9]', password):
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
        
        # Check length
        if len(password) < self.min_length:
            errors.append(f"Password must be at least {self.min_length} characters long")
        
        # Check character requirements
        if self.require_lowercase and not re.search(r'[a-z]', password):
            errors.append("Password must contain at least one lowercase letter")
        
        if self.require_uppercase and not re.search(r'[A-Z]', password):
            errors.append("Password must contain at least one uppercase letter")
        
        if self.require_digit and not re.search(r'\d', password):
            errors.append("Password must contain at least one digit")
        
        if self.require_special and not re.search(r'[^a-zA-Z0-9]', password):
            errors.append("Password must contain at least one special character")
        
        return errors
    
    def check_history(self, password: str, password_history: List[str]) -> bool:
        """
        Checks if a password is in the password history.
        
        Args:
            password: The password to check
            password_history: List of previous password hashes
            
        Returns:
            True if password is in history, False otherwise
        """
        return check_password_history(password, password_history)
    
    def update_history(self, new_password_hash: str, password_history: List[str]) -> List[str]:
        """
        Updates password history with a new password hash.
        
        Args:
            new_password_hash: The new password hash to add
            password_history: Current password history list
            
        Returns:
            Updated password history list
        """
        return update_password_history(new_password_hash, password_history)