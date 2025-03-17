"""
Core JWT module that implements JSON Web Token generation, validation, and management
for the Molecular Data Management and CRO Integration Platform.

This module provides secure token-based authentication without external dependencies,
supporting the platform's requirement for fully local deployment.
"""

import jwt  # PyJWT 2.7+
from datetime import datetime, timedelta  # standard library
import uuid  # standard library
from typing import Dict, Any, Optional, Union  # standard library

from .config import get_settings
from ..exceptions import AuthenticationException
from ..schemas.token import TokenData
from ..constants import TOKEN_EXPIRY_MINUTES, REFRESH_TOKEN_EXPIRY_DAYS, ALGORITHM


def create_access_token(data: Dict[str, Any], expires_delta: Optional[timedelta] = None) -> str:
    """
    Creates a JWT access token with the provided data and expiration time.
    
    Args:
        data: Dictionary containing claims to be encoded in the token
        expires_delta: Optional custom expiration time, defaults to TOKEN_EXPIRY_MINUTES
        
    Returns:
        Encoded JWT access token
    """
    settings = get_settings()
    to_encode = data.copy()
    
    # Set token type
    to_encode.update({"token_type": "access"})
    
    # Add a unique token ID
    jti = str(uuid.uuid4())
    to_encode.update({"jti": jti})
    
    # Calculate expiration time
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(minutes=TOKEN_EXPIRY_MINUTES)
    
    # Add exp and iat claims
    to_encode.update({"exp": expire})
    to_encode.update({"iat": datetime.utcnow()})
    
    # Encode the token using RS256 algorithm and private key
    encoded_jwt = jwt.encode(
        to_encode, 
        settings.get_jwt_private_key, 
        algorithm=ALGORITHM
    )
    
    return encoded_jwt


def create_refresh_token(data: Dict[str, Any], expires_delta: Optional[timedelta] = None) -> str:
    """
    Creates a JWT refresh token with the provided data and longer expiration time.
    
    Args:
        data: Dictionary containing claims to be encoded in the token
        expires_delta: Optional custom expiration time, defaults to REFRESH_TOKEN_EXPIRY_DAYS
        
    Returns:
        Encoded JWT refresh token
    """
    settings = get_settings()
    to_encode = data.copy()
    
    # Set token type
    to_encode.update({"token_type": "refresh"})
    
    # Add a unique token ID
    jti = str(uuid.uuid4())
    to_encode.update({"jti": jti})
    
    # Calculate expiration time
    if expires_delta:
        expire = datetime.utcnow() + expires_delta
    else:
        expire = datetime.utcnow() + timedelta(days=REFRESH_TOKEN_EXPIRY_DAYS)
    
    # Add exp and iat claims
    to_encode.update({"exp": expire})
    to_encode.update({"iat": datetime.utcnow()})
    
    # Encode the token using RS256 algorithm and private key
    encoded_jwt = jwt.encode(
        to_encode, 
        settings.get_jwt_private_key, 
        algorithm=ALGORITHM
    )
    
    return encoded_jwt


def verify_token(token: str) -> Dict[str, Any]:
    """
    Verifies and decodes a JWT token, returning the payload if valid.
    
    Args:
        token: JWT token string to verify
        
    Returns:
        Decoded token payload
        
    Raises:
        AuthenticationException: If token is invalid, expired, or has signature issues
    """
    settings = get_settings()
    
    try:
        # Decode the token using the public key
        payload = jwt.decode(
            token, 
            settings.get_jwt_public_key, 
            algorithms=[ALGORITHM]
        )
        return payload
    except jwt.ExpiredSignatureError:
        raise AuthenticationException("Token expired")
    except jwt.InvalidSignatureError:
        raise AuthenticationException("Invalid token")
    except jwt.PyJWTError as e:
        raise AuthenticationException(f"Invalid token: {str(e)}")


def decode_token(token: str) -> Dict[str, Any]:
    """
    Decodes a JWT token without verification, for debugging purposes.
    
    Args:
        token: JWT token string to decode
        
    Returns:
        Decoded token payload without verification
    """
    decoded = jwt.decode_complete(token, options={"verify_signature": False})
    return decoded["payload"]


def get_token_data(token: str) -> TokenData:
    """
    Extracts and validates token data using the TokenData schema.
    
    Args:
        token: JWT token string
        
    Returns:
        Validated token data
        
    Raises:
        AuthenticationException: If token is invalid
    """
    payload = verify_token(token)
    token_data = TokenData(**payload)
    return token_data


def create_token_for_user(
    user_id: int, 
    email: str, 
    role: str, 
    expires_delta: Optional[timedelta] = None
) -> str:
    """
    Creates an access token for a user with standard claims.
    
    Args:
        user_id: User ID to encode in the token
        email: User email to encode in the token
        role: User role to encode in the token
        expires_delta: Optional custom expiration time
        
    Returns:
        Encoded JWT access token
    """
    payload = {
        "user_id": user_id,
        "email": email,
        "role": role
    }
    return create_access_token(payload, expires_delta)


def create_refresh_token_for_user(
    user_id: int, 
    email: str, 
    role: str, 
    expires_delta: Optional[timedelta] = None
) -> str:
    """
    Creates a refresh token for a user with standard claims.
    
    Args:
        user_id: User ID to encode in the token
        email: User email to encode in the token
        role: User role to encode in the token
        expires_delta: Optional custom expiration time
        
    Returns:
        Encoded JWT refresh token
    """
    payload = {
        "user_id": user_id,
        "email": email,
        "role": role
    }
    return create_refresh_token(payload, expires_delta)


def get_token_expiration(payload: Dict[str, Any]) -> int:
    """
    Calculates token expiration time in seconds from now.
    
    Args:
        payload: Decoded token payload containing expiration time
        
    Returns:
        Seconds until token expiration (or 0 if expired/no expiration)
    """
    exp = payload.get("exp")
    if not exp:
        return 0
    
    now = datetime.utcnow().timestamp()
    return max(int(exp - now), 0)


def is_token_valid(token: str) -> bool:
    """
    Checks if a token is valid without raising exceptions.
    
    Args:
        token: JWT token string to validate
        
    Returns:
        True if token is valid, False otherwise
    """
    try:
        verify_token(token)
        return True
    except Exception:
        return False


def create_token_with_custom_claims(
    claims: Dict[str, Any], 
    token_type: str, 
    expires_delta: Optional[timedelta] = None
) -> str:
    """
    Creates a JWT token with custom claims and expiration.
    
    Args:
        claims: Dictionary containing custom claims
        token_type: Token type ("access" or "refresh")
        expires_delta: Optional custom expiration time
        
    Returns:
        Encoded JWT token
        
    Raises:
        ValueError: If token_type is not "access" or "refresh"
    """
    payload = claims.copy()
    payload.update({"token_type": token_type})
    
    if token_type == "access":
        return create_access_token(payload, expires_delta)
    elif token_type == "refresh":
        return create_refresh_token(payload, expires_delta)
    else:
        raise ValueError("Invalid token type. Must be 'access' or 'refresh'")