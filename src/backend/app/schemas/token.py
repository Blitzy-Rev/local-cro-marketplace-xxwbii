from typing import Optional, Dict, Any, List
from pydantic import BaseModel, Field

class Token(BaseModel):
    """
    Schema for a single JWT token.
    
    This model represents a single token response that would be returned
    when a token is generated.
    """
    token: str = Field(..., description="The JWT token string")
    token_type: str = Field(..., description="The type of token (e.g., 'bearer')")


class TokenData(BaseModel):
    """
    Schema for decoded JWT token data.
    
    This model represents the data extracted from a decoded JWT token
    and is used for internal token validation and processing.
    """
    user_id: Optional[int] = Field(None, description="ID of the user associated with the token")
    email: Optional[str] = Field(None, description="Email of the user associated with the token")
    role: Optional[str] = Field(None, description="Role of the user (pharma, cro, admin)")
    token_type: Optional[str] = Field(None, description="Type of token (access or refresh)")
    jti: Optional[str] = Field(None, description="Unique identifier for the token (JWT ID)")
    exp: Optional[int] = Field(None, description="Expiration timestamp (Unix time)")
    iat: Optional[int] = Field(None, description="Issued at timestamp (Unix time)")
    permissions: Optional[Dict[str, Any]] = Field(None, description="User permissions associated with the token")


class TokenResponse(BaseModel):
    """
    Schema for token response containing access and refresh tokens.
    
    This model represents the response returned to clients during
    authentication, containing both access and refresh tokens.
    """
    access_token: str = Field(..., description="JWT access token for short-term authentication")
    refresh_token: str = Field(..., description="JWT refresh token for obtaining new access tokens")
    token_type: str = Field(default="bearer", description="The type of token (typically 'bearer')")


class TokenPayload(BaseModel):
    """
    Schema for JWT token payload.
    
    This model represents the payload structure that will be encoded
    into a JWT token. It includes standard JWT claims and custom claims.
    """
    user_id: Optional[int] = Field(None, description="ID of the user associated with the token")
    email: Optional[str] = Field(None, description="Email of the user associated with the token")
    role: Optional[str] = Field(None, description="Role of the user (pharma, cro, admin)")
    token_type: Optional[str] = Field(None, description="Type of token (access or refresh)")
    jti: Optional[str] = Field(None, description="Unique identifier for the token (JWT ID)")
    exp: Optional[int] = Field(None, description="Expiration timestamp (Unix time)")
    iat: Optional[int] = Field(None, description="Issued at timestamp (Unix time)")
    custom_claims: Optional[Dict[str, Any]] = Field(None, description="Additional custom claims to include in the token")