from typing import Dict, Any
from fastapi import APIRouter, Depends, HTTPException, status, Request
from sqlalchemy.orm import Session

from ..deps import get_db_session, get_current_user_from_token, get_current_active_user_from_token
from ...services.auth_service import AuthService
from ..schemas.auth import (
    LoginRequest, LoginResponse, RefreshTokenRequest, RegistrationRequest,
    RegistrationResponse, EmailVerificationRequest, PasswordResetRequest,
    PasswordResetConfirmRequest, ChangePasswordRequest, MessageResponse
)
from ...schemas.token import TokenResponse
from ...schemas.user import UserResponse
from ...exceptions import AuthenticationException
from ...logging_config import logger

router = APIRouter(prefix="/auth", tags=["authentication"])


@router.post("/login", response_model=LoginResponse, status_code=status.HTTP_200_OK)
async def login(
    login_data: LoginRequest,
    db: Session = Depends(get_db_session),
    request: Request = None
) -> LoginResponse:
    """
    Authenticate a user and return access and refresh tokens
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Authenticate user with credentials
        user, access_token, refresh_token = auth_service.authenticate(
            login_data.email, login_data.password
        )
        
        # Create response
        tokens = TokenResponse(
            access_token=access_token,
            refresh_token=refresh_token,
            token_type="bearer"
        )
        
        user_response = UserResponse(
            id=user.id,
            email=user.email,
            role=user.role.name
        )
        
        # Log successful login
        logger.info(f"User {user.email} successfully logged in", {"user_id": user.id})
        
        return LoginResponse(tokens=tokens, user=user_response)
    
    except AuthenticationException as e:
        # Log failed login attempt
        logger.warning(f"Failed login attempt for {login_data.email}: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=str(e)
        )


@router.post("/refresh-token", status_code=status.HTTP_200_OK)
async def refresh_token(
    refresh_data: RefreshTokenRequest,
    db: Session = Depends(get_db_session)
) -> Dict[str, Any]:
    """
    Generate a new access token using a valid refresh token
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Refresh token
        access_token = auth_service.refresh_token(refresh_data.refresh_token)
        
        # Log successful token refresh
        logger.info("Token successfully refreshed")
        
        return {"access_token": access_token, "token_type": "bearer"}
    
    except AuthenticationException as e:
        # Log failed token refresh
        logger.warning(f"Failed token refresh: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail=str(e)
        )


@router.post("/logout", response_model=MessageResponse, status_code=status.HTTP_200_OK)
async def logout(
    refresh_data: RefreshTokenRequest,
    db: Session = Depends(get_db_session),
    current_user: Dict[str, Any] = Depends(get_current_user_from_token)
) -> MessageResponse:
    """
    Invalidate a refresh token on user logout
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Invalidate refresh token
        auth_service.invalidate_refresh_token(refresh_data.refresh_token)
        
        # Log successful logout
        logger.info(f"User {current_user['email']} successfully logged out", 
                  {"user_id": current_user["user_id"]})
        
        return MessageResponse(success=True, message="Successfully logged out")
    
    except AuthenticationException as e:
        # Log failed logout
        logger.warning(f"Failed logout: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.post("/register", response_model=RegistrationResponse, status_code=status.HTTP_201_CREATED)
async def register(
    registration_data: RegistrationRequest,
    db: Session = Depends(get_db_session),
    request: Request = None
) -> RegistrationResponse:
    """
    Register a new user in the system
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Validate password match
        if registration_data.password != registration_data.confirm_password:
            raise AuthenticationException("Passwords do not match")
        
        # Get base URL for verification email
        base_url = str(request.base_url).rstrip("/")
        
        # Register user
        new_user = auth_service.register_user(
            db,
            registration_data,
            base_url
        )
        
        # Send verification email
        auth_service.send_verification_email(
            new_user,
            base_url
        )
        
        # Log successful registration
        logger.info(f"User {new_user.email} successfully registered", {"user_id": new_user.id})
        
        # Create response
        user_response = UserResponse(
            id=new_user.id,
            email=new_user.email,
            role=new_user.role.name
        )
        
        return RegistrationResponse(
            success=True,
            user=user_response,
            message="Registration successful. Please check your email for verification."
        )
    
    except AuthenticationException as e:
        # Log failed registration
        logger.warning(f"Failed registration for {registration_data.email}: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.post("/verify-email", response_model=MessageResponse, status_code=status.HTTP_200_OK)
async def verify_email(
    verification_data: EmailVerificationRequest,
    db: Session = Depends(get_db_session)
) -> MessageResponse:
    """
    Verify a user's email using a verification token
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Verify email
        verified_user = auth_service.verify_email(
            db,
            verification_data.token
        )
        
        # Log successful verification
        logger.info(f"Email successfully verified for user {verified_user.email}", 
                  {"user_id": verified_user.id})
        
        return MessageResponse(
            success=True,
            message="Email successfully verified. You can now log in."
        )
    
    except AuthenticationException as e:
        # Log failed verification
        logger.warning(f"Failed email verification: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.post("/password-reset/request", response_model=MessageResponse, status_code=status.HTTP_200_OK)
async def request_password_reset(
    reset_data: PasswordResetRequest,
    db: Session = Depends(get_db_session),
    request: Request = None
) -> MessageResponse:
    """
    Request a password reset for a user
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Get base URL for reset email
        base_url = str(request.base_url).rstrip("/")
        
        # Request password reset
        auth_service.request_password_reset(
            db,
            reset_data.email,
            base_url
        )
        
        # Log password reset request
        logger.info(f"Password reset requested for {reset_data.email}")
        
        # Always return success (even if email doesn't exist) for security
        return MessageResponse(
            success=True,
            message="If your email is registered, you will receive a password reset link."
        )
    
    except Exception as e:
        # Log error but don't expose it to the client
        logger.error(f"Error in password reset request for {reset_data.email}: {str(e)}")
        
        # Still return success for security
        return MessageResponse(
            success=True,
            message="If your email is registered, you will receive a password reset link."
        )


@router.post("/password-reset/confirm", response_model=MessageResponse, status_code=status.HTTP_200_OK)
async def confirm_password_reset(
    reset_data: PasswordResetConfirmRequest,
    db: Session = Depends(get_db_session)
) -> MessageResponse:
    """
    Reset a user's password using a reset token
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Validate password match
        if reset_data.new_password != reset_data.confirm_password:
            raise AuthenticationException("Passwords do not match")
        
        # Confirm password reset
        auth_service.confirm_password_reset(
            db,
            reset_data.token,
            reset_data.new_password
        )
        
        # Log successful password reset
        logger.info("Password successfully reset")
        
        return MessageResponse(
            success=True,
            message="Password has been reset successfully. You can now log in with your new password."
        )
    
    except AuthenticationException as e:
        # Log failed password reset
        logger.warning(f"Failed password reset: {str(e)}")
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.post("/change-password", response_model=MessageResponse, status_code=status.HTTP_200_OK)
async def change_password(
    password_data: ChangePasswordRequest,
    db: Session = Depends(get_db_session),
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
) -> MessageResponse:
    """
    Change a user's password
    """
    try:
        # Create AuthService instance
        auth_service = AuthService(db)
        
        # Validate password match
        if password_data.new_password != password_data.confirm_password:
            raise AuthenticationException("Passwords do not match")
        
        # Change password
        auth_service.change_password(
            db,
            current_user["user_id"],
            password_data.current_password,
            password_data.new_password
        )
        
        # Log successful password change
        logger.info(f"Password successfully changed for user {current_user['email']}", 
                  {"user_id": current_user["user_id"]})
        
        return MessageResponse(
            success=True,
            message="Password has been changed successfully."
        )
    
    except AuthenticationException as e:
        # Log failed password change
        logger.warning(f"Failed password change for user {current_user['email']}: {str(e)}", 
                     {"user_id": current_user["user_id"]})
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=str(e)
        )


@router.get("/me", response_model=UserResponse, status_code=status.HTTP_200_OK)
async def get_current_user(
    current_user: Dict[str, Any] = Depends(get_current_active_user_from_token)
) -> UserResponse:
    """
    Get the current authenticated user's information
    """
    return UserResponse(
        id=current_user["user_id"],
        email=current_user["email"],
        role=current_user["role"]
    )