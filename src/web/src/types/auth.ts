/**
 * Authentication-related TypeScript definitions for the Molecular Data Management and CRO Integration Platform.
 * These types support the JWT-based authentication system that operates without external dependencies.
 * @packageDocumentation
 */

import { UserRole, UserStatus } from './user';

/**
 * Interface representing a user in the system.
 * Contains all standard user information required for authentication.
 */
export interface User {
  /** Unique identifier for the user */
  id: number;
  /** User's email address */
  email: string;
  /** User's role in the system */
  role: UserRole;
  /** Current account status */
  status: UserStatus;
  /** Whether the user account is currently active */
  is_active: boolean;
  /** Whether the user's email has been verified */
  email_verified: boolean;
  /** Timestamp when the user account was created */
  created_at: string;
  /** Timestamp of the user's last successful login */
  last_login: string | null;
}

/**
 * Interface for user login credentials.
 * Contains the required fields for authenticating a user.
 */
export interface LoginCredentials {
  /** User's email address */
  email: string;
  /** User's password */
  password: string;
}

/**
 * Interface for the response from a successful login attempt.
 * Contains authentication tokens and user information.
 */
export interface LoginResponse {
  /** JWT access token for API authorization */
  access_token: string;
  /** JWT refresh token for obtaining new access tokens */
  refresh_token: string;
  /** The authenticated user's information */
  user: User;
}

/**
 * Interface for user registration data.
 * Contains all fields needed to register a new user.
 */
export interface RegisterData {
  /** New user's email address */
  email: string;
  /** New user's password */
  password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
  /** User's role in the system */
  role: UserRole;
}

/**
 * Interface for the response from a user registration attempt.
 * Indicates success/failure and includes created user information.
 */
export interface RegistrationResponse {
  /** Whether the registration was successful */
  success: boolean;
  /** Message describing the result of the registration */
  message: string;
  /** The newly registered user's information (if successful) */
  user?: User;
}

/**
 * Interface for a refresh token request.
 * Used to obtain a new access token without re-authentication.
 */
export interface RefreshTokenRequest {
  /** The JWT refresh token */
  refresh_token: string;
}

/**
 * Interface for an email verification request.
 * Used to verify a user's email address.
 */
export interface EmailVerificationRequest {
  /** The verification token sent to the user's email */
  token: string;
}

/**
 * Interface for a password reset request.
 * Used to initiate the password reset process.
 */
export interface PasswordResetRequest {
  /** The email address of the account to reset */
  email: string;
}

/**
 * Interface for a password reset confirmation request.
 * Used to complete the password reset process.
 */
export interface PasswordResetConfirmRequest {
  /** The password reset token sent to the user's email */
  token: string;
  /** The new password to set */
  new_password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
}

/**
 * Interface for a password change request.
 * Used when a user wants to change their password.
 */
export interface ChangePasswordRequest {
  /** The user's current password for verification */
  current_password: string;
  /** The new password to set */
  new_password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
}

/**
 * Interface for decoded JWT token data.
 * Contains the standard JWT claims and user-specific data.
 */
export interface TokenData {
  /** Subject identifier (usually user ID) */
  sub: string;
  /** Expiration timestamp */
  exp: number;
  /** Issued at timestamp */
  iat: number;
  /** User ID */
  user_id: number;
  /** User email */
  email: string;
  /** User role */
  role: UserRole;
}

/**
 * Interface for the authentication state in the Redux store.
 * Tracks the current authentication status and user information.
 */
export interface AuthState {
  /** Whether the user is currently authenticated */
  isAuthenticated: boolean;
  /** The authenticated user's information */
  user: User | null;
  /** The current JWT access token */
  accessToken: string | null;
  /** The current JWT refresh token */
  refreshToken: string | null;
  /** Whether authentication operations are in progress */
  loading: boolean;
  /** Any error that occurred during authentication */
  error: string | null;
}

/**
 * Type for permission strings in the system.
 * Represents specific actions users can perform.
 */
export type Permission = string;

/**
 * Interface for generic message responses from the API.
 * Used for operations that return simple success/failure messages.
 */
export interface MessageResponse {
  /** Whether the operation was successful */
  success: boolean;
  /** Message describing the result of the operation */
  message: string;
}