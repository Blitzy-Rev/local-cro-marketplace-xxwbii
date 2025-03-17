/**
 * API module that provides authentication-related API functions for the
 * Molecular Data Management and CRO Integration Platform frontend.
 * This module handles user login, registration, token refresh, password management,
 * and logout operations without external dependencies.
 * 
 * @module api/auth
 */

import { apiClient, handleApiError } from './client'; // axios ^1.4.0
import {
  LoginCredentials,
  LoginResponse,
  RegisterData,
  RegistrationResponse,
  RefreshTokenRequest,
  EmailVerificationRequest,
  PasswordResetRequest,
  PasswordResetConfirmRequest,
  ChangePasswordRequest,
  MessageResponse,
  User
} from '../types/auth';

/**
 * API endpoints for authentication operations
 */
const AUTH_API_ENDPOINTS = {
  LOGIN: '/auth/login',
  REGISTER: '/auth/register',
  REFRESH_TOKEN: '/auth/refresh-token',
  LOGOUT: '/auth/logout',
  VERIFY_EMAIL: '/auth/verify-email',
  RESET_PASSWORD: '/auth/reset-password',
  RESET_PASSWORD_CONFIRM: '/auth/reset-password-confirm',
  CHANGE_PASSWORD: '/auth/change-password'
};

/**
 * Authenticates a user with email and password
 * 
 * @param credentials User login credentials
 * @returns Promise resolving to login response with tokens and user data
 */
export const login = async (credentials: LoginCredentials): Promise<LoginResponse> => {
  try {
    const response = await apiClient.post<LoginResponse>(
      AUTH_API_ENDPOINTS.LOGIN,
      credentials
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Registers a new user with the system
 * 
 * @param registrationData User registration data
 * @returns Promise resolving to registration response with user data
 */
export const register = async (registrationData: RegisterData): Promise<RegistrationResponse> => {
  try {
    const response = await apiClient.post<RegistrationResponse>(
      AUTH_API_ENDPOINTS.REGISTER,
      registrationData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Refreshes the authentication token using a refresh token
 * 
 * @param refreshToken The refresh token
 * @returns Promise resolving to login response with new tokens
 */
export const refreshToken = async (refreshToken: string): Promise<LoginResponse> => {
  try {
    const request: RefreshTokenRequest = { refresh_token: refreshToken };
    const response = await apiClient.post<LoginResponse>(
      AUTH_API_ENDPOINTS.REFRESH_TOKEN,
      request
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Logs out a user by invalidating their refresh token
 * 
 * @param refreshToken The refresh token to invalidate
 * @returns Promise resolving to message response
 */
export const logout = async (refreshToken: string): Promise<MessageResponse> => {
  try {
    const request: RefreshTokenRequest = { refresh_token: refreshToken };
    const response = await apiClient.post<MessageResponse>(
      AUTH_API_ENDPOINTS.LOGOUT,
      request
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Verifies a user's email address using a verification token
 * 
 * @param token Email verification token
 * @returns Promise resolving to message response
 */
export const verifyEmail = async (token: string): Promise<MessageResponse> => {
  try {
    const request: EmailVerificationRequest = { token };
    const response = await apiClient.post<MessageResponse>(
      AUTH_API_ENDPOINTS.VERIFY_EMAIL,
      request
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Initiates a password reset process for a user
 * 
 * @param email The email address of the user
 * @returns Promise resolving to message response
 */
export const resetPassword = async (email: string): Promise<MessageResponse> => {
  try {
    const request: PasswordResetRequest = { email };
    const response = await apiClient.post<MessageResponse>(
      AUTH_API_ENDPOINTS.RESET_PASSWORD,
      request
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Completes a password reset process with a new password
 * 
 * @param resetData Password reset confirmation data
 * @returns Promise resolving to message response
 */
export const resetPasswordConfirm = async (
  resetData: PasswordResetConfirmRequest
): Promise<MessageResponse> => {
  try {
    const response = await apiClient.post<MessageResponse>(
      AUTH_API_ENDPOINTS.RESET_PASSWORD_CONFIRM,
      resetData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Changes a user's password while authenticated
 * 
 * @param passwordData Password change data
 * @returns Promise resolving to message response
 */
export const changePassword = async (
  passwordData: ChangePasswordRequest
): Promise<MessageResponse> => {
  try {
    const response = await apiClient.post<MessageResponse>(
      AUTH_API_ENDPOINTS.CHANGE_PASSWORD,
      passwordData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};