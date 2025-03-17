/**
 * Custom React hook that provides authentication functionality for the Molecular Data Management and CRO Integration Platform.
 * This hook encapsulates authentication logic, user management, and role-based access control without external dependencies.
 * 
 * @module features/auth/hooks/useAuth
 */

import { useState, useEffect, useCallback } from 'react'; // react 18.2+
import { useAppDispatch, useAppSelector } from '../../../store';
import { 
  loginUser, 
  registerUser, 
  logoutUser, 
  refreshUserToken,
  authActions
} from '../../../store/auth/authSlice';
import { LoginCredentials, RegisterData, User, Permission } from '../../../types/auth';
import { UserRole } from '../../../types/user';
import { hasPermission, hasRole } from '../../../utils/auth';

/**
 * Custom hook that provides authentication functionality
 * 
 * @returns Authentication methods and state
 */
export function useAuth() {
  // Get the Redux dispatch function
  const dispatch = useAppDispatch();
  
  // Select authentication state from Redux store
  const { isAuthenticated, user, loading, error } = useAppSelector(state => state.auth);
  
  /**
   * Login with email and password
   * 
   * @param credentials - Login credentials (email, password)
   * @returns Promise that resolves when login completes
   */
  const login = useCallback(
    (credentials: LoginCredentials) => {
      return dispatch(loginUser(credentials)).unwrap();
    },
    [dispatch]
  );
  
  /**
   * Register a new user
   * 
   * @param data - Registration data (email, password, confirm_password, role)
   * @returns Promise that resolves when registration completes
   */
  const register = useCallback(
    (data: RegisterData) => {
      return dispatch(registerUser(data)).unwrap();
    },
    [dispatch]
  );
  
  /**
   * Logout the current user
   * 
   * @returns Promise that resolves when logout completes
   */
  const logout = useCallback(() => {
    return dispatch(logoutUser()).unwrap();
  }, [dispatch]);
  
  /**
   * Refresh the authentication token
   * 
   * @returns Promise that resolves when token refresh completes
   */
  const refreshToken = useCallback(() => {
    return dispatch(refreshUserToken()).unwrap();
  }, [dispatch]);
  
  /**
   * Clear authentication error
   */
  const clearError = useCallback(() => {
    dispatch(authActions.clearAuthError());
  }, [dispatch]);
  
  /**
   * Check if user has the specified permission
   * 
   * @param permission - Permission to check
   * @returns True if user has permission, false otherwise
   */
  const checkPermission = useCallback(
    (permission: Permission): boolean => {
      return hasPermission(user, permission);
    },
    [user]
  );
  
  /**
   * Check if user has the specified role
   * 
   * @param role - Role to check
   * @returns True if user has role, false otherwise
   */
  const checkRole = useCallback(
    (role: UserRole): boolean => {
      return hasRole(user, role);
    },
    [user]
  );
  
  // Return authentication state and methods
  return {
    // Authentication state
    isAuthenticated,
    user,
    loading,
    error,
    
    // Authentication methods
    login,
    register,
    logout,
    refreshToken,
    clearError,
    
    // Permission and role checking
    hasPermission: checkPermission,
    hasRole: checkRole
  };
}

export default useAuth;