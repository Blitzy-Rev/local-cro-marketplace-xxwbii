/**
 * Core authentication utility module that provides functions for JWT token management,
 * user data persistence, permission checking, and authentication state management.
 * This module serves as the foundation for the application's authentication system
 * without external dependencies.
 */

import { User, UserRole, TokenData, Permission, AuthState } from '../types/auth';
import { setLocalStorageItem, getLocalStorageItem, removeLocalStorageItem } from './storage';
import jwtDecode from 'jwt-decode'; // jwt-decode ^3.1.2

// Storage keys
const AUTH_TOKEN_KEY = 'auth_token';
const REFRESH_TOKEN_KEY = 'refresh_token';
const USER_KEY = 'user_data';

// Define permissions for each role
const pharmaPermissions: Permission[] = [
  'molecules:view',
  'molecules:create',
  'molecules:edit',
  'molecules:delete',
  'libraries:view',
  'libraries:create',
  'libraries:edit',
  'libraries:delete',
  'experiments:view',
  'experiments:create',
  'experiments:edit',
  'submissions:view',
  'submissions:create',
  'submissions:approve',
  'results:view',
];

const croPermissions: Permission[] = [
  'submissions:view',
  'submissions:quote',
  'submissions:process',
  'results:view',
  'results:create',
  'results:edit',
  'communications:view',
  'communications:create',
];

const adminPermissions: Permission[] = [
  // Admin has all permissions
  'molecules:view',
  'molecules:create',
  'molecules:edit',
  'molecules:delete',
  'libraries:view',
  'libraries:create',
  'libraries:edit',
  'libraries:delete',
  'experiments:view',
  'experiments:create',
  'experiments:edit',
  'submissions:view',
  'submissions:create',
  'submissions:approve',
  'submissions:quote',
  'submissions:process',
  'results:view',
  'results:create',
  'results:edit',
  'communications:view',
  'communications:create',
  'users:view',
  'users:create',
  'users:edit',
  'users:delete',
  'system:configure',
  'system:monitor',
  'logs:view',
];

// Role permissions mapping
const ROLE_PERMISSIONS = {
  [UserRole.PHARMA]: [...pharmaPermissions],
  [UserRole.CRO]: [...croPermissions],
  [UserRole.ADMIN]: [...adminPermissions]
};

/**
 * Stores the authentication token in localStorage
 * @param token JWT token string
 */
export const setAuthToken = (token: string): void => {
  if (!token) {
    removeAuthToken();
    return;
  }
  setLocalStorageItem(AUTH_TOKEN_KEY, token);
};

/**
 * Retrieves the authentication token from localStorage
 * @returns The stored authentication token or null if not found
 */
export const getAuthToken = (): string | null => {
  return getLocalStorageItem<string>(AUTH_TOKEN_KEY);
};

/**
 * Removes the authentication token from localStorage
 */
export const removeAuthToken = (): void => {
  removeLocalStorageItem(AUTH_TOKEN_KEY);
};

/**
 * Stores the refresh token in localStorage
 * @param token JWT refresh token string
 */
export const setRefreshToken = (token: string): void => {
  if (!token) {
    removeRefreshToken();
    return;
  }
  setLocalStorageItem(REFRESH_TOKEN_KEY, token);
};

/**
 * Retrieves the refresh token from localStorage
 * @returns The stored refresh token or null if not found
 */
export const getRefreshToken = (): string | null => {
  return getLocalStorageItem<string>(REFRESH_TOKEN_KEY);
};

/**
 * Removes the refresh token from localStorage
 */
export const removeRefreshToken = (): void => {
  removeLocalStorageItem(REFRESH_TOKEN_KEY);
};

/**
 * Stores the user data in localStorage
 * @param user User object
 */
export const setUser = (user: User): void => {
  if (!user) {
    removeUser();
    return;
  }
  setLocalStorageItem(USER_KEY, user);
};

/**
 * Retrieves the user data from localStorage
 * @returns The stored user data or null if not found
 */
export const getUser = (): User | null => {
  return getLocalStorageItem<User>(USER_KEY);
};

/**
 * Removes the user data from localStorage
 */
export const removeUser = (): void => {
  removeLocalStorageItem(USER_KEY);
};

/**
 * Clears all authentication data from localStorage
 */
export const clearAuth = (): void => {
  removeAuthToken();
  removeRefreshToken();
  removeUser();
};

/**
 * Decodes a JWT token to extract its payload
 * @param token JWT token string
 * @returns Decoded token data or null if token is invalid
 */
export const decodeToken = (token: string): TokenData | null => {
  if (!token) {
    return null;
  }
  
  try {
    return jwtDecode<TokenData>(token);
  } catch (error) {
    console.error('Failed to decode token:', error);
    return null;
  }
};

/**
 * Checks if a JWT token is expired
 * @param token JWT token string
 * @returns True if the token is expired or invalid, false otherwise
 */
export const isTokenExpired = (token: string): boolean => {
  const decoded = decodeToken(token);
  
  if (!decoded) {
    return true;
  }
  
  const currentTime = Math.floor(Date.now() / 1000);
  
  return decoded.exp < currentTime;
};

/**
 * Initializes the authentication state from localStorage
 * @returns The initial authentication state
 */
export const getInitialAuthState = (): AuthState => {
  const token = getAuthToken();
  const refreshToken = getRefreshToken();
  const user = getUser();
  
  // Check if the token is valid
  let isAuthenticated = false;
  
  if (token && !isTokenExpired(token) && user) {
    isAuthenticated = true;
  } else if (token && isTokenExpired(token) && refreshToken) {
    // Token is expired but we have a refresh token
    isAuthenticated = false;
  } else {
    isAuthenticated = false;
  }
  
  return {
    isAuthenticated,
    user,
    accessToken: token,
    refreshToken
  };
};

/**
 * Gets the permissions associated with a user role
 * @param role User role
 * @returns Array of permissions for the specified role
 */
export const getRolePermissions = (role: UserRole): Permission[] => {
  if (!role) {
    return [];
  }
  
  return ROLE_PERMISSIONS[role] || [];
};

/**
 * Checks if a user has a specific permission
 * @param user User object or null
 * @param permission Permission to check
 * @returns True if the user has the permission, false otherwise
 */
export const hasPermission = (user: User | null, permission: Permission): boolean => {
  if (!user) {
    return false;
  }
  
  // Admins have all permissions
  if (user.role === UserRole.ADMIN) {
    return true;
  }
  
  const permissions = getRolePermissions(user.role);
  
  return permissions.includes(permission);
};

/**
 * Checks if a user has a specific role
 * @param user User object or null
 * @param role Role to check
 * @returns True if the user has the role, false otherwise
 */
export const hasRole = (user: User | null, role: UserRole): boolean => {
  if (!user) {
    return false;
  }
  
  return user.role === role;
};