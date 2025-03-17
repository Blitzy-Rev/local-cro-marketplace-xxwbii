/**
 * User-related TypeScript definitions for the Molecular Data Management and CRO Integration Platform.
 * These types support the role-based authentication system and user management functionality.
 * @packageDocumentation
 */

/**
 * Enumeration of user roles for role-based access control.
 * Used to determine permissions and access to features within the application.
 */
export enum UserRole {
  /** Pharmaceutical company researcher role with molecule management capabilities */
  PHARMA = 'pharma',
  /** Contract Research Organization role with experiment processing capabilities */
  CRO = 'cro',
  /** Administrator role with full system access */
  ADMIN = 'admin'
}

/**
 * Enumeration of possible user account statuses.
 * Controls whether users can access the system and what they can do.
 */
export enum UserStatus {
  /** Account created but not yet verified */
  PENDING = 'pending',
  /** Fully functional account */
  ACTIVE = 'active',
  /** Deactivated account with no system access */
  INACTIVE = 'inactive',
  /** Account locked due to security concerns */
  LOCKED = 'locked'
}

/**
 * Base interface for user data with common fields.
 * Used as a foundation for more specific user interfaces.
 */
export interface UserBase {
  /** User's email address, serves as unique identifier and login credential */
  email: string;
}

/**
 * Interface for user creation data.
 * Contains all fields needed to register a new user.
 */
export interface UserCreate extends UserBase {
  /** New user's password */
  password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
  /** User's role in the system */
  role: UserRole;
}

/**
 * Interface for updating user information.
 * All fields are optional as updates may modify only a subset of user data.
 */
export interface UserUpdate {
  /** User's email address */
  email?: string;
  /** User's password - only provided when changing password */
  password?: string;
  /** User's role in the system */
  role?: UserRole;
  /** User's account status */
  status?: UserStatus;
}

/**
 * Interface representing a user in the system.
 * Contains all standard user information.
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
 * Interface for detailed user information including all fields.
 * Used for admin views and detailed user profiles.
 */
export interface UserDetail extends User {
  /** Timestamp when the user was last updated */
  updated_at: string;
}

/**
 * Interface for user list items with essential information.
 * Used in user management tables and lists.
 */
export interface UserListItem {
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
}

/**
 * Interface for paginated list of users response.
 * Used for API responses that return lists of users.
 */
export interface UserListResponse {
  /** Array of user items */
  items: UserListItem[];
  /** Total number of users matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
}

/**
 * Interface for admin user update request with additional fields.
 * Allows administrators to modify more user properties than regular users.
 */
export interface UserAdminUpdateRequest {
  /** User's email address */
  email?: string;
  /** User's role in the system */
  role?: UserRole;
  /** User's account status */
  status?: UserStatus;
  /** Whether the user's email has been verified */
  email_verified?: boolean;
  /** Whether the user account is active */
  is_active?: boolean;
}

/**
 * Interface for password reset request.
 * Used when resetting password via a reset token.
 */
export interface UserPasswordResetRequest {
  /** New password to set */
  new_password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
}

/**
 * Interface for changing user password.
 * Requires current password for security verification.
 */
export interface UserPasswordChangeRequest {
  /** Current password for verification */
  current_password: string;
  /** New password to set */
  new_password: string;
  /** Password confirmation to prevent typing errors */
  confirm_password: string;
}

/**
 * Interface for user statistics data.
 * Used for admin dashboards and system monitoring.
 */
export interface UserStats {
  /** Total number of users in the system */
  total: number;
  /** Number of active users */
  active: number;
  /** User counts broken down by role */
  by_role: Record<UserRole, number>;
}

/**
 * Type for permission strings in the system.
 * Represents specific actions users can perform.
 */
export type Permission = string;

/**
 * Type for mapping roles to their permissions.
 * Used for role-based access control.
 */
export type RolePermissions = Record<UserRole, Permission[]>;