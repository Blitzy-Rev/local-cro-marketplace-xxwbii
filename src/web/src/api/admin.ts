/**
 * API client module for admin-related operations in the Molecular Data Management and CRO Integration Platform.
 * This module provides functions to interact with the admin endpoints for user management, system monitoring, 
 * activity logging, and system alerts management.
 * 
 * @module api/admin
 */

import { AxiosResponse } from 'axios'; // ^1.4.0
import { apiClient, handleApiError } from './client';
import { ApiResponse, User, UserDetail, UserListResponse, UserAdminUpdateRequest, UserPasswordResetRequest } from '../types';

/**
 * Interface for system statistics response
 */
interface SystemStatsResponse {
  /** Total number of users in the system */
  user_count: number;
  /** Number of active users */
  active_users: number;
  /** User counts broken down by role */
  users_by_role: {
    pharma: number;
    cro: number;
    admin: number;
  };
  /** Total number of molecules in the system */
  molecule_count: number;
  /** Total number of libraries in the system */
  library_count: number;
  /** Total number of experiments in the system */
  experiment_count: number;
  /** Total number of submissions in the system */
  submission_count: number;
  /** System uptime in seconds */
  uptime: number;
}

/**
 * Interface for system resources response
 */
interface SystemResourcesResponse {
  /** CPU usage as a percentage */
  cpu_usage: number;
  /** Memory usage statistics */
  memory: {
    /** Total memory in bytes */
    total: number;
    /** Used memory in bytes */
    used: number;
    /** Free memory in bytes */
    free: number;
    /** Memory usage as a percentage */
    usage_percent: number;
  };
  /** Disk usage statistics */
  disk: {
    /** Total disk space in bytes */
    total: number;
    /** Used disk space in bytes */
    used: number;
    /** Free disk space in bytes */
    free: number;
    /** Disk usage as a percentage */
    usage_percent: number;
  };
  /** Database statistics */
  database: {
    /** Number of active connections */
    connections: number;
    /** Database size in bytes */
    size: number;
    /** Database usage as a percentage */
    usage_percent: number;
  };
}

/**
 * Interface for activity log entry
 */
interface ActivityLog {
  /** Unique identifier for the log entry */
  id: number;
  /** User who performed the action */
  user_id: number | null;
  /** User's email address */
  user_email: string | null;
  /** Type of action performed */
  action: string;
  /** Target resource type */
  resource_type: string;
  /** Target resource identifier */
  resource_id: string | null;
  /** Additional details about the action */
  details: Record<string, any> | null;
  /** IP address where the action originated */
  ip_address: string;
  /** Timestamp when the action occurred */
  timestamp: string;
}

/**
 * Interface for paginated activity log list response
 */
interface ActivityLogListResponse {
  /** List of activity log entries */
  items: ActivityLog[];
  /** Total number of log entries matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
}

/**
 * Enum for system alert severity levels
 */
enum AlertSeverity {
  INFO = 'info',
  WARNING = 'warning',
  ERROR = 'error',
  CRITICAL = 'critical'
}

/**
 * Interface for system alert
 */
interface SystemAlert {
  /** Unique identifier for the alert */
  id: number;
  /** Alert severity level */
  severity: AlertSeverity;
  /** Alert type category */
  alert_type: string;
  /** Alert message */
  message: string;
  /** Whether the alert has been resolved */
  resolved: boolean;
  /** Timestamp when the alert was created */
  created_at: string;
  /** Timestamp when the alert was resolved, if applicable */
  resolved_at: string | null;
  /** Additional details about the alert */
  details: Record<string, any> | null;
}

/**
 * Interface for paginated system alert list response
 */
interface SystemAlertListResponse {
  /** List of system alerts */
  items: SystemAlert[];
  /** Total number of alerts matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
}

/**
 * Interface for system alert response
 */
interface SystemAlertResponse {
  /** The system alert data */
  alert: SystemAlert;
}

/**
 * Interface for system alert update request
 */
interface SystemAlertUpdateRequest {
  /** Whether the alert is resolved */
  resolved: boolean;
  /** Notes about the resolution */
  resolution_notes?: string;
}

/**
 * Get system statistics including user counts, data counts, and resource usage
 * 
 * @returns Promise resolving to system statistics
 */
export const getSystemStats = async (): Promise<ApiResponse<SystemStatsResponse>> => {
  try {
    const response = await apiClient.get<ApiResponse<SystemStatsResponse>>('/admin/stats');
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Get detailed system resource usage information
 * 
 * @returns Promise resolving to system resources
 */
export const getSystemResources = async (): Promise<ApiResponse<SystemResourcesResponse>> => {
  try {
    const response = await apiClient.get<ApiResponse<SystemResourcesResponse>>('/admin/resources');
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Perform a comprehensive system health check
 * 
 * @returns Promise resolving to health check results
 */
export const checkSystemHealth = async (): Promise<ApiResponse<Record<string, any>>> => {
  try {
    const response = await apiClient.get<ApiResponse<Record<string, any>>>('/admin/health');
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Get list of all users with optional filtering
 * 
 * @param params - Filter parameters for the user list
 * @returns Promise resolving to paginated list of users
 */
export const getUsers = async (params: {
  role?: string;
  status?: string;
  page?: number;
  limit?: number;
}): Promise<ApiResponse<UserListResponse>> => {
  try {
    const queryParams = new URLSearchParams();
    
    if (params.role) queryParams.append('role', params.role);
    if (params.status) queryParams.append('status', params.status);
    if (params.page) queryParams.append('page', params.page.toString());
    if (params.limit) queryParams.append('limit', params.limit.toString());
    
    const response = await apiClient.get<ApiResponse<UserListResponse>>(
      `/admin/users?${queryParams.toString()}`
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Get user by ID
 * 
 * @param userId - ID of the user to retrieve
 * @returns Promise resolving to user details
 */
export const getUser = async (userId: number): Promise<ApiResponse<UserDetail>> => {
  try {
    const response = await apiClient.get<ApiResponse<UserDetail>>(`/admin/users/${userId}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Create a new user
 * 
 * @param userData - User data for creation
 * @returns Promise resolving to created user
 */
export const createUser = async (userData: any): Promise<ApiResponse<User>> => {
  try {
    const response = await apiClient.post<ApiResponse<User>>('/admin/users', userData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Update user information
 * 
 * @param userId - ID of the user to update
 * @param userData - User data for update
 * @returns Promise resolving to updated user
 */
export const updateUser = async (
  userId: number,
  userData: UserAdminUpdateRequest
): Promise<ApiResponse<User>> => {
  try {
    const response = await apiClient.put<ApiResponse<User>>(`/admin/users/${userId}`, userData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Delete a user
 * 
 * @param userId - ID of the user to delete
 * @returns Promise resolving to success message
 */
export const deleteUser = async (userId: number): Promise<ApiResponse<{ message: string }>> => {
  try {
    const response = await apiClient.delete<ApiResponse<{ message: string }>>(`/admin/users/${userId}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Reset a user's password
 * 
 * @param userId - ID of the user to reset password for
 * @param passwordData - Password reset data
 * @returns Promise resolving to success message
 */
export const resetUserPassword = async (
  userId: number,
  passwordData: UserPasswordResetRequest
): Promise<ApiResponse<{ message: string }>> => {
  try {
    const response = await apiClient.post<ApiResponse<{ message: string }>>(
      `/admin/users/${userId}/reset-password`,
      passwordData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Get activity logs with optional filtering
 * 
 * @param params - Filter parameters for activity logs
 * @returns Promise resolving to paginated list of activity logs
 */
export const getActivityLogs = async (params: {
  user_id?: number;
  action?: string;
  resource_type?: string;
  resource_id?: string;
  from_date?: string;
  to_date?: string;
  page?: number;
  size?: number;
}): Promise<ApiResponse<ActivityLogListResponse>> => {
  try {
    const queryParams = new URLSearchParams();
    
    if (params.user_id) queryParams.append('user_id', params.user_id.toString());
    if (params.action) queryParams.append('action', params.action);
    if (params.resource_type) queryParams.append('resource_type', params.resource_type);
    if (params.resource_id) queryParams.append('resource_id', params.resource_id);
    if (params.from_date) queryParams.append('from_date', params.from_date);
    if (params.to_date) queryParams.append('to_date', params.to_date);
    if (params.page) queryParams.append('page', params.page.toString());
    if (params.size) queryParams.append('size', params.size.toString());
    
    const response = await apiClient.get<ApiResponse<ActivityLogListResponse>>(
      `/admin/activity-logs?${queryParams.toString()}`
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Get system alerts with optional filtering
 * 
 * @param params - Filter parameters for system alerts
 * @returns Promise resolving to list of system alerts
 */
export const getSystemAlerts = async (params: {
  severity?: AlertSeverity;
  alert_type?: string;
  resolved?: boolean;
  page?: number;
  size?: number;
}): Promise<ApiResponse<SystemAlertListResponse>> => {
  try {
    const queryParams = new URLSearchParams();
    
    if (params.severity) queryParams.append('severity', params.severity);
    if (params.alert_type) queryParams.append('alert_type', params.alert_type);
    if (params.resolved !== undefined) queryParams.append('resolved', params.resolved.toString());
    if (params.page) queryParams.append('page', params.page.toString());
    if (params.size) queryParams.append('size', params.size.toString());
    
    const response = await apiClient.get<ApiResponse<SystemAlertListResponse>>(
      `/admin/alerts?${queryParams.toString()}`
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Update system alert status
 * 
 * @param alertId - ID of the alert to update
 * @param alertData - Alert update data
 * @returns Promise resolving to updated system alert
 */
export const updateSystemAlert = async (
  alertId: number,
  alertData: SystemAlertUpdateRequest
): Promise<ApiResponse<SystemAlertResponse>> => {
  try {
    const response = await apiClient.put<ApiResponse<SystemAlertResponse>>(
      `/admin/alerts/${alertId}`,
      alertData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};