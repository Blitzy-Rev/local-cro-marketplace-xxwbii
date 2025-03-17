import { useState, useEffect, useCallback, useMemo } from 'react'; // react ^18.2.0
import { useMutation, useQuery } from 'react-query'; // react-query ^3.39.0
import { useAppDispatch, useAppSelector } from '../../../store';
import { useToast } from '../../../hooks/useToast';
import { usePermissions } from '../../../hooks/usePermissions';
import {
  getUsers,
  getUser,
  createUser,
  updateUser,
  deleteUser,
  resetUserPassword,
  getSystemStats,
  getSystemResources,
  checkSystemHealth,
  getActivityLogs,
  getSystemAlerts,
  updateSystemAlert,
} from '../../../api/admin';
import {
  User,
  UserRole,
  UserStatus,
  UserAdminUpdateRequest,
  UserPasswordResetRequest,
  UserListResponse,
} from '../../../types/user';
import { ApiResponse } from '../../../types';

/**
 * Interface for user filtering options
 */
interface UserFilters {
  role: UserRole | null;
  status: UserStatus | null;
  email: string | null;
}

/**
 * Interface for pagination state
 */
interface Pagination {
  page: number;
  pageSize: number;
  totalPages: number;
  totalItems: number;
}

/**
 * Interface for system statistics response
 */
interface SystemStatsResponse {
  userStats: any; // Replace 'any' with a more specific type if available
  moleculeStats: any; // Replace 'any' with a more specific type if available
  experimentStats: any; // Replace 'any' with a more specific type if available
  resourceUsage: any; // Replace 'any' with a more specific type if available
}

/**
 * Interface for system resources response
 */
interface SystemResourcesResponse {
  cpu: any; // Replace 'any' with a more specific type if available
  memory: any; // Replace 'any' with a more specific type if available
  disk: any; // Replace 'any' with a more specific type if available
  network: any; // Replace 'any' with a more specific type if available
}

/**
 * Interface for activity log item
 */
interface ActivityLogItem {
  id: number;
  userId: number;
  userEmail: string;
  action: string;
  resourceType: string;
  resourceId: string;
  timestamp: string;
  details: any; // Replace 'any' with a more specific type if available
}

/**
 * Interface for paginated activity logs response
 */
interface ActivityLogListResponse {
  items: ActivityLogItem[];
  total: number;
  page: number;
  size: number;
}

/**
 * Enum for system alert severity levels
 */
enum SystemAlertSeverity {
  INFO = 'INFO',
  WARNING = 'WARNING',
  CRITICAL = 'CRITICAL',
}

/**
 * Interface for system alert item
 */
interface SystemAlertResponse {
  id: number;
  severity: SystemAlertSeverity;
  alertType: string;
  message: string;
  details: any; // Replace 'any' with a more specific type if available
  timestamp: string;
  resolved: boolean;
  resolvedBy: number | null;
  resolvedAt: string | null;
}

/**
 * Interface for paginated system alerts response
 */
interface SystemAlertListResponse {
  items: SystemAlertResponse[];
  total: number;
  page: number;
  size: number;
}

/**
 * Interface for updating system alert
 */
interface SystemAlertUpdateRequest {
  resolved: boolean;
  notes: string;
}

/**
 * Custom hook that provides admin functionality for user management, system monitoring, and activity logs
 * @returns Object containing admin state and functions
 */
export const useAdmin = () => {
  // Check if user has admin permissions using usePermissions hook
  const { isAdmin } = usePermissions();

  // Initialize state for users, filters, pagination, system stats, resources, health, activity logs, and alerts
  const [filters, setFilters] = useState<UserFilters>({
    role: null,
    status: null,
    email: null,
  });
  const [pagination, setPagination] = useState<Pagination>({
    page: 1,
    pageSize: 10,
    totalPages: 1,
    totalItems: 0,
  });

  // Set up query for fetching users with useQuery
  const {
    data: usersData,
    isLoading: isLoadingUsers,
    refetch: refetchUsers,
  } = useQuery<ApiResponse<UserListResponse>, Error>(
    ['users', filters, pagination.page, pagination.pageSize],
    () =>
      getUsers({
        role: filters.role || undefined,
        status: filters.status || undefined,
        page: pagination.page,
        limit: pagination.pageSize,
      }),
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
      select: (data) => {
        if (data?.success && data.data) {
          setPagination({
            ...pagination,
            totalPages: Math.ceil(data.data.total / pagination.pageSize),
            totalItems: data.data.total,
          });
          return data.data.items;
        }
        return [];
      },
    }
  );

  // Set up query for fetching system stats with useQuery
  const {
    data: systemStats,
    isLoading: isLoadingSystemStats,
    refetch: refetchSystemStats,
  } = useQuery<ApiResponse<SystemStatsResponse>, Error>(
    'systemStats',
    getSystemStats,
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
    }
  );

  // Set up query for fetching system resources with useQuery
  const {
    data: systemResources,
    isLoading: isLoadingSystemResources,
    refetch: refetchSystemResources,
  } = useQuery<ApiResponse<SystemResourcesResponse>, Error>(
    'systemResources',
    getSystemResources,
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
    }
  );

  // Set up query for fetching system health with useQuery
  const {
    data: systemHealth,
    isLoading: isLoadingSystemHealth,
    refetch: refetchSystemHealth,
  } = useQuery<ApiResponse<Record<string, any>>, Error>(
    'systemHealth',
    checkSystemHealth,
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
    }
  );

  // Set up query for fetching activity logs with useQuery
  const [activityLogPage, setActivityLogPage] = useState(1);
  const [activityLogPageSize, setActivityLogPageSize] = useState(10);

  const {
    data: activityLogsData,
    isLoading: isLoadingActivityLogs,
    refetch: refetchActivityLogs,
  } = useQuery<ApiResponse<ActivityLogListResponse>, Error>(
    ['activityLogs', activityLogPage, activityLogPageSize],
    () =>
      getActivityLogs({
        page: activityLogPage,
        size: activityLogPageSize,
      }),
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
    }
  );

  // Set up query for fetching system alerts with useQuery
  const [systemAlertPage, setSystemAlertPage] = useState(1);
  const [systemAlertPageSize, setSystemAlertPageSize] = useState(10);

  const {
    data: systemAlertsData,
    isLoading: isLoadingSystemAlerts,
    refetch: refetchSystemAlerts,
  } = useQuery<ApiResponse<SystemAlertListResponse>, Error>(
    ['systemAlerts', systemAlertPage, systemAlertPageSize],
    () =>
      getSystemAlerts({
        page: systemAlertPage,
        size: systemAlertPageSize,
      }),
    {
      enabled: isAdmin(), // Only fetch if the user is an admin
    }
  );

  // Set up mutation for creating users with useMutation
  const createUserMutation = useMutation(createUser, {
    onSuccess: () => {
      refetchUsers();
    },
  });

  // Set up mutation for updating users with useMutation
  const updateUserMutation = useMutation(
    ({ userId, userData }: { userId: number; userData: UserAdminUpdateRequest }) =>
      updateUser(userId, userData),
    {
      onSuccess: () => {
        refetchUsers();
      },
    }
  );

  // Set up mutation for deleting users with useMutation
  const deleteUserMutation = useMutation(deleteUser, {
    onSuccess: () => {
      refetchUsers();
    },
  });

  // Set up mutation for resetting user passwords with useMutation
  const resetPasswordMutation = useMutation(
    ({ userId, passwordData }: { userId: number; passwordData: UserPasswordResetRequest }) =>
      resetUserPassword(userId, passwordData),
    {
      onSuccess: () => {
        // Optionally show a success message
      },
    }
  );

  // Set up mutation for updating system alerts with useMutation
  const updateAlertMutation = useMutation(
    ({ alertId, alertData }: { alertId: number; alertData: SystemAlertUpdateRequest }) =>
      updateSystemAlert(alertId, alertData),
    {
      onSuccess: () => {
        refetchSystemAlerts();
      },
    }
  );

  // Create callback function for handling filter changes
  const handleFilterChange = useCallback((newFilters: UserFilters) => {
    setFilters(newFilters);
    setPagination({ ...pagination, page: 1 }); // Reset pagination on filter change
  }, [pagination]);

  // Create callback function for handling pagination changes
  const handlePageChange = useCallback((newPage: number) => {
    setPagination({ ...pagination, page: newPage });
  }, [pagination]);

  // Return object with all state and functions
  return {
    users: usersData || [],
    isLoadingUsers,
    filters,
    pagination,
    systemStats,
    isLoadingSystemStats,
    systemResources,
    isLoadingSystemResources,
    systemHealth,
    isLoadingSystemHealth,
    activityLogs: activityLogsData?.items || [],
    isLoadingActivityLogs,
    systemAlerts: systemAlertsData?.items || [],
    isLoadingSystemAlerts,
    handleFilterChange,
    handlePageChange,
    createUserMutation,
    updateUserMutation,
    deleteUserMutation,
    resetPasswordMutation,
    updateAlertMutation,
    refetchUsers,
    refetchSystemStats,
    refetchSystemResources,
    refetchSystemHealth,
    refetchActivityLogs,
    refetchSystemAlerts,
  };
};