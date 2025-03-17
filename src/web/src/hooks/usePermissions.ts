import { useCallback } from 'react'; // react ^18.2.0
import { useAppSelector } from '../store';
import { User, Permission } from '../types/auth';
import { UserRole } from '../types/user';
import { hasPermission, hasRole } from '../utils/auth';

/**
 * Custom hook that provides permission checking functionality for the Molecular
 * Data Management and CRO Integration Platform. This hook encapsulates
 * role-based access control to ensure consistent permission verification
 * throughout the application.
 * 
 * @returns Object containing permission checking functions
 */
export const usePermissions = () => {
  // Get the current user from Redux store
  const user = useAppSelector((state) => state.auth.user);

  /**
   * Checks if the current user has a specific permission
   * @param permission - The permission to check
   * @returns Boolean indicating if user has the permission
   */
  const checkPermission = useCallback(
    (permission: Permission): boolean => {
      return hasPermission(user, permission);
    },
    [user]
  );

  /**
   * Checks if the current user has a specific role
   * @param role - The role to check
   * @returns Boolean indicating if user has the role
   */
  const checkRole = useCallback(
    (role: UserRole): boolean => {
      return hasRole(user, role);
    },
    [user]
  );

  /**
   * Checks if the current user has admin role
   * @returns Boolean indicating if user is an admin
   */
  const isAdmin = useCallback((): boolean => {
    return checkRole(UserRole.ADMIN);
  }, [checkRole]);

  /**
   * Checks if the current user has pharma role
   * @returns Boolean indicating if user is a pharma user
   */
  const isPharma = useCallback((): boolean => {
    return checkRole(UserRole.PHARMA);
  }, [checkRole]);

  /**
   * Checks if the current user has CRO role
   * @returns Boolean indicating if user is a CRO user
   */
  const isCRO = useCallback((): boolean => {
    return checkRole(UserRole.CRO);
  }, [checkRole]);

  /**
   * Checks if the current user can access molecule data
   * @returns Boolean indicating if user can access molecules
   */
  const canAccessMolecules = useCallback((): boolean => {
    return checkPermission('molecules:view');
  }, [checkPermission]);

  /**
   * Checks if the current user can manage libraries
   * @returns Boolean indicating if user can view and create libraries
   */
  const canManageLibraries = useCallback((): boolean => {
    return checkPermission('libraries:view') && checkPermission('libraries:create');
  }, [checkPermission]);

  /**
   * Checks if the current user can submit experiments to CROs
   * @returns Boolean indicating if user can create and submit experiments
   */
  const canSubmitExperiments = useCallback((): boolean => {
    return checkPermission('experiments:create') && checkPermission('submissions:create');
  }, [checkPermission]);

  /**
   * Checks if the current user can review submissions (typical CRO action)
   * @returns Boolean indicating if user can review submissions
   */
  const canReviewSubmissions = useCallback((): boolean => {
    return checkPermission('submissions:view') && checkRole(UserRole.CRO);
  }, [checkPermission, checkRole]);

  /**
   * Checks if the current user can upload experimental results
   * @returns Boolean indicating if user can create result records
   */
  const canUploadResults = useCallback((): boolean => {
    return checkPermission('results:create') && checkRole(UserRole.CRO);
  }, [checkPermission, checkRole]);

  /**
   * Checks if the current user can manage system users
   * @returns Boolean indicating if user can view and edit users
   */
  const canManageUsers = useCallback((): boolean => {
    return checkPermission('users:view') && checkPermission('users:edit');
  }, [checkPermission]);

  return {
    checkPermission,
    checkRole,
    isAdmin,
    isPharma,
    isCRO,
    canAccessMolecules,
    canManageLibraries,
    canSubmitExperiments,
    canReviewSubmissions,
    canUploadResults,
    canManageUsers
  };
};