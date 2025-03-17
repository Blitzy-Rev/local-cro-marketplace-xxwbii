/**
 * TypeScript interfaces for notification-related data structures in the
 * Molecular Data Management and CRO Integration Platform.
 */

/**
 * Enumeration of notification types
 * Corresponds to backend NotificationType enum
 */
export enum NotificationType {
  EXPERIMENT_STATUS_CHANGE = 'EXPERIMENT_STATUS_CHANGE',
  SUBMISSION_CREATED = 'SUBMISSION_CREATED',
  QUOTE_PROVIDED = 'QUOTE_PROVIDED',
  RESULTS_UPLOADED = 'RESULTS_UPLOADED',
  SYSTEM_ALERT = 'SYSTEM_ALERT',
  USER_MENTION = 'USER_MENTION'
}

/**
 * Interface for notification data payload with type-specific properties
 * Contains optional fields that may be present depending on notification type
 */
export interface NotificationData {
  experiment_id?: number;
  experiment_name?: string;
  submission_id?: number;
  result_id?: number;
  status?: string;
  price?: number;
  cro_name?: string;
  pharma_company?: string;
}

/**
 * Interface for notification data structure
 * Corresponds to backend NotificationRead schema
 */
export interface Notification {
  id: number;
  user_id: number;
  type: NotificationType;
  message: string;
  data?: NotificationData;
  read_status: boolean;
  created_at: string;
  read_at?: string;
}

/**
 * Interface for notification filtering parameters
 * Corresponds to backend NotificationFilter schema
 */
export interface NotificationFilter {
  read_status?: boolean | null;
  type?: NotificationType;
  from_date?: string;
  to_date?: string;
}

/**
 * Interface for paginated notification list response
 * Corresponds to backend NotificationList schema
 */
export interface NotificationList {
  items: Notification[];
  total: number;
  page: number;
  size: number;
  pages: number;
}

/**
 * Interface for notification count response
 * Corresponds to backend NotificationCount schema
 */
export interface NotificationCount {
  total: number;
  unread: number;
}

/**
 * Interface for notification update request
 * Corresponds to backend NotificationUpdate schema
 */
export interface NotificationUpdate {
  read_status: boolean;
}

/**
 * Interface for bulk notification update request
 * Corresponds to backend NotificationBulkUpdate schema
 */
export interface NotificationBulkUpdate {
  notification_ids: number[];
  read_status: boolean;
}

/**
 * Interface for notification Redux state structure
 */
export interface NotificationState {
  notifications: Notification[];
  total: number;
  unreadCount: number;
  loading: boolean;
  error: string | null;
  page: number;
  pageSize: number;
}