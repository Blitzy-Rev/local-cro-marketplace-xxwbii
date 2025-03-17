/**
 * API client functions for notification operations in the Molecular Data Management and CRO Integration Platform.
 * This module provides methods for fetching, updating, and managing user notifications through the backend API.
 * 
 * @module api/notifications
 */

import { apiClient, handleApiError } from './client';
import {
  Notification,
  NotificationFilter,
  NotificationList,
  NotificationCount,
  NotificationUpdate,
  NotificationBulkUpdate
} from '../types/notification';

/**
 * Base URL for notification endpoints.
 * Points to the current user's notifications.
 */
const BASE_URL = '/users/me/notifications';

/**
 * Fetches paginated notifications for the current user with optional filtering.
 * 
 * @param page - Page number (1-based indexing)
 * @param pageSize - Number of items per page
 * @param filter - Optional filter criteria
 * @returns Promise resolving to paginated notification data
 */
export const getNotifications = async (
  page = 1,
  pageSize = 10,
  filter: NotificationFilter = {}
): Promise<{ items: Notification[], total: number, page: number, size: number, pages: number, unread_count: number }> => {
  try {
    // Construct query parameters
    const params = {
      page,
      size: pageSize,
      ...filter
    };

    // Make the API request
    const response = await apiClient.get<NotificationList>(BASE_URL, { params });
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches total and unread notification counts for the current user.
 * 
 * @returns Promise resolving to notification count data
 */
export const getNotificationCount = async (): Promise<NotificationCount> => {
  try {
    const response = await apiClient.get<NotificationCount>(`${BASE_URL}/count`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Marks a specific notification as read.
 * 
 * @param id - ID of the notification to mark as read
 * @returns Promise resolving to the updated notification
 */
export const markNotificationAsRead = async (id: number): Promise<Notification> => {
  try {
    const update: NotificationUpdate = { read_status: true };
    const response = await apiClient.put<Notification>(`${BASE_URL}/${id}`, update);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Marks all notifications for the current user as read.
 * 
 * @returns Promise resolving to success status and count of updated notifications
 */
export const markAllNotificationsAsRead = async (): Promise<{ success: boolean, count: number }> => {
  try {
    const response = await apiClient.put<{ success: boolean, count: number }>(`${BASE_URL}/read-all`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes a specific notification.
 * 
 * @param id - ID of the notification to delete
 * @returns Promise resolving to success status
 */
export const deleteNotification = async (id: number): Promise<{ success: boolean }> => {
  try {
    const response = await apiClient.delete<{ success: boolean }>(`${BASE_URL}/${id}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes all notifications for the current user.
 * 
 * @returns Promise resolving to success status and count of deleted notifications
 */
export const deleteAllNotifications = async (): Promise<{ success: boolean, count: number }> => {
  try {
    const response = await apiClient.delete<{ success: boolean, count: number }>(`${BASE_URL}/delete-all`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates multiple notifications with the same status.
 * 
 * @param updateData - Object containing notification IDs and desired status
 * @returns Promise resolving to success status and count of updated notifications
 */
export const bulkUpdateNotifications = async (
  updateData: NotificationBulkUpdate
): Promise<{ success: boolean, count: number }> => {
  try {
    const response = await apiClient.put<{ success: boolean, count: number }>(
      `${BASE_URL}/bulk-update`,
      updateData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};