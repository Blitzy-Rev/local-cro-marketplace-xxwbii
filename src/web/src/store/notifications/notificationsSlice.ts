/**
 * Redux slice for managing notification state in the Molecular Data Management and CRO Integration Platform.
 * This slice handles fetching, updating, and managing user notifications with Redux Toolkit.
 * 
 * @packageDocumentation
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import {
  Notification,
  NotificationState,
  NotificationFilter,
  NotificationUpdate,
  NotificationBulkUpdate
} from '../../types/notification';
import {
  getNotifications,
  getNotificationCount,
  markNotificationAsRead,
  markAllNotificationsAsRead,
  deleteNotification,
  deleteAllNotifications,
  bulkUpdateNotifications
} from '../../api/notifications';

// Define the initial state
const initialState: NotificationState = {
  notifications: [],
  total: 0,
  unreadCount: 0,
  loading: false,
  error: null,
  page: 1,
  pageSize: 10
};

// Define the interface for the root state used by selectors
export interface NotificationsRootState {
  notifications: NotificationState;
}

/**
 * Async thunk that fetches notifications from the API
 */
export const fetchNotifications = createAsyncThunk(
  'notifications/fetchNotifications',
  async ({ page, pageSize, filter }: { 
    page?: number, 
    pageSize?: number, 
    filter?: NotificationFilter 
  }, thunkAPI) => {
    try {
      const response = await getNotifications(
        page || 1, 
        pageSize || 10, 
        filter || {}
      );
      return response;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch notifications');
    }
  }
);

/**
 * Async thunk that fetches notification counts from the API
 */
export const fetchNotificationCount = createAsyncThunk(
  'notifications/fetchNotificationCount',
  async (_, thunkAPI) => {
    try {
      const response = await getNotificationCount();
      return response;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch notification count');
    }
  }
);

/**
 * Async thunk that marks a notification as read
 */
export const markAsRead = createAsyncThunk(
  'notifications/markAsRead',
  async (id: number, thunkAPI) => {
    try {
      const response = await markNotificationAsRead(id);
      return response;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to mark notification as read');
    }
  }
);

/**
 * Async thunk that marks all notifications as read
 */
export const markAllAsRead = createAsyncThunk(
  'notifications/markAllAsRead',
  async (_, thunkAPI) => {
    try {
      const response = await markAllNotificationsAsRead();
      return response;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to mark all notifications as read');
    }
  }
);

/**
 * Async thunk that deletes a notification
 */
export const deleteNotificationById = createAsyncThunk(
  'notifications/deleteNotification',
  async (id: number, thunkAPI) => {
    try {
      const response = await deleteNotification(id);
      return { success: response.success, id };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to delete notification');
    }
  }
);

/**
 * Async thunk that deletes all notifications
 */
export const clearAllNotifications = createAsyncThunk(
  'notifications/clearAllNotifications',
  async (_, thunkAPI) => {
    try {
      const response = await deleteAllNotifications();
      return response;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to clear notifications');
    }
  }
);

/**
 * Async thunk that updates multiple notifications with the same status
 */
export const bulkUpdateNotificationStatus = createAsyncThunk(
  'notifications/bulkUpdateStatus',
  async (updateData: NotificationBulkUpdate, thunkAPI) => {
    try {
      const response = await bulkUpdateNotifications(updateData);
      return { ...response, updateData };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to update notifications');
    }
  }
);

// Create the notifications slice
const notificationsSlice = createSlice({
  name: 'notifications',
  initialState,
  reducers: {
    // Set the current page
    setPage(state, action: PayloadAction<number>) {
      state.page = action.payload;
    },
    // Reset state to initial values
    resetState() {
      return initialState;
    }
  },
  extraReducers: (builder) => {
    // fetchNotifications
    builder
      .addCase(fetchNotifications.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(fetchNotifications.fulfilled, (state, action) => {
        state.loading = false;
        state.notifications = action.payload.items;
        state.total = action.payload.total;
        state.unreadCount = action.payload.unread_count || 0;
        state.page = action.payload.page;
        state.pageSize = action.payload.size;
      })
      .addCase(fetchNotifications.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to fetch notifications';
      });

    // fetchNotificationCount
    builder
      .addCase(fetchNotificationCount.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(fetchNotificationCount.fulfilled, (state, action) => {
        state.loading = false;
        state.total = action.payload.total;
        state.unreadCount = action.payload.unread;
      })
      .addCase(fetchNotificationCount.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to fetch notification count';
      });

    // markAsRead
    builder
      .addCase(markAsRead.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(markAsRead.fulfilled, (state, action) => {
        state.loading = false;
        // Update the notification in the list
        const index = state.notifications.findIndex(n => n.id === action.payload.id);
        if (index !== -1) {
          const wasUnread = !state.notifications[index].read_status;
          state.notifications[index] = action.payload;
          
          // If it was unread and is now read, decrement the unread count
          if (wasUnread && action.payload.read_status) {
            state.unreadCount = Math.max(0, state.unreadCount - 1);
          }
        }
      })
      .addCase(markAsRead.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to mark notification as read';
      });

    // markAllAsRead
    builder
      .addCase(markAllAsRead.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(markAllAsRead.fulfilled, (state, action) => {
        state.loading = false;
        // Mark all notifications as read in the list
        state.notifications = state.notifications.map(notification => ({
          ...notification,
          read_status: true,
          read_at: new Date().toISOString()
        }));
        // Reset unread count
        state.unreadCount = 0;
      })
      .addCase(markAllAsRead.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to mark all notifications as read';
      });

    // deleteNotificationById
    builder
      .addCase(deleteNotificationById.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(deleteNotificationById.fulfilled, (state, action) => {
        state.loading = false;
        // Remove the notification from the list
        const deletedNotification = state.notifications.find(n => n.id === action.payload.id);
        state.notifications = state.notifications.filter(n => n.id !== action.payload.id);
        
        // Decrement total
        state.total = Math.max(0, state.total - 1);
        
        // Decrement unread count if the notification was unread
        if (deletedNotification && !deletedNotification.read_status) {
          state.unreadCount = Math.max(0, state.unreadCount - 1);
        }
      })
      .addCase(deleteNotificationById.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to delete notification';
      });

    // clearAllNotifications
    builder
      .addCase(clearAllNotifications.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(clearAllNotifications.fulfilled, (state) => {
        state.loading = false;
        // Clear all notifications
        state.notifications = [];
        state.total = 0;
        state.unreadCount = 0;
      })
      .addCase(clearAllNotifications.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to clear notifications';
      });

    // bulkUpdateNotificationStatus
    builder
      .addCase(bulkUpdateNotificationStatus.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(bulkUpdateNotificationStatus.fulfilled, (state, action) => {
        state.loading = false;
        
        const { notification_ids, read_status } = action.meta.arg;
        
        // Update notifications in the list
        state.notifications = state.notifications.map(notification => {
          if (notification_ids.includes(notification.id)) {
            return {
              ...notification,
              read_status,
              read_at: read_status ? new Date().toISOString() : notification.read_at
            };
          }
          return notification;
        });
        
        // Update unread count
        if (read_status) {
          // If marking as read, count how many were previously unread
          const unreadCount = state.notifications.filter(
            n => notification_ids.includes(n.id) && !n.read_status
          ).length;
          state.unreadCount = Math.max(0, state.unreadCount - unreadCount);
        } else {
          // If marking as unread, count how many were previously read
          const readCount = state.notifications.filter(
            n => notification_ids.includes(n.id) && n.read_status
          ).length;
          state.unreadCount += readCount;
        }
      })
      .addCase(bulkUpdateNotificationStatus.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to update notifications';
      });
  }
});

// Export actions
export const notificationsActions = notificationsSlice.actions;

// Export selectors
export const selectNotifications = (state: NotificationsRootState) => state.notifications.notifications;
export const selectNotificationsLoading = (state: NotificationsRootState) => state.notifications.loading;
export const selectNotificationsError = (state: NotificationsRootState) => state.notifications.error;
export const selectUnreadCount = (state: NotificationsRootState) => state.notifications.unreadCount;
export const selectTotalCount = (state: NotificationsRootState) => state.notifications.total;
export const selectCurrentPage = (state: NotificationsRootState) => state.notifications.page;

// Export reducer as default
export default notificationsSlice.reducer;