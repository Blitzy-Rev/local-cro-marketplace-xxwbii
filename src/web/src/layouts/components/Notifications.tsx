import React, { useState, useEffect, useCallback } from 'react'; // react v18.2+
import { useSelector, useDispatch } from 'react-redux'; // react-redux v8.0+
import {
  Menu,
  MenuItem,
  Typography,
  Box,
  IconButton,
  Divider,
  Button,
  CircularProgress,
  Tooltip,
  ListItemIcon,
  ListItemText,
} from '@mui/material'; // @mui/material v5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13+
import {
  MarkEmailReadIcon,
  DeleteIcon,
  NotificationsActiveIcon,
  NotificationsOffIcon,
  CheckCircleIcon,
  ErrorIcon,
  WarningIcon,
  InfoIcon,
  ExperimentIcon,
  LabIcon,
  AssignmentIcon,
} from '@mui/icons-material'; // @mui/icons-material v5.13+

import {
  selectNotifications,
  selectNotificationsLoading,
  selectUnreadCount,
  fetchNotifications,
  markAsRead,
  markAllAsRead,
  deleteNotificationById,
} from '../../store/notifications/notificationsSlice';
import { Notification, NotificationType } from '../../types/notification';
import { formatRelativeTime } from '../../utils/formatters';
import Badge from '../../components/common/Badge';
import theme from '../../theme';

/**
 * Props interface for the Notifications component
 */
interface NotificationsProps {
  anchorEl: HTMLElement | null;
  onClose: () => void;
}

/**
 * Styled components for the notification UI elements
 */
const NotificationMenu = styled(Menu)(({ theme }) => ({
  maxWidth: '400px',
  maxHeight: '500px',
  overflow: 'auto',
  marginTop: '8px',
  boxShadow: theme.shadows[3],
  borderRadius: '8px',
  padding: '0',
}));

const NotificationHeader = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  padding: '12px 16px',
  borderBottom: '1px solid rgba(0, 0, 0, 0.12)',
}));

const NotificationItem = styled(MenuItem, {
  shouldForwardProp: (prop) => prop !== 'read',
})<{ read: boolean }>(({ theme, read }) => ({
  padding: '12px 16px',
  borderBottom: '1px solid rgba(0, 0, 0, 0.12)',
  opacity: read ? 0.7 : 1,
  backgroundColor: read ? 'transparent' : 'rgba(0, 0, 0, 0.04)',
  '&:hover': {
    backgroundColor: 'rgba(0, 0, 0, 0.04)',
  },
}));

const NotificationContent = styled(Box)(({ theme }) => ({
  flex: '1',
  marginRight: '8px',
}));

const NotificationActions = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
}));

const NotificationFooter = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'center',
  padding: '12px 16px',
}));

const EmptyNotifications = styled(Box)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '32px 16px',
  color: 'rgba(0, 0, 0, 0.6)',
}));

/**
 * A dropdown component that displays user notifications with real-time updates,
 * read/unread status, and interactive actions. This component is used in the
 * application header to provide users with timely information about system
 * events, experiment status changes, and CRO interactions.
 */
const Notifications: React.FC<NotificationsProps> = ({ anchorEl, onClose }) => {
  const dispatch = useDispatch();
  const notifications = useSelector(selectNotifications);
  const loading = useSelector(selectNotificationsLoading);
  const unreadCount = useSelector(selectUnreadCount);
  
  // Check if menu is open
  const open = Boolean(anchorEl);
  
  // Fetch notifications when menu opens
  useEffect(() => {
    if (open) {
      dispatch(fetchNotifications({ page: 1, pageSize: 10 }));
    }
  }, [open, dispatch]);
  
  /**
   * Handles click on a notification item
   */
  const handleNotificationClick = useCallback((notification: Notification) => {
    // Mark notification as read if not already read
    if (!notification.read_status) {
      dispatch(markAsRead(notification.id));
    }
    
    // Navigate based on notification type and data
    // Implementation would depend on the routing system
    // For example:
    // if (notification.type === NotificationType.EXPERIMENT_STATUS_CHANGE) {
    //   history.push(`/experiments/${notification.data?.experiment_id}`);
    // } else if (notification.type === NotificationType.RESULTS_UPLOADED) {
    //   history.push(`/results/${notification.data?.result_id}`);
    // }
    
    onClose();
  }, [dispatch, onClose]);
  
  /**
   * Marks a notification as read
   */
  const handleMarkAsRead = useCallback((event: React.MouseEvent<HTMLButtonElement>, id: number) => {
    event.stopPropagation(); // Prevent menu from closing
    dispatch(markAsRead(id));
  }, [dispatch]);
  
  /**
   * Deletes a notification
   */
  const handleDeleteNotification = useCallback((event: React.MouseEvent<HTMLButtonElement>, id: number) => {
    event.stopPropagation(); // Prevent menu from closing
    dispatch(deleteNotificationById(id));
  }, [dispatch]);
  
  /**
   * Marks all notifications as read
   */
  const handleMarkAllAsRead = useCallback(() => {
    dispatch(markAllAsRead());
  }, [dispatch]);
  
  /**
   * Returns the appropriate icon for a notification type
   */
  const getNotificationIcon = (type: NotificationType): React.ReactNode => {
    switch (type) {
      case NotificationType.EXPERIMENT_STATUS_CHANGE:
        return <ExperimentIcon color="primary" />;
      case NotificationType.SUBMISSION_CREATED:
        return <AssignmentIcon color="primary" />;
      case NotificationType.QUOTE_PROVIDED:
        return <InfoIcon color="info" />;
      case NotificationType.RESULTS_UPLOADED:
        return <LabIcon color="success" />;
      case NotificationType.SYSTEM_ALERT:
        return <WarningIcon color="warning" />;
      case NotificationType.USER_MENTION:
        return <InfoIcon color="secondary" />;
      default:
        return <InfoIcon color="info" />;
    }
  };
  
  /**
   * Returns the appropriate color for a notification type
   */
  const getNotificationColor = (type: NotificationType): string => {
    switch (type) {
      case NotificationType.EXPERIMENT_STATUS_CHANGE:
        return 'primary';
      case NotificationType.SUBMISSION_CREATED:
        return 'primary';
      case NotificationType.QUOTE_PROVIDED:
        return 'info';
      case NotificationType.RESULTS_UPLOADED:
        return 'success';
      case NotificationType.SYSTEM_ALERT:
        return 'warning';
      case NotificationType.USER_MENTION:
        return 'secondary';
      default:
        return 'info';
    }
  };
  
  /**
   * Returns a human-readable label for a notification type
   */
  const getNotificationLabel = (type: NotificationType): string => {
    switch (type) {
      case NotificationType.EXPERIMENT_STATUS_CHANGE:
        return 'Experiment Update';
      case NotificationType.SUBMISSION_CREATED:
        return 'New Submission';
      case NotificationType.QUOTE_PROVIDED:
        return 'Quote Received';
      case NotificationType.RESULTS_UPLOADED:
        return 'Results Available';
      case NotificationType.SYSTEM_ALERT:
        return 'System Alert';
      case NotificationType.USER_MENTION:
        return 'Mention';
      default:
        return 'Notification';
    }
  };
  
  return (
    <NotificationMenu
      anchorEl={anchorEl}
      open={open}
      onClose={onClose}
      anchorOrigin={{
        vertical: 'bottom',
        horizontal: 'right',
      }}
      transformOrigin={{
        vertical: 'top',
        horizontal: 'right',
      }}
    >
      <NotificationHeader>
        <Typography variant="subtitle1" fontWeight={600}>
          Notifications {unreadCount > 0 && `(${unreadCount} unread)`}
        </Typography>
        {unreadCount > 0 && (
          <Tooltip title="Mark all as read">
            <IconButton size="small" onClick={handleMarkAllAsRead}>
              <MarkEmailReadIcon fontSize="small" />
            </IconButton>
          </Tooltip>
        )}
      </NotificationHeader>
      
      {loading ? (
        <Box sx={{ display: 'flex', justifyContent: 'center', padding: '24px' }}>
          <CircularProgress size={24} />
        </Box>
      ) : notifications.length === 0 ? (
        <EmptyNotifications>
          <NotificationsOffIcon sx={{ fontSize: 48, mb: 2, opacity: 0.6 }} />
          <Typography variant="body2">No notifications yet</Typography>
        </EmptyNotifications>
      ) : (
        <>
          {notifications.map((notification) => (
            <NotificationItem
              key={notification.id}
              read={notification.read_status}
              onClick={() => handleNotificationClick(notification)}
            >
              <ListItemIcon>
                {getNotificationIcon(notification.type)}
              </ListItemIcon>
              <NotificationContent>
                <Box sx={{ display: 'flex', alignItems: 'center', mb: 0.5 }}>
                  <Badge 
                    label={getNotificationLabel(notification.type)} 
                    color={getNotificationColor(notification.type) as any}
                    size="small"
                    style={{ marginRight: '8px' }}
                  />
                  <Typography variant="caption" color="text.secondary">
                    {formatRelativeTime(notification.created_at)}
                  </Typography>
                </Box>
                <Typography variant="body2">{notification.message}</Typography>
                {notification.data && notification.data.experiment_name && (
                  <Typography variant="caption" color="text.secondary">
                    {notification.data.experiment_name}
                  </Typography>
                )}
              </NotificationContent>
              <NotificationActions>
                {!notification.read_status && (
                  <Tooltip title="Mark as read">
                    <IconButton
                      size="small"
                      onClick={(e) => handleMarkAsRead(e, notification.id)}
                    >
                      <MarkEmailReadIcon fontSize="small" />
                    </IconButton>
                  </Tooltip>
                )}
                <Tooltip title="Delete">
                  <IconButton
                    size="small"
                    onClick={(e) => handleDeleteNotification(e, notification.id)}
                  >
                    <DeleteIcon fontSize="small" />
                  </IconButton>
                </Tooltip>
              </NotificationActions>
            </NotificationItem>
          ))}
          
          {unreadCount > 0 && (
            <NotificationFooter>
              <Button
                size="small"
                startIcon={<MarkEmailReadIcon />}
                onClick={handleMarkAllAsRead}
              >
                Mark all as read
              </Button>
            </NotificationFooter>
          )}
        </>
      )}
    </NotificationMenu>
  );
};

export default Notifications;