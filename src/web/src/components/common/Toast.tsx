import React, { useEffect, useState } from 'react'; // react v18.2+
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13+
import { Alert, AlertProps } from '@mui/material'; // @mui/material v5.13+
import { Snackbar, SnackbarProps } from '@mui/material'; // @mui/material v5.13+
import { IconButton } from '@mui/material'; // @mui/material v5.13+
import { Close as CloseIcon } from '@mui/icons-material'; // @mui/icons-material v5.13+
import theme from '../../theme';

/**
 * Interface defining the props for the Toast component
 */
interface ToastProps {
  /** Unique identifier for the toast */
  id: string;
  /** Type of toast notification that determines styling */
  type: 'success' | 'error' | 'warning' | 'info';
  /** Content message to display in the toast */
  message: string;
  /** Time in milliseconds before auto-dismissal */
  duration?: number;
  /** Callback function when toast is closed */
  onClose: (id: string) => void;
}

/**
 * Interface defining the props for the ToastContainer component
 */
interface ToastContainerProps {
  /** Array of toast notifications to display */
  toasts: Array<{
    id: string;
    type: 'success' | 'error' | 'warning' | 'info';
    message: string;
    duration?: number;
  }>;
  /** Callback function when a toast is closed */
  onClose: (id: string) => void;
}

/**
 * Helper function to get the appropriate background color based on toast type
 * @param severity Alert severity type
 * @param theme Current theme object
 * @returns CSS color value
 */
const getBackgroundColor = (severity: string | undefined, theme: typeof import('../../theme').default) => {
  switch (severity) {
    case 'success':
      return theme.palette.success.lighter;
    case 'error':
      return theme.palette.error.lighter;
    case 'warning':
      return theme.palette.warning.lighter;
    case 'info':
      return theme.palette.info.lighter;
    default:
      return theme.palette.grey[100];
  }
};

/**
 * Helper function to get the appropriate text color based on toast type
 * @param severity Alert severity type
 * @param theme Current theme object
 * @returns CSS color value
 */
const getTextColor = (severity: string | undefined, theme: typeof import('../../theme').default) => {
  switch (severity) {
    case 'success':
      return theme.palette.success.dark;
    case 'error':
      return theme.palette.error.dark;
    case 'warning':
      return theme.palette.warning.dark;
    case 'info':
      return theme.palette.info.dark;
    default:
      return theme.palette.text.primary;
  }
};

/**
 * Styled Alert component with custom styling based on severity
 */
const StyledAlert = styled(Alert)<AlertProps>(({ severity }) => ({
  width: '100%',
  boxShadow: theme.shadows[3],
  borderRadius: '4px',
  display: 'flex',
  alignItems: 'center',
  padding: '6px 16px',
  backgroundColor: getBackgroundColor(severity, theme),
  color: getTextColor(severity, theme),
  fontWeight: 500,
  fontSize: '0.875rem',
  lineHeight: '1.43',
  letterSpacing: '0.01071em',
}));

/**
 * Container for positioning multiple toasts
 */
const ToastWrapper = styled('div')({
  position: 'fixed',
  bottom: '24px',
  right: '24px',
  zIndex: 2000,
  display: 'flex',
  flexDirection: 'column-reverse',
  gap: '8px',
  maxWidth: 'calc(100% - 48px)',
  maxHeight: 'calc(100vh - 48px)',
  overflowY: 'auto',
  padding: '8px',
  paddingBottom: '0',
  '@media (max-width: 600px)': {
    bottom: '16px',
    right: '16px',
    maxWidth: 'calc(100% - 32px)',
  },
});

/**
 * A component that displays a toast notification with auto-dismissal
 * @param props Component props
 * @returns Toast component
 */
export const Toast: React.FC<ToastProps> = (props) => {
  const { id, type, message, duration = 5000, onClose } = props;
  const [open, setOpen] = useState<boolean>(true);

  // Auto-dismiss the toast after duration
  useEffect(() => {
    if (duration !== 0) {
      const timer = setTimeout(() => {
        setOpen(false);
      }, duration);
      return () => clearTimeout(timer);
    }
    return undefined;
  }, [duration]);

  /**
   * Handle closing the toast
   * @param event Close event
   * @param reason Reason for closing
   */
  const handleClose = (event?: React.SyntheticEvent | Event, reason?: string) => {
    if (reason === 'clickaway') {
      return;
    }
    setOpen(false);
    onClose(id);
  };

  return (
    <Snackbar
      open={open}
      onClose={handleClose}
      anchorOrigin={{ vertical: 'bottom', horizontal: 'right' }}
      key={id}
      autoHideDuration={duration}
      sx={{ position: 'static', mb: 1 }}
    >
      <StyledAlert
        severity={type}
        variant="filled"
        action={
          <IconButton
            size="small"
            aria-label="close"
            color="inherit"
            onClick={handleClose}
          >
            <CloseIcon fontSize="small" />
          </IconButton>
        }
      >
        {message}
      </StyledAlert>
    </Snackbar>
  );
};

/**
 * A container component that manages multiple toast notifications
 * @param props Component props
 * @returns ToastContainer component
 */
const ToastContainer: React.FC<ToastContainerProps> = (props) => {
  const { toasts, onClose } = props;

  if (!toasts.length) {
    return null;
  }

  return (
    <ToastWrapper>
      {toasts.map((toast) => (
        <Toast
          key={toast.id}
          id={toast.id}
          type={toast.type}
          message={toast.message}
          duration={toast.duration}
          onClose={onClose}
        />
      ))}
    </ToastWrapper>
  );
};

export default ToastContainer;