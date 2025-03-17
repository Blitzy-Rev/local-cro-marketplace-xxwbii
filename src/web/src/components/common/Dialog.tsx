import React from 'react'; // React 18.2+
import {
  Dialog as MuiDialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  IconButton,
  Typography,
  Slide,
} from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import CloseIcon from '@mui/icons-material/Close'; // v5.13+
import { TransitionProps } from '@mui/material'; // v5.13+
import theme from '../../theme'; // Import theme for consistent styling
import Button from './Button'; // Reusable button component for dialog actions

// Custom slide transition for dialogs
const Transition = React.forwardRef(function Transition(
  props: React.PropsWithChildren<TransitionProps>,
  ref: React.Ref<unknown>,
) {
  return <Slide direction="up" ref={ref} {...props} />;
});

// Props interface for the Dialog component
interface CustomDialogProps {
  open: boolean;
  onClose: () => void;
  title: string;
  children: React.ReactNode;
  contentText?: string;
  maxWidth?: 'xs' | 'sm' | 'md' | 'lg' | 'xl' | false;
  fullWidth?: boolean;
  fullScreen?: boolean;
  disableBackdropClick?: boolean;
  disableEscapeKeyDown?: boolean;
  showCloseButton?: boolean;
  actions?: React.ReactNode;
  confirmButtonText?: string;
  cancelButtonText?: string;
  onConfirm?: () => void;
  onCancel?: () => void;
  loading?: boolean;
  className?: string;
  style?: React.CSSProperties;
  hideActions?: boolean;
}

// Styled components for dialog elements
const StyledDialog = styled(MuiDialog)(({ theme }) => ({
  '& .MuiDialog-paper': {
    borderRadius: '8px',
    boxShadow: theme.shadows[5],
  }
}));

const StyledDialogTitle = styled(DialogTitle)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  padding: '16px 24px',
  backgroundColor: theme.palette.background.paper,
  borderBottom: '1px solid rgba(0, 0, 0, 0.12)'
}));

const StyledDialogContent = styled(DialogContent)({
  padding: '24px',
  minHeight: '120px'
});

const StyledDialogActions = styled(DialogActions)(({ theme }) => ({
  padding: '16px 24px',
  borderTop: '1px solid rgba(0, 0, 0, 0.12)'
}));

const CloseButton = styled(IconButton)(({ theme }) => ({
  position: 'absolute',
  right: '8px',
  top: '8px',
  color: theme.palette.grey[500]
}));

/**
 * A customizable dialog component extending Material-UI Dialog with
 * additional styling options and features.
 */
const CustomDialog: React.FC<CustomDialogProps> = ({
  open,
  onClose,
  title,
  children,
  contentText,
  maxWidth = 'sm',
  fullWidth = true,
  fullScreen = false,
  disableBackdropClick = false,
  disableEscapeKeyDown = false,
  showCloseButton = true,
  actions,
  confirmButtonText = 'Confirm',
  cancelButtonText = 'Cancel',
  onConfirm,
  onCancel,
  loading = false,
  className,
  style,
  hideActions = false,
}) => {
  // Handle backdrop click based on disableBackdropClick prop
  const handleBackdropClick = (event: React.MouseEvent<HTMLDivElement>) => {
    if (disableBackdropClick) {
      event.stopPropagation();
    } else {
      onClose();
    }
  };

  return (
    <StyledDialog
      open={open}
      onClose={onClose}
      maxWidth={maxWidth}
      fullWidth={fullWidth}
      fullScreen={fullScreen}
      TransitionComponent={Transition}
      disableEscapeKeyDown={disableEscapeKeyDown}
      onClick={handleBackdropClick}
      className={className}
      style={style}
      aria-labelledby="dialog-title"
      aria-describedby={contentText ? "dialog-description" : undefined}
    >
      <StyledDialogTitle id="dialog-title">
        <Typography variant="h6" component="div">
          {title}
        </Typography>
        {showCloseButton && (
          <CloseButton aria-label="close" onClick={onClose}>
            <CloseIcon />
          </CloseButton>
        )}
      </StyledDialogTitle>
      
      <StyledDialogContent>
        {contentText && (
          <DialogContentText id="dialog-description" gutterBottom>
            {contentText}
          </DialogContentText>
        )}
        {children}
      </StyledDialogContent>
      
      {!hideActions && (
        <StyledDialogActions>
          {actions ? (
            actions
          ) : (
            <>
              {onCancel && (
                <Button 
                  variant="outlined" 
                  color="primary" 
                  onClick={onCancel}
                  disabled={loading}
                >
                  {cancelButtonText}
                </Button>
              )}
              {onConfirm && (
                <Button 
                  variant="contained" 
                  color="primary" 
                  onClick={onConfirm}
                  loading={loading}
                  autoFocus
                >
                  {confirmButtonText}
                </Button>
              )}
            </>
          )}
        </StyledDialogActions>
      )}
    </StyledDialog>
  );
};

export default CustomDialog;