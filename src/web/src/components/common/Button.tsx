import React from 'react'; // React 18.2+
import { Button as MuiButton, ButtonProps, CircularProgress } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import theme from '../../theme'; // Import theme for consistent styling

/**
 * Props interface for the Button component extending Material-UI ButtonProps
 */
interface CustomButtonProps extends ButtonProps {
  loading?: boolean;  // Whether to show a loading spinner
  rounded?: boolean;  // Whether to use rounded corners
}

/**
 * A styled version of the Material-UI Button with custom styling
 */
const StyledButton = styled(MuiButton, {
  shouldForwardProp: (prop) => prop !== 'loading' && prop !== 'rounded',
})<CustomButtonProps>(({ theme, variant, size, loading, rounded }) => ({
  borderRadius: rounded ? '24px' : '8px',
  textTransform: 'none',
  fontWeight: 500,
  boxShadow: variant === 'contained' ? theme.shadows[2] : 'none',
  padding: size === 'small' ? '6px 16px' : size === 'medium' ? '8px 20px' : '10px 24px',
  position: 'relative',
  '&:hover': {
    boxShadow: variant === 'contained' ? theme.shadows[4] : 'none',
  },
  '&:disabled': {
    opacity: 0.7,
    pointerEvents: loading ? 'none' : 'auto',
  },
}));

/**
 * Props for the LoadingIndicator component
 */
interface LoadingIndicatorProps {
  variant?: 'text' | 'outlined' | 'contained';
  color?: 'inherit' | 'primary' | 'secondary' | 'success' | 'error' | 'info' | 'warning' | 'default';
  size?: number;
}

/**
 * A styled CircularProgress component for the loading indicator
 */
const LoadingIndicator = styled(CircularProgress)<LoadingIndicatorProps>(({ theme, variant, color = 'primary' }) => ({
  position: 'absolute',
  top: '50%',
  left: '50%',
  marginTop: '-12px',
  marginLeft: '-12px',
  color: variant === 'contained' ? 'white' : (color === 'inherit' ? 'currentColor' : theme.palette[color === 'default' ? 'primary' : color].main),
}));

/**
 * A customizable button component extending Material-UI Button with additional
 * styling options and features like loading state.
 */
const Button: React.FC<CustomButtonProps> = ({
  children,
  variant = 'contained',
  color = 'primary',
  size = 'medium',
  fullWidth = false,
  disabled = false,
  type = 'button',
  onClick,
  className,
  style,
  loading = false,
  startIcon,
  endIcon,
  href,
  target,
  rounded = false,
  ...props
}) => {
  return (
    <StyledButton
      variant={variant}
      color={color}
      size={size}
      fullWidth={fullWidth}
      disabled={disabled || loading}
      type={type}
      onClick={onClick}
      className={className}
      style={style}
      loading={loading}
      rounded={rounded}
      href={href}
      target={target}
      startIcon={startIcon}
      endIcon={endIcon}
      {...props}
    >
      {loading && <LoadingIndicator size={24} color={color} variant={variant} />}
      <span style={{ visibility: loading ? 'hidden' : 'visible' }}>
        {children}
      </span>
    </StyledButton>
  );
};

export default Button;