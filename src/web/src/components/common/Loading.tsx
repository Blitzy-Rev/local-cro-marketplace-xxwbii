import React from 'react';
import { CircularProgress, Box, Typography, styled, LinearProgress } from '@mui/material'; // v5.13+
import theme from '../../theme';

/**
 * Props interface for the Loading component
 */
interface LoadingProps {
  message?: string;
  size?: 'small' | 'medium' | 'large';
  color?: 'primary' | 'secondary' | 'success' | 'error' | 'info' | 'warning';
  fullScreen?: boolean;
  overlay?: boolean;
  variant?: 'circular' | 'linear';
  className?: string;
  style?: React.CSSProperties;
  thickness?: number;
  disableShrink?: boolean;
}

const LoadingContainer = styled(Box, {
  shouldForwardProp: (prop) => 
    prop !== 'fullScreen' && prop !== 'overlay',
})<{ fullScreen?: boolean; overlay?: boolean }>(({ theme, fullScreen, overlay }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  position: fullScreen || overlay ? 'fixed' : 'relative',
  top: fullScreen || overlay ? 0 : 'auto',
  left: fullScreen || overlay ? 0 : 'auto',
  right: fullScreen || overlay ? 0 : 'auto',
  bottom: fullScreen || overlay ? 0 : 'auto',
  width: fullScreen || overlay ? '100%' : '100%',
  height: fullScreen || overlay ? '100vh' : '100%',
  backgroundColor: overlay ? 'rgba(255, 255, 255, 0.8)' : 'transparent',
  zIndex: fullScreen || overlay ? 9999 : 1,
  padding: '16px',
}));

const StyledCircularProgress = styled(CircularProgress)(({ theme }) => ({
  // Base styling uses MUI defaults with props applied for customization
}));

const StyledLinearProgress = styled(LinearProgress)(({ theme }) => ({
  width: '100%',
  maxWidth: '300px',
  borderRadius: '4px',
}));

const LoadingMessage = styled(Typography)<{ customSize?: 'small' | 'medium' | 'large' }>(
  ({ theme, customSize = 'medium' }) => ({
    marginTop: '16px',
    color: theme.palette.text.secondary,
    fontSize: customSize === 'small' ? '0.875rem' : customSize === 'medium' ? '1rem' : '1.25rem',
    fontWeight: 500,
    textAlign: 'center',
  })
);

/**
 * A customizable loading component that displays a circular or linear progress 
 * indicator with optional text
 */
const Loading: React.FC<LoadingProps> = ({
  message = 'Loading...',
  size = 'medium',
  color = 'primary',
  fullScreen = false,
  overlay = false,
  variant = 'circular',
  className,
  style,
  thickness = 3.6,
  disableShrink = false,
}) => {
  // Calculate the numeric size for the circular progress
  const circularSize = size === 'small' ? 24 : size === 'medium' ? 40 : 56;
  
  return (
    <LoadingContainer 
      fullScreen={fullScreen} 
      overlay={overlay} 
      className={className} 
      style={style}
      role="status"
      aria-live="polite"
      aria-label={message || "Loading"}
    >
      {variant === 'circular' ? (
        <StyledCircularProgress 
          color={color}
          size={circularSize}
          thickness={thickness}
          disableShrink={disableShrink}
        />
      ) : (
        <StyledLinearProgress 
          color={color}
          sx={{
            height: size === 'small' ? 4 : size === 'medium' ? 6 : 8,
          }}
        />
      )}
      {message && (
        <LoadingMessage customSize={size} variant="body2">
          {message}
        </LoadingMessage>
      )}
    </LoadingContainer>
  );
};

export default Loading;