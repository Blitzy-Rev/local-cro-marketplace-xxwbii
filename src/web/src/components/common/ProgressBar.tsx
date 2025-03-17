import React from 'react'; // React 18.2+
import { LinearProgress, Box, Typography } from '@mui/material'; // Material-UI 5.13+
import { styled } from '@mui/material/styles'; // Material-UI 5.13+
import theme from '../../theme'; // Import theme for consistent styling

/**
 * Props interface for the ProgressBar component
 */
interface ProgressBarProps {
  /** Current progress value (0-100) */
  value?: number;
  /** Type of progress bar to display */
  variant?: 'determinate' | 'indeterminate' | 'buffer';
  /** Color of the progress bar */
  color?: 'primary' | 'secondary' | 'success' | 'error' | 'warning' | 'info';
  /** Size of the progress bar */
  size?: 'small' | 'medium' | 'large';
  /** Whether to show percentage text */
  showPercentage?: boolean;
  /** Optional label to display above the progress bar */
  label?: string;
  /** Additional CSS class name */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
  /** Buffer value for buffer variant (0-100) */
  bufferValue?: number;
  /** Whether to use rounded corners */
  rounded?: boolean;
}

// Styled LinearProgress component with customizations
const StyledLinearProgress = styled(LinearProgress, {
  shouldForwardProp: (prop) => prop !== 'size' && prop !== 'rounded',
})<{ size?: 'small' | 'medium' | 'large'; rounded?: boolean; color: string }>(
  ({ size, rounded, color }) => ({
    height: size === 'small' ? 4 : size === 'medium' ? 8 : 12,
    borderRadius: rounded
      ? size === 'small'
        ? 2
        : size === 'medium'
        ? 4
        : 6
      : 0,
    backgroundColor:
      color === 'primary'
        ? theme.palette.primary.light
        : color === 'secondary'
        ? theme.palette.secondary.light
        : theme.palette[color].light,
    '.MuiLinearProgress-bar': {
      backgroundColor:
        color === 'primary'
          ? theme.palette.primary.main
          : color === 'secondary'
          ? theme.palette.secondary.main
          : theme.palette[color].main,
    },
  })
);

// Container for the progress bar and label
const ProgressContainer = styled(Box)<{ label?: string }>(({ label }) => ({
  width: '100%',
  marginBottom: label ? 8 : 0,
  marginTop: label ? 8 : 0,
}));

// Typography component for the label
const LabelText = styled(Typography)<{ size?: 'small' | 'medium' | 'large' }>(({ size }) => ({
  fontSize: size === 'small' ? '0.75rem' : size === 'medium' ? '0.875rem' : '1rem',
  fontWeight: 500,
  marginBottom: '4px',
  color: theme.palette.text.secondary,
}));

// Typography component for the percentage text
const PercentageText = styled(Typography)<{ size?: 'small' | 'medium' | 'large' }>(({ size }) => ({
  fontSize: size === 'small' ? '0.75rem' : size === 'medium' ? '0.875rem' : '1rem',
  fontWeight: 500,
  marginLeft: '8px',
  color: theme.palette.text.secondary,
}));

// Wrapper for the progress bar and percentage text
const ProgressWrapper = styled(Box)({
  display: 'flex',
  alignItems: 'center',
});

/**
 * A customizable progress bar component for visualizing operation progress
 * used for file uploads, data processing, and experiment status tracking
 */
const ProgressBar: React.FC<ProgressBarProps> = ({
  value = 0,
  variant = 'determinate',
  color = 'primary',
  size = 'medium',
  showPercentage = false,
  label,
  className,
  style,
  bufferValue = 0,
  rounded = true,
}) => {
  // Ensure value is between 0 and 100
  const normalizedValue = Math.min(Math.max(value, 0), 100);
  const normalizedBufferValue = Math.min(Math.max(bufferValue, 0), 100);

  return (
    <ProgressContainer className={className} style={style} label={label}>
      {label && <LabelText size={size}>{label}</LabelText>}
      <ProgressWrapper>
        <StyledLinearProgress
          variant={variant}
          value={normalizedValue}
          valueBuffer={normalizedBufferValue}
          size={size}
          rounded={rounded}
          color={color}
          sx={{ flexGrow: 1 }}
        />
        {showPercentage && variant === 'determinate' && (
          <PercentageText size={size}>{`${Math.round(normalizedValue)}%`}</PercentageText>
        )}
      </ProgressWrapper>
    </ProgressContainer>
  );
};

export default ProgressBar;