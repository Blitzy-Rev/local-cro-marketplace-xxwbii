import React from 'react'; // React v18.2+
import { styled } from '@mui/material/styles'; // v5.13+
import { Chip, Typography } from '@mui/material'; // v5.13+
import theme from '../../theme';

/**
 * Props interface for the Badge component
 */
interface BadgeProps {
  /** Text to display in the badge */
  label: string;
  /** Color scheme for the badge */
  color?: 'default' | 'primary' | 'secondary' | 'success' | 'error' | 'warning' | 'info';
  /** Size variation of the badge */
  size?: 'small' | 'medium' | 'large';
  /** Style variation of the badge */
  variant?: 'filled' | 'outlined';
  /** Custom inline styles */
  style?: React.CSSProperties;
  /** CSS class name */
  className?: string;
  /** Whether the badge has pill shape (rounded corners) */
  pill?: boolean;
  /** Click handler */
  onClick?: () => void;
}

/**
 * A styled Chip component for the Badge
 */
const StyledChip = styled(Chip, {
  shouldForwardProp: (prop) => 
    !['badgeColor', 'badgeSize', 'badgeVariant', 'pill'].includes(String(prop)),
})(({ 
  badgeColor, 
  badgeSize, 
  badgeVariant, 
  pill, 
  onClick 
}) => ({
  height: badgeSize === 'small' ? '24px' : badgeSize === 'medium' ? '32px' : '40px',
  fontSize: badgeSize === 'small' ? '0.75rem' : badgeSize === 'medium' ? '0.875rem' : '1rem',
  fontWeight: 500,
  borderRadius: pill ? '16px' : '4px',
  textTransform: 'none',
  cursor: onClick ? 'pointer' : 'default',
  backgroundColor: badgeVariant === 'filled' ? theme.palette[badgeColor].main : 'transparent',
  color: badgeVariant === 'filled' ? theme.palette[badgeColor].contrastText : theme.palette[badgeColor].main,
  border: badgeVariant === 'outlined' ? `1px solid ${theme.palette[badgeColor].main}` : 'none',
  '&:hover': {
    backgroundColor: onClick 
      ? (badgeVariant === 'filled' ? theme.palette[badgeColor].dark : theme.palette[badgeColor].lighter) 
      : undefined,
    opacity: onClick ? 0.9 : undefined,
  },
}));

/**
 * Badge component for displaying status, categories, or properties
 * 
 * This component provides a customizable badge with different colors, sizes,
 * variants, and shapes for indicating status, categories, or properties.
 * It's used throughout the application for visual labeling and status indication.
 *
 * @example
 * // Basic usage
 * <Badge label="Active" />
 * 
 * // With custom color and variant
 * <Badge label="Warning" color="warning" variant="outlined" />
 * 
 * // With click handler
 * <Badge label="Click Me" color="primary" onClick={() => console.log('Clicked!')} />
 */
const Badge: React.FC<BadgeProps> = ({
  label,
  color = 'default',
  size = 'small',
  variant = 'filled',
  style,
  className,
  pill = true,
  onClick,
}) => {
  // Convert size to Chip-compatible size
  const chipSize = size === 'large' ? 'medium' : size;

  return (
    <StyledChip
      label={label}
      size={chipSize}
      style={style}
      className={className}
      onClick={onClick}
      badgeColor={color}
      badgeSize={size}
      badgeVariant={variant}
      pill={pill}
    />
  );
};

export default Badge;