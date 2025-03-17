import React from 'react'; // react v18.2+
import { Card as MuiCard, CardProps, CardContent, CardActions, Box } from '@mui/material'; // @mui/material v5.13+
import { styled, alpha } from '@mui/material/styles'; // @mui/material/styles v5.13+
import theme from '../../theme';

/**
 * Props interface for the Card component extending Material-UI CardProps
 */
interface CustomCardProps extends CardProps {
  children?: React.ReactNode;
  selected?: boolean;
  interactive?: boolean;
  elevation?: number;
  noPadding?: boolean;
  outlined?: boolean;
  variant?: 'elevation' | 'outlined';
  className?: string;
  style?: React.CSSProperties;
  onClick?: () => void;
  header?: React.ReactNode;
  footer?: React.ReactNode;
  fullHeight?: boolean;
  disableHoverEffect?: boolean;
}

// Styled card component with customizable appearance
const StyledCard = styled(MuiCard, {
  shouldForwardProp: (prop) => 
    prop !== 'selected' && 
    prop !== 'interactive' && 
    prop !== 'noPadding' && 
    prop !== 'outlined' && 
    prop !== 'fullHeight' && 
    prop !== 'disableHoverEffect'
})<{
  selected?: boolean;
  interactive?: boolean;
  elevation?: number;
  outlined?: boolean;
  variant?: 'elevation' | 'outlined';
  fullHeight?: boolean;
  disableHoverEffect?: boolean;
}>(({ selected, interactive, elevation = 1, outlined, variant, fullHeight, disableHoverEffect }) => ({
  position: 'relative',
  transition: 'all 0.2s ease-in-out',
  height: fullHeight ? '100%' : 'auto',
  display: 'flex',
  flexDirection: 'column',
  border: selected 
    ? `2px solid ${theme.palette.primary.main}` 
    : outlined 
      ? `1px solid ${theme.palette.divider}` 
      : 'none',
  backgroundColor: selected ? alpha(theme.palette.primary.main, 0.08) : 'inherit',
  boxShadow: variant === 'elevation' ? theme.shadows[elevation] : 'none',
  cursor: interactive ? 'pointer' : 'default',
  '&:hover': {
    transform: interactive && !disableHoverEffect ? 'translateY(-4px)' : 'none',
    boxShadow: interactive && !disableHoverEffect && variant === 'elevation' 
      ? theme.shadows[elevation + 2] 
      : variant === 'elevation' 
        ? theme.shadows[elevation] 
        : 'none',
  }
}));

// Custom card content with adjustable padding
const CardContentWrapper = styled(CardContent, {
  shouldForwardProp: (prop) => prop !== 'noPadding'
})<{ noPadding?: boolean }>(({ noPadding }) => ({
  padding: noPadding ? 0 : 16,
  flex: '1',
  '&:last-child': {
    paddingBottom: noPadding ? 0 : 16
  }
}));

// Custom header component with optional border
const CardHeaderWrapper = styled(Box, {
  shouldForwardProp: (prop) => prop !== 'noPadding' && prop !== 'header'
})<{ noPadding?: boolean; header?: React.ReactNode }>(({ noPadding, header }) => ({
  padding: noPadding ? 0 : '16px 16px 0 16px',
  borderBottom: header ? `1px solid ${theme.palette.divider}` : 'none'
}));

// Custom footer component with optional border
const CardFooterWrapper = styled(CardActions, {
  shouldForwardProp: (prop) => prop !== 'noPadding' && prop !== 'footer'
})<{ noPadding?: boolean; footer?: React.ReactNode }>(({ noPadding, footer }) => ({
  padding: noPadding ? 0 : 16,
  borderTop: footer ? `1px solid ${theme.palette.divider}` : 'none'
}));

/**
 * A customizable card component extending Material-UI Card with additional styling options and features.
 * Supports interactive states, custom header/footer sections, and responsive layout options.
 */
const Card: React.FC<CustomCardProps> = ({
  children,
  selected = false,
  interactive = false,
  elevation = 1,
  noPadding = false,
  outlined = false,
  variant = 'elevation',
  className,
  style,
  onClick,
  header,
  footer,
  fullHeight = false,
  disableHoverEffect = false,
  ...rest
}) => {
  return (
    <StyledCard
      selected={selected}
      interactive={interactive}
      elevation={elevation}
      outlined={outlined}
      variant={variant}
      className={className}
      style={style}
      onClick={interactive ? onClick : undefined}
      fullHeight={fullHeight}
      disableHoverEffect={disableHoverEffect}
      {...rest}
    >
      {header && (
        <CardHeaderWrapper noPadding={noPadding} header={header}>
          {header}
        </CardHeaderWrapper>
      )}
      <CardContentWrapper noPadding={noPadding}>
        {children}
      </CardContentWrapper>
      {footer && (
        <CardFooterWrapper noPadding={noPadding} footer={footer}>
          {footer}
        </CardFooterWrapper>
      )}
    </StyledCard>
  );
};

export default Card;