import React, { Component, ErrorInfo, ReactNode } from 'react'; // react 18.2+
import { styled } from '@mui/material/styles'; // v5.13+
import { Box, Typography, Paper } from '@mui/material'; // v5.13+
import { ErrorOutline } from '@mui/icons-material'; // v5.13+
import { useToast } from '../../hooks/useToast';
import Button from './Button';

/**
 * Props interface for the ErrorBoundary component
 */
interface ErrorBoundaryProps {
  children: ReactNode;
  FallbackComponent?: React.ComponentType<FallbackProps>;
  onError?: (error: Error, errorInfo: ErrorInfo) => void;
  onReset?: () => void;
  resetKeys?: Array<unknown>;
  logError?: (error: Error, errorInfo: ErrorInfo) => void;
  showNotification?: boolean;
}

/**
 * Props interface for the fallback component
 */
interface FallbackProps {
  error: Error;
  resetErrorBoundary: () => void;
}

/**
 * State interface for the ErrorBoundary component
 */
interface ErrorBoundaryState {
  hasError: boolean;
  error: Error | null;
}

/**
 * Styled component for the error container
 */
const ErrorContainer = styled(Paper)(({ theme }) => ({
  padding: '2rem',
  margin: '1rem',
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  textAlign: 'center',
  borderRadius: '8px',
  backgroundColor: theme.palette.error.light,
  color: theme.palette.error.contrastText,
}));

/**
 * Styled component for the error icon
 */
const ErrorIcon = styled(ErrorOutline)(({ theme }) => ({
  fontSize: '4rem',
  color: theme.palette.error.main,
  marginBottom: '1rem',
}));

/**
 * Styled component for the error message
 */
const ErrorMessage = styled(Typography)(({ theme }) => ({
  marginBottom: '1rem',
  color: theme.palette.error.dark,
  fontWeight: '500',
}));

/**
 * Styled component for the error details
 */
const ErrorDetails = styled(Box)(({ theme }) => ({
  backgroundColor: 'rgba(0, 0, 0, 0.05)',
  padding: '1rem',
  borderRadius: '4px',
  marginBottom: '1rem',
  maxWidth: '100%',
  overflow: 'auto',
  fontFamily: 'monospace',
  fontSize: '0.875rem',
  color: theme.palette.text.secondary,
}));

/**
 * Default fallback UI component displayed when an error occurs
 */
const DefaultFallbackComponent: React.FC<FallbackProps> = ({ error, resetErrorBoundary }) => {
  return (
    <ErrorContainer>
      <ErrorIcon />
      <ErrorMessage variant="h6">
        Oops! Something went wrong.
      </ErrorMessage>
      <ErrorMessage variant="body2">
        {error.name}: {error.message}
      </ErrorMessage>
      {process.env.NODE_ENV === 'development' && (
        <ErrorDetails>
          {error.stack}
        </ErrorDetails>
      )}
      <Button onClick={resetErrorBoundary}>Try again</Button>
    </ErrorContainer>
  );
};

/**
 * A class component that catches errors in its child component tree and displays a fallback UI
 */
class ErrorBoundary extends Component<ErrorBoundaryProps, ErrorBoundaryState> {
  static getDerivedStateFromError(error: Error) {
    // Update state so the next render will show the fallback UI.
    return { hasError: true, error: error };
  }

  constructor(props: ErrorBoundaryProps) {
    super(props);
    this.state = { hasError: false, error: null };
    this.resetErrorBoundary = this.resetErrorBoundary.bind(this);
  }

  state: ErrorBoundaryState = {
    hasError: false,
    error: null,
  };

  componentDidCatch(error: Error, errorInfo: ErrorInfo) {
    // You can also log the error to an error reporting service
    console.error('Caught error: ', error, errorInfo);
    if (this.props.onError) {
      this.props.onError(error, errorInfo);
    }
    if (this.props.logError) {
      this.props.logError(error, errorInfo);
    }
  }

  resetErrorBoundary() {
    this.setState({ hasError: false, error: null });
    if (this.props.onReset) {
      this.props.onReset();
    }
  }

  render() {
    if (this.state.hasError) {
      // You can render any custom fallback UI
      const FallbackComponent = this.props.FallbackComponent || DefaultFallbackComponent;
      return <FallbackComponent error={this.state.error as Error} resetErrorBoundary={this.resetErrorBoundary} />;
    }

    return this.props.children;
  }
}

export default ErrorBoundary;