import React from 'react'; // react 18.2+
import { Navigate, Outlet, useLocation } from 'react-router-dom'; // react-router-dom 6.4+
import { Box, Container, CssBaseline, Paper, Typography, styled } from '@mui/material'; // @mui/material 5.13+
import theme from '../theme'; // src/web/src/theme/index.ts
import useAuth from '../features/auth/hooks/useAuth'; // src/web/src/features/auth/hooks/useAuth.ts

/**
 * Interface for AuthLayout props (currently empty)
 */
interface AuthLayoutProps {}

/**
 * Styled component for the main authentication container
 */
const AuthContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  minHeight: '100vh',
  backgroundColor: theme.palette.background.default,
});

/**
 * Styled component for the authentication content container
 */
const AuthContent = styled(Container)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  flex: '1',
  padding: theme.spacing(3),
});

/**
 * Styled component for the authentication card
 */
const AuthCard = styled(Paper)({
  padding: theme.spacing(4),
  width: '100%',
  maxWidth: '450px',
  borderRadius: theme.shape.borderRadius,
  boxShadow: theme.shadows[3],
});

/**
 * Styled component for the logo container
 */
const LogoContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'center',
  marginBottom: theme.spacing(4),
});

/**
 * Layout component for authentication-related pages (login, register, etc.)
 *
 * This component provides a consistent layout for authentication pages,
 * including a background, logo, and centered content area.
 */
const AuthLayout: React.FC<AuthLayoutProps> = () => {
  // Use the useAuth hook to get authentication status
  const { isAuthenticated, user } = useAuth();

  // Use the useLocation hook to get the current route location
  const location = useLocation();

  // If the user is already authenticated, redirect them to the appropriate dashboard
  if (isAuthenticated) {
    let redirectPath = '/'; // Default redirect path
    if (user?.role === 'pharma') {
      redirectPath = '/molecules';
    } else if (user?.role === 'cro') {
      redirectPath = '/submissions';
    } else if (user?.role === 'admin') {
      redirectPath = '/admin';
    }
    return <Navigate to={redirectPath} state={{ from: location }} replace />;
  }

  // Render the authentication layout
  return (
    <AuthContainer>
      <CssBaseline />
      <AuthContent maxWidth="sm">
        <AuthCard>
          <LogoContainer>
            <Typography variant="h4" component="h1" gutterBottom>
              Molecular Data Platform
            </Typography>
          </LogoContainer>
          <Outlet />
        </AuthCard>
      </AuthContent>
    </AuthContainer>
  );
};

export default AuthLayout;