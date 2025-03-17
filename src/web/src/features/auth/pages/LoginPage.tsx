import React from 'react'; // React 18.2+
import { Navigate, useNavigate } from 'react-router-dom'; // react-router-dom ^6.10.0
import { Box, Typography } from '@mui/material'; // @mui/material ^5.13.0
import LoginForm from '../components/LoginForm';
import useAuth from '../hooks/useAuth';
import { UserRole } from '../../../types/user';

/**
 * Component that renders the login page and handles redirects based on authentication status
 * @returns {JSX.Element} The rendered login page component or a redirect if already authenticated
 */
const LoginPage: React.FC = () => {
  // LD1: Get authentication state and user information using useAuth hook
  const { isAuthenticated, user } = useAuth();

  // LD1: Get navigate function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Check if user is already authenticated
  if (isAuthenticated) {
    // LD1: If authenticated, redirect to the appropriate dashboard based on user role
    let redirectPath = '/app/dashboard'; // Default path

    // LD1: Check user role from authentication state
    if (user) {
      // LD1: Navigate to the appropriate dashboard based on role
      switch (user.role) {
        // LD1: For PHARMA users, navigate to /app/dashboard
        case UserRole.PHARMA:
          redirectPath = '/app/dashboard';
          break;
        // LD1: For CRO users, navigate to /cro/dashboard
        case UserRole.CRO:
          redirectPath = '/cro/dashboard';
          break;
        // LD1: For ADMIN users, navigate to /admin/dashboard
        case UserRole.ADMIN:
          redirectPath = '/admin/dashboard';
          break;
        // LD1: For unknown roles, default to /app/dashboard
        default:
          redirectPath = '/app/dashboard';
          break;
      }
    }

    // LD1: Return a Navigate component to redirect
    return <Navigate to={redirectPath} replace />;
  }

  // LD1: If not authenticated, render the login form
  return (
    <Box
      sx={{
        display: 'flex',
        flexDirection: 'column',
        alignItems: 'center',
        justifyContent: 'center',
        height: '100vh',
        bgcolor: 'background.default',
      }}
    >
      <Typography variant="h4" component="h1" gutterBottom>
        Welcome to Molecular Data Platform
      </Typography>
      {/* LD1: Render the login form with the handleLoginSuccess callback */}
      <LoginForm onSuccess={() => {
        // LD1: Check user role from authentication state
        if (user) {
          // LD1: Navigate to the appropriate dashboard based on role
          switch (user.role) {
            // LD1: For PHARMA users, navigate to /app/dashboard
            case UserRole.PHARMA:
              navigate('/app/dashboard');
              break;
            // LD1: For CRO users, navigate to /cro/dashboard
            case UserRole.CRO:
              navigate('/cro/dashboard');
              break;
            // LD1: For ADMIN users, navigate to /admin/dashboard
            case UserRole.ADMIN:
              navigate('/admin/dashboard');
              break;
            // LD1: For unknown roles, default to /app/dashboard
            default:
              navigate('/app/dashboard');
              break;
          }
        }
      }} />
    </Box>
  );
};

// IE3: Export the LoginPage component for use in routing
export default LoginPage;