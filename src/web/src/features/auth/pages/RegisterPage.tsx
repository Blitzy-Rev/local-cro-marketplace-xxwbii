import React from 'react'; // React v18.2+
import { Navigate, useNavigate } from 'react-router-dom'; // React Router v6.10+
import { Box, Typography } from '@mui/material'; // Material-UI v5.13+

// Internal imports
import RegisterForm from '../components/RegisterForm';
import useAuth from '../hooks/useAuth';
import { UserRole } from '../../../types/user';

/**
 * Component that renders the registration page and handles redirects based on authentication status
 * @returns The rendered registration page component or a redirect if already authenticated
 */
const RegisterPage: React.FC = () => {
  // Get authentication state and user information using useAuth hook
  const { isAuthenticated, user } = useAuth();

  // Get navigate function from useNavigate hook
  const navigate = useNavigate();

  // Check if user is already authenticated
  if (isAuthenticated) {
    // If authenticated, redirect to the appropriate dashboard based on user role
    let redirectPath: string;
    switch (user?.role) {
      case UserRole.PHARMA:
        redirectPath = '/molecules';
        break;
      case UserRole.CRO:
        redirectPath = '/submissions';
        break;
      case UserRole.ADMIN:
        redirectPath = '/admin';
        break;
      default:
        redirectPath = '/'; // Default redirect
    }
    return <Navigate to={redirectPath} replace />;
  }

  /**
   * Handles successful registration by redirecting to the login page
   */
  const handleRegistrationSuccess = () => {
    // Navigate to the login page
    navigate('/login');
  };

  /**
   * Handles click on the login link
   */
  const handleLoginClick = () => {
    // Navigate to the login page
    navigate('/login');
  };

  // If not authenticated, render the registration form
  return (
    <Box
      display="flex"
      flexDirection="column"
      alignItems="center"
      justifyContent="center"
      minHeight="100vh"
      bgcolor="background.default"
      color="text.primary"
      padding={3}
    >
      <Typography variant="h4" component="h1" gutterBottom>
        Create an Account
      </Typography>
      <RegisterForm onSuccess={handleRegistrationSuccess} onLoginClick={handleLoginClick} />
    </Box>
  );
};

export default RegisterPage;