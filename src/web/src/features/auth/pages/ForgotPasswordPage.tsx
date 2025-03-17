import React from 'react'; // react 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom ^6.10.0
import { Box, Typography } from '@mui/material'; // @mui/material ^5.13.0

import AuthLayout from '../../../layouts/AuthLayout';
import ForgotPasswordForm from '../components/ForgotPasswordForm';
import useAuth from '../hooks/useAuth';

/**
 * Component that renders the forgot password page
 * @returns {JSX.Element} The rendered forgot password page
 */
const ForgotPasswordPage: React.FC = () => {
  // LD1: Get navigation function using useNavigate hook
  const navigate = useNavigate();

  // LD1: Get authentication state using useAuth hook
  const {  } = useAuth();

  // LD1: Define handleSuccess function to navigate to login page after successful submission
  const handleSuccess = () => {
    // LD1: Navigate to the login page after successful submission
    navigate('/login');
  };

  // LD1: Render AuthLayout as the container component
  return (
    <AuthLayout>
      {/* LD1: Render Typography component for the page title */}
      <Typography variant="h4" component="h1" align="center" gutterBottom>
        Forgot Password
      </Typography>

      {/* LD1: Render Typography component for the page description */}
      <Typography variant="body1" align="center" sx={{ mb: 3 }}>
        Enter your email address below and we'll send you a link to reset your password.
      </Typography>

      {/* LD1: Render ForgotPasswordForm component with handleSuccess callback */}
      <ForgotPasswordForm onSuccess={handleSuccess} />
    </AuthLayout>
  );
  // LD1: Return the JSX structure
};

// IE3: Export the ForgotPasswordPage component for use in routing
export default ForgotPasswordPage;