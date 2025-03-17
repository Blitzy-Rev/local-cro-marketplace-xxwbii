import React from 'react'; // React 18.2+
import { useNavigate, useParams } from 'react-router-dom'; // react-router-dom 6.10+
import { Box, Typography } from '@mui/material'; // @mui/material 5.13+
import AuthLayout from '../../../layouts/AuthLayout'; // Layout component for authentication pages
import ResetPasswordForm from '../components/ResetPasswordForm'; // Form component for resetting password
import useToast from '../../../hooks/useToast'; // Hook for displaying toast notifications

/**
 * Page component for password reset functionality
 */
const ResetPasswordPage: React.FC = () => {
  // Extract token parameter from URL using useParams hook
  const { token } = useParams<{ token: string }>();

  // Initialize navigate function using useNavigate hook
  const navigate = useNavigate();

  // Initialize showToast function from useToast hook
  const { showToast } = useToast();

  /**
   * Define handleSuccess function to handle successful password reset
   */
  const handleSuccess = () => {
    showToast({ type: 'success', message: 'Password reset successful!' });
    setTimeout(() => {
      navigate('/login');
    }, 3000);
  };

  // Render AuthLayout component as the container
  return (
    <AuthLayout>
      {/* Render Typography component for the page title */}
      <Typography variant="h4" align="center" gutterBottom>
        Reset Password
      </Typography>
      {/* Render ResetPasswordForm component with token and onSuccess props */}
      {token ? (
        <ResetPasswordForm token={token} onSuccess={handleSuccess} />
      ) : (
        // If token is missing, display an error message
        <Box textAlign="center" mt={3}>
          <Typography color="error">
            Invalid or missing reset token. Please check your email.
          </Typography>
        </Box>
      )}
    </AuthLayout>
  );
};

// Export the ResetPasswordPage component for use in routing
export default ResetPasswordPage;