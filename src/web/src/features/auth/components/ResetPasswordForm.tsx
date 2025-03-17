import React, { useState } from 'react'; // React 18.2+
import { useNavigate, useParams } from 'react-router-dom'; // react-router-dom 6.10+
import { Box, Typography, Paper } from '@mui/material'; // @mui/material 5.13+
import Input from '../../../components/common/Input'; // Input component for form fields
import Button from '../../../components/common/Button'; // Button component for form submission
import { useToast } from '../../../hooks/useToast'; // Hook for displaying toast notifications
import { resetPasswordConfirm } from '../../../api/auth'; // API function to confirm password reset
import { PasswordResetConfirmRequest } from '../../../types/auth'; // Type definition for password reset confirmation request

/**
 * Props interface for the ResetPasswordForm component
 */
interface ResetPasswordFormProps {
  /** The reset token */
  token?: string;
  /** Callback function to execute on successful password reset */
  onSuccess?: () => void;
}

/**
 * A form component that allows users to reset their password using a reset token
 */
const ResetPasswordForm: React.FC<ResetPasswordFormProps> = ({
  token = '',
  onSuccess = () => {},
}) => {
  // Component state
  const [newPassword, setNewPassword] = useState(''); // State for the new password input
  const [confirmPassword, setConfirmPassword] = useState(''); // State for the confirm password input
  const [errors, setErrors] = useState<{ newPassword: string; confirmPassword: string }>({ newPassword: '', confirmPassword: '' }); // State for form validation errors
  const [loading, setLoading] = useState(false); // State for tracking form submission status
  const [success, setSuccess] = useState(false); // State for tracking successful password reset

  // Hooks
  const navigate = useNavigate(); // Hook for programmatic navigation
  const { id: paramId, token: paramToken } = useParams<{ id: string; token: string }>(); // Hook to access URL parameters
  const { showToast } = useToast(); // Hook for displaying toast notifications

  /**
   * Validates the form inputs
   * @returns Whether the form is valid
   */
  const validateForm = (): boolean => {
    let valid = true; // Initialize a valid flag as true
    const newErrors: { newPassword: string; confirmPassword: string } = { newPassword: '', confirmPassword: '' }; // Create a new errors object

    if (!newPassword) { // If newPassword is empty, set error and set valid to false
      newErrors.newPassword = 'New password is required';
      valid = false;
    } else if (newPassword.length < 10) { // If newPassword is less than 10 characters, set error and set valid to false
      newErrors.newPassword = 'Password must be at least 10 characters';
      valid = false;
    } else if (!/^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[!@#$%^&*()_+{}\[\]:;<>,.?~\\/-]).*$/.test(newPassword)) { // If newPassword doesn't meet complexity requirements, set error and set valid to false
      newErrors.newPassword = 'Password must contain at least one uppercase letter, one lowercase letter, one number, and one special character';
      valid = false;
    }

    if (!confirmPassword) { // If confirmPassword is empty, set error and set valid to false
      newErrors.confirmPassword = 'Confirm password is required';
      valid = false;
    } else if (newPassword !== confirmPassword) { // If confirmPassword doesn't match newPassword, set error and set valid to false
      newErrors.confirmPassword = 'Passwords do not match';
      valid = false;
    }

    setErrors(newErrors); // Update the errors state
    return valid; // Return the valid flag
  };

  /**
   * Handles form submission
   * @param event React.FormEvent
   */
  const handleSubmit = async (event: React.FormEvent): Promise<void> => {
    event.preventDefault(); // Prevent default form submission

    if (!validateForm()) { // If form is not valid, return early
      return;
    }

    setLoading(true); // Set loading state to true

    const tokenValue = token || paramToken; // Get token from props or URL params

    const resetData: PasswordResetConfirmRequest = { // Create reset data object with token, new_password, and confirm_password
      token: tokenValue,
      new_password: newPassword,
      confirm_password: confirmPassword,
    };

    try {
      await resetPasswordConfirm(resetData); // Try to call resetPasswordConfirm API function
      setSuccess(true); // If successful, set success state to true
      showToast({ type: 'success', message: 'Password reset successfully!' }); // Show success toast notification
      onSuccess(); // Call onSuccess callback if provided

      setTimeout(() => { // Navigate to login page after a delay
        navigate('/login');
      }, 3000);
    } catch (error: any) {
      showToast({ type: 'error', message: error.error || 'Failed to reset password' }); // If error occurs, show error toast notification
    } finally {
      setLoading(false); // Finally, set loading state to false
    }
  };

  /**
   * Handles changes to the password input
   * @param e React.ChangeEvent<HTMLInputElement>
   */
  const handlePasswordChange = (e: React.ChangeEvent<HTMLInputElement>): void => {
    setNewPassword(e.target.value); // Update newPassword state with input value
    setErrors(prevErrors => ({ ...prevErrors, newPassword: '' })); // Clear the newPassword error if it exists
  };

  /**
   * Handles changes to the confirm password input
   * @param e React.ChangeEvent<HTMLInputElement>
   */
  const handleConfirmPasswordChange = (e: React.ChangeEvent<HTMLInputElement>): void => {
    setConfirmPassword(e.target.value); // Update confirmPassword state with input value
    setErrors(prevErrors => ({ ...prevErrors, confirmPassword: '' })); // Clear the confirmPassword error if it exists
  };

  return (
    <Box
      display="flex"
      justifyContent="center"
      alignItems="center"
      minHeight="80vh"
    >
      <Paper elevation={3} sx={{ padding: 4, width: '100%', maxWidth: 400 }}>
        {success ? ( // If success is true, render success message with link to login
          <Box textAlign="center">
            <Typography variant="h6" gutterBottom>
              Password reset successful!
            </Typography>
            <Typography>
              Redirecting to <a href="/login">login page</a>...
            </Typography>
          </Box>
        ) : ( // Otherwise, render the password reset form
          <form onSubmit={handleSubmit}> {/* Form includes title and instructions */}
            <Typography variant="h5" align="center" gutterBottom>
              Reset Password
            </Typography>
            <Typography variant="body2" align="center" paragraph>
              Enter your new password below.
            </Typography>
            <Input // Form includes Input component for new password with error state
              label="New Password"
              type="password"
              fullWidth
              margin="normal"
              value={newPassword}
              onChange={handlePasswordChange}
              error={errors.newPassword}
              helperText={errors.newPassword}
              required
            />
            <Input // Form includes Input component for confirm password with error state
              label="Confirm New Password"
              type="password"
              fullWidth
              margin="normal"
              value={confirmPassword}
              onChange={handleConfirmPasswordChange}
              error={errors.confirmPassword}
              helperText={errors.confirmPassword}
              required
            />
            <Box mt={3}>
              <Button // Form includes Button component for submission with loading state
                type="submit"
                fullWidth
                variant="contained"
                color="primary"
                loading={loading}
              >
                Reset Password
              </Button>
            </Box>
          </form> // Form has onSubmit handler attached
        )}
      </Paper>
    </Box>
  );
};

export default ResetPasswordForm; // Export the ResetPasswordForm component for use in password reset flow