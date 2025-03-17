import React, { useState, useEffect } from 'react';
import { Link } from 'react-router-dom'; // react-router-dom ^6.10.0
import { Box, Typography, Alert } from '@mui/material'; // @mui/material ^5.13.0

import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import { validateEmail } from '../../../utils/validation';
import { resetPassword } from '../../../api/auth';

/**
 * Props interface for the ForgotPasswordForm component
 */
interface ForgotPasswordFormProps {
  /**
   * Callback function to be called after a successful password reset request
   */
  onSuccess: () => void;
}

/**
 * A form component that handles password reset requests with validation and error handling
 */
const ForgotPasswordForm: React.FC<ForgotPasswordFormProps> = ({ onSuccess }) => {
  // State for storing user input email
  const [email, setEmail] = useState('');
  // State for storing API error message
  const [error, setError] = useState<string | null>(null);
  // State for storing email validation error
  const [validationError, setValidationError] = useState<string | null>(null);
  // State for tracking form submission status
  const [loading, setLoading] = useState(false);
  // State for tracking successful submission
  const [success, setSuccess] = useState(false);
  // State for tracking if email field has been touched
  const [touched, setTouched] = useState(false);

  // Reset error state when email changes
  useEffect(() => {
    if (error) {
      setError(null);
    }
  }, [email]);

  // Handles email input field changes and updates state
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    setEmail(value);
    
    // Only validate if field has been touched
    if (touched) {
      const { isValid, error } = validateEmail(value);
      setValidationError(isValid ? null : error);
    }
  };

  // Handles email input field blur event for validation
  const handleBlur = () => {
    setTouched(true);
    const { isValid, error } = validateEmail(email);
    setValidationError(isValid ? null : error);
  };

  // Handles form submission for password reset request
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    // Validate email before submission
    const { isValid, error: emailError } = validateEmail(email);
    if (!isValid) {
      setValidationError(emailError);
      return;
    }

    setLoading(true);
    
    try {
      // Call the resetPassword API function with email
      const response = await resetPassword(email);
      
      // Set success state and call onSuccess callback
      setSuccess(true);
      onSuccess();
    } catch (err: any) {
      // Handle API errors by setting error state
      setError(err.error || 'An error occurred while requesting a password reset.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <Box
      component="form"
      onSubmit={handleSubmit}
      sx={{
        display: 'flex',
        flexDirection: 'column',
        gap: 2,
        width: '100%',
        maxWidth: '400px',
        margin: '0 auto',
      }}
      noValidate
      aria-labelledby="forgot-password-title"
    >
      {success ? (
        // Show success message if request was successful
        <Box sx={{ textAlign: 'center' }}>
          <Alert severity="success" sx={{ mb: 2 }}>
            Password reset link sent!
          </Alert>
          <Typography variant="body1" sx={{ mb: 2 }}>
            If an account exists with the email you provided, you will receive instructions 
            to reset your password shortly. Please check your email and follow the instructions.
          </Typography>
          <Typography variant="body2" sx={{ mb: 3, color: 'text.secondary' }}>
            If you don't see the email, please check your spam folder.
          </Typography>
          <Link to="/login" style={{ textDecoration: 'none' }}>
            <Button variant="outlined" fullWidth>
              Back to Login
            </Button>
          </Link>
        </Box>
      ) : (
        // Show the form if request has not been submitted yet
        <>
          <Typography id="forgot-password-title" variant="h5" component="h1" align="center" gutterBottom>
            Forgot Password
          </Typography>
          
          <Typography variant="body1" align="center" sx={{ mb: 2 }}>
            Enter your email address below and we'll send you a link to reset your password.
          </Typography>
          
          {error && (
            <Alert severity="error" sx={{ mb: 2 }}>
              {error}
            </Alert>
          )}
          
          <Input
            label="Email Address"
            type="email"
            name="email"
            value={email}
            onChange={handleChange}
            onBlur={handleBlur}
            error={validationError}
            required
            disabled={loading}
            fullWidth
            aria-required="true"
            placeholder="your-email@example.com"
          />
          
          <Button
            type="submit"
            variant="contained"
            color="primary"
            fullWidth
            loading={loading}
            disabled={loading || !!validationError}
            sx={{ mt: 1 }}
            aria-label="Request password reset"
          >
            Reset Password
          </Button>
          
          <Box textAlign="center" mt={2}>
            <Link to="/login" style={{ textDecoration: 'none' }}>
              Back to Login
            </Link>
          </Box>
        </>
      )}
    </Box>
  );
};

export default ForgotPasswordForm;