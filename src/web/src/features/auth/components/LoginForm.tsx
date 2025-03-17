import React, { useState, useEffect } from 'react'; // React 18.2+
import { Link } from 'react-router-dom'; // react-router-dom ^6.10.0
import { Box, Typography, Checkbox, FormControlLabel, Alert } from '@mui/material'; // @mui/material ^5.13.0
import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import useAuth from '../hooks/useAuth';
import { validateLoginCredentials, validateEmail, validatePassword } from '../../../utils/validation';
import { LoginCredentials } from '../../../types/auth';

/**
 * Props interface for the LoginForm component
 */
interface LoginFormProps {
  onSuccess: () => void;
  redirectPath: string;
}

/**
 * A form component that handles user login with validation and error handling
 */
const LoginForm: React.FC<LoginFormProps> = ({ onSuccess, redirectPath = '/dashboard' }) => {
  // Authentication hook for login functionality
  const { login, error, loading, clearError } = useAuth();

  // State for storing user input credentials
  const [credentials, setCredentials] = useState<LoginCredentials>({ email: '', password: '' });

  // State for storing validation errors
  const [errors, setErrors] = useState<Record<string, string>>({});

  // State for remember me checkbox
  const [rememberMe, setRememberMe] = useState<boolean>(false);

  // State for tracking which fields have been touched
  const [touched, setTouched] = useState<Record<string, boolean>>({ email: false, password: false });

  /**
   * Handles input field changes and updates state
   * @param e - React change event
   */
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setCredentials(prev => ({ ...prev, [name]: value }));
    setTouched(prev => ({ ...prev, [name]: true }));
    validateField(name, value);
  };

  /**
   * Handles input field blur events for validation
   * @param field - The name of the field that was blurred
   */
  const handleBlur = (field: string) => {
    setTouched(prev => ({ ...prev, [field]: true }));
    const value = credentials[field as keyof LoginCredentials];
    validateField(field, value);
  };

  /**
   * Validates a specific form field
   * @param field - The name of the field to validate
   * @param value - The value of the field
   * @returns Error message or empty string
   */
  const validateField = (field: string, value: string): string => {
    let errorMessage = '';

    if (field === 'email') {
      const emailValidation = validateEmail(value);
      if (!emailValidation.isValid) {
        errorMessage = emailValidation.error || '';
      }
    } else if (field === 'password') {
      const passwordValidation = validatePassword(value);
      if (!passwordValidation.isValid) {
        errorMessage = passwordValidation.error || '';
      }
    }

    setErrors(prev => ({ ...prev, [field]: errorMessage }));
    return errorMessage;
  };

  /**
   * Handles form submission
   * @param e - React form event
   */
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    // Validate all form fields
    const validationResult = validateLoginCredentials(credentials);
    if (!validationResult.isValid) {
      setErrors(validationResult.errors);
      return;
    }

    try {
      await login(credentials);
      if (onSuccess) {
        onSuccess();
      }

      // Store remember me preference
      localStorage.setItem('rememberMe', rememberMe.toString());
    } catch (loginError: any) {
      // Handle login error (already handled by useAuth hook)
    }
  };

  /**
   * Clears authentication errors when component unmounts
   */
  useEffect(() => {
    return () => {
      clearError();
    };
  }, [clearError]);

  return (
    <Box
      component="form"
      onSubmit={handleSubmit}
      sx={{
        maxWidth: 400,
        mx: 'auto',
        p: 3,
        borderRadius: 2,
        bgcolor: 'background.paper',
        boxShadow: 1,
      }}
    >
      <Typography variant="h5" align="center" mb={3}>
        Login
      </Typography>
      {error && (
        <Alert severity="error" sx={{ mb: 2 }}>
          {error}
        </Alert>
      )}
      <Input
        label="Email"
        name="email"
        type="email"
        value={credentials.email}
        onChange={handleChange}
        onBlur={() => handleBlur('email')}
        error={touched.email && errors.email}
        helperText={touched.email && errors.email}
        required
        fullWidth
        margin="normal"
      />
      <Input
        label="Password"
        name="password"
        type="password"
        value={credentials.password}
        onChange={handleChange}
        onBlur={() => handleBlur('password')}
        error={touched.password && errors.password}
        helperText={touched.password && errors.password}
        required
        fullWidth
        margin="normal"
      />
      <FormControlLabel
        control={
          <Checkbox
            checked={rememberMe}
            onChange={(e) => setRememberMe(e.target.checked)}
            name="rememberMe"
            color="primary"
          />
        }
        label="Remember me"
        sx={{ display: 'block', mt: 1, mb: 2 }}
      />
      <Button
        type="submit"
        fullWidth
        variant="contained"
        color="primary"
        loading={loading}
        disabled={loading}
      >
        Log In
      </Button>
      <Box mt={2} textAlign="center">
        <Link to="/forgot-password" style={{ display: 'block', marginBottom: '0.5rem' }}>
          Forgot password?
        </Link>
        <Link to="/register">
          Don't have an account? Register
        </Link>
      </Box>
    </Box>
  );
};

export default LoginForm;