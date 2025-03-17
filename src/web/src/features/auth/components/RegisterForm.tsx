import React, { useState, useEffect } from 'react'; // React 18.2+
import {
  FormControl,
  FormLabel,
  RadioGroup,
  Radio,
  FormControlLabel,
  Box,
  Typography,
  Link,
  Paper,
} from '@mui/material'; // @mui/material v5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13.0

// Internal imports
import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import { useAuth } from '../hooks/useAuth';
import { useToast } from '../../../hooks/useToast';
import { validateRegistrationData } from '../../../utils/validation';
import { RegisterData } from '../../../types/auth';
import { UserRole } from '../../../types/user';

/**
 * Props interface for the RegisterForm component
 */
interface RegisterFormProps {
  /** Callback function to execute after successful registration */
  onSuccess: () => void;
  /** Callback function to execute when the login link is clicked */
  onLoginClick: () => void;
}

/**
 * A styled container for the registration form
 */
const FormContainer = styled(Paper)(({ theme }) => ({
  padding: theme.spacing(3), // 24px
  maxWidth: '500px',
  width: '100%',
  margin: '0 auto',
  borderRadius: theme.shape.borderRadius, // 8px
}));

/**
 * A form component for user registration with email, password, and role selection
 */
const RegisterForm: React.FC<RegisterFormProps> = ({ onSuccess, onLoginClick }) => {
  // State for form input values
  const [formData, setFormData] = useState<RegisterData>({
    email: '',
    password: '',
    confirm_password: '',
    role: UserRole.PHARMA,
  });

  // State for validation error messages
  const [errors, setErrors] = useState<Record<string, string>>({});

  // State for tracking form submission status
  const [isSubmitting, setIsSubmitting] = useState<boolean>(false);

  // State for tracking successful registration
  const [registrationSuccess, setRegistrationSuccess] = useState<boolean>(false);

  // Use the custom authentication hook
  const { register, loading, error } = useAuth();

  // Use the custom toast hook
  const { showToast } = useToast();

  /**
   * Handles input field changes and updates form state
   * @param e React.ChangeEvent<HTMLInputElement>
   */
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setFormData({
      ...formData,
      [name]: value,
    });
    setErrors({
      ...errors,
      [name]: '', // Clear error for the changed field
    });
  };

  /**
   * Handles role selection change
   * @param e React.ChangeEvent<HTMLInputElement>
   */
  const handleRoleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { value } = e.target;
    setFormData({
      ...formData,
      role: value as UserRole,
    });
    setErrors({
      ...errors,
      role: '', // Clear error for the role field
    });
  };

  /**
   * Handles form submission
   * @param e React.FormEvent<HTMLFormElement>
   */
  const handleSubmit = async (e: React.FormEvent<HTMLFormElement>) => {
    e.preventDefault();
    setIsSubmitting(true);

    // Validate form data
    const validationResult = validateRegistrationData(formData);
    if (!validationResult.isValid) {
      setErrors(validationResult.errors);
      setIsSubmitting(false);
      return;
    }

    try {
      // Call the register function from the useAuth hook
      await register(formData);

      // Set registration success to true
      setRegistrationSuccess(true);

      // Show success toast
      showToast({
        type: 'success',
        message: 'Registration successful! Please check your email to verify your account.',
      });

      // Call the onSuccess prop if it exists
      if (onSuccess) {
        onSuccess();
      }
    } catch (err: any) {
      // Show error toast
      showToast({
        type: 'error',
        message: err.message || 'Registration failed. Please try again.',
      });
    } finally {
      setIsSubmitting(false);
    }
  };

  /**
   * Handle authentication errors
   */
  useEffect(() => {
    if (error) {
      setErrors({ ...errors, auth: error });
    }
  }, [error]);

  /**
   * Handle successful registration
   */
  useEffect(() => {
    if (registrationSuccess) {
      // Reset form data
      setFormData({
        email: '',
        password: '',
        confirm_password: '',
        role: UserRole.PHARMA,
      });

      // Clear errors
      setErrors({});
    }
  }, [registrationSuccess]);

  return (
    <FormContainer component="form" onSubmit={handleSubmit}>
      <Typography variant="h5" gutterBottom>
        Create an Account
      </Typography>
      <Typography variant="body2" color="textSecondary" paragraph>
        Join our platform to manage molecular data and streamline CRO interactions.
      </Typography>

      <Input
        label="Email"
        name="email"
        type="email"
        value={formData.email}
        onChange={handleChange}
        error={errors.email}
        required
      />

      <Input
        label="Password"
        name="password"
        type="password"
        value={formData.password}
        onChange={handleChange}
        error={errors.password}
        required
      />

      <Input
        label="Confirm Password"
        name="confirm_password"
        type="password"
        value={formData.confirm_password}
        onChange={handleChange}
        error={errors.confirm_password}
        required
      />

      <FormControl component="fieldset" margin="normal">
        <FormLabel component="legend">Role</FormLabel>
        <RadioGroup
          name="role"
          value={formData.role}
          onChange={handleRoleChange}
          row
        >
          <FormControlLabel
            value={UserRole.PHARMA}
            control={<Radio />}
            label="Pharma User"
          />
          <FormControlLabel
            value={UserRole.CRO}
            control={<Radio />}
            label="CRO User"
          />
        </RadioGroup>
        {errors.role && (
          <Typography variant="caption" color="error">
            {errors.role}
          </Typography>
        )}
      </FormControl>

      <Box mt={2}>
        <Button
          type="submit"
          fullWidth
          disabled={isSubmitting}
          loading={loading}
        >
          Register
        </Button>
      </Box>

      {registrationSuccess && (
        <Box mt={2} color="success.main">
          <Typography variant="body2">
            Registration successful! Please check your email to verify your account.
          </Typography>
        </Box>
      )}

      <Box mt={2} textAlign="center">
        <Typography variant="body2">
          Already have an account?{' '}
          <Link component="button" variant="body2" onClick={onLoginClick}>
            Log In
          </Link>
        </Typography>
      </Box>
    </FormContainer>
  );
};

export default RegisterForm;