import React, { useState, useEffect } from 'react'; // react ^18.2.0
import { Box, Typography, Grid, Divider, Switch, FormControlLabel, Dialog, DialogTitle, DialogContent, DialogActions } from '@mui/material'; // @mui/material ^5.13.0
import { format } from 'date-fns'; // date-fns ^2.29.0
import { LockReset, Edit, Save, Cancel } from '@mui/icons-material'; // @mui/icons-material ^5.11.16

import { useAdmin } from '../hooks/useAdmin';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Input from '../../../components/common/Input';
import Dropdown from '../../../components/common/Dropdown';
import { useToast } from '../../../hooks/useToast';
import { User, UserRole, UserStatus, UserAdminUpdateRequest } from '../../../types/user';

/**
 * Interface defining the props for the UserDetail component.
 */
interface UserDetailProps {
  user: User;
  onClose: () => void;
  onUpdate: () => void;
}

/**
 * Interface defining the structure for user form data.
 */
interface UserFormData {
  email: string;
  role: UserRole;
  status: UserStatus;
  is_active: boolean;
  email_verified: boolean;
}

/**
 * Interface defining the structure for password reset form data.
 */
interface PasswordResetFormData {
  new_password: string;
  confirm_password: string;
}

/**
 * Formats a date string into a readable format.
 * @param dateString The date string to format.
 * @returns The formatted date string or 'Never' if date is null or undefined.
 */
const formatDate = (dateString: string | null | undefined): string => {
  if (!dateString) {
    return 'Never';
  }
  return format(new Date(dateString), 'MMM dd, yyyy hh:mm a');
};

/**
 * Component for displaying and editing user details in the admin interface.
 * @param user The user object containing user details.
 * @param onClose Function to call when closing the detail view.
 * @param onUpdate Function to call after successful user update to refresh data.
 * @returns Rendered component.
 */
const UserDetail: React.FC<UserDetailProps> = ({ user, onClose, onUpdate }) => {
  // State variables for managing edit mode, form data, and dialog visibility
  const [editMode, setEditMode] = useState(false);
  const [formData, setFormData] = useState<UserFormData>({
    email: user.email,
    role: user.role,
    status: user.status,
    is_active: user.is_active,
    email_verified: user.email_verified,
  });
  const [formErrors, setFormErrors] = useState<Record<string, string>>({});
  const [isPasswordResetDialogOpen, setIsPasswordResetDialogOpen] = useState(false);
  const [passwordResetData, setPasswordResetData] = useState<PasswordResetFormData>({
    new_password: '',
    confirm_password: '',
  });
  const [passwordResetErrors, setPasswordResetErrors] = useState<Record<string, string>>({});

  // Access admin functions and toast notification function
  const { updateUserMutation, resetPasswordMutation } = useAdmin();
  const { showToast } = useToast();

  // Initialize form data when user prop changes
  useEffect(() => {
    setFormData({
      email: user.email,
      role: user.role,
      status: user.status,
      is_active: user.is_active,
      email_verified: user.email_verified,
    });
  }, [user]);

  /**
   * Toggles edit mode on and off.
   */
  const handleEditToggle = () => {
    if (editMode) {
      // Reset form data to original user values when canceling edit
      setFormData({
        email: user.email,
        role: user.role,
        status: user.status,
        is_active: user.is_active,
        email_verified: user.email_verified,
      });
      setFormErrors({});
    }
    setEditMode(!editMode);
  };

  /**
   * Updates form data when form fields change.
   * @param event The change event from the input field.
   */
  const handleFormChange = (event: React.ChangeEvent<HTMLInputElement | { name?: string; value: unknown }>) => {
    const { name, value } = event.target;
    setFormData({
      ...formData,
      [name as string]: value,
    });
  };

  /**
   * Updates boolean form fields when switches change.
   * @param name The name of the field to update.
   * @param checked The new checked value.
   */
  const handleSwitchChange = (name: string, checked: boolean) => {
    setFormData({
      ...formData,
      [name]: checked,
    });
  };

  /**
   * Validates the user edit form data.
   * @returns True if the form data is valid, false otherwise.
   */
  const validateForm = (): boolean => {
    let errors: Record<string, string> = {};

    if (!formData.email) {
      errors.email = 'Email is required';
    } else if (!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(formData.email)) {
      errors.email = 'Invalid email format';
    }

    if (!formData.role) {
      errors.role = 'Role is required';
    }

    if (!formData.status) {
      errors.status = 'Status is required';
    }

    setFormErrors(errors);
    return Object.keys(errors).length === 0;
  };

  /**
   * Handles submission of the user edit form.
   * @param event The form submit event.
   */
  const handleSubmit = (event: React.FormEvent) => {
    event.preventDefault();
    if (validateForm()) {
      const updateData: UserAdminUpdateRequest = {
        email: formData.email,
        role: formData.role,
        status: formData.status,
        is_active: formData.is_active,
        email_verified: formData.email_verified,
      };

      updateUserMutation.mutate(
        { userId: user.id, userData: updateData },
        {
          onSuccess: () => {
            showToast({ type: 'success', message: 'User updated successfully' });
            setEditMode(false);
            onUpdate();
          },
          onError: (error: any) => {
            showToast({ type: 'error', message: error.message || 'Failed to update user' });
          },
        }
      );
    }
  };

  /**
   * Opens the password reset dialog.
   */
  const handlePasswordResetDialogOpen = () => {
    setPasswordResetData({ new_password: '', confirm_password: '' });
    setPasswordResetErrors({});
    setIsPasswordResetDialogOpen(true);
  };

  /**
   * Closes the password reset dialog.
   */
  const handlePasswordResetDialogClose = () => {
    setIsPasswordResetDialogOpen(false);
  };

  /**
   * Updates password reset form data when form fields change.
   * @param event The change event from the input field.
   */
  const handlePasswordResetFormChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = event.target;
    setPasswordResetData({
      ...passwordResetData,
      [name]: value,
    });
  };

  /**
   * Validates the password reset form data.
   * @returns True if the form data is valid, false otherwise.
   */
  const validatePasswordResetForm = (): boolean => {
    let errors: Record<string, string> = {};

    if (passwordResetData.new_password.length < 8) {
      errors.new_password = 'Password must be at least 8 characters';
    }

    if (passwordResetData.new_password !== passwordResetData.confirm_password) {
      errors.confirm_password = 'Passwords do not match';
    }

    setPasswordResetErrors(errors);
    return Object.keys(errors).length === 0;
  };

  /**
   * Handles submission of the password reset form.
   * @param event The form submit event.
   */
  const handlePasswordResetSubmit = (event: React.FormEvent) => {
    event.preventDefault();
    if (validatePasswordResetForm()) {
      resetPasswordMutation.mutate(
        { userId: user.id, passwordData: passwordResetData },
        {
          onSuccess: () => {
            showToast({ type: 'success', message: 'Password reset successfully' });
            setIsPasswordResetDialogOpen(false);
          },
          onError: (error: any) => {
            showToast({ type: 'error', message: error.message || 'Failed to reset password' });
          },
        }
      );
    }
  };

  return (
    <Box>
      <Card>
        <Typography variant="h6">Account Information</Typography>
        <Divider />
        <Grid container spacing={2} sx={{ padding: 2 }}>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">User ID:</Typography>
            <Typography>{user.id}</Typography>
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Email:</Typography>
            {editMode ? (
              <Input
                name="email"
                value={formData.email}
                onChange={handleFormChange}
                error={formErrors.email}
              />
            ) : (
              <Typography>{user.email}</Typography>
            )}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Role:</Typography>
            {editMode ? (
              <Dropdown
                name="role"
                label="Role"
                value={formData.role}
                onChange={handleFormChange}
                options={[
                  { value: 'pharma', label: 'Pharma' },
                  { value: 'cro', label: 'CRO' },
                  { value: 'admin', label: 'Admin' },
                ]}
                error={formErrors.role}
              />
            ) : (
              <Typography>{user.role}</Typography>
            )}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Status:</Typography>
            {editMode ? (
              <Dropdown
                name="status"
                label="Status"
                value={formData.status}
                onChange={handleFormChange}
                options={[
                  { value: 'pending', label: 'Pending' },
                  { value: 'active', label: 'Active' },
                  { value: 'inactive', label: 'Inactive' },
                  { value: 'locked', label: 'Locked' },
                ]}
                error={formErrors.status}
              />
            ) : (
              <Typography>{user.status}</Typography>
            )}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Active:</Typography>
            {editMode ? (
              <FormControlLabel
                control={
                  <Switch
                    checked={formData.is_active}
                    onChange={(e) => handleSwitchChange('is_active', e.target.checked)}
                    name="is_active"
                  />
                }
                label="Is Active"
              />
            ) : (
              <Typography>{formData.is_active ? 'Yes' : 'No'}</Typography>
            )}
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Email Verified:</Typography>
            {editMode ? (
              <FormControlLabel
                control={
                  <Switch
                    checked={formData.email_verified}
                    onChange={(e) => handleSwitchChange('email_verified', e.target.checked)}
                    name="email_verified"
                  />
                }
                label="Email Verified"
              />
            ) : (
              <Typography>{formData.email_verified ? 'Yes' : 'No'}</Typography>
            )}
          </Grid>
        </Grid>
      </Card>

      <Card sx={{ mt: 2 }}>
        <Typography variant="h6">Activity Information</Typography>
        <Divider />
        <Grid container spacing={2} sx={{ padding: 2 }}>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Created At:</Typography>
            <Typography>{formatDate(user.created_at)}</Typography>
          </Grid>
          <Grid item xs={12} md={6}>
            <Typography variant="subtitle2">Last Login:</Typography>
            <Typography>{formatDate(user.last_login)}</Typography>
          </Grid>
        </Grid>
      </Card>

      <Box sx={{ mt: 2, display: 'flex', justifyContent: 'space-between' }}>
        {editMode ? (
          <Box>
            <Button variant="contained" color="primary" onClick={handleSubmit}>
              <Save />
              Save
            </Button>
            <Button variant="outlined" onClick={handleEditToggle} sx={{ ml: 1 }}>
              <Cancel />
              Cancel
            </Button>
          </Box>
        ) : (
          <Button variant="outlined" onClick={handleEditToggle}>
            <Edit />
            Edit
          </Button>
        )}
        <Button variant="text" onClick={handlePasswordResetDialogOpen}>
          <LockReset />
          Reset Password
        </Button>
      </Box>

      <Dialog open={isPasswordResetDialogOpen} onClose={handlePasswordResetDialogClose}>
        <DialogTitle>Reset Password for {user.email}</DialogTitle>
        <DialogContent>
          <Grid container spacing={2} sx={{ mt: 1 }}>
            <Grid item xs={12}>
              <Input
                label="New Password"
                type="password"
                name="new_password"
                value={passwordResetData.new_password}
                onChange={handlePasswordResetFormChange}
                error={passwordResetErrors.new_password}
                fullWidth
              />
            </Grid>
            <Grid item xs={12}>
              <Input
                label="Confirm Password"
                type="password"
                name="confirm_password"
                value={passwordResetData.confirm_password}
                onChange={handlePasswordResetFormChange}
                error={passwordResetErrors.confirm_password}
                fullWidth
              />
            </Grid>
          </Grid>
        </DialogContent>
        <DialogActions>
          <Button onClick={handlePasswordResetDialogClose}>Cancel</Button>
          <Button onClick={handlePasswordResetSubmit} variant="contained" color="primary">
            Reset Password
          </Button>
        </DialogActions>
      </Dialog>
    </Box>
  );
};

export default UserDetail;