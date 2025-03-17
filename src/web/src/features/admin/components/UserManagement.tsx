import React, { useState, useEffect, useCallback, useMemo } from 'react'; // react ^18.2.0
import { Box, Typography, Grid, Paper, Divider, Chip, IconButton, Tooltip } from '@mui/material'; // @mui/material ^5.13.0
import { Add, Edit, Delete, Refresh, FilterList, Search, PersonAdd } from '@mui/icons-material'; // @mui/icons-material ^5.11.16
import { format } from 'date-fns'; // date-fns ^2.29.0

import { useAdmin } from '../hooks/useAdmin';
import UserDetail from './UserDetail';
import Table from '../../../components/common/Table';
import Button from '../../../components/common/Button';
import Input from '../../../components/common/Input';
import Dropdown from '../../../components/common/Dropdown';
import Dialog from '../../../components/common/Dialog';
import { useToast } from '../../../hooks/useToast';
import { usePermissions } from '../../../hooks/usePermissions';
import { User, UserRole, UserStatus, UserCreate } from '../../../types/user';

/**
 * Interface for form data used in the create user dialog.
 */
interface UserFormData {
  email: string;
  password: string;
  confirm_password: string;
  role: UserRole;
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
 * Component for managing users in the admin interface.
 * @returns Rendered component.
 */
const UserManagement: React.FC = () => {
  // Initialize state for selected user, create user dialog, delete confirmation dialog, and form data
  const [selectedUser, setSelectedUser] = useState<User | null>(null);
  const [isCreateDialogOpen, setIsCreateDialogOpen] = useState(false);
  const [isDeleteDialogOpen, setIsDeleteDialogOpen] = useState(false);
  const [userToDelete, setUserToDelete] = useState<User | null>(null);
  const [formData, setFormData] = useState<UserFormData>({
    email: '',
    password: '',
    confirm_password: '',
    role: UserRole.PHARMA,
  });
  const [formErrors, setFormErrors] = useState<Record<string, string>>({});
  const [searchTerm, setSearchTerm] = useState('');
  const [roleFilter, setRoleFilter] = useState<UserRole | null>(null);
  const [statusFilter, setStatusFilter] = useState<UserStatus | null>(null);

  // Get admin functions from useAdmin hook
  const {
    users,
    isLoadingUsers,
    filters,
    pagination,
    handleFilterChange,
    handlePageChange,
    createUserMutation,
    deleteUserMutation,
    refetchUsers,
  } = useAdmin();

  // Get toast notification function from useToast hook
  const { showToast } = useToast();

  // Get permission check functions from usePermissions hook
  const { isAdmin, canManageUsers } = usePermissions();

  // Define table columns for user list
  const columns = useMemo(
    () => [
      { field: 'id', headerName: 'ID', width: 70 },
      { field: 'email', headerName: 'Email', width: 250 },
      {
        field: 'role',
        headerName: 'Role',
        width: 150,
        renderCell: (row: User) => (
          <Chip label={row.role} color="primary" size="small" />
        ),
      },
      {
        field: 'status',
        headerName: 'Status',
        width: 150,
        renderCell: (row: User) => (
          <Chip label={row.status} color="secondary" size="small" />
        ),
      },
      {
        field: 'actions',
        headerName: 'Actions',
        width: 150,
        renderCell: (row: User) => (
          <>
            <Tooltip title="Edit User">
              <IconButton
                aria-label="edit"
                onClick={() => handleUserSelect(row)}
                disabled={!canManageUsers()}
              >
                <Edit />
              </IconButton>
            </Tooltip>
            <Tooltip title="Delete User">
              <IconButton
                aria-label="delete"
                onClick={() => handleDeleteDialogOpen(row)}
                disabled={!canManageUsers()}
              >
                <Delete />
              </IconButton>
            </Tooltip>
          </>
        ),
      },
    ],
    [canManageUsers, handleUserSelect, handleDeleteDialogOpen]
  );

  // Handle filter changes for role, status, and email search
  const handleSearchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setSearchTerm(event.target.value);
    handleFilterChange({ ...filters, email: event.target.value });
  };

  const handleRoleFilterChange = (role: UserRole | null) => {
    setRoleFilter(role);
    handleFilterChange({ ...filters, role: role });
  };

  const handleStatusFilterChange = (status: UserStatus | null) => {
    setStatusFilter(status);
    handleFilterChange({ ...filters, status: status });
  };

  // Handle pagination changes
  const handlePageChangeLocal = (page: number) => {
    handlePageChange(page);
  };

  // Handle user selection for viewing/editing details
  const handleUserSelect = (user: User) => {
    setSelectedUser(user);
  };

  // Handle closing the user detail view
  const handleUserDetailClose = () => {
    setSelectedUser(null);
  };

  // Handle create user dialog open/close
  const handleCreateDialogOpen = () => {
    setFormData({ email: '', password: '', confirm_password: '', role: UserRole.PHARMA });
    setFormErrors({});
    setIsCreateDialogOpen(true);
  };

  const handleCreateDialogClose = () => {
    setIsCreateDialogOpen(false);
  };

  // Handle delete confirmation dialog open/close
  const handleDeleteDialogOpen = (user: User) => {
    setUserToDelete(user);
    setIsDeleteDialogOpen(true);
  };

  const handleDeleteDialogClose = () => {
    setIsDeleteDialogOpen(false);
    setUserToDelete(null);
  };

  // Handle form field changes for create user form
  const handleFormChange = (event: React.ChangeEvent<HTMLInputElement | { name?: string; value: unknown }>) => {
    const { name, value } = event.target;
    setFormData({
      ...formData,
      [name as string]: value,
    });
  };

  // Validate form data before submission
  const validateForm = (): boolean => {
    let errors: Record<string, string> = {};

    if (!formData.email) {
      errors.email = 'Email is required';
    } else if (!/^[^s@]+@[^s@]+.[^s@]+$/.test(formData.email)) {
      errors.email = 'Invalid email format';
    }

    if (formData.password.length < 8) {
      errors.password = 'Password must be at least 8 characters';
    }

    if (formData.password !== formData.confirm_password) {
      errors.confirm_password = 'Passwords do not match';
    }

    if (!formData.role) {
      errors.role = 'Role is required';
    }

    setFormErrors(errors);
    return Object.keys(errors).length === 0;
  };

  // Handle form submission for creating new user
  const handleCreateUser = (event: React.FormEvent) => {
    event.preventDefault();
    if (validateForm()) {
      const { confirm_password, ...createRequest } = formData;
      createUserMutation.mutate(
        createRequest,
        {
          onSuccess: () => {
            showToast({ type: 'success', message: 'User created successfully' });
            handleCreateDialogClose();
            refetchUsers();
          },
          onError: (error: any) => {
            showToast({ type: 'error', message: error.message || 'Failed to create user' });
          },
        }
      );
    }
  };

  // Handle user deletion
  const handleDeleteUser = () => {
    if (userToDelete) {
      deleteUserMutation.mutate(
        userToDelete.id,
        {
          onSuccess: () => {
            showToast({ type: 'success', message: 'User deleted successfully' });
            handleDeleteDialogClose();
            refetchUsers();
          },
          onError: (error: any) => {
            showToast({ type: 'error', message: error.message || 'Failed to delete user' });
          },
        }
      );
    }
  };

  return (
    <Box>
      {/* Header section with title and action buttons */}
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
        <Typography variant="h6">User Management</Typography>
        {canManageUsers() && (
          <Box>
            <Tooltip title="Refresh Users">
              <IconButton aria-label="refresh" onClick={refetchUsers}>
                <Refresh />
              </IconButton>
            </Tooltip>
            <Button variant="contained" startIcon={<PersonAdd />} onClick={handleCreateDialogOpen}>
              Add User
            </Button>
          </Box>
        )}
      </Box>

      {/* Filter section with search input, role filter, and status filter */}
      <Paper sx={{ p: 2, mb: 2 }}>
        <Grid container spacing={2} alignItems="center">
          <Grid item xs={12} md={4}>
            <Input
              label="Search by Email"
              value={searchTerm}
              onChange={handleSearchChange}
              fullWidth
              InputProps={{
                startAdornment: (
                  <Search color="action" />
                ),
              }}
            />
          </Grid>
          <Grid item xs={12} md={4}>
            <Dropdown
              id="role-filter"
              name="role"
              label="Filter by Role"
              value={roleFilter || ''}
              onChange={(e) => handleRoleFilterChange(e.target.value as UserRole || null)}
              options={[
                { value: '', label: 'All Roles' },
                { value: UserRole.PHARMA, label: 'Pharma' },
                { value: UserRole.CRO, label: 'CRO' },
                { value: UserRole.ADMIN, label: 'Admin' },
              ]}
              fullWidth
            />
          </Grid>
          <Grid item xs={12} md={4}>
            <Dropdown
              id="status-filter"
              name="status"
              label="Filter by Status"
              value={statusFilter || ''}
              onChange={(e) => handleStatusFilterChange(e.target.value as UserStatus || null)}
              options={[
                { value: '', label: 'All Statuses' },
                { value: UserStatus.ACTIVE, label: 'Active' },
                { value: UserStatus.INACTIVE, label: 'Inactive' },
                { value: UserStatus.PENDING, label: 'Pending' },
                { value: UserStatus.LOCKED, label: 'Locked' },
              ]}
              fullWidth
            />
          </Grid>
        </Grid>
      </Paper>

      {/* User table */}
      <Table
        data={users}
        columns={columns}
        loading={isLoadingUsers}
        error={createUserMutation.isError ? createUserMutation.error?.message : null}
        pagination={{
          page: pagination.page,
          pageSize: pagination.pageSize,
          totalItems: pagination.totalItems,
          onPageChange: handlePageChangeLocal,
          onPageSizeChange: (pageSize) => handleFilterChange({ ...filters, pageSize }),
          pageSizeOptions: [5, 10, 20, 50],
          showPageSizeSelector: true,
        }}
      />

      {/* User detail component */}
      {selectedUser && (
        <UserDetail user={selectedUser} onClose={handleUserDetailClose} onUpdate={refetchUsers} />
      )}

      {/* Create user dialog */}
      <Dialog open={isCreateDialogOpen} onClose={handleCreateDialogClose} title="Create New User" onConfirm={handleCreateUser} confirmButtonText="Create" loading={createUserMutation.isLoading} onCancel={handleCreateDialogClose}>
        <Grid container spacing={2}>
          <Grid item xs={12}>
            <Input
              label="Email"
              name="email"
              value={formData.email}
              onChange={handleFormChange}
              error={formErrors.email}
              fullWidth
            />
          </Grid>
          <Grid item xs={12}>
            <Input
              label="Password"
              type="password"
              name="password"
              value={formData.password}
              onChange={handleFormChange}
              error={formErrors.password}
              fullWidth
            />
          </Grid>
          <Grid item xs={12}>
            <Input
              label="Confirm Password"
              type="password"
              name="confirm_password"
              value={formData.confirm_password}
              onChange={handleFormChange}
              error={formErrors.confirm_password}
              fullWidth
            />
          </Grid>
          <Grid item xs={12}>
            <Dropdown
              id="role"
              name="role"
              label="Role"
              value={formData.role}
              onChange={handleFormChange}
              options={[
                { value: UserRole.PHARMA, label: 'Pharma' },
                { value: UserRole.CRO, label: 'CRO' },
                { value: UserRole.ADMIN, label: 'Admin' },
              ]}
              fullWidth
            />
          </Grid>
        </Grid>
      </Dialog>

      {/* Delete confirmation dialog */}
      <Dialog
        open={isDeleteDialogOpen}
        onClose={handleDeleteDialogClose}
        title="Confirm Delete"
        contentText={`Are you sure you want to delete user ${userToDelete?.email}?`}
        onConfirm={handleDeleteUser}
        confirmButtonText="Delete"
        onCancel={handleDeleteDialogClose}
        loading={deleteUserMutation.isLoading}
      />
    </Box>
  );
};

export default UserManagement;