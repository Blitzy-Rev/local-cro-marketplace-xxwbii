import React from 'react'; // React v18.2+
import { Box, Typography } from '@mui/material'; // Material-UI components v5.13+

import AdminLayout from '../../../layouts/AdminLayout';
import UserManagement from '../components/UserManagement';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Page component that renders the user management interface within the admin layout
 * @returns Rendered component
 */
const UserManagementPage: React.FC = () => {
  // Check if user has permission to manage users using usePermissions hook
  const { canManageUsers } = usePermissions();

  // Render AdminLayout component to ensure only admin users can access this page
  return (
    <AdminLayout>
      <Box sx={{ padding: 3 }}>
        {/* Render page title and description */}
        <Typography variant="h4" component="h1" gutterBottom>
          User Management
        </Typography>
        <Typography variant="body1" paragraph>
          Manage user accounts and roles within the system.
        </Typography>

        {/* Render UserManagement component that provides the user management functionality */}
        {canManageUsers() ? (
          <UserManagement />
        ) : (
          <Typography variant="body1" color="error">
            You do not have permission to manage users.
          </Typography>
        )}
      </Box>
    </AdminLayout>
  );
};

// Export the UserManagementPage component for use in routing
export default UserManagementPage;