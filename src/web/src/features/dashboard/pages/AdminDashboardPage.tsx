import React, { useEffect } from 'react'; // React 18.2+
import { Grid, Box, Typography, Container, useTheme, useMediaQuery } from '@mui/material'; // @mui/material v5.13+
import Stats from '../../dashboard/components/Stats';
import QuickActions from '../../dashboard/components/QuickActions';
import SystemStatus from '../../admin/components/SystemStatus';
import SystemAlerts from '../../admin/components/SystemAlerts';
import ActivityLog from '../../admin/components/ActivityLog';
import { useAdmin } from '../../admin/hooks/useAdmin';

/**
 * Main component for the admin dashboard page
 * @returns Rendered admin dashboard page
 */
const AdminDashboardPage: React.FC = () => {
  // Initialize theme using useTheme hook
  const theme = useTheme();

  // Check if screen is mobile using useMediaQuery hook
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // Initialize admin functionality using useAdmin hook
  const { refetchSystemHealth, refetchSystemResources } = useAdmin();

  // Set up effect to refresh system data at regular intervals
  useEffect(() => {
    const intervalId = setInterval(() => {
      refetchSystemHealth();
      refetchSystemResources();
    }, 60000); // Refresh every 60 seconds

    return () => clearInterval(intervalId); // Clear interval on unmount
  }, [refetchSystemHealth, refetchSystemResources]);

  // Render page title and container
  return (
    <Container maxWidth="xl">
      <Box sx={{ mb: 4 }}>
        <Typography variant="h4" component="h1" gutterBottom>
          Admin Dashboard
        </Typography>
      </Box>

      {/* Render Stats component for key system statistics */}
      <Box sx={{ mb: 3 }}>
        <Stats />
      </Box>

      {/* Render Grid layout with SystemStatus and SystemAlerts components */}
      <Grid container spacing={3} sx={{ mb: 3 }}>
        <Grid item xs={12} md={6}>
          <SystemStatus />
        </Grid>
        <Grid item xs={12} md={6}>
          <SystemAlerts />
        </Grid>
      </Grid>

      {/* Render QuickActions component for admin quick actions */}
      <Box sx={{ mb: 3 }}>
        <QuickActions />
      </Box>

      {/* Render ActivityLog component for system activity monitoring */}
      <Box>
        <ActivityLog />
      </Box>
    </Container>
  );
};

export default AdminDashboardPage;