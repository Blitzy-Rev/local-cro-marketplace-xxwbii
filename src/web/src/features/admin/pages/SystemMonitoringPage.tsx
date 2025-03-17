import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { Box, Typography, Container, Grid, Paper } from '@mui/material'; // @mui/material v5.13+
import { Dashboard as DashboardIcon, Warning as AlertIcon, History as ActivityIcon } from '@mui/icons-material'; // @mui/icons-material v5.13+
import AdminLayout from '../../../layouts/AdminLayout';
import { useAdmin } from '../hooks/useAdmin';
import SystemStatus from '../components/SystemStatus';
import SystemAlerts from '../components/SystemAlerts';
import ActivityLog from '../components/ActivityLog';
import { Tabs, TabPanel } from '../../../components/common/Tabs';

/**
 * @file A React component that serves as the main system monitoring page for administrators.
 * It provides a comprehensive dashboard for monitoring system health, resource usage, activity logs, and system alerts in the Molecular Data Management and CRO Integration Platform.
 */

/**
 * @interface TabConfig
 * @description Configuration for a tab in the monitoring page
 */
interface TabConfig {
  /**
   * @property label
   * @type string
   * @description Display label for the tab
   */
  label: string;
  /**
   * @property icon
   * @type React.ReactNode
   * @description Icon to display with the tab label
   */
  icon: React.ReactNode;
}

/**
 * @function SystemMonitoringPage
 * @description Main component for the system monitoring page
 * @returns {JSX.Element} Rendered component
 */
const SystemMonitoringPage: React.FC = () => {
  // LD1: Initialize state for the currently selected tab
  const [selectedTab, setSelectedTab] = useState(0);

  // LD1: Use the useAdmin hook to get system monitoring data
  const {
    systemStats,
    systemResources,
    systemHealth,
    activityLogs,
    systemAlerts,
    refetchSystemStats,
    refetchSystemHealth,
    refetchSystemResources,
    refetchActivityLogs,
    refetchSystemAlerts,
  } = useAdmin();

  // LD1: Set up auto-refresh interval for system data
  useEffect(() => {
    const intervalId = setInterval(() => {
      refetchSystemStats();
      refetchSystemHealth();
      refetchSystemResources();
      refetchActivityLogs();
      refetchSystemAlerts();
    }, 60000); // Refresh every 60 seconds

    // LD1: Clean up interval on component unmount
    return () => clearInterval(intervalId);
  }, [refetchSystemStats, refetchSystemHealth, refetchSystemResources, refetchActivityLogs, refetchSystemAlerts]);

  // LD1: Handle tab change events
  const handleTabChange = useCallback((event: React.SyntheticEvent, newValue: number) => {
    setSelectedTab(newValue);
  }, []);

  // LD1: Define tab configurations
  const tabsConfig: TabConfig[] = [
    { label: 'Dashboard', icon: <DashboardIcon /> },
    { label: 'Alerts', icon: <AlertIcon /> },
    { label: 'Activity', icon: <ActivityIcon /> },
  ];

  // LD1: Render the AdminLayout component as the page container
  return (
    <AdminLayout>
      {/* LD1: Render page header with title and description */}
      <Container maxWidth="xl">
        <Box sx={{ mb: 4 }}>
          <Typography variant="h4" component="h1" gutterBottom>
            System Monitoring
          </Typography>
          <Typography variant="subtitle1">
            Monitor system health, resource usage, and recent activity.
          </Typography>
        </Box>

        {/* LD1: Render Tabs component with Dashboard, Alerts, and Activity tabs */}
        <Tabs tabs={tabsConfig} value={selectedTab} onChange={handleTabChange}>
          {/* LD1: Render TabPanel for System Status showing SystemStatus component */}
          <TabPanel value={0} index={0}>
            <SystemStatus
              systemStats={systemStats}
              systemResources={systemResources}
              systemHealth={systemHealth}
            />
          </TabPanel>

          {/* LD1: Render TabPanel for System Alerts showing SystemAlerts component */}
          <TabPanel value={1} index={1}>
            <SystemAlerts systemAlerts={systemAlerts} />
          </TabPanel>

          {/* LD1: Render TabPanel for Activity Log showing ActivityLog component */}
          <TabPanel value={2} index={2}>
            <ActivityLog activityLogs={activityLogs} />
          </TabPanel>
        </Tabs>
      </Container>
    </AdminLayout>
  );
};

// IE3: Export the SystemMonitoringPage component for use in routing
export default SystemMonitoringPage;