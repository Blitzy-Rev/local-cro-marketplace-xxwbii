import React, { useEffect, useState } from 'react'; // react ^18.2.0
import { Box, Typography, Grid, Chip, IconButton, Tooltip } from '@mui/material'; // @mui/material ^5.13.0
import { Refresh as RefreshIcon, CheckCircle as CheckIcon, Error as ErrorIcon, Warning as WarningIcon } from '@mui/icons-material'; // @mui/icons-material ^5.13.0
import { useAdmin } from '../hooks/useAdmin';
import Card from '../../../components/common/Card';
import ProgressBar from '../../../components/common/ProgressBar';

/**
 * @packageDocumentation
 * @module src/web/src/features/admin/components/SystemStatus
 */

/**
 * Interface for system health status data
 */
interface SystemHealthStatus {
  database: string;
  storage: string;
  services: string;
  lastChecked: string;
}

/**
 * Interface for system resource usage data
 */
interface SystemResourceUsage {
  memory: number;
  cpu: number;
  disk: number;
  network: object;
}

/**
 * Props interface for the SystemStatus component
 */
interface SystemStatusProps {
  refreshInterval?: number;
  showHeader?: boolean;
  title?: string;
  className?: string;
  style?: React.CSSProperties;
}

/**
 * Determines the color to use for status indicators based on the status value
 * @param status 
 * @returns Color to use for the status indicator (success, error, warning, or info)
 */
const getStatusColor = (status: string): string => {
  if (status === 'OK' || status === 'Healthy') {
    return 'success';
  }
  if (status === 'Error' || status === 'Critical') {
    return 'error';
  }
  if (status === 'Warning') {
    return 'warning';
  }
  return 'info';
};

/**
 * Determines the color to use for progress bars based on the percentage value
 * @param percentage 
 * @returns Color to use for the progress bar (success, error, warning, or info)
 */
const getProgressBarColor = (percentage: number): string => {
  if (percentage < 70) {
    return 'success';
  }
  if (percentage < 85) {
    return 'warning';
  }
  return 'error';
};

/**
 * Component that displays system status information
 * @param props 
 * @returns Rendered component
 */
export const SystemStatus: React.FC<SystemStatusProps> = ({ refreshInterval = 60000, showHeader = true, title = 'System Status', className, style }) => {
  // Use the useAdmin hook to get system health and resources data
  const { systemHealth, systemResources, refetchSystemHealth, refetchSystemResources } = useAdmin();

  // Set up state for auto-refresh interval
  const [intervalId, setIntervalId] = useState<NodeJS.Timeout | null>(null);

  // Set up useEffect to periodically refresh data based on refreshInterval
  useEffect(() => {
    if (refreshInterval > 0) {
      const id = setInterval(() => {
        refetchSystemHealth();
        refetchSystemResources();
      }, refreshInterval);
      setIntervalId(id);
      return () => clearInterval(id);
    }
    return () => {};
  }, [refreshInterval, refetchSystemHealth, refetchSystemResources]);

  // Handle manual refresh button click
  const handleRefresh = () => {
    refetchSystemHealth();
    refetchSystemResources();
  };

  // Render a Card component containing the system status information
  return (
    <Card className={className} style={style}>
      {/* Display header with title and refresh button */}
      {showHeader && (
        <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
          <Typography variant="h6">{title}</Typography>
          <Tooltip title="Refresh Status">
            <IconButton onClick={handleRefresh}>
              <RefreshIcon />
            </IconButton>
          </Tooltip>
        </Box>
      )}

      {/* Create a Grid layout for organizing status items */}
      <Grid container spacing={2}>
        {/* For each system component (Database, Storage, Services), display status with appropriate icon */}
        <Grid item xs={12} sm={6} md={4}>
          <Typography variant="subtitle1">Database:</Typography>
          <Chip
            icon={
              systemHealth?.database === 'OK' ? (
                <CheckIcon />
              ) : systemHealth?.database === 'Error' ? (
                <ErrorIcon />
              ) : (
                <WarningIcon />
              )
            }
            label={systemHealth?.database || 'Loading...'}
            color={getStatusColor(systemHealth?.database || '')}
          />
        </Grid>

        <Grid item xs={12} sm={6} md={4}>
          <Typography variant="subtitle1">Storage:</Typography>
          <Chip
            icon={
              systemHealth?.storage === 'OK' ? (
                <CheckIcon />
              ) : systemHealth?.storage === 'Error' ? (
                <ErrorIcon />
              ) : (
                <WarningIcon />
              )
            }
            label={systemHealth?.storage || 'Loading...'}
            color={getStatusColor(systemHealth?.storage || '')}
          />
        </Grid>

        <Grid item xs={12} sm={6} md={4}>
          <Typography variant="subtitle1">Services:</Typography>
          <Chip
            icon={
              systemHealth?.services === 'OK' ? (
                <CheckIcon />
              ) : systemHealth?.services === 'Error' ? (
                <ErrorIcon />
              ) : (
                <WarningIcon />
              )
            }
            label={systemHealth?.services || 'Loading...'}
            color={getStatusColor(systemHealth?.services || '')}
          />
        </Grid>

        {/* For resource usage (Memory, CPU, Disk), display ProgressBar with appropriate color based on usage percentage */}
        <Grid item xs={12} sm={6}>
          <Typography variant="subtitle1">Memory Usage:</Typography>
          <ProgressBar
            value={systemResources?.memory}
            color={getProgressBarColor(systemResources?.memory || 0)}
            showPercentage
          />
        </Grid>

        <Grid item xs={12} sm={6}>
          <Typography variant="subtitle1">CPU Usage:</Typography>
          <ProgressBar
            value={systemResources?.cpu}
            color={getProgressBarColor(systemResources?.cpu || 0)}
            showPercentage
          />
        </Grid>

        <Grid item xs={12} sm={6}>
          <Typography variant="subtitle1">Disk Usage:</Typography>
          <ProgressBar
            value={systemResources?.disk}
            color={getProgressBarColor(systemResources?.disk || 0)}
            showPercentage
          />
        </Grid>
      </Grid>
    </Card>
  );
};