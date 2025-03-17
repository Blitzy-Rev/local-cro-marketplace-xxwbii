import React, { useState, useEffect, useMemo, useCallback } from 'react'; // React v18.2+
import {
  Box,
  Typography,
  Button,
  Chip,
  Dialog,
  DialogTitle,
  DialogContent,
  DialogActions,
  TextField,
  FormControl,
  InputLabel,
  Select,
  MenuItem,
  IconButton,
  Tooltip,
} from '@mui/material'; // Material-UI components v5.13+
import { FilterList, Refresh, CheckCircle, Error, Warning, Info, Close } from '@mui/icons-material'; // Material-UI icons v5.13+
import { format } from 'date-fns'; // Date formatting utility v2.30+
import { useAdmin } from '../hooks/useAdmin'; // Hook for accessing and managing system alerts
import Card from '../../../components/common/Card'; // Card component for containing the alerts display
import Badge from '../../../components/common/Badge'; // Badge component for displaying alert severity
import Table from '../../../components/common/Table'; // Table component for displaying alerts in a structured format
import { useToast } from '../../../hooks/useToast'; // Hook for displaying toast notifications
import { SystemAlertSeverity } from '../../../types'; // Enum for alert severity levels
import { SystemAlertResponse } from '../../../types';

/**
 * Props interface for the SystemAlerts component
 */
interface SystemAlertsProps {
  /** Optional CSS class name for styling */
  className?: string;
}

/**
 * Interface for alert filtering state
 */
interface AlertFilterState {
  /** Filter by alert severity */
  severity: SystemAlertSeverity | 'ALL';
  /** Filter by resolution status */
  resolved: boolean | null;
  /** Filter by alert type */
  alertType: string | null;
}

/**
 * Interface for the alert detail dialog props
 */
interface AlertDetailDialogProps {
  /** Alert data to display */
  alert: SystemAlertResponse;
  /** Whether the dialog is open */
  open: boolean;
  /** Handler for dialog close */
  onClose: () => void;
  /** Handler for resolving the alert */
  onResolve: (alertId: number, notes: string) => void;
}

/**
 * Maps alert severity to corresponding color for visual indication
 * @param severity Alert severity
 * @returns Color name for the severity level
 */
const getSeverityColor = (severity: SystemAlertSeverity): string => {
  switch (severity) {
    case SystemAlertSeverity.CRITICAL:
      return 'error';
    case SystemAlertSeverity.WARNING:
      return 'warning';
    case SystemAlertSeverity.INFO:
      return 'info';
    default:
      return 'default';
  }
};

/**
 * Returns the appropriate icon component for a given severity level
 * @param severity Alert severity
 * @returns Icon component for the severity level
 */
const getSeverityIcon = (severity: SystemAlertSeverity): React.ReactNode => {
  switch (severity) {
    case SystemAlertSeverity.CRITICAL:
      return <Error />;
    case SystemAlertSeverity.WARNING:
      return <Warning />;
    case SystemAlertSeverity.INFO:
      return <Info />;
    default:
      return null;
  }
};

/**
 * Component that displays system alerts with filtering and management capabilities
 */
const SystemAlerts: React.FC<SystemAlertsProps> = ({ className }) => {
  // Access system alerts and management functions from the useAdmin hook
  const { systemAlerts, refetchSystemAlerts, updateAlertMutation } = useAdmin();
  // Access the showToast function from the useToast hook
  const { showToast } = useToast();

  // State for managing filters
  const [filters, setFilters] = useState<AlertFilterState>({ severity: 'ALL', resolved: false, alertType: null });
  // State for managing filter dialog visibility
  const [showFilterDialog, setShowFilterDialog] = useState<boolean>(false);
  // State for managing selected alert for detail view
  const [selectedAlert, setSelectedAlert] = useState<SystemAlertResponse | null>(null);
  // State for managing alert detail dialog visibility
  const [showAlertDetail, setShowAlertDetail] = useState<boolean>(false);

  /**
   * Updates the filter state and applies filters
   * @param newFilters New filter state
   */
  const handleFilterChange = useCallback((newFilters: AlertFilterState) => {
    setFilters(newFilters);
    refetchSystemAlerts();
  }, [refetchSystemAlerts]);

  /**
   * Resolves an alert with optional notes
   * @param alertId ID of the alert to resolve
   * @param notes Optional notes for the resolution
   */
  const handleResolveAlert = (alertId: number, notes: string) => {
    updateAlertMutation.mutate({ alertId, alertData: { resolved: true, notes } }, {
      onSuccess: () => {
        refetchSystemAlerts();
        showToast({ type: 'success', message: 'Alert resolved successfully!' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: `Failed to resolve alert: ${error.message || 'Unknown error'}` });
      },
    });
  };

  /**
   * Opens the alert detail dialog for a specific alert
   * @param alert Alert data
   */
  const handleViewAlertDetail = (alert: SystemAlertResponse) => {
    setSelectedAlert(alert);
    setShowAlertDetail(true);
  };

  /**
   * Closes the alert detail dialog
   */
  const handleCloseAlertDetail = () => {
    setShowAlertDetail(false);
    setSelectedAlert(null);
  };

  // Define table columns
  const columns = useMemo(() => [
    {
      field: 'timestamp',
      headerName: 'Timestamp',
      width: 150,
      valueGetter: (row: SystemAlertResponse) => format(new Date(row.timestamp), 'yyyy-MM-dd HH:mm:ss'),
    },
    {
      field: 'severity',
      headerName: 'Severity',
      width: 120,
      renderCell: (row: SystemAlertResponse) => (
        <Badge label={row.severity} color={getSeverityColor(row.severity)} icon={getSeverityIcon(row.severity)} />
      ),
    },
    {
      field: 'alertType',
      headerName: 'Type',
      width: 150,
    },
    {
      field: 'message',
      headerName: 'Message',
      width: 300,
    },
    {
      field: 'actions',
      headerName: 'Actions',
      width: 150,
      renderCell: (row: SystemAlertResponse) => (
        <>
          <Tooltip title="View Details">
            <IconButton onClick={() => handleViewAlertDetail(row)} size="small">
              <Info />
            </IconButton>
          </Tooltip>
          {!row.resolved && (
            <Tooltip title="Resolve Alert">
              <IconButton onClick={() => handleViewAlertDetail(row)} size="small">
                <CheckCircle />
              </IconButton>
            </Tooltip>
          )}
        </>
      ),
    },
  ], [handleViewAlertDetail]);

  return (
    <Card className={className} header={
      <Box display="flex" justifyContent="space-between" alignItems="center">
        <Typography variant="h6">System Alerts</Typography>
        <Box>
          <IconButton onClick={() => setShowFilterDialog(true)} aria-label="filter">
            <FilterList />
          </IconButton>
          <IconButton onClick={() => refetchSystemAlerts()} aria-label="refresh">
            <Refresh />
          </IconButton>
        </Box>
      </Box>
    }>
      <Table
        data={systemAlerts}
        columns={columns}
        loading={updateAlertMutation.isLoading}
        emptyMessage="No system alerts to display."
      />
      <FilterDialog
        open={showFilterDialog}
        onClose={() => setShowFilterDialog(false)}
        filters={filters}
        onApplyFilters={handleFilterChange}
      />
      <AlertDetailDialog
        alert={selectedAlert!}
        open={showAlertDetail}
        onClose={handleCloseAlertDetail}
        onResolve={handleResolveAlert}
      />
    </Card>
  );
};

/**
 * Props interface for the alert filter dialog
 */
interface FilterDialogProps {
  /** Whether the dialog is open */
  open: boolean;
  /** Handler for dialog close */
  onClose: () => void;
  /** Current filter state */
  filters: AlertFilterState;
  /** Handler for applying filters */
  onApplyFilters: (filters: AlertFilterState) => void;
}

/**
 * Dialog component for filtering alerts by various criteria
 */
const FilterDialog: React.FC<FilterDialogProps> = ({ open, onClose, filters, onApplyFilters }) => {
  // Local state for managing filters within the dialog
  const [localFilters, setLocalFilters] = useState<AlertFilterState>(filters);

  /**
   * Updates local filter state when a filter option changes
   * @param event Change event
   * @param field Field to update
   */
  const handleChange = (event: React.ChangeEvent<HTMLInputElement | { value: unknown }>, field: string) => {
    setLocalFilters({ ...localFilters, [field]: event.target.value });
  };

  /**
   * Applies the current filters and closes the dialog
   */
  const handleApply = () => {
    onApplyFilters(localFilters);
    onClose();
  };

  /**
   * Resets filters to default values
   */
  const handleReset = () => {
    setLocalFilters({ severity: 'ALL', resolved: false, alertType: null });
  };

  return (
    <Dialog open={open} onClose={onClose} aria-labelledby="filter-dialog-title">
      <DialogTitle id="filter-dialog-title">Filter Alerts</DialogTitle>
      <DialogContent>
        <FormControl fullWidth margin="dense">
          <InputLabel id="severity-select-label">Severity</InputLabel>
          <Select
            labelId="severity-select-label"
            id="severity-select"
            value={localFilters.severity}
            label="Severity"
            onChange={(e) => handleChange(e, 'severity')}
          >
            <MenuItem value="ALL">All</MenuItem>
            <MenuItem value={SystemAlertSeverity.INFO}>Info</MenuItem>
            <MenuItem value={SystemAlertSeverity.WARNING}>Warning</MenuItem>
            <MenuItem value={SystemAlertSeverity.CRITICAL}>Critical</MenuItem>
          </Select>
        </FormControl>
        <FormControl fullWidth margin="dense">
          <InputLabel id="resolved-select-label">Resolution Status</InputLabel>
          <Select
            labelId="resolved-select-label"
            id="resolved-select"
            value={localFilters.resolved === null ? '' : localFilters.resolved}
            label="Resolution Status"
            onChange={(e) => handleChange(e, 'resolved')}
          >
            <MenuItem value="">All</MenuItem>
            <MenuItem value={true}>Resolved</MenuItem>
            <MenuItem value={false}>Unresolved</MenuItem>
          </Select>
        </FormControl>
        <TextField
          fullWidth
          margin="dense"
          id="alert-type"
          label="Alert Type"
          type="text"
          value={localFilters.alertType || ''}
          onChange={(e) => handleChange(e, 'alertType')}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={handleReset}>Reset</Button>
        <Button onClick={onClose}>Cancel</Button>
        <Button onClick={handleApply} variant="contained">Apply</Button>
      </DialogActions>
    </Dialog>
  );
};

/**
 * Props interface for the alert detail dialog
 */
interface AlertDetailDialogProps {
  /** Alert data to display */
  alert: SystemAlertResponse;
  /** Whether the dialog is open */
  open: boolean;
  /** Handler for dialog close */
  onClose: () => void;
  /** Handler for resolving the alert */
  onResolve: (alertId: number, notes: string) => void;
}

/**
 * Dialog component for displaying detailed information about an alert and allowing resolution
 */
const AlertDetailDialog: React.FC<AlertDetailDialogProps> = ({ alert, open, onClose, onResolve }) => {
  // Local state for managing resolution notes
  const [resolutionNotes, setResolutionNotes] = useState<string>('');

  /**
   * Resolves the alert with optional notes
   */
  const handleResolve = () => {
    onResolve(alert.id, resolutionNotes);
    onClose();
  };

  return (
    <Dialog open={open} onClose={onClose} aria-labelledby="alert-detail-dialog-title">
      <DialogTitle id="alert-detail-dialog-title">
        Alert Details
      </DialogTitle>
      <DialogContent>
        <Box mb={2}>
          <Badge label={alert.severity} color={getSeverityColor(alert.severity)} icon={getSeverityIcon(alert.severity)} />
        </Box>
        <Typography variant="subtitle2">Timestamp:</Typography>
        <Typography variant="body2">{format(new Date(alert.timestamp), 'yyyy-MM-dd HH:mm:ss')}</Typography>
        <Typography variant="subtitle2">Type:</Typography>
        <Typography variant="body2">{alert.alertType}</Typography>
        <Typography variant="subtitle2">Message:</Typography>
        <Typography variant="body2">{alert.message}</Typography>
        <Typography variant="subtitle2">Details:</Typography>
        <Typography variant="body2">{JSON.stringify(alert.details, null, 2)}</Typography>
        {!alert.resolved && (
          <TextField
            fullWidth
            margin="dense"
            id="resolution-notes"
            label="Resolution Notes"
            multiline
            rows={4}
            value={resolutionNotes}
            onChange={(e) => setResolutionNotes(e.target.value)}
          />
        )}
      </DialogContent>
      <DialogActions>
        {!alert.resolved && (
          <Button onClick={handleResolve} variant="contained" color="primary">
            Resolve
          </Button>
        )}
        <Button onClick={onClose}>Close</Button>
      </DialogActions>
    </Dialog>
  );
};

export default SystemAlerts;