# src/web/src/features/admin/components/ActivityLog.tsx
```typescript
import React, { useState, useEffect, useMemo } from 'react'; // react ^18.2.0
import {
  Box,
  Typography,
  TextField,
  MenuItem,
  Select,
  FormControl,
  InputLabel,
  Button,
  IconButton,
  Tooltip,
} from '@mui/material'; // @mui/material ^5.13.0
import { FilterList, Refresh, CalendarToday } from '@mui/icons-material'; // @mui/icons-material ^5.13.0
import { DatePicker } from '@mui/x-date-pickers'; // @mui/x-date-pickers ^6.0.0
import dayjs from 'dayjs'; // dayjs ^1.11.7
import { useAdmin } from '../hooks/useAdmin';
import Card from '../../../components/common/Card';
import Table from '../../../components/common/Table';
import { ActivityLogItem } from '../../../types';

/**
 * Interface for filter options for activity logs
 */
interface ActivityLogFilters {
  user: string;
  action: string;
  resourceType: string;
  startDate: Date | null;
  endDate: Date | null;
}

/**
 * Props for the ActivityLog component
 */
interface ActivityLogProps {
  className?: string;
  style?: React.CSSProperties;
}

/**
 * Component that displays system activity logs with filtering and pagination
 * @param props - Props for the ActivityLog component
 * @returns Rendered ActivityLog component
 */
const ActivityLog: React.FC<ActivityLogProps> = (props) => {
  // Get activity logs and related functions from useAdmin hook
  const { activityLogs, isLoadingActivityLogs, refetchActivityLogs } = useAdmin();

  // Initialize state for filters (user, action, resource type, date range)
  const [filters, setFilters] = useState<ActivityLogFilters>({
    user: '',
    action: '',
    resourceType: '',
    startDate: null,
    endDate: null,
  });

  // Initialize state for pagination (page, pageSize)
  const [page, setPage] = useState(1);
  const [pageSize, setPageSize] = useState(10);

  // Define table columns for activity log display
  const columns = useMemo(
    () => [
      {
        field: 'timestamp',
        headerName: 'Timestamp',
        width: 150,
        renderCell: (row: ActivityLogItem) => dayjs(row.timestamp).format('YYYY-MM-DD HH:mm:ss'),
      },
      {
        field: 'userEmail',
        headerName: 'User',
        width: 200,
        renderCell: (row: ActivityLogItem) => row.userEmail || 'System',
      },
      {
        field: 'action',
        headerName: 'Action',
        width: 150,
      },
      {
        field: 'resourceType',
        headerName: 'Resource Type',
        width: 150,
      },
      {
        field: 'resourceId',
        headerName: 'Resource ID',
        width: 100,
      },
      {
        field: 'details',
        headerName: 'Details',
        width: 300,
        renderCell: (row: ActivityLogItem) => (
          <Tooltip title={JSON.stringify(row.details, null, 2)} placement="right">
            <span>View Details</span>
          </Tooltip>
        ),
      },
    ],
    []
  );

  // Create memoized filtered logs based on current filters
  const filteredLogs = useMemo(() => {
    return activityLogs.filter((log) => {
      if (filters.user && !log.userEmail?.toLowerCase().includes(filters.user.toLowerCase())) {
        return false;
      }
      if (filters.action && !log.action.toLowerCase().includes(filters.action.toLowerCase())) {
        return false;
      }
      if (filters.resourceType && !log.resourceType.toLowerCase().includes(filters.resourceType.toLowerCase())) {
        return false;
      }
      if (filters.startDate && dayjs(log.timestamp) < dayjs(filters.startDate)) {
        return false;
      }
      if (filters.endDate && dayjs(log.timestamp) > dayjs(filters.endDate)) {
        return false;
      }
      return true;
    });
  }, [activityLogs, filters]);

  // Handle filter changes and update state
  const handleFilterChange = (field: string, value: string | Date | null) => {
    setFilters((prevFilters) => ({
      ...prevFilters,
      [field]: value,
    }));
  };

  // Handle pagination changes
  const handlePageChange = (newPage: number) => {
    setPage(newPage);
  };

  // Implement refresh function to reload activity logs
  const handleRefresh = () => {
    refetchActivityLogs();
  };

  return (
    <Card className={props.className} style={props.style}>
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2, p: 2 }}>
        <Typography variant="h6">System Activity Logs</Typography>
        <Box>
          <IconButton onClick={handleRefresh} disabled={isLoadingActivityLogs}>
            <Refresh />
          </IconButton>
        </Box>
      </Box>

      {/* Render filter section with text inputs, selects, and date pickers */}
      <Box sx={{ display: 'flex', gap: 2, p: 2, alignItems: 'center', flexWrap: 'wrap' }}>
        <TextField
          label="User"
          value={filters.user}
          onChange={(e) => handleFilterChange('user', e.target.value)}
          size="small"
        />
        <TextField
          label="Action"
          value={filters.action}
          onChange={(e) => handleFilterChange('action', e.target.value)}
          size="small"
        />
        <TextField
          label="Resource Type"
          value={filters.resourceType}
          onChange={(e) => handleFilterChange('resourceType', e.target.value)}
          size="small"
        />
        <DatePicker
          label="Start Date"
          value={filters.startDate}
          onChange={(date) => handleFilterChange('startDate', date)}
          slotProps={{ textField: { size: 'small' } }}
        />
        <DatePicker
          label="End Date"
          value={filters.endDate}
          onChange={(date) => handleFilterChange('endDate', date)}
          slotProps={{ textField: { size: 'small' } }}
        />
      </Box>

      {/* Render Table component with activity log data */}
      <Table
        data={filteredLogs.slice((page - 1) * pageSize, page * pageSize)}
        columns={columns}
        loading={isLoadingActivityLogs}
        emptyMessage="No activity logs found"
        pagination={{
          page,
          pageSize,
          totalItems: filteredLogs.length,
          onPageChange: handlePageChange,
          onPageSizeChange: (newPageSize) => setPageSize(newPageSize),
          pageSizeOptions: [5, 10, 20, 50],
          showPageSizeSelector: true,
        }}
      />
    </Card>
  );
};

export default ActivityLog;