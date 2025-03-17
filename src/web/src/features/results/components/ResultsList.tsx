import React, { useState, useEffect, useMemo, useCallback } from 'react'; // react ^18.2.0
import {
  Box,
  Typography,
  Button,
  Chip,
  TextField,
  MenuItem,
  FormControl,
  InputLabel,
  Select,
  IconButton,
  Tooltip
} from '@mui/material'; // @mui/material ^5.13.0
import {
  CheckCircle,
  Cancel,
  FilterList,
  Refresh,
  CheckBox,
  CheckBoxOutlineBlank
} from '@mui/icons-material'; // @mui/icons-material ^5.11.16
import { format, parseISO } from 'date-fns'; // date-fns ^2.30.0

import useResults from '../hooks/useResults';
import { Result, ResultStatus, ResultFilter } from '../../../types/result';
import Table from '../../../components/common/Table';
import useToast from '../../../hooks/useToast';
import useDialog from '../../../hooks/useDialog';
import usePermissions from '../../../hooks/usePermissions';
import { UserRole } from '../../../types/user';

interface ResultsListProps {
  submissionId?: string;
  onResultSelect: (resultId: string) => void;
  initialFilter?: Partial<ResultFilter>;
}

/**
 * Formats a date string into a readable format
 * @param dateString 
 * @returns Formatted date string
 */
const formatDate = (dateString: string): string => {
  const date = parseISO(dateString);
  return format(date, 'MMM d, yyyy');
};

/**
 * Renders a status chip with appropriate color and icon based on result status
 * @param status 
 * @returns A Material-UI Chip component with status information
 */
const getStatusChip = (status: ResultStatus): JSX.Element => {
  let color: 'success' | 'error' | 'info' | 'default' = 'default';
  let icon = null;

  switch (status) {
    case ResultStatus.APPROVED:
      color = 'success';
      icon = <CheckCircle />;
      break;
    case ResultStatus.REJECTED:
      color = 'error';
      icon = <Cancel />;
      break;
    case ResultStatus.UPLOADED:
      color = 'info';
      break;
    default:
      break;
  }

  return (
    <Chip
      label={status}
      color={color}
      icon={icon}
      size="small"
    />
  );
};

/**
 * A component that displays a paginated, filterable list of experimental results
 * @param props 
 * @returns The rendered component
 */
const ResultsList: React.FC<ResultsListProps> = ({ submissionId, onResultSelect, initialFilter }) => {
  // State for selected results
  const [selectedResults, setSelectedResults] = useState<string[]>([]);

  // Initialize useResults hook
  const {
    results,
    loading,
    error,
    pagination,
    filter,
    setFilter,
    setPage,
    setPageSize,
    fetchResultDetail,
    batchApproveResults,
    batchRejectResults
  } = useResults({ initialFilter, submissionId });

  // Initialize useToast hook
  const { showToast } = useToast();

  // Initialize useDialog hook
  const { showDialog } = useDialog();

  // Initialize usePermissions hook
  const { canApproveResults } = usePermissions();

  // State for status filter
  const [statusFilter, setStatusFilter] = useState<ResultStatus | ''>('');

  // State for date range filter
  const [dateRangeFilter, setDateRangeFilter] = useState<{ startDate: string | null; endDate: string | null }>({
    startDate: null,
    endDate: null
  });

  // State for notes
  const [notes, setNotes] = useState<string>('');

  /**
   * Handles selecting a result to view details
   * @param resultId 
   */
  const handleResultSelect = (resultId: string) => {
    onResultSelect(resultId);
  };

  /**
   * Handles approving multiple selected results
   */
  const handleBatchApprove = async () => {
    showDialog({
      title: 'Approve Results',
      content: 'Are you sure you want to approve the selected results?',
      onConfirm: async () => {
        try {
          await batchApproveResults(selectedResults, notes);
          showToast({ message: 'Results approved successfully', type: 'success' });
          setSelectedResults([]);
        } catch (err: any) {
          showToast({ message: err.message || 'Failed to approve results', type: 'error' });
        }
      }
    });
  };

  /**
   * Handles rejecting multiple selected results
   */
  const handleBatchReject = async () => {
    showDialog({
      title: 'Reject Results',
      content: 'Are you sure you want to reject the selected results?',
      onConfirm: async () => {
        try {
          await batchRejectResults(selectedResults, notes);
          showToast({ message: 'Results rejected successfully', type: 'success' });
          setSelectedResults([]);
        } catch (err: any) {
          showToast({ message: err.message || 'Failed to reject results', type: 'error' });
        }
      }
    });
  };

  /**
   * Handles changes to the status filter
   * @param event 
   */
  const handleStatusFilterChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setStatusFilter(event.target.value as ResultStatus | '');
    setFilter({ ...filter, status: event.target.value as ResultStatus | '' });
  };

  /**
   * Handles changes to the date range filter
   * @param type 
   * @param value 
   */
  const handleDateRangeChange = (type: 'startDate' | 'endDate', value: string) => {
    setDateRangeFilter({ ...dateRangeFilter, [type]: value });
    setFilter({ ...filter, [`uploaded_${type}`]: value });
  };

  /**
   * Handles clearing all filters
   */
  const handleClearFilters = () => {
    setStatusFilter('');
    setDateRangeFilter({ startDate: null, endDate: null });
    setFilter({});
  };

  // Define table columns
  const columns = useMemo(() => [
    {
      field: 'status',
      headerName: 'Status',
      width: '120px',
      renderCell: (row: Result) => getStatusChip(row.status)
    },
    {
      field: 'uploaded_at',
      headerName: 'Date',
      width: '150px',
      renderCell: (row: Result) => formatDate(row.uploaded_at)
    },
    {
      field: 'submission.experiment_id',
      headerName: 'Experiment',
      width: '200px',
      renderCell: (row: Result) => row.submission?.experiment_id
    },
    {
      field: 'file_count',
      headerName: 'Files',
      width: '100px',
      renderCell: (row: Result) => row.file_count
    },
    {
      field: 'actions',
      headerName: 'Actions',
      width: '120px',
      renderCell: (row: Result) => (
        <Button size="small" onClick={() => handleResultSelect(row.id)}>
          View
        </Button>
      )
    }
  ], [handleResultSelect]);

  return (
    <Box sx={{ width: '100%', mt: 2, display: 'flex', flexDirection: 'column', gap: 2 }}>
      {/* Filter Section */}
      <Box sx={{ display: 'flex', gap: 2, alignItems: 'center' }}>
        <FormControl sx={{ minWidth: 120 }} size="small">
          <InputLabel id="status-filter-label">Status</InputLabel>
          <Select
            labelId="status-filter-label"
            id="status-filter"
            value={statusFilter}
            label="Status"
            onChange={handleStatusFilterChange}
          >
            <MenuItem value="">All</MenuItem>
            <MenuItem value={ResultStatus.PENDING}>Pending</MenuItem>
            <MenuItem value={ResultStatus.UPLOADED}>Uploaded</MenuItem>
            <MenuItem value={ResultStatus.APPROVED}>Approved</MenuItem>
            <MenuItem value={ResultStatus.REJECTED}>Rejected</MenuItem>
          </Select>
        </FormControl>

        <TextField
          id="start-date-filter"
          label="Start Date"
          type="date"
          size="small"
          value={dateRangeFilter.startDate || ''}
          onChange={(e) => handleDateRangeChange('startDate', e.target.value)}
          InputLabelProps={{ shrink: true }}
        />

        <TextField
          id="end-date-filter"
          label="End Date"
          type="date"
          size="small"
          value={dateRangeFilter.endDate || ''}
          onChange={(e) => handleDateRangeChange('endDate', e.target.value)}
          InputLabelProps={{ shrink: true }}
        />

        <Button variant="outlined" size="small" onClick={handleClearFilters}>
          Clear Filters
        </Button>
      </Box>

      {/* Batch Actions Section */}
      {selectedResults.length > 0 && (
        <Box sx={{ display: 'flex', gap: 2, alignItems: 'center' }}>
          <TextField
            label="Notes"
            multiline
            rows={2}
            value={notes}
            onChange={(e) => setNotes(e.target.value)}
            size="small"
            sx={{ flexGrow: 1 }}
          />
          {canApproveResults() && (
            <Button variant="contained" color="primary" onClick={handleBatchApprove}>
              Approve
            </Button>
          )}
          {canApproveResults() && (
            <Button variant="contained" color="error" onClick={handleBatchReject}>
              Reject
            </Button>
          )}
        </Box>
      )}

      {/* Table Component */}
      <Table
        data={results || []}
        columns={columns}
        loading={loading}
        error={error}
        pagination={{
          page: pagination?.currentPage || 1,
          pageSize: pagination?.pageSize || 10,
          totalItems: pagination?.totalResults || 0,
          onPageChange: setPage,
          onPageSizeChange: setPageSize,
          pageSizeOptions: [10, 25, 50, 100],
          showPageSizeSelector: true,
        }}
        selectable
        selectedRows={selectedResults}
        onSelectionChange={setSelectedResults}
        idField="id"
      />
    </Box>
  );
};

export default ResultsList;