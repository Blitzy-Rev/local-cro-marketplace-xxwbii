# src/web/src/features/experiments/components/ExperimentList.tsx
```typescript
import React, { useState, useEffect, useMemo, useCallback } from 'react'; // react 18.2+
import { useNavigate, Link } from 'react-router-dom'; // react-router-dom ^6.11.0
import {
  Box,
  Typography,
  Paper,
  Grid,
  IconButton,
  Tooltip,
  Chip,
  TextField,
  MenuItem,
  Select,
  FormControl,
  InputLabel
} from '@mui/material'; // @mui/material ^5.13.0
import {
  Add as AddIcon,
  Visibility as ViewIcon,
  Edit as EditIcon,
  Delete as DeleteIcon,
  PlayArrow as QueueIcon,
  Send as SubmitIcon,
  Cancel as CancelIcon
} from '@mui/icons-material'; // @mui/icons-material ^5.13.0
import { format as formatDate } from 'date-fns'; // date-fns ^2.30.0
import Table from '../../../components/common/Table'; // Reusable table component for displaying experiment data
import Button from '../../../components/common/Button'; // Reusable button component for actions
import ExperimentStatus from './ExperimentStatus'; // Component to display experiment status with appropriate styling
import { useExperiments } from '../hooks/useExperiments'; // Hook for managing experiment data and operations
import { usePermissions } from '../../../hooks/usePermissions'; // Hook for checking user permissions
import { Experiment, ExperimentStatus as ExperimentStatusEnum, ExperimentFilter } from '../../../types/experiment'; // Type definitions for experiment data

/**
 * Interface defining the props for the ExperimentList component
 */
interface ExperimentListProps {
  /** Whether to show filter controls */
  showFilters?: boolean;
  /** Whether to show action buttons */
  showActions?: boolean;
  /** Whether to show only the current user's experiments */
  onlyMyExperiments?: boolean;
  /** Initial filter values */
  initialFilter?: Partial<ExperimentFilter>;
  /** Callback when an experiment is selected */
  onExperimentSelect?: (experiment: Experiment) => void;
  /** Additional CSS class */
  className?: string;
  /** Inline styles */
  style?: React.CSSProperties;
}

/**
 * A component that displays a list of experiments with filtering, sorting, and pagination
 * @param props - The component props
 * @returns The rendered experiment list component
 */
const ExperimentList: React.FC<ExperimentListProps> = (props) => {
  // 1. Extract props with default values
  const {
    showFilters = true,
    showActions = true,
    onlyMyExperiments = false,
    initialFilter,
    onExperimentSelect,
    className,
    style
  } = props;

  // 2. Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // 3. Get experiment data and functions from useExperiments hook
  const {
    experiments,
    loading,
    error,
    totalExperiments,
    currentPage,
    pageSize,
    filter,
    fetchExperiments,
    fetchMyExperiments,
    setFilter,
    setPage,
    setSize
  } = useExperiments();

  // 4. Get permission checking functions from usePermissions hook
  const { isPharma, isCRO } = usePermissions();

  // 5. Set up state for filter values
  const [filterValues, setFilterValues] = useState({
    name: '',
    status: '',
    type_id: ''
  });

  // 6. Create memoized table columns configuration with appropriate cell renderers
  const columns = useMemo(() => {
    const baseColumns = [
      {
        field: 'name',
        headerName: 'Name',
        width: 200,
        renderCell: (row: Experiment) => (
          <Link to={`/experiments/${row.id}`} style={{ textDecoration: 'none' }}>
            {row.name}
          </Link>
        )
      },
      {
        field: 'status',
        headerName: 'Status',
        width: 150,
        renderCell: (row: Experiment) => (
          <ExperimentStatus status={row.status as ExperimentStatusEnum} />
        )
      },
      {
        field: 'type',
        headerName: 'Type',
        width: 150,
        valueGetter: (row: Experiment) => row.type?.name || 'Unknown'
      },
      {
        field: 'created_at',
        headerName: 'Date',
        width: 150,
        valueGetter: (row: Experiment) => formatDate(new Date(row.created_at), 'yyyy-MM-dd')
      },
      {
        field: 'molecule_count',
        headerName: 'Molecules',
        width: 100,
        align: 'right'
      }
    ];

    if (showActions) {
      baseColumns.push({
        field: 'actions',
        headerName: 'Actions',
        width: 200,
        renderCell: (row: Experiment) => (
          <Box display="flex" gap={1}>
            <Tooltip title="View Details">
              <IconButton size="small" onClick={() => handleView(row)}>
                <ViewIcon />
              </IconButton>
            </Tooltip>
            {isPharma() && (
              <>
                <Tooltip title="Edit Experiment">
                  <IconButton size="small" onClick={() => handleEdit(row)}>
                    <EditIcon />
                  </IconButton>
                </Tooltip>
                {row.status === ExperimentStatusEnum.DRAFT && (
                  <Tooltip title="Add to Queue">
                    <IconButton size="small" onClick={() => handleQueue(row)}>
                      <QueueIcon />
                    </IconButton>
                  </Tooltip>
                )}
                {row.status === ExperimentStatusEnum.QUEUED && (
                  <Tooltip title="Submit to CRO">
                    <IconButton size="small" onClick={() => handleSubmit(row)}>
                      <SubmitIcon />
                    </IconButton>
                  </Tooltip>
                )}
                {row.status === ExperimentStatusEnum.QUOTE_PENDING && (
                  <Tooltip title="Cancel Experiment">
                    <IconButton size="small" onClick={() => handleCancel(row)}>
                      <CancelIcon />
                    </IconButton>
                  </Tooltip>
                )}
              </>
            )}
            {isCRO() && row.status === ExperimentStatusEnum.SUBMITTED && (
              <Tooltip title="Review Submission">
                <IconButton size="small" onClick={() => handleView(row)}>
                  <ViewIcon />
                </IconButton>
              </Tooltip>
            )}
          </Box>
        )
      });
    }

    return baseColumns;
  }, [showActions, isPharma, isCRO, navigate]);

  // 7. Create memoized action handlers for view, edit, queue, submit, and cancel operations
  const handleView = useCallback((experiment: Experiment) => {
    navigate(`/experiments/${experiment.id}`);
  }, [navigate]);

  const handleEdit = useCallback((experiment: Experiment) => {
    navigate(`/experiments/${experiment.id}/edit`);
  }, [navigate]);

  const handleQueue = useCallback((experiment: Experiment) => {
    // Dispatch queueExperiment action
    console.log('Queue experiment:', experiment.id);
  }, []);

  const handleSubmit = useCallback((experiment: Experiment) => {
    // Dispatch submitExperiment action
    console.log('Submit experiment:', experiment.id);
  }, []);

  const handleCancel = useCallback((experiment: Experiment) => {
    // Dispatch cancelExperiment action
    console.log('Cancel experiment:', experiment.id);
  }, []);

  // 8. Set up effect to fetch experiments on component mount and when filters change
  useEffect(() => {
    if (onlyMyExperiments) {
      fetchMyExperiments();
    } else {
      fetchExperiments();
    }
  }, [fetchExperiments, fetchMyExperiments, onlyMyExperiments, filter]);

  // 9. Set up effect to apply initial filter if provided
  useEffect(() => {
    if (initialFilter) {
      setFilter(initialFilter);
      setFilterValues(initialFilter);
    }
  }, [initialFilter, setFilter]);

  // 10. Create memoized filter change handlers
  const handleFilterChange = useCallback((event: React.ChangeEvent<HTMLInputElement | HTMLSelectElement>) => {
    const { name, value } = event.target;
    setFilterValues(prev => ({ ...prev, [name]: value }));
  }, []);

  const handleApplyFilters = useCallback(() => {
    setFilter(filterValues);
  }, [filterValues, setFilter]);

  // 11. Create memoized pagination handlers
  const handlePageChange = useCallback((newPage: number) => {
    setPage(newPage);
  }, [setPage]);

  const handlePageSizeChange = useCallback((newPageSize: number) => {
    setSize(newPageSize);
  }, [setSize]);

  // 12. Render filter section if showFilters is true
  const renderFilters = () => (
    <Paper sx={{ p: 2, mb: 2 }}>
      <Grid container spacing={2} alignItems="center">
        <Grid item xs={12} sm={4}>
          <TextField
            fullWidth
            label="Name"
            name="name"
            value={filterValues.name || ''}
            onChange={handleFilterChange}
          />
        </Grid>
        <Grid item xs={12} sm={4}>
          <FormControl fullWidth>
            <InputLabel id="status-select-label">Status</InputLabel>
            <Select
              labelId="status-select-label"
              name="status"
              value={filterValues.status || ''}
              label="Status"
              onChange={handleFilterChange}
            >
              <MenuItem value="">All</MenuItem>
              {Object.values(ExperimentStatusEnum).map(status => (
                <MenuItem key={status} value={status}>{status}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Grid>
        <Grid item xs={12} sm={4}>
          <FormControl fullWidth>
            <InputLabel id="type-select-label">Type</InputLabel>
            <Select
              labelId="type-select-label"
              name="type_id"
              value={filterValues.type_id || ''}
              label="Type"
              onChange={handleFilterChange}
            >
              <MenuItem value="">All</MenuItem>
              {/* Assuming experimentTypes is an array of { id: string, name: string } */}
              {experimentTypes.map(type => (
                <MenuItem key={type.id} value={type.id}>{type.name}</MenuItem>
              ))}
            </Select>
          </FormControl>
        </Grid>
        <Grid item xs={12}>
          <Button variant="contained" onClick={handleApplyFilters}>
            Apply Filters
          </Button>
        </Grid>
      </Grid>
    </Paper>
  );

  // 13. Render table with experiment data, columns, and pagination
  return (
    <Box className={className} style={style}>
      {showFilters && renderFilters()}
      <Table
        data={experiments}
        columns={columns}
        loading={loading}
        error={error}
        pagination={{
          page: currentPage,
          pageSize: pageSize,
          totalItems: totalExperiments,
          onPageChange: handlePageChange,
          onPageSizeChange: handlePageSizeChange
        }}
        emptyMessage="No experiments found."
      />
    </Box>
  );
};

ExperimentList.defaultProps = {
  showFilters: true,
  showActions: true,
  onlyMyExperiments: false
};

export default ExperimentList;