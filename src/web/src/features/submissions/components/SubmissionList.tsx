import React, { useState, useEffect, useMemo, useCallback } from 'react'; // react 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.10.0
import { Box, Typography, Grid, Paper, Stack, Divider } from '@mui/material'; // @mui/material v5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13.0
import { Add, FilterList, Refresh } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import Table from '../../../components/common/Table';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Dropdown from '../../../components/common/Dropdown';
import Badge from '../../../components/common/Badge';
import useSubmissions from '../hooks/useSubmissions';
import { Submission, SubmissionStatus } from '../../../types/submission';
import usePermissions from '../../../hooks/usePermissions';
import { formatDate } from '../../../utils/formatters';

/**
 * Interface for the SubmissionList component props
 */
interface SubmissionListProps {
  onSubmissionSelect?: (submissionId: string) => void;
  onCreateSubmission?: () => void;
  experimentId?: string;
  croId?: string;
  status?: SubmissionStatus;
  className?: string;
}

/**
 * Styled component for the header of the submission list
 */
const Header = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  padding: theme.spacing(2),
}));

/**
 * Styled component for the filter section
 */
const FilterSection = styled(Stack)(({ theme }) => ({
  padding: theme.spacing(2),
  backgroundColor: theme.palette.background.paper,
  borderRadius: theme.shape.borderRadius,
}));

/**
 * Styled component for the submission status badge
 */
const StatusBadge = styled(Badge)(({ theme }) => ({
  marginLeft: theme.spacing(1),
}));

/**
 * Determines the appropriate color for a submission status badge
 * @param status The submission status
 * @returns Color name for the badge (success, warning, error, info, etc.)
 */
const getSubmissionStatusColor = (status: SubmissionStatus): string => {
  switch (status) {
    case SubmissionStatus.APPROVED:
    case SubmissionStatus.COMPLETED:
      return 'success';
    case SubmissionStatus.REJECTED:
    case SubmissionStatus.CANCELLED:
      return 'error';
    case SubmissionStatus.QUOTE_PROVIDED:
    case SubmissionStatus.QUOTE_REJECTED:
      return 'warning';
    case SubmissionStatus.IN_PROGRESS:
      return 'info';
    default:
      return 'default';
  }
};

/**
 * Converts submission status enum value to a user-friendly display label
 * @param status Human-readable status label
 */
const getSubmissionStatusLabel = (status: SubmissionStatus): string => {
  switch (status) {
    case SubmissionStatus.PENDING:
      return 'Pending';
    case SubmissionStatus.REJECTED:
      return 'Rejected';
    case SubmissionStatus.QUOTE_PROVIDED:
      return 'Quote Received';
    case SubmissionStatus.QUOTE_REJECTED:
      return 'Quote Rejected';
    case SubmissionStatus.APPROVED:
      return 'Approved';
    case SubmissionStatus.IN_PROGRESS:
      return 'In Progress';
    case SubmissionStatus.COMPLETED:
      return 'Completed';
    case SubmissionStatus.CANCELLED:
      return 'Cancelled';
    default:
      return 'Unknown';
  }
};

/**
 * Component that displays a list of submissions with filtering and pagination
 * @param props The component props
 * @returns The rendered submission list component
 */
const SubmissionList: React.FC<SubmissionListProps> = ({
  onSubmissionSelect,
  onCreateSubmission,
  experimentId,
  croId,
  status,
  className,
}) => {
  // Initialize useSubmissions hook with initial filter based on props
  const {
    submissions,
    loading,
    error,
    filter,
    totalSubmissions,
    pagination,
    setFilter,
    refreshSubmissions,
  } = useSubmissions({
    initialFilter: {
      experiment_id: experimentId,
      cro_id: croId ? parseInt(croId, 10) : undefined,
      status: status,
    },
    enableAutoFetch: true,
  });

  // Initialize usePermissions hook to check user role
  const { isPharma, isCRO } = usePermissions();

  // Initialize useNavigate hook for navigation
  const navigate = useNavigate();

  // Define table columns with appropriate headers and cell renderers
  const columns = useMemo(
    () => [
      {
        field: 'id',
        headerName: 'ID',
        width: 80,
      },
      {
        field: 'experiment.name',
        headerName: 'Experiment',
        width: 200,
        valueGetter: (row: Submission) => row.experiment?.name,
      },
      {
        field: 'cro.email',
        headerName: 'CRO',
        width: 200,
        valueGetter: (row: Submission) => row.cro?.email,
      },
      {
        field: 'status',
        headerName: 'Status',
        width: 150,
        renderCell: (row: Submission) => (
          <StatusBadge
            label={getSubmissionStatusLabel(row.status)}
            color={getSubmissionStatusColor(row.status)}
            size="small"
          />
        ),
      },
      {
        field: 'submitted_at',
        headerName: 'Date',
        width: 150,
        valueGetter: (row: Submission) => formatDate(row.submitted_at),
      },
    ],
    []
  );

  // Define status filter options based on SubmissionStatus enum
  const statusFilterOptions = useMemo(
    () =>
      Object.values(SubmissionStatus).map((status) => ({
        value: status,
        label: getSubmissionStatusLabel(status),
      })),
    []
  );

  // Handle filter changes for experiment, CRO, and status
  const handleFilterChange = useCallback(
    (newFilter: Partial<typeof filter>) => {
      setFilter({ ...filter, ...newFilter });
    },
    [filter, setFilter]
  );

  // Handle row click to navigate to submission details
  const handleRowClick = useCallback(
    (submissionId: string) => {
      if (onSubmissionSelect) {
        onSubmissionSelect(submissionId);
      } else {
        navigate(`/submissions/${submissionId}`);
      }
    },
    [navigate, onSubmissionSelect]
  );

  // Handle create submission button click
  const handleCreateSubmissionClick = useCallback(() => {
    if (onCreateSubmission) {
      onCreateSubmission();
    } else {
      navigate('/submissions/new');
    }
  }, [navigate, onCreateSubmission]);

  // Handle refresh button click
  const handleRefreshClick = useCallback(() => {
    refreshSubmissions();
  }, [refreshSubmissions]);

  return (
    <Card className={className} fullHeight>
      <Header>
        <Typography variant="h6">Submissions</Typography>
        <Box>
          {/* Render action buttons for creating submissions and refreshing */}
          {isPharma() && (
            <Button
              variant="contained"
              color="primary"
              startIcon={<Add />}
              onClick={handleCreateSubmissionClick}
            >
              Create Submission
            </Button>
          )}
          <Button
            variant="outlined"
            color="primary"
            startIcon={<Refresh />}
            onClick={handleRefreshClick}
            style={{ marginLeft: '8px' }}
          >
            Refresh
          </Button>
        </Box>
      </Header>
      <Divider />
      {/* Render filter controls for status, experiment, and CRO */}
      <FilterSection direction="row" spacing={2}>
        <Dropdown
          id="status-filter"
          name="status"
          label="Status"
          value={filter.status || ''}
          onChange={(e) => handleFilterChange({ status: e.target.value as SubmissionStatus })}
          options={statusFilterOptions}
          placeholder="All Statuses"
          clearable
        />
      </FilterSection>
      {/* Render Table component with submissions data */}
      <Table
        data={submissions}
        columns={columns}
        loading={loading}
        error={error}
        pagination={{
          ...pagination,
          onPageChange: pagination.goToPage,
          onPageSizeChange: pagination.setPageSize,
        }}
        sortable
        onSortChange={(field, direction) =>
          setFilter({ sort_by: field, sort_desc: direction === 'desc' })
        }
        emptyMessage="No submissions found."
        onClick={(row) => handleRowClick(row.id)}
        style={{ cursor: 'pointer' }}
      />
    </Card>
  );
};

export default SubmissionList;