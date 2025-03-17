import React, { useState, useEffect, useCallback } from 'react'; // React, { useState, useEffect, useCallback } ^18.2.0
import {
  Box,
  Typography,
  Grid,
  Tabs,
  Tab,
  TextField,
  MenuItem,
  Select,
  FormControl,
  InputLabel,
  Chip
} from '@mui/material'; // @mui/material ^5.13.0
import { Search, FilterList, Refresh } from '@mui/icons-material'; // @mui/icons-material ^5.13.0
import { useNavigate, useParams } from 'react-router-dom'; // { useNavigate, useParams } ^6.4.0

import useCROInterface from '../hooks/useCROInterface';
import SubmissionReview from '../components/SubmissionReview';
import CROLayout from '../../../layouts/CROLayout';
import Table from '../../../components/common/Table';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import useToast from '../../../hooks/useToast';
import { SubmissionStatus } from '../../../types/submission';

/**
 * Main component for displaying and managing CRO submissions
 * @returns {JSX.Element} Rendered component
 */
const CROSubmissionsPage: React.FC = () => {
  // LD1: Initialize state for selected submission ID
  const [selectedSubmissionId, setSelectedSubmissionId] = useState<string | null>(null);
  // LD1: Initialize state for current tab index
  const [tabIndex, setTabIndex] = useState<number>(0);
  // LD1: Initialize state for search term
  const [searchTerm, setSearchTerm] = useState<string>('');

  // LD1: Get navigation function from useNavigate hook
  const navigate = useNavigate();
  // LD1: Get URL parameters from useParams hook
  const { submissionId } = useParams<{ submissionId?: string }>();
  // LD1: Get toast notification function from useToast hook
  const { showToast } = useToast();

  // LD1: Get CRO interface functionality from useCROInterface hook
  const {
    submissions,
    loading,
    error,
    filters,
    setFilters,
    refreshSubmissions,
    handleProvideQuote,
    handleUpdateStatus
  } = useCROInterface();

  /**
   * Define handleTabChange function to switch between submission status tabs
   * @param {React.SyntheticEvent} event - The event object
   * @param {number} newValue - The new tab index
   */
  const handleTabChange = (event: React.SyntheticEvent, newValue: number) => {
    setTabIndex(newValue);
    // Update filters based on the selected tab
    switch (newValue) {
      case 0: // All
        setFilters({});
        break;
      case 1: // Pending
        setFilters({ status: SubmissionStatus.PENDING });
        break;
      case 2: // Approved
        setFilters({ status: SubmissionStatus.APPROVED });
        break;
      case 3: // In Progress
        setFilters({ status: SubmissionStatus.IN_PROGRESS });
        break;
      case 4: // Completed
        setFilters({ status: SubmissionStatus.COMPLETED });
        break;
      default:
        setFilters({});
        break;
    }
  };

  /**
   * Define handleSearchChange function to update search term
   * @param {React.ChangeEvent<HTMLInputElement>} event - The event object
   */
  const handleSearchChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setSearchTerm(event.target.value);
    // Implement search logic here (e.g., update filters)
  };

  /**
   * Define handleFilterChange function to update submission filters
   * @param {string} filterName - The name of the filter to update
   * @param {any} value - The new value for the filter
   */
  const handleFilterChange = (filterName: string, value: any) => {
    setFilters({ ...filters, [filterName]: value });
  };

  /**
   * Define handleRefresh function to refresh submission data
   */
  const handleRefresh = () => {
    refreshSubmissions();
  };

  /**
   * Define handleViewSubmission function to view submission details
   * @param {string} submissionId - The ID of the submission to view
   */
  const handleViewSubmission = (submissionId: string) => {
    setSelectedSubmissionId(submissionId);
    navigate(`/cro/submissions/${submissionId}`);
  };

  /**
   * Define handleBackToList function to return to submission list view
   */
  const handleBackToList = () => {
    setSelectedSubmissionId(null);
    navigate('/cro/submissions');
  };

  /**
   * Define handleQuoteSubmitted function to handle successful quote submission
   */
  const handleQuoteSubmitted = () => {
    showToast({ type: 'success', message: 'Quote submitted successfully!' });
    refreshSubmissions();
  };

  /**
   * Define handleDeclineSubmission function to decline a submission
   */
  const handleDeclineSubmission = () => {
    showToast({ type: 'info', message: 'Submission declined.' });
    refreshSubmissions();
  };

  /**
   * Define getStatusChip function to render status chips with appropriate colors
   * @param {string} status - The submission status
   */
  const getStatusChip = (status: string) => {
    switch (status) {
      case SubmissionStatus.PENDING:
        return <Chip label={status} color="primary" />;
      case SubmissionStatus.APPROVED:
        return <Chip label={status} color="success" />;
      case SubmissionStatus.REJECTED:
        return <Chip label={status} color="error" />;
      case SubmissionStatus.IN_PROGRESS:
        return <Chip label={status} color="info" />;
      case SubmissionStatus.COMPLETED:
        return <Chip label={status} color="success" />;
      default:
        return <Chip label={status} />;
    }
  };

  /**
   * Define table columns configuration for submission list
   */
  const columns = React.useMemo(
    () => [
      {
        field: 'id',
        headerName: 'ID',
        width: 100,
      },
      {
        field: 'experiment.name',
        headerName: 'Experiment',
        width: 200,
        valueGetter: (params: any) => params.row.experiment?.name,
      },
      {
        field: 'experiment.type.name',
        headerName: 'Type',
        width: 150,
        valueGetter: (params: any) => params.row.experiment?.type?.name,
      },
      {
        field: 'submitted_at',
        headerName: 'Date',
        width: 150,
      },
      {
        field: 'status',
        headerName: 'Status',
        width: 150,
        renderCell: (params: any) => getStatusChip(params.row.status),
      },
      {
        field: 'actions',
        headerName: 'Actions',
        width: 150,
        renderCell: (params: any) => (
          <Button size="small" onClick={() => handleViewSubmission(params.row.id)}>
            View
          </Button>
        ),
      },
    ],
    [handleViewSubmission, getStatusChip]
  );

  // LD1: Set up effect to update filters when tab changes
  useEffect(() => {
    if (submissionId) {
      setSelectedSubmissionId(submissionId);
    }
  }, [submissionId]);

  // LD1: Set up effect to initialize selected submission from URL parameters
  useEffect(() => {
    if (submissionId) {
      setSelectedSubmissionId(submissionId);
    }
  }, [submissionId]);

  // LD1: Render CROLayout as the page container
  return (
    <CROLayout>
      {/* LD1: If selectedSubmissionId is set, render SubmissionReview component */}
      {selectedSubmissionId ? (
        <SubmissionReview
          submissionId={selectedSubmissionId}
          onQuoteSubmitted={handleQuoteSubmitted}
          onDecline={handleDeclineSubmission}
          onBack={handleBackToList}
        />
      ) : (
        // LD1: Otherwise, render submission list view with:
        <Box>
          {/* - Page header with title and refresh button */}
          <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
            <Typography variant="h4">CRO Submissions</Typography>
            <Button variant="outlined" startIcon={<Refresh />} onClick={handleRefresh}>
              Refresh
            </Button>
          </Box>

          {/* - Filter section with tabs for different submission statuses */}
          <Card>
            <Box sx={{ borderBottom: 1, borderColor: 'divider' }}>
              <Tabs value={tabIndex} onChange={handleTabChange} aria-label="submission status tabs">
                <Tab label="All" />
                <Tab label="Pending" />
                <Tab label="Approved" />
                <Tab label="In Progress" />
                <Tab label="Completed" />
              </Tabs>
            </Box>

            {/* - Search and filter controls */}
            <Box p={2} display="flex" alignItems="center" justifyContent="space-between">
              <TextField
                size="small"
                label="Search"
                variant="outlined"
                value={searchTerm}
                onChange={handleSearchChange}
                InputProps={{
                  startAdornment: <Search />,
                }}
              />
              <FormControl variant="outlined" size="small">
                <InputLabel id="filter-sort-label">Sort By</InputLabel>
                <Select
                  labelId="filter-sort-label"
                  id="filter-sort"
                  value={filters.sortBy || ''}
                  onChange={(e) => handleFilterChange('sortBy', e.target.value)}
                  label="Sort By"
                >
                  <MenuItem value="">None</MenuItem>
                  <MenuItem value="date">Date</MenuItem>
                  <MenuItem value="experiment">Experiment</MenuItem>
                </Select>
              </FormControl>
            </Box>

            {/* - Table component displaying submissions with defined columns */}
            <Table
              columns={columns}
              data={submissions || []}
              loading={loading}
              error={error}
              rowClick={handleViewSubmission}
            />
          </Card>

          {/* - Loading state and error handling */}
          {loading && <Typography>Loading submissions...</Typography>}
          {error && <Typography color="error">Error: {error}</Typography>}
        </Box>
      )}
    </CROLayout>
  );
};

// LD1: Export the CRO submissions page component for use in routing
export default CROSubmissionsPage;