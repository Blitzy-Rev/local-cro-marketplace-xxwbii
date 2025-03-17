import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.14.0
import { Box, Typography, Container, Paper, Grid, Divider } from '@mui/material'; // @mui/material v5.13.0
import { Add } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import SubmissionList from '../components/SubmissionList';
import SubmissionDetail from '../components/SubmissionDetail';
import useSubmissions from '../hooks/useSubmissions';
import useAuth from '../../auth/hooks/useAuth';
import usePermissions from '../../../hooks/usePermissions';
import useToast from '../../../hooks/useToast';

/**
 * Main page component for viewing and managing submissions
 * @returns Rendered submissions page component
 */
const SubmissionsPage: React.FC = () => {
  // LD1: Initialize state for selected submission ID
  const [selectedSubmissionId, setSelectedSubmissionId] = useState<string | null>(null);

  // LD1: Initialize state for view mode (list or detail)
  const [isDetailView, setIsDetailView] = useState(false);

  // LD1: Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Get user authentication state and role from useAuth and usePermissions hooks
  const { user, isAuthenticated } = useAuth();
  const { isPharma, isCRO } = usePermissions();

  // LD1: Initialize useSubmissions hook for submission data and operations
  const {
    submissions,
    loading,
    error,
    filter,
    totalSubmissions,
    pagination,
    setFilter,
    refreshSubmissions,
    fetchSubmissionDetails,
    createSubmission,
    updateSubmission,
    updateSubmissionStatus,
    provideQuote,
    respondToQuote,
    cancelSubmission,
    fetchByExperiment,
    fetchByCRO,
    fetchByStatus,
    clearError
  } = useSubmissions();

  // LD1: Initialize useToast hook for notifications
  const { showToast } = useToast();

  /**
   * LD1: Define handler for selecting a submission
   * @param submissionId - The ID of the selected submission
   */
  const handleSubmissionSelect = useCallback((submissionId: string) => {
    setSelectedSubmissionId(submissionId);
    setIsDetailView(true);
  }, []);

  /**
   * LD1: Define handler for returning to list view
   */
  const handleBackToList = useCallback(() => {
    setSelectedSubmissionId(null);
    setIsDetailView(false);
  }, []);

  /**
   * LD1: Define handler for creating a new submission
   */
  const handleCreateSubmission = useCallback(() => {
    navigate('/submissions/new');
  }, [navigate]);

  /**
   * LD1: Define handler for submission status changes
   */
  const handleStatusChange = useCallback(() => {
    refreshSubmissions();
  }, [refreshSubmissions]);

  // LD1: Render page title and description
  return (
    <Container maxWidth="lg">
      <Box sx={{ mt: 4, mb: 2 }}>
        <Typography variant="h4" component="h1" gutterBottom>
          Submissions
        </Typography>
        <Typography variant="body1">
          View and manage your submissions to Contract Research Organizations (CROs).
        </Typography>
      </Box>

      {/* LD1: Render SubmissionList component when in list view */}
      {!isDetailView && (
        <SubmissionList
          onSubmissionSelect={handleSubmissionSelect}
          onCreateSubmission={handleCreateSubmission}
        />
      )}

      {/* LD1: Render SubmissionDetail component when in detail view */}
      {isDetailView && selectedSubmissionId && (
        <SubmissionDetail
          submissionId={selectedSubmissionId}
          onBack={handleBackToList}
          onStatusChange={handleStatusChange}
        />
      )}

      {/* LD1: Handle loading and error states appropriately */}
      {loading && <Typography>Loading submissions...</Typography>}
      {error && <Typography color="error">Error: {error}</Typography>}
    </Container>
  );
};

export default SubmissionsPage;