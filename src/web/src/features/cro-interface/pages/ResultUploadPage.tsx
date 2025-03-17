import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useParams, useNavigate } from 'react-router-dom'; // react-router-dom v6.10.0
import { Box, Typography, Paper, CircularProgress, Alert } from '@mui/material'; // @mui/material v5.13.0

import CROLayout from '../../../../layouts/CROLayout';
import ResultUploadForm from '../../components/ResultUploadForm';
import { useCROInterface } from '../../hooks/useCROInterface';
import { useToast } from '../../../../hooks/useToast';
import { getSubmission } from '../../../../api/submissions';
import { SubmissionDetailed, SubmissionStatus } from '../../../../types/submission';

/**
 * ResultUploadPage component for CRO users to upload experimental results
 * @returns {JSX.Element} - The rendered component
 */
const ResultUploadPage: React.FC = () => {
  // Extract submissionId from URL parameters using useParams hook
  const { submissionId } = useParams<{ submissionId: string }>();

  // Initialize state for submission data, loading status, and error handling
  const [submission, setSubmission] = useState<SubmissionDetailed | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);

  // Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // Custom hook for CRO interface functionality
  const {  } = useCROInterface();

  // Custom hook for displaying toast notifications
  const { showToast } = useToast();

  /**
   * Fetches the submission data from the API
   */
  const fetchSubmission = useCallback(async () => {
    // Set loading state to true
    setLoading(true);
    // Clear any previous errors
    setError(null);

    try {
      // Call getSubmission API with submissionId
      const response = await getSubmission(submissionId as string);

      // If successful, store submission data in state
      setSubmission(response);

      // If submission is not in IN_PROGRESS status, set error message
      if (response.status !== SubmissionStatus.IN_PROGRESS) {
        setError('Submission is not in progress and cannot be updated.');
      }
    } catch (err: any) {
      // If API call fails, set error message
      setError(err.message || 'Failed to load submission details.');
    } finally {
      // Set loading state to false regardless of outcome
      setLoading(false);
    }
  }, [submissionId]);

  /**
   * Handler for successful result upload
   */
  const handleSuccess = useCallback(() => {
    // Show success toast notification
    showToast({ type: 'success', message: 'Results uploaded successfully!' });

    // Navigate back to submission details page
    navigate(`/cro/submissions/${submissionId}`);
  }, [navigate, submissionId, showToast]);

  /**
   * Handler for cancelling result upload
   */
  const handleCancel = useCallback(() => {
    // Navigate back to submission details page without saving
    navigate(`/cro/submissions/${submissionId}`);
  }, [navigate, submissionId]);

  // Fetch submission details when component mounts or submissionId changes
  useEffect(() => {
    if (submissionId) {
      fetchSubmission();
    }
  }, [submissionId, fetchSubmission]);

  return (
    <CROLayout>
      <Paper elevation={1} sx={{ p: 2, mb: 3 }}>
        <Typography variant="h6" gutterBottom>
          Upload Results
        </Typography>
        {loading && (
          <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '200px' }}>
            <CircularProgress />
          </Box>
        )}
        {error && (
          <Alert severity="error">
            {error}
          </Alert>
        )}
        {submission && submission.status === SubmissionStatus.IN_PROGRESS && (
          <ResultUploadForm submission={submission} onSuccess={handleSuccess} onCancel={handleCancel} />
        )}
      </Paper>
    </CROLayout>
  );
};

export default ResultUploadPage;