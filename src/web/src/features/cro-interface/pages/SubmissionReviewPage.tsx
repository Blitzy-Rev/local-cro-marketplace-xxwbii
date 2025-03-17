import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useParams, useNavigate } from 'react-router-dom'; // react-router-dom v6.4+
import { Container, Paper, Typography, Box } from '@mui/material'; // @mui/material v5.13.0

import CROLayout from '../../../layouts/CROLayout';
import SubmissionReview from '../components/SubmissionReview';
import useCROInterface from '../hooks/useCROInterface';
import Loading from '../../../components/common/Loading';
import useToast from '../../../hooks/useToast';

/**
 * Page component for reviewing submission details and providing quotes
 */
const SubmissionReviewPage: React.FC = () => {
  // LD1: Get submissionId from URL parameters using useParams hook
  const { submissionId } = useParams<{ submissionId: string }>();

  // LD1: Get navigation function using useNavigate hook
  const navigate = useNavigate();

  // LD1: Get refreshSubmissions function from useCROInterface hook
  const { refreshSubmissions } = useCROInterface();

  // LD1: Get showToast function from useToast hook
  const { showToast } = useToast();

  // LD1: Initialize loading state
  const [loading, setLoading] = useState<boolean>(false);

  /**
   * Define handleQuoteSubmitted function to handle successful quote submission
   */
  const handleQuoteSubmitted = useCallback(() => {
    // Show success toast notification
    showToast({
      type: 'success',
      message: 'Quote submitted successfully!',
    });

    // Refresh submissions list to update data
    refreshSubmissions();

    // Navigate back to submissions list page
    navigate('/cro/submissions');
  }, [showToast, refreshSubmissions, navigate]);

  /**
   * Define handleDeclineSubmission function to handle submission decline
   */
  const handleDeclineSubmission = useCallback(() => {
    // Show info toast notification about declined submission
    showToast({
      type: 'info',
      message: 'Submission declined.',
    });

    // Refresh submissions list to update data
    refreshSubmissions();

    // Navigate back to submissions list page
    navigate('/cro/submissions');
  }, [showToast, refreshSubmissions, navigate]);

  /**
   * Define handleBack function to navigate back to submissions list
   */
  const handleBack = useCallback(() => {
    // Navigate to the CRO submissions list page
    navigate('/cro/submissions');
  }, [navigate]);

  // LD1: Render CROLayout as the base layout
  return (
    <CROLayout>
      {/* LD1: Render Container component with appropriate sizing */}
      <Container maxWidth="lg">
        {/* LD1: Render Paper component for the content container */}
        <Paper elevation={3} style={{ padding: '20px' }}>
          {/* LD1: If submissionId is not available, display appropriate message */}
          {!submissionId ? (
            <Typography variant="h6" color="error" align="center">
              {/* LD1: Loading component if in loading state */}
              {loading ? <Loading message="Loading..." /> : 'Submission ID is missing.'}
            </Typography>
          ) : (
            // LD1: If submissionId is available, render SubmissionReview component with:
            //   - submissionId prop
            //   - onQuoteSubmitted callback
            //   - onDecline callback
            //   - onBack callback
            <SubmissionReview
              submissionId={submissionId}
              onQuoteSubmitted={handleQuoteSubmitted}
              onDecline={handleDeclineSubmission}
              onBack={handleBack}
            />
          )}
        </Paper>
      </Container>
    </CROLayout>
  );
};

export default SubmissionReviewPage;