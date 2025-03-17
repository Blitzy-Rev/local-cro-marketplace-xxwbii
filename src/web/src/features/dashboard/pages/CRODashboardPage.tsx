import React, { useEffect } from 'react'; // React 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.10+
import { Grid, Box, Typography, useTheme, useMediaQuery } from '@mui/material'; // @mui/material v5.13+
import Stats from '../components/Stats';
import QuickActions from '../components/QuickActions';
import SubmissionReview from '../../cro-interface/components/SubmissionReview';
import ExperimentManagement from '../../cro-interface/components/ExperimentManagement';
import ResultUploadForm from '../../cro-interface/components/ResultUploadForm';
import CROCommunications from '../components/CROCommunications';
import useCROInterface from '../../cro-interface/hooks/useCROInterface';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Main dashboard page component that displays various dashboard widgets for CRO users
 * @returns Rendered CRO dashboard page
 */
const CRODashboardPage: React.FC = () => {
  // LD1: Initialize navigate function using useNavigate hook
  const navigate = useNavigate();

  // LD1: Initialize theme using useTheme hook
  const theme = useTheme();

  // LD1: Check if screen is mobile using useMediaQuery hook
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // LD1: Check user role using usePermissions hook
  const { isPharma, isCRO, isAdmin } = usePermissions();

  // LD1: Initialize state for new submissions count
  const [newSubmissionsCount, setNewSubmissionsCount] = React.useState(0);

  // LD1: Get submissions data and loading state from useCROInterface hook
  const { submissions, loading: submissionsLoading, fetchCROSubmissions } = useCROInterface();

  // LD1: Redirect pharma users to pharma dashboard if they somehow access this page
  useEffect(() => {
    if (isPharma()) {
      navigate('/app');
    }
  }, [isPharma, navigate]);

  // LD1: Redirect admin users to admin dashboard if they somehow access this page
  useEffect(() => {
    if (isAdmin()) {
      navigate('/admin');
    }
  }, [isAdmin, navigate]);

  // LD1: Use useEffect to fetch CRO submissions on component mount
  useEffect(() => {
    fetchCROSubmissions();
  }, [fetchCROSubmissions]);

  // LD1: Define navigation handlers for each dashboard widget
  const handleNavigateToSubmissions = () => {
    navigate('/cro/submissions');
  };

  const handleNavigateToInProgress = () => {
    navigate('/cro/experiments');
  };

  const handleNavigateToCompleted = () => {
    navigate('/cro/results');
  };

  const handleNavigateToCommunications = () => {
    navigate('/cro/communications');
  };

  const handleNavigateToReviewSubmissions = () => {
    navigate('/cro/submissions');
  };

  const handleNavigateToUploadResults = () => {
    navigate('/cro/results/upload');
  };

  const handleNavigateToUpdateStatus = () => {
    navigate('/cro/experiments/status');
  };

  // LD1: Render page title and greeting
  return (
    <Box sx={{ flexGrow: 1, padding: 3 }}>
      <Typography variant="h4" gutterBottom>
        CRO Dashboard
      </Typography>
      <Typography variant="subtitle1">
        Welcome, CRO User! Here's an overview of your tasks.
      </Typography>

      {/* LD1: Render Stats component with CRO-specific statistics */}
      <Stats />

      {/* LD1: Render QuickActions component with CRO-specific actions */}
      <Grid container spacing={3} mt={3}>
        <Grid item xs={12} md={4}>
          <QuickActions />
        </Grid>
        <Grid item xs={12} md={8}>
          <Grid container spacing={3}>
            <Grid item xs={12} md={6}>
              <SubmissionReview
                submissionId="123"
                onQuoteSubmitted={() => console.log('Quote Submitted')}
                onDecline={() => console.log('Submission Declined')}
                onBack={() => console.log('Back to Submissions')}
              />
            </Grid>
            <Grid item xs={12} md={6}>
              <ExperimentManagement />
            </Grid>
            <Grid item xs={12} md={6}>
              <ResultUploadForm
                submission={{}}
                onSuccess={() => console.log('Results Uploaded')}
                onCancel={() => console.log('Upload Cancelled')}
              />
            </Grid>
            <Grid item xs={12} md={6}>
              <CROCommunications />
            </Grid>
          </Grid>
        </Grid>
      </Grid>
    </Box>
  );
};

export default CRODashboardPage;