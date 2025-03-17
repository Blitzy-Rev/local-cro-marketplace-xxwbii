import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useNavigate, useParams } from 'react-router-dom'; // react-router-dom v6.4+
import { Container, Typography, Paper, Box, Breadcrumbs, Link } from '@mui/material'; // @mui/material v5.13+

import MainLayout from '../../../layouts/MainLayout';
import SubmissionForm from '../components/SubmissionForm';
import useSubmissions from '../hooks/useSubmissions';
import useExperiments from '../../experiments/hooks/useExperiments';
import useToast from '../../../hooks/useToast';
import Loading from '../../../components/common/Loading';
import { Submission } from '../../../types/submission';

/**
 * Page component for creating new submissions to CROs.
 * This component allows users to select an experiment, choose a CRO, add submission details, and submit the form.
 * It handles navigation, form submission, and error handling.
 */
const CreateSubmissionPage: React.FC = () => {
  // LD1: Get navigation function for redirecting after submission
  const navigate = useNavigate();

  // LD1: Get URL parameters including experimentId if provided
  const { experimentId: experimentIdParam } = useParams<{ experimentId?: string }>();

  // LD1: Use state to store the experiment ID, which can come from URL params
  const [experimentId, setExperimentId] = useState<number | null>(null);

  // LD1: Access submission-related functionality using the useSubmissions hook
  const { createSubmission, loading: submissionLoading } = useSubmissions();

  // LD1: Access experiment data and functionality using the useExperiments hook
  const { fetchExperimentById, currentExperiment, loading: experimentLoading } = useExperiments();

  // LD1: Access toast notification functionality using the useToast hook
  const { showToast } = useToast();

  // LD1: Extract experimentId from URL parameters
  useEffect(() => {
    if (experimentIdParam) {
      setExperimentId(Number(experimentIdParam));
    }
  }, [experimentIdParam]);

  // LD1: Fetch experiment details when experimentId changes
  useEffect(() => {
    if (experimentId) {
      fetchExperimentById(experimentId.toString());
    }
  }, [experimentId, fetchExperimentById]);

  /**
   * Handles the submission of the form
   * @param submission - The submission data
   */
  const handleSubmit = async (submission: Submission) => {
    // LD1: Call createSubmission function from useSubmissions hook with the submission data
    await createSubmission(submission);

    // LD1: Show success toast notification
    showToast({ type: 'success', message: 'Submission created successfully!' });

    // LD1: Navigate to the submissions list page
    navigate('/app/submissions');
  };

  /**
   * Handles cancellation of the form
   */
  const handleCancel = () => {
    // LD1: Navigate back to the submissions list page
    navigate('/app/submissions');
  };

  // LD1: Render the component
  return (
    <MainLayout>
      <Container maxWidth="md">
        <Box sx={{ mt: 4, mb: 2 }}>
          <Breadcrumbs aria-label="breadcrumb">
            <Link underline="hover" color="inherit" href="/app/dashboard">
              Dashboard
            </Link>
            <Link underline="hover" color="inherit" href="/app/experiments">
              Experiments
            </Link>
            <Typography color="text.primary">Create Submission</Typography>
          </Breadcrumbs>
        </Box>
        <Typography variant="h4" gutterBottom>
          Create Submission
        </Typography>
        {(submissionLoading || experimentLoading) ? (
          <Loading />
        ) : (
          <Paper elevation={3} style={{ padding: '20px', width: '100%' }}>
            <SubmissionForm
              experimentId={experimentId}
              onSubmit={handleSubmit}
              onCancel={handleCancel}
            />
          </Paper>
        )}
      </Container>
    </MainLayout>
  );
};

export default CreateSubmissionPage;