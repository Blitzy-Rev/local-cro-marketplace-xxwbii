import React from 'react'; // React v18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom ^6.11.2
import { Box, Typography, Paper, Breadcrumbs, Link } from '@mui/material'; // @mui/material ^5.13.0

import CreateExperimentForm from '../components/CreateExperimentForm';
import useToast from '../../../hooks/useToast';
import MainLayout from '../../../layouts/MainLayout';

/**
 * Page component for creating a new experiment
 */
const CreateExperimentPage: React.FC = () => {
  // LD1: Get navigate function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Get showToast function from useToast hook
  const { showToast } = useToast();

  // LD1: Define handleSuccess function to handle successful experiment creation
  const handleSuccess = (experimentId: string) => {
    // Show success toast notification
    showToast({
      type: 'success',
      message: 'Experiment created successfully!',
    });

    // Navigate to the newly created experiment's detail page
    navigate(`/experiments/${experimentId}`);
  };

  // LD1: Render page layout with breadcrumbs navigation
  return (
    <MainLayout>
      <Box sx={{ padding: 3 }}>
        {/* LD1: Render Breadcrumbs navigation with links to Dashboard, Experiments, and Create Experiment */}
        <Breadcrumbs aria-label="breadcrumb" sx={{ mb: 2 }}>
          <Link underline="hover" color="inherit" href="/">
            Dashboard
          </Link>
          <Link underline="hover" color="inherit" href="/experiments">
            Experiments
          </Link>
          <Typography color="text.primary">Create Experiment</Typography>
        </Breadcrumbs>

        {/* LD1: Render Typography component for page title */}
        <Typography variant="h4" component="h1" gutterBottom>
          Create New Experiment
        </Typography>

        {/* LD1: Render Typography component for page description */}
        <Typography variant="body1" paragraph>
          Define the parameters and molecules for your new experiment.
        </Typography>

        {/* LD1: Render Paper component with padding to contain the form */}
        <Paper elevation={3} sx={{ padding: 3 }}>
          {/* LD1: Render CreateExperimentForm component with onSuccess handler */}
          <CreateExperimentForm onSuccess={handleSuccess} />
        </Paper>
      </Box>
    </MainLayout>
  );
};

export default CreateExperimentPage;