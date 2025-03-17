import React, { useEffect } from 'react'; // React 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.10+
import { Grid, Box, Typography, useTheme, useMediaQuery } from '@mui/material'; // @mui/material v5.13+

// Internal imports
import Stats from '../components/Stats';
import QuickActions from '../components/QuickActions';
import ActiveExperiments from '../components/ActiveExperiments';
import RecentResults from '../components/RecentResults';
import MoleculeLibraries from '../components/MoleculeLibraries';
import CROCommunications from '../components/CROCommunications';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Main dashboard page component that displays various dashboard widgets for Pharma users
 * @returns Rendered dashboard page
 */
const DashboardPage: React.FC = () => {
  // LD1: Initialize navigate function using useNavigate hook
  const navigate = useNavigate();

  // LD1: Initialize theme using useTheme hook
  const theme = useTheme();

  // LD1: Check if screen is mobile using useMediaQuery hook
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // LD1: Check user role using usePermissions hook
  const { isPharma, isCRO, isAdmin } = usePermissions();

  // LD1: Redirect admin users to admin dashboard if they somehow access this page
  useEffect(() => {
    if (isAdmin()) {
      navigate('/admin');
    }
  }, [isAdmin, navigate]);

  // LD1: Redirect CRO users to CRO dashboard if they somehow access this page
  useEffect(() => {
    if (isCRO()) {
      navigate('/cro');
    }
  }, [isCRO, navigate]);

  // LD1: Define navigation handlers for each dashboard widget
  const handleNavigateToExperiments = () => {
    navigate('/app/experiments');
  };

  const handleNavigateToResults = () => {
    navigate('/app/results');
  };

  const handleNavigateToLibraries = () => {
    navigate('/app/libraries');
  };

  const handleNavigateToCommunications = () => {
    navigate('/app/communications');
  };

  const handleNavigateToCSVUpload = () => {
    navigate('/app/molecules/upload');
  };

  const handleNavigateToCreateLibrary = () => {
    navigate('/app/libraries/create');
  };

  const handleNavigateToCreateExperiment = () => {
    navigate('/app/experiments/create');
  };

  // LD1: Render page title and greeting
  return (
    <Box sx={{ flexGrow: 1, padding: theme.spacing(3) }}>
      <Typography variant="h4" component="h1" gutterBottom>
        Dashboard
      </Typography>
      <Typography variant="subtitle1" paragraph>
        Welcome to the Molecular Data Management Platform!
      </Typography>

      {/* LD1: Render Stats component with Pharma-specific statistics */}
      {isPharma() && <Stats />}

      {/* LD1: Render QuickActions component with Pharma-specific actions */}
      {isPharma() && (
        <QuickActions
          onCreateLibrary={handleNavigateToCreateLibrary}
          onUploadCSV={handleNavigateToCSVUpload}
          onCreateExperiment={handleNavigateToCreateExperiment}
        />
      )}

      {/* LD1: Render grid layout with ActiveExperiments, RecentResults, MoleculeLibraries, and CROCommunications components */}
      <Grid container spacing={isMobile ? 2 : 3}>
        <Grid item xs={12} md={6}>
          <ActiveExperiments onViewAll={handleNavigateToExperiments} />
        </Grid>
        <Grid item xs={12} md={6}>
          <RecentResults onViewAll={handleNavigateToResults} />
        </Grid>
        <Grid item xs={12} md={6}>
          <MoleculeLibraries
            onViewAll={handleNavigateToLibraries}
            onCreateLibrary={handleNavigateToCreateLibrary}
            maxItems={3}
          />
        </Grid>
        <Grid item xs={12} md={6}>
          <CROCommunications onOpenMessages={handleNavigateToCommunications} />
        </Grid>
      </Grid>
    </Box>
  );
};

export default DashboardPage;