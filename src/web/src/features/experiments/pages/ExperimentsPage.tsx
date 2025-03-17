# src/web/src/features/experiments/pages/ExperimentsPage.tsx
```typescript
import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useNavigate, useParams, Routes, Route } from 'react-router-dom'; // react-router-dom ^6.11.0
import { Box, Typography, Paper, Grid, Tabs, Tab } from '@mui/material'; // @mui/material ^5.13.0
import { Add as AddIcon, FilterList as FilterIcon } from '@mui/icons-material'; // @mui/icons-material ^5.13.0

import MainLayout from '../../../layouts/MainLayout'; // Main layout wrapper for authenticated pages
import ExperimentList from '../components/ExperimentList'; // Component for displaying and filtering experiments
import ExperimentDetail from '../components/ExperimentDetail'; // Component for displaying experiment details
import Button from '../../../components/common/Button'; // Reusable button component for actions
import usePermissions from '../../../hooks/usePermissions'; // Hook for checking user permissions
import useExperiments from '../hooks/useExperiments'; // Hook for managing experiment data and operations

/**
 * Main page component for experiments management
 */
const ExperimentsPage: React.FC = () => {
  // LD1: Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Get URL parameters using useParams hook
  const { id } = useParams<{ id: string }>();

  // LD1: Get permission checking functions from usePermissions hook
  const { isPharma, can } = usePermissions();

  // LD1: Set up state for view mode (list or detail)
  const [viewMode, setViewMode] = useState<'list' | 'detail'>('list');

  // LD1: Set up state for selected experiment ID
  const [selectedExperimentId, setSelectedExperimentId] = useState<string | null>(null);

  // LD1: Set up state for filter mode (all or my experiments)
  const [showOnlyMyExperiments, setShowOnlyMyExperiments] = useState(false);

  // LD1: Create handler for navigating to experiment creation page
  const handleCreateExperiment = useCallback(() => {
    navigate('/experiments/create');
  }, [navigate]);

  // LD1: Create handler for selecting an experiment to view details
  const handleExperimentSelect = useCallback((experiment: any) => {
    setSelectedExperimentId(experiment.id);
    setViewMode('detail');
  }, []);

  // LD1: Create handler for returning to experiment list view
  const handleBackToList = useCallback(() => {
    setViewMode('list');
    setSelectedExperimentId(null);
    navigate('/experiments');
  }, [navigate]);

  // LD1: Create handler for refreshing experiment list after status changes
  const handleStatusChange = useCallback(() => {
    // Refresh experiment list by calling fetchExperiments or fetchMyExperiments based on showOnlyMyExperiments state
    // Implementation depends on how you fetch data in your component
  }, []);

  // LD1: Create handler for toggling between all experiments and my experiments
  const handleToggleMyExperiments = useCallback(() => {
    setShowOnlyMyExperiments((prev) => !prev);
  }, []);

  // LD1: Set viewMode to 'detail' and selectedExperimentId to id if id is provided in URL params
  useEffect(() => {
    if (id) {
      setViewMode('detail');
      setSelectedExperimentId(id);
    }
  }, [id]);

  // LD1: Define styled components for layout and visual enhancements
  const PageContainer = styled(Box)({
    padding: theme.spacing(3),
    maxWidth: '1200px',
    margin: '0 auto',
  });

  const PageHeader = styled(Box)({
    display: 'flex',
    justifyContent: 'space-between',
    alignItems: 'center',
    marginBottom: theme.spacing(3),
  });

  const FilterContainer = styled(Box)({
    display: 'flex',
    alignItems: 'center',
    marginBottom: theme.spacing(2),
  });

  // LD1: Render page header with title and action buttons
  return (
    <MainLayout>
      <PageContainer>
        <PageHeader>
          <Typography variant="h5" component="h2">
            Experiments
          </Typography>
          <Box>
            {isPharma() && can('create_experiment') && (
              <Button variant="contained" startIcon={<AddIcon />} onClick={handleCreateExperiment}>
                Create Experiment
              </Button>
            )}
          </Box>
        </PageHeader>

        {/* LD1: Render filter toggle for switching between all and my experiments */}
        <FilterContainer>
          <Button
            variant="outlined"
            startIcon={<FilterIcon />}
            onClick={handleToggleMyExperiments}
          >
            {showOnlyMyExperiments ? 'Show All Experiments' : 'Show My Experiments'}
          </Button>
        </FilterContainer>

        {/* LD1: Conditionally render based on viewMode */}
        {viewMode === 'list' && (
          <ExperimentList
            showFilters={true}
            showActions={true}
            onlyMyExperiments={showOnlyMyExperiments}
            onExperimentSelect={handleExperimentSelect}
          />
        )}
        {viewMode === 'detail' && selectedExperimentId && (
          <ExperimentDetail
            experimentId={selectedExperimentId}
            onBack={handleBackToList}
            onStatusChange={handleStatusChange}
          />
        )}
      </PageContainer>
    </MainLayout>
  );
};

export default ExperimentsPage;