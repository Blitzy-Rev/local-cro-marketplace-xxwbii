# src/web/src/features/results/pages/ResultsPage.tsx
```typescript
import React, { useCallback } from 'react'; // React v18.2+
import { useNavigate, useParams } from 'react-router-dom'; // react-router-dom v6.11.2
import { Box, Typography, Container, Paper } from '@mui/material'; // @mui/material v5.13.0

import MainLayout from '../../../layouts/MainLayout';
import ResultsList from '../components/ResultsList';
import ResultDetail from '../components/ResultDetail';
import useResults from '../hooks/useResults';
import { ResultFilter } from '../../../types/result';

/**
 * Main page component for viewing and managing experimental results
 * @returns The rendered component
 */
const ResultsPage: React.FC = () => {
  // Get URL parameters using useParams hook
  const { id } = useParams<{ id: string }>();

  // Initialize navigation using useNavigate hook
  const navigate = useNavigate();

  // Initialize useResults hook for results data and operations
  const { clearError } = useResults();

  /**
   * Handles selecting a result to view details
   * @param resultId 
   */
  const handleResultSelect = useCallback((resultId: string) => {
    clearError();
    navigate(`/app/results/${resultId}`);
  }, [navigate, clearError]);

  /**
   * Handles navigating back to the results list
   */
  const handleBackToList = useCallback(() => {
    clearError();
    navigate('/app/results');
  }, [navigate, clearError]);

  return (
    <MainLayout>
      <Container maxWidth="lg">
        <Paper elevation={3} sx={{ p: 3, display: 'flex', flexDirection: 'column', gap: 2 }}>
          <Typography variant="h4" component="h1" gutterBottom>
            Experimental Results
          </Typography>

          {/* Conditional rendering based on URL parameter */}
          {id ? (
            <ResultDetail resultId={id} onBack={handleBackToList} />
          ) : (
            <ResultsList onResultSelect={handleResultSelect} />
          )}
        </Paper>
      </Container>
    </MainLayout>
  );
};

export default ResultsPage;