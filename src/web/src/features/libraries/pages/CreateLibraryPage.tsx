# src/web/src/features/libraries/pages/CreateLibraryPage.tsx
```typescript
import React, { useCallback } from 'react'; // React v18.2+
import { useNavigate, useLocation } from 'react-router-dom'; // react-router-dom v6.4+
import { Box, Typography, Container } from '@mui/material'; // @mui/material v5.13+

import MainLayout from '../../../layouts/MainLayout';
import CreateLibraryForm from '../components/CreateLibraryForm';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import useToast from '../../../hooks/useToast';
import useMolecules from '../../molecules/hooks/useMolecules';

/**
 * Page component for creating a new molecule library
 */
const CreateLibraryPage: React.FC = () => {
  // HOOK: Initialize useNavigate for navigation
  const navigate = useNavigate();

  // HOOK: Initialize useLocation to get current location
  const location = useLocation();

  // HOOK: Initialize useToast for displaying toast notifications
  const { showToast } = useToast();

  // HOOK: Initialize useMolecules to access selected molecules
  const { selectedMolecules, deselectAllMolecules } = useMolecules();

  /**
   * Handles successful library creation
   * @param libraryId The ID of the created library
   */
  const handleLibraryCreated = useCallback((libraryId: string) => {
    // STEP: Clear any selected molecules using deselectAllMolecules
    deselectAllMolecules();

    // STEP: Show success toast notification
    showToast({
      type: 'success',
      message: 'Library created successfully!',
    });

    // STEP: Navigate to the library detail page using the created library ID
    navigate(`/app/libraries/${libraryId}`);
  }, [deselectAllMolecules, navigate, showToast]);

  /**
   * Handles cancellation of library creation
   */
  const handleCancel = useCallback(() => {
    // STEP: Navigate to the libraries list page
    navigate('/app/libraries');
  }, [navigate]);

  return (
    <MainLayout>
      <Container maxWidth="md">
        <Box sx={{ mt: 3, mb: 2 }}>
          <Typography variant="h4" component="h1" gutterBottom>
            Create New Library
          </Typography>
        </Box>
        <Card>
          <CreateLibraryForm
            onSuccess={handleLibraryCreated}
            onCancel={handleCancel}
            initialMolecules={selectedMolecules}
          />
        </Card>
      </Container>
    </MainLayout>
  );
};

export default CreateLibraryPage;