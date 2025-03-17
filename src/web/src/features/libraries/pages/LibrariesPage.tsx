import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useNavigate, useParams } from 'react-router-dom'; // react-router-dom v6.4+
import { Box, Grid, Typography, CircularProgress, Divider } from '@mui/material'; // @mui/material v5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13+

import MainLayout from '../../../layouts/MainLayout';
import LibraryList from '../components/LibraryList';
import LibraryDetail from '../components/LibraryDetail';
import CreateLibraryForm from '../components/CreateLibraryForm';
import Dialog from '../../../components/common/Dialog';
import useLibraries from '../hooks/useLibraries';
import { Library } from '../../../types/library';
import { Molecule } from '../../../types/molecule';
import useToast from '../../../hooks/useToast';

/**
 * @file A page component that provides a comprehensive interface for managing molecule libraries in the Molecular Data Management and CRO Integration Platform. It allows users to view, create, edit, and interact with libraries, and displays detailed information about selected libraries with drag-and-drop functionality for organizing molecules.
 */

/**
 * @component PageContainer
 * @description Styled Box component for the main page container
 */
const PageContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: '100%',
  width: '100%',
});

/**
 * @component ContentContainer
 * @description Styled Box component for the main content area
 */
const ContentContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  flex: '1',
  overflow: 'hidden'
});

/**
 * @component LibrariesPage
 * @description Main page component for managing molecule libraries
 */
const LibrariesPage: React.FC = () => {
  // LD1: Initialize state variables for managing dialogs and selected library
  const [isCreateDialogOpen, setIsCreateDialogOpen] = useState<boolean>(false);
  const [isEditDialogOpen, setIsEditDialogOpen] = useState<boolean>(false);
  const [libraryToEdit, setLibraryToEdit] = useState<Library | null>(null);

  // LD1: Get libraryId from URL parameters using useParams
  const { libraryId } = useParams<{ libraryId: string }>();

  // LD1: Use useLibraries hook to get library data and functions
  const { 
    currentLibrary,
    loading,
    error,
    fetchLibrary,
    clearCurrentLibrary
  } = useLibraries();

  // LD1: Use useToast hook for notifications
  const { showToast } = useToast();

  // LD1: Use useNavigate hook for navigation
  const navigate = useNavigate();

  // LD1: useEffect hook to fetch library details when libraryId changes
  useEffect(() => {
    if (libraryId) {
      fetchLibrary(libraryId);
    }
    return () => {
      clearCurrentLibrary();
    };
  }, [libraryId, fetchLibrary, clearCurrentLibrary]);

  /**
   * @function handleLibrarySelect
   * @description Handles the selection of a library from the list
   * @param {Library} library - The selected library
   * @returns {void} - No return value
   */
  const handleLibrarySelect = useCallback((library: Library) => {
    // LD1: Navigate to the library detail page using the library ID
    navigate(`/app/libraries/${library.id}`);
  }, [navigate]);

  /**
   * @function handleCreateLibrary
   * @description Opens the create library dialog
   * @returns {void} - No return value
   */
  const handleCreateLibrary = useCallback(() => {
    // LD1: Set isCreateDialogOpen state to true
    setIsCreateDialogOpen(true);
  }, []);

  /**
   * @function handleCloseCreateDialog
   * @description Closes the create library dialog
   * @returns {void} - No return value
   */
  const handleCloseCreateDialog = useCallback(() => {
    // LD1: Set isCreateDialogOpen state to false
    setIsCreateDialogOpen(false);
  }, []);

  /**
   * @function handleLibraryCreated
   * @description Handles the successful creation of a library
   * @param {string} libraryId - The ID of the newly created library
   * @returns {void} - No return value
   */
  const handleLibraryCreated = useCallback((libraryId: string) => {
    // LD1: Close the create library dialog
    handleCloseCreateDialog();
    // LD1: Navigate to the library detail page using the new library ID
    navigate(`/app/libraries/${libraryId}`);
    // LD1: Show a success toast notification
    showToast({
      type: 'success',
      message: 'Library created successfully!'
    });
  }, [navigate, handleCloseCreateDialog, showToast]);

  /**
   * @function handleEditLibrary
   * @description Handles editing a library
   * @param {Library} library - The library to edit
   * @returns {void} - No return value
   */
  const handleEditLibrary = useCallback((library: Library) => {
    // LD1: Set the library to edit in state
    setLibraryToEdit(library);
    // LD1: Open the edit library dialog
    setIsEditDialogOpen(true);
  }, []);

  /**
   * @function handleBackToList
   * @description Handles navigation back to the library list
   * @returns {void} - No return value
   */
  const handleBackToList = useCallback(() => {
    // LD1: Navigate to the main libraries page
    navigate('/app/libraries');
    // LD1: Clear the selected library
    clearCurrentLibrary();
  }, [navigate, clearCurrentLibrary]);

  /**
   * @function handleMoleculeSelect
   * @description Handles the selection of a molecule from a library
   * @param {Molecule} molecule - The selected molecule
   * @returns {void} - No return value
   */
  const handleMoleculeSelect = useCallback((molecule: Molecule) => {
    // LD1: Navigate to the molecule detail page using the molecule ID
    navigate(`/app/molecules/${molecule.id}`);
  }, [navigate]);

  // LD1: Render the component
  return (
    <MainLayout>
      <PageContainer>
        <ContentContainer>
          {loading && <CircularProgress />}
          {error && <Typography color="error">{error}</Typography>}
          {libraryId && currentLibrary ? (
            <LibraryDetail
              libraryId={libraryId}
              onEdit={() => handleEditLibrary(currentLibrary)}
              onBack={handleBackToList}
              onMoleculeSelect={handleMoleculeSelect}
              className="library-detail"
            />
          ) : (
            <LibraryList
              onLibrarySelect={handleLibrarySelect}
              onCreateLibrary={handleCreateLibrary}
              onEditLibrary={handleEditLibrary}
              className="library-list"
            />
          )}
        </ContentContainer>
      </PageContainer>

      {/* LD1: Render CreateLibraryForm in a Dialog when isCreateDialogOpen is true */}
      <Dialog
        open={isCreateDialogOpen}
        onClose={handleCloseCreateDialog}
        title="Create New Library"
        onConfirm={handleCloseCreateDialog}
        onCancel={handleCloseCreateDialog}
      >
        <CreateLibraryForm onSuccess={handleLibraryCreated} onCancel={handleCloseCreateDialog} />
      </Dialog>

      {/* LD1: Render EditLibraryForm in a Dialog when isEditDialogOpen is true */}
      {isEditDialogOpen && libraryToEdit && (
        <Dialog
          open={isEditDialogOpen}
          onClose={() => setIsEditDialogOpen(false)}
          title="Edit Library"
          onConfirm={() => setIsEditDialogOpen(false)}
          onCancel={() => setIsEditDialogOpen(false)}
        >
          {/* TODO: Implement EditLibraryForm */}
          <Typography>Edit Library Form</Typography>
        </Dialog>
      )}
    </MainLayout>
  );
};

export default LibrariesPage;