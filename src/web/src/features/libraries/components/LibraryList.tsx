import React, { useState, useCallback } from 'react'; // react v18.2+
import { Box, Typography, Grid, IconButton, Tooltip } from '@mui/material'; // @mui/material v5.13+
import { Add, Delete, Edit, Visibility } from '@mui/icons-material'; // @mui/icons-material v5.13+
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Pagination from '../../../components/common/Pagination';
import { Library } from '../../../types/library';
import useLibraries from '../hooks/useLibraries';

/**
 * Props interface for the LibraryList component
 */
interface LibraryListProps {
  /** Callback function called when a library is selected */
  onLibrarySelect: (library: Library) => void;
  /** Callback function called when the create library button is clicked */
  onCreateLibrary: () => void;
  /** Callback function called when the edit library button is clicked */
  onEditLibrary: (library: Library) => void;
  /** Filter criteria for libraries */
  filter?: any; // TODO: Replace 'any' with the actual LibraryFilter type
  /** Whether to show action buttons (view, edit, delete) */
  showActions?: boolean;
  /** Additional CSS class for the container */
  className?: string;
}

/**
 * Component that displays a list of molecule libraries with pagination and interactive functionality
 */
const LibraryList: React.FC<LibraryListProps> = ({
  onLibrarySelect,
  onCreateLibrary,
  onEditLibrary,
  filter,
  showActions = true,
  className,
}) => {
  // LD1: Use the useLibraries hook to manage library data and operations
  const {
    libraries,
    loading,
    totalLibraries,
    pagination,
    setCurrentLibrary,
    deleteLibrary,
    refreshLibraries,
  } = useLibraries({ filter });

  // LD1: Use useState to manage the ID of the library to confirm deletion
  const [confirmDeleteId, setConfirmDeleteId] = useState<string | null>(null);

  // LD1: useCallback hook to handle library selection
  const handleLibrarySelect = useCallback((library: Library) => {
    setCurrentLibrary(library);
    onLibrarySelect(library);
  }, [setCurrentLibrary, onLibrarySelect]);

  // LD1: useCallback hook to handle library deletion
  const handleDeleteLibrary = useCallback((libraryId: string) => {
    setConfirmDeleteId(libraryId);
  }, []);

  // LD1: useCallback hook to confirm library deletion
  const handleConfirmDelete = useCallback(() => {
    if (confirmDeleteId) {
      deleteLibrary(confirmDeleteId);
      setConfirmDeleteId(null);
    }
  }, [deleteLibrary, confirmDeleteId]);

  // LD1: useCallback hook to cancel library deletion
  const handleCancelDelete = useCallback(() => {
    setConfirmDeleteId(null);
  }, []);

  // LD1: useCallback hook to handle library editing
  const handleEditLibrary = useCallback((library: Library) => {
    onEditLibrary(library);
  }, [onEditLibrary]);

  // LD1: Render a container Box with className prop
  return (
    <Box className={className}>
      {/* LD1: Render header with title and create button if showActions is true */}
      {showActions && (
        <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
          <Typography variant="h6">Molecule Libraries</Typography>
          <Button variant="contained" startIcon={<Add />} onClick={onCreateLibrary}>
            Create Library
          </Button>
        </Box>
      )}

      {/* LD1: Render loading state if loading is true */}
      {loading ? (
        <Typography>Loading libraries...</Typography>
      ) : (
        <>
          {/* LD1: Render empty state if no libraries are available */}
          {libraries.length === 0 ? (
            <Typography>No libraries found.</Typography>
          ) : (
            <>
              {/* LD1: Render Grid container with library cards */}
              <Grid container spacing={2}>
                {libraries.map((library) => (
                  <Grid item xs={12} sm={6} md={4} key={library.id}>
                    {/* LD1: For each library, render a Card with library information */}
                    <Card interactive onClick={() => handleLibrarySelect(library)}>
                      {/* LD1: Include library name, description, molecule count, and creation date */}
                      <Typography variant="h6">{library.name}</Typography>
                      <Typography variant="body2" color="textSecondary">
                        {library.description}
                      </Typography>
                      <Typography variant="caption" color="textSecondary">
                        {library.molecule_count} molecules
                      </Typography>
                      <Typography variant="caption" color="textSecondary">
                        Created: {new Date(library.created_at).toLocaleDateString()}
                      </Typography>

                      {/* LD1: If showActions is true, render action buttons (view, edit, delete) */}
                      {showActions && (
                        <Box mt={2} display="flex" justifyContent="flex-end">
                          <Tooltip title="View Library">
                            <IconButton onClick={(e) => {
                                e.stopPropagation();
                                handleLibrarySelect(library);
                              }} aria-label="view">
                              <Visibility />
                            </IconButton>
                          </Tooltip>
                          <Tooltip title="Edit Library">
                            <IconButton onClick={(e) => {
                                e.stopPropagation();
                                handleEditLibrary(library);
                              }} aria-label="edit">
                              <Edit />
                            </IconButton>
                          </Tooltip>
                          <Tooltip title="Delete Library">
                            <IconButton onClick={(e) => {
                                e.stopPropagation();
                                handleDeleteLibrary(library.id);
                              }} aria-label="delete">
                              <Delete />
                            </IconButton>
                          </Tooltip>
                        </Box>
                      )}
                    </Card>
                  </Grid>
                ))}
              </Grid>

              {/* LD1: Render confirmation dialog when delete is clicked */}
              {confirmDeleteId && (
                <Box>
                  <Typography>Are you sure you want to delete this library?</Typography>
                  <Button onClick={handleConfirmDelete}>Confirm</Button>
                  <Button onClick={handleCancelDelete}>Cancel</Button>
                </Box>
              )}

              {/* LD1: Render Pagination component at the bottom for navigation */}
              <Pagination
                page={pagination.page}
                pageSize={pagination.pageSize}
                totalItems={totalLibraries}
                onPageChange={pagination.goToPage}
                onPageSizeChange={pagination.setPageSize}
              />
            </>
          )}
        </>
      )}
    </Box>
  );
};

export default LibraryList;