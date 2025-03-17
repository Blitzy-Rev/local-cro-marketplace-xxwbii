import React, { useCallback, useEffect, useMemo, useState } from 'react'; // React v18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.4+
import {
  Box,
  Typography,
  Paper,
  Grid,
  Divider,
  Menu,
  MenuItem,
  IconButton,
  Tooltip,
} from '@mui/material'; // @mui/material v5.13+
import {
  Add,
  FilterList,
  Refresh,
  Upload,
  LibraryAdd,
  Science,
  Flag,
  MoreVert,
} from '@mui/icons-material'; // @mui/icons-material v5.13+

import MainLayout from '../../../layouts/MainLayout';
import MoleculeTable from '../components/MoleculeTable';
import MoleculeFilter from '../components/MoleculeFilter';
import Button from '../../../components/common/Button';
import useMolecules from '../hooks/useMolecules';
import useToast from '../../../hooks/useToast';
import useDialog from '../../../hooks/useDialog';
import {
  Molecule,
  MoleculeFilter as MoleculeFilterType,
} from '../../../types/molecule';
import { DEFAULT_MOLECULE_PROPERTIES } from '../../../utils/molecularUtils';

/**
 * @file MoleculesListPage.tsx
 * @src_subfolder web
 * @description A page component that displays a list of molecules with filtering, sorting, and selection capabilities. It serves as the main interface for users to browse, search, and interact with molecular data in the platform.
 * @requirements_addressed
 * - {@link Technical Specifications/2.2.3 Molecule Sorting & Organization (F-003) | Molecule Sorting & Organization}: Implements interactive interface for sorting, filtering, and organizing molecules by properties
 * - {@link Technical Specifications/2.2.4 Molecule Management & Experiment Queuing (F-004) | Molecule Management & Experiment Queuing}: Provides functionality for flagging molecules and adding them to experimental queues
 * - {@link Technical Specifications/7.4.4 Molecule List View | User Interface Design}: Implements the molecule list view as specified in the wireframes
 * - {@link Technical Specifications/7.8 Interaction Patterns | Real-time Filtering}: Implements debounced filtering to update results as user types
 */

/**
 * @function MoleculesListPage
 * @description Main page component for displaying and interacting with the list of molecules
 * @parameters None
 * @returns {JSX.Element} The rendered page component
 */
const MoleculesListPage: React.FC = () => {
  // IE1: Initialize navigation with useNavigate hook
  const navigate = useNavigate();

  // IE1: Initialize toast notifications with useToast hook
  const { showToast } = useToast();

  // IE1: Initialize dialog management with useDialog hook
  const { handleOpen: openLibraryDialog } = useDialog();
  const { handleOpen: openExperimentDialog } = useDialog();

  // IE1: Set up molecule data management with useMolecules hook
  const {
    molecules,
    loading,
    error,
    filter,
    selectedMolecules,
    totalMolecules,
    pagination,
    setFilter,
    refreshMolecules,
    selectMolecule,
    deselectMolecule,
    selectAllMolecules,
    deselectAllMolecules,
    updateMoleculeFlag,
    addToLibrary,
    addToExperiment,
  } = useMolecules();

  // LD1: Set up state for action menu anchor element
  const [anchorEl, setAnchorEl] = useState<null | HTMLElement>(null);

  // LD1: Define handler for viewing molecule details
  const handleViewDetails = useCallback(
    (molecule: Molecule) => {
      // IE1: Navigate to the molecule detail page using the molecule ID
      navigate(`/app/molecules/${molecule.id}`);
    },
    [navigate]
  );

  // LD1: Define handler for adding molecule to queue
  const handleAddToQueue = useCallback(
    (molecule: Molecule) => {
      // IE1: Navigate to the experiment creation page with molecule ID as parameter
      navigate(`/app/experiments/create?moleculeId=${molecule.id}`);
    },
    [navigate]
  );

  // LD1: Define handler for adding selected molecules to library
  const handleAddToLibrary = useCallback(() => {
    // LD1: Check if any molecules are selected
    if (selectedMolecules.length === 0) {
      // IE1: If no molecules selected, show warning toast
      showToast({
        type: 'warning',
        message: 'Please select molecules to add to library',
      });
      return;
    }

    // IE1: If molecules selected, open library selection dialog
    openLibraryDialog();

    // IE1: When library selected, call addToLibrary function from useMolecules hook
    // IE1: Show success or error toast based on result
  }, [selectedMolecules, showToast, openLibraryDialog]);

  // LD1: Define handler for adding selected molecules to experiment
  const handleAddToExperiment = useCallback(() => {
    // LD1: Check if any molecules are selected
    if (selectedMolecules.length === 0) {
      // IE1: If no molecules selected, show warning toast
      showToast({
        type: 'warning',
        message: 'Please select molecules to add to experiment',
      });
      return;
    }

    // IE1: If molecules selected, open experiment selection dialog
    openExperimentDialog();

    // IE1: When experiment selected, call addToExperiment function from useMolecules hook
    // IE1: Show success or error toast based on result
  }, [selectedMolecules, showToast, openExperimentDialog]);

  // LD1: Define handler for CSV upload navigation
  const handleCSVUpload = useCallback(() => {
    // IE1: Navigate to the CSV upload page
    navigate('/app/molecules/upload');
  }, [navigate]);

  // LD1: Define handler for opening and closing action menu
  const handleMenuOpen = useCallback((event: React.MouseEvent<HTMLElement>) => {
    // LD1: Set the anchor element to the current target of the event
    setAnchorEl(event.currentTarget);
  }, []);

  const handleMenuClose = useCallback(() => {
    // LD1: Set the anchor element to null
    setAnchorEl(null);
  }, []);

  // LD1: Define handler for bulk actions on selected molecules
  const handleRefresh = useCallback(() => {
    // IE1: Call refreshMolecules function from useMolecules hook
    refreshMolecules();
    // IE1: Show toast notification that data is refreshing
    showToast({
      type: 'info',
      message: 'Refreshing molecule data...',
    });
  }, [refreshMolecules, showToast]);

  const handleBulkActions = () => {
    // TODO: Implement bulk actions
  };

  return (
    <MainLayout>
      {/* LD1: Render page header with title and action buttons */}
      <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
        <Typography variant="h4" component="h1">
          Molecules
        </Typography>
        <Box>
          <Button variant="contained" startIcon={<Add />} onClick={() => navigate('/app/molecules/create')}>
            Add Molecule
          </Button>
          <Button variant="outlined" startIcon={<Upload />} onClick={handleCSVUpload} sx={{ ml: 2 }}>
            Upload CSV
          </Button>
          {selectedMolecules.length > 0 && (
            <IconButton onClick={handleMenuOpen} aria-label="Actions" sx={{ ml: 2 }}>
              <MoreVert />
            </IconButton>
          )}
          <IconButton onClick={handleRefresh} aria-label="Refresh" sx={{ ml: 2 }}>
            <Refresh />
          </IconButton>
        </Box>
      </Box>

      {/* LD1: Render MoleculeFilter component for filtering molecules */}
      <MoleculeFilter initialFilter={filter} onFilterChange={setFilter} />

      {/* LD1: Render MoleculeTable component for displaying molecules */}
      <MoleculeTable
        molecules={molecules}
        loading={loading}
        error={error}
        selectedMolecules={selectedMolecules}
        onSelectionChange={selectAllMolecules}
        onViewDetails={handleViewDetails}
        onAddToQueue={handleAddToQueue}
        onFlagChange={updateMoleculeFlag}
        pagination={{
          ...pagination,
          onPageChange: pagination.goToPage,
          onPageSizeChange: pagination.setPageSize,
        }}
        displayProperties={DEFAULT_MOLECULE_PROPERTIES}
        selectable={true}
        showActions={true}
      />

      {/* LD1: Render bulk action controls when molecules are selected */}
      <Menu
        anchorEl={anchorEl}
        open={Boolean(anchorEl)}
        onClose={handleMenuClose}
      >
        <MenuItem onClick={handleAddToLibrary}>
          <ListItemIcon>
            <LibraryAdd />
          </ListItemIcon>
          Add to Library
        </MenuItem>
        <MenuItem onClick={handleAddToExperiment}>
          <ListItemIcon>
            <Science />
          </ListItemIcon>
          Add to Experiment
        </MenuItem>
        <MenuItem onClick={handleBulkActions}>
          <ListItemIcon>
            <Flag />
          </ListItemIcon>
          Flag as Important
        </MenuItem>
      </Menu>

      {/* LD1: Handle loading and error states appropriately */}
      {loading && (
        <Box display="flex" justifyContent="center" alignItems="center">
          <Typography variant="body1">Loading molecules...</Typography>
        </Box>
      )}
      {error && (
        <Box display="flex" justifyContent="center" alignItems="center">
          <Typography variant="body1" color="error">
            Error: {error}
          </Typography>
        </Box>
      )}
    </MainLayout>
  );
};

export default MoleculesListPage;