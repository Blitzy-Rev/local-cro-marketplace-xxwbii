import React, { useState, useEffect, useCallback, useMemo } from 'react'; // React 18.2+
import { DndProvider } from 'react-dnd'; // react-dnd v16.0.1
import { HTML5Backend } from 'react-dnd-html5-backend'; // react-dnd-html5-backend v16.0.1
import { Box, Typography, Grid, Divider, Paper, CircularProgress, Tooltip, IconButton } from '@mui/material'; // @mui/material v5.13+
import { Edit, Delete, Add, Download, Share } from '@mui/icons-material'; // @mui/icons-material v5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13+

import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import MoleculeDragItem from './MoleculeDragItem';
import { LibraryDetail as LibraryDetailType, LibraryMoleculeOperation, LibraryOperationType } from '../../../types/library';
import { Molecule } from '../../../types/molecule';
import useLibraries from '../hooks/useLibraries';
import { formatDate } from '../../../utils/formatters';
import { useToast } from '../../../hooks/useToast';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Interface for the props of the LibraryDetail component
 */
interface LibraryDetailProps {
  libraryId: string;
  onEdit: () => void;
  onBack: () => void;
  onMoleculeSelect: (molecule: Molecule) => void;
  readOnly: boolean;
  className: string;
}

// Styled components for layout and appearance
const LibraryContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: '100%',
  width: '100%',
  overflow: 'hidden'
});

const LibraryHeader = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  marginBottom: '16px'
});

const LibraryInfo = styled(Card)({
  marginBottom: '16px'
});

const LibraryContent = styled(Box)({
  flex: '1',
  overflow: 'auto'
});

const MoleculesGrid = styled(Grid)({
  container: true,
  spacing: 2
});

const MoleculeItem = styled(Grid)({
  item: true,
  xs: 12,
  sm: 6,
  md: 4,
  lg: 3
});

const DropZone = styled(Paper)(({ isOver }: { isOver: boolean }) => ({
  padding: '16px',
  textAlign: 'center',
  backgroundColor: isOver ? 'rgba(0, 0, 0, 0.05)' : 'transparent',
  border: isOver ? '2px dashed #1976d2' : '2px dashed #ccc',
  borderRadius: '4px',
  marginTop: '16px',
  transition: 'all 0.2s'
}));

const EmptyState = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '32px',
  textAlign: 'center'
});

const LoadingContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'center',
  alignItems: 'center',
  height: '200px'
});

/**
 * A component that displays detailed information about a molecule library and its contents
 */
const LibraryDetail: React.FC<LibraryDetailProps> = ({
  libraryId,
  onEdit,
  onBack,
  onMoleculeSelect,
  readOnly = false,
  className
}) => {
  // IE1: Use the useLibraries hook to manage library data and operations
  const {
    currentLibrary,
    loading,
    error,
    dragDrop,
    fetchLibrary,
  } = useLibraries();

  // IE1: Fetch the library details when the component mounts or the libraryId changes
  useEffect(() => {
    fetchLibrary(libraryId);
  }, [libraryId, fetchLibrary]);

  // IE1: Handle molecule drag start
  const handleMoleculeDragStart = (molecule: Molecule, sourceLibraryId: string | null) => {
    dragDrop.handleDragStart(molecule, sourceLibraryId);
  };

  // IE1: Handle molecule drag end
  const handleMoleculeDragEnd = () => {
    dragDrop.handleDragEnd();
  };

  // IE1: Handle molecule drop
  const handleMoleculeDrop = () => {
    dragDrop.handleDrop(libraryId, currentLibrary?.name || '');
  };

  // IE1: Handle molecule remove
  const handleMoleculeRemove = (molecule: Molecule) => {
    // IE1: Call removeMolecules function from useLibraries hook with the library ID and molecule ID
    // IE1: Show a success toast notification
  };

  // IE1: Handle molecule click
  const handleMoleculeClick = (molecule: Molecule) => {
    // IE1: If onMoleculeSelect prop is provided, call it with the selected molecule
    if (onMoleculeSelect) {
      onMoleculeSelect(molecule);
    }
  };

  // IE1: Handle export
  const handleExport = (format: string) => {
    // IE1: Call the API to export the library in the specified format
    // IE1: Download the exported file
    // IE1: Show a success toast notification
  };

  // LD1: Render loading state
  if (loading) {
    return (
      <LibraryContainer className={className}>
        <LoadingContainer>
          <CircularProgress />
        </LoadingContainer>
      </LibraryContainer>
    );
  }

  // LD1: Render error state
  if (error) {
    return (
      <LibraryContainer className={className}>
        <Typography color="error">{error}</Typography>
      </LibraryContainer>
    );
  }

  // LD1: Render library details
  return (
    <LibraryContainer className={className}>
      <LibraryHeader>
        <Typography variant="h6">{currentLibrary?.name}</Typography>
        <div>
          <Button variant="outlined" onClick={onBack}>Back</Button>
          {!readOnly && (
            <Button startIcon={<Edit />} onClick={onEdit}>
              Edit
            </Button>
          )}
        </div>
      </LibraryHeader>
      <LibraryInfo>
        <Typography variant="body2">Description: {currentLibrary?.description}</Typography>
        <Typography variant="body2">Created: {formatDate(currentLibrary?.created_at || '')}</Typography>
        <Typography variant="body2">Number of Molecules: {currentLibrary?.molecule_count}</Typography>
      </LibraryInfo>
      <LibraryContent ref={dragDrop.libraryDropRef}>
        {currentLibrary?.molecules && currentLibrary.molecules.length > 0 ? (
          <MoleculesGrid container spacing={2}>
            {currentLibrary.molecules.map((molecule) => (
              <MoleculeItem item xs={12} sm={6} md={4} lg={3} key={molecule.id}>
                <MoleculeDragItem
                  molecule={molecule}
                  sourceLibraryId={libraryId}
                  onDragStart={handleMoleculeDragStart}
                  onDragEnd={handleMoleculeDragEnd}
                  onClick={handleMoleculeClick}
                />
              </MoleculeItem>
            ))}
          </MoleculesGrid>
        ) : (
          <EmptyState>
            <Typography variant="body1">No molecules in this library.</Typography>
            <Typography variant="body2">Drag molecules here to add them.</Typography>
          </EmptyState>
        )}
        <DropZone isOver={dragDrop.isOver}>
          Drop molecules here to add them to the library
        </DropZone>
      </LibraryContent>
    </LibraryContainer>
  );
};

export default LibraryDetail;