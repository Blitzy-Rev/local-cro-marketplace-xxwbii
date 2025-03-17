import React, { CSSProperties } from 'react'; // react 18.2+
import { Box, Paper, Typography, styled } from '@mui/material'; // @mui/material 5.13+
import useDragDrop from '../hooks/useDragDrop';
import { Library } from '../../../types/library';
import { Molecule } from '../../../types/molecule';
import { useToast } from '../../../hooks/useToast';

/**
 * Props interface for the DragDropLibrary component
 */
interface DragDropLibraryProps {
  /** The library object that serves as the drop target */
  library: Library;
  /** Callback function triggered when a molecule is dropped on the library */
  onDrop: (libraryId: string, libraryName: string) => void;
  /** Child components to render inside the drop zone */
  children?: React.ReactNode;
  /** Optional CSS class name for styling */
  className?: string;
  /** Optional inline styles */
  style?: CSSProperties;
  /** Message to display when the library is empty */
  emptyStateMessage?: string;
  /** Whether to show a visual indicator for the drop zone */
  showDropIndicator?: boolean;
}

/**
 * Styled Paper component to act as the drop container
 */
const DropContainer = styled(Paper, {
  shouldForwardProp: (prop) => prop !== 'isOver',
})<{ isOver: boolean }>(({ theme, isOver }) => ({
  padding: theme.spacing(2),
  minHeight: '150px',
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  border: isOver ? `2px dashed ${theme.palette.primary.main}` : `2px dashed ${theme.palette.divider}`,
  backgroundColor: isOver ? 'rgba(25, 118, 210, 0.08)' : 'transparent',
  transition: 'all 0.2s ease-in-out',
  borderRadius: theme.shape.borderRadius,
  position: 'relative',
  overflow: 'hidden'
}));

/**
 * Styled Box component to act as the drop indicator
 */
const DropIndicator = styled(Box, {
  shouldForwardProp: (prop) => prop !== 'isOver',
})<{ isOver: boolean }>(({ isOver }) => ({
  position: 'absolute',
  top: 0,
  left: 0,
  right: 0,
  bottom: 0,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  backgroundColor: 'rgba(255, 255, 255, 0.8)',
  zIndex: 1,
  opacity: isOver ? 1 : 0,
  pointerEvents: 'none',
  transition: 'opacity 0.2s ease-in-out'
}));

/**
 * Styled Typography component for the empty state message
 */
const EmptyState = styled(Typography)(({ theme }) => ({
  color: theme.palette.text.secondary,
  textAlign: 'center',
  fontStyle: 'italic'
}));

/**
 * Styled Box component for the content container
 */
const ContentContainer = styled(Box)({
  width: '100%',
  height: '100%',
  position: 'relative',
  zIndex: 0
});

/**
 * A component that provides a drop target for dragging molecules into a library
 */
const DragDropLibrary: React.FC<DragDropLibraryProps> = ({
  library,
  onDrop,
  children,
  className = '',
  style = {},
  emptyStateMessage = 'Drag molecules here to add to library',
  showDropIndicator = true,
}) => {
  // Initialize the useDragDrop hook
  const { libraryDropRef, isOver } = useDragDrop({
    onDrop: (molecule: Molecule, libraryId: string, sourceLibraryId: string | null) => {
      // Call the onDrop callback with the library ID and name
      onDrop(library.id, library.name);
    }
  });

  return (
    <DropContainer
      className={className}
      style={style}
      isOver={isOver}
      ref={libraryDropRef as React.Ref<Paper>}
    >
      <ContentContainer>
        {children ? (
          children
        ) : (
          <EmptyState>{emptyStateMessage}</EmptyState>
        )}
      </ContentContainer>
      {showDropIndicator && (
        <DropIndicator isOver={isOver}>
          <Typography variant="h6">Drop here to add</Typography>
        </DropIndicator>
      )}
    </DropContainer>
  );
};

export default DragDropLibrary;