import React from 'react'; // React 18.2+
import { useState, useEffect } from 'react'; // React 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom 6.10+
import { Box, Grid, Typography, Skeleton, Divider } from '@mui/material'; // @mui/material v5.13+
import { LibraryBooks as LibraryIcon, Add as AddIcon } from '@mui/icons-material'; // @mui/icons-material v5.13+
import { useTheme, useMediaQuery } from '@mui/material'; // @mui/material v5.13+

// Internal imports
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import useLibraries from '../../libraries/hooks/useLibraries';
import { Library } from '../../../types/library';
import { formatDate } from '../../../utils/formatters';

/**
 * Interface for the MoleculeLibraries component props
 */
interface MoleculeLibrariesProps {
  /** Callback function to navigate to the full libraries list */
  onViewAll: () => void;
  /** Callback function to create a new library */
  onCreateLibrary: () => void;
  /** Maximum number of library items to display */
  maxItems: number;
  /** Optional CSS class name for styling */
  className?: string;
  /** Optional inline styles for the component */
  style?: React.CSSProperties;
}

/**
 * Dashboard widget component that displays a summary of the user's molecule libraries
 */
const MoleculeLibraries: React.FC<MoleculeLibrariesProps> = ({
  onViewAll,
  onCreateLibrary,
  maxItems,
  className,
  style
}) => {
  // LD1: Initialize navigation hook for routing
  const navigate = useNavigate();

  // LD1: Initialize theme and media query for responsive design
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // LD1: Initialize libraries hook with limited number of items
  const { libraries, loading, totalLibraries } = useLibraries({
    initialPageSize: maxItems,
    fetchMyLibrariesOnly: true
  });

  // LD1: Handle navigation to a specific library detail page
  const handleNavigateToLibrary = (libraryId: string) => {
    navigate(`/libraries/${libraryId}`);
  };

  // LD1: Render a card for a single library with its details
  const renderLibraryCard = (library: Library) => (
    <Grid item xs={12} sm={6} md={4} key={library.id}>
      <Card
        interactive
        onClick={() => handleNavigateToLibrary(library.id)}
        style={{ height: '100%' }}
      >
        <Typography variant="h6" gutterBottom>
          {library.name}
        </Typography>
        <Typography variant="body2" color="textSecondary" paragraph>
          {library.description}
        </Typography>
        <Box display="flex" alignItems="center">
          <LibraryIcon color="primary" sx={{ mr: 1 }} />
          <Typography variant="body2">
            {library.molecule_count} Molecules
          </Typography>
        </Box>
      </Card>
    </Grid>
  );

  // LD1: Render a skeleton placeholder card for loading state
  const renderSkeletonCard = () => (
    <Grid item xs={12} sm={6} md={4} key={`skeleton-${Math.random()}`}>
      <Card style={{ height: '100%' }}>
        <Skeleton variant="text" height={24} width="70%" />
        <Skeleton variant="text" height={20} width="90%" />
        <Skeleton variant="circular" width={24} height={24} />
      </Card>
    </Grid>
  );

  // LD1: Apply responsive layout based on screen size
  const containerProps = isMobile ? { flexDirection: 'column' } : {};

  return (
    <Card className={className} style={style}>
      <Box display="flex" justifyContent="space-between" alignItems="center" mb={2}>
        <Typography variant="h5" component="h3">
          Molecule Libraries
        </Typography>
        <Typography variant="body2" color="textSecondary">
          {totalLibraries} Total
        </Typography>
      </Box>
      <Divider sx={{ mb: 2 }} />
      {loading ? (
        <Grid container spacing={2} {...containerProps}>
          {Array.from({ length: maxItems }).map((_, index) => renderSkeletonCard())}
        </Grid>
      ) : (
        <Grid container spacing={2} {...containerProps}>
          {libraries.map(library => renderLibraryCard(library))}
        </Grid>
      )}
      <Box mt={2} display="flex" justifyContent="space-between" alignItems="center">
        <Button size="small" onClick={onViewAll}>
          View All
        </Button>
        <Button size="small" variant="outlined" startIcon={<AddIcon />} onClick={onCreateLibrary}>
          Create Library
        </Button>
      </Box>
    </Card>
  );
};

export default MoleculeLibraries;