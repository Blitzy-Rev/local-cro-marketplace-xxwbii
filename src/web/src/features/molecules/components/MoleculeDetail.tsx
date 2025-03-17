import React, { useMemo } from 'react';
import { 
  Box, 
  Typography, 
  Divider, 
  Chip, 
  Grid, 
  Tooltip, 
  IconButton, 
  Paper, 
  Tab, 
  Tabs,
  Table,
  TableBody,
  TableCell,
  TableContainer,
  TableHead,
  TableRow 
} from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { ArrowBack, Flag, FlagOutlined, Star, StarOutlined, Add, Science, LibraryBooks } from '@mui/icons-material'; // v5.13+

import MoleculeStructure from './MoleculeStructure';
import PropertyTable from './PropertyTable';
import Button from '../../../components/common/Button';
import Card from '../../../components/common/Card';
import Badge from '../../../components/common/Badge';
import { Molecule, FlagStatus } from '../../../types/molecule';
import { Library } from '../../../types/library';
import { Experiment, ExperimentStatus } from '../../../types/experiment';
import theme from '../../../theme';
import { formatDate } from '../../../utils/formatters';
import { getMoleculeNameFromSmiles } from '../../../utils/molecularUtils';

/**
 * Props interface for the MoleculeDetail component
 */
interface MoleculeDetailProps {
  /** Molecule object to display details for */
  molecule: Molecule;
  /** Libraries that contain this molecule */
  libraries: Library[];
  /** Experiments that include this molecule */
  experiments: Experiment[];
  /** Callback when back button is clicked */
  onBack: () => void;
  /** Callback when flag status is changed */
  onFlagChange: (molecule: Molecule, flagStatus: FlagStatus | null) => void;
  /** Callback when add to library button is clicked */
  onAddToLibrary: (molecule: Molecule) => void;
  /** Callback when add to experiment button is clicked */
  onAddToExperiment: (molecule: Molecule) => void;
  /** Additional CSS class for the component */
  className?: string;
  /** Additional inline styles for the component */
  style?: React.CSSProperties;
}

const DetailContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  gap: theme.spacing(3),
  width: '100%',
  padding: theme.spacing(3),
});

const HeaderSection = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  width: '100%',
});

const ActionButtons = styled(Box)({
  display: 'flex',
  gap: theme.spacing(1),
});

const SectionGrid = styled(Grid)({
  container: true,
  spacing: 3,
  marginTop: theme.spacing(2),
});

const ExperimentStatusChip = styled(Chip)({
  size: 'small',
  marginLeft: theme.spacing(1),
});

/**
 * A component that displays detailed information about a specific molecule,
 * including its structure, properties, libraries, and experiments. This component
 * serves as the main content for the molecule detail page and provides actions
 * for flagging molecules and adding them to libraries or experiments.
 */
const MoleculeDetail: React.FC<MoleculeDetailProps> = ({
  molecule,
  libraries = [],
  experiments = [],
  onBack,
  onFlagChange,
  onAddToLibrary,
  onAddToExperiment,
  className,
  style,
}) => {
  // Handle case where molecule data isn't loaded yet
  if (!molecule) {
    return (
      <Box p={3} textAlign="center">
        <Typography variant="body1">Loading molecule details...</Typography>
      </Box>
    );
  }

  // Handler for flag change
  const handleFlagChange = () => {
    const newFlagStatus = molecule.flag_status === FlagStatus.IMPORTANT ? null : FlagStatus.IMPORTANT;
    onFlagChange(molecule, newFlagStatus);
  };

  // Determine molecule name
  const moleculeName = useMemo(() => {
    const nameProperty = molecule.properties?.find(prop => prop.property_name.toLowerCase() === 'name');
    if (nameProperty) return nameProperty.property_value.toString();
    return getMoleculeNameFromSmiles(molecule.smiles);
  }, [molecule]);

  // Get status color based on experiment status
  const getStatusColor = (status: ExperimentStatus): 'default' | 'primary' | 'secondary' | 'success' | 'error' | 'info' | 'warning' => {
    switch (status) {
      case ExperimentStatus.COMPLETED:
        return 'success';
      case ExperimentStatus.IN_PROGRESS:
        return 'info';
      case ExperimentStatus.RESULTS_AVAILABLE:
        return 'primary';
      case ExperimentStatus.REJECTED:
      case ExperimentStatus.QUOTE_REJECTED:
      case ExperimentStatus.RESULTS_REJECTED:
        return 'error';
      case ExperimentStatus.QUOTE_PENDING:
        return 'warning';
      default:
        return 'default';
    }
  };

  return (
    <DetailContainer className={className} style={style}>
      {/* Header with title */}
      <HeaderSection>
        <Box display="flex" alignItems="center">
          <IconButton onClick={onBack} aria-label="back" size="small" sx={{ mr: 1 }}>
            <ArrowBack />
          </IconButton>
          <Typography variant="h5">MOLECULE DETAILS</Typography>
        </Box>
      </HeaderSection>

      {/* Main content grid */}
      <Grid container spacing={3}>
        {/* Structure section */}
        <Grid item xs={12} md={6}>
          <Card>
            <Box p={2}>
              <Typography variant="h6" gutterBottom>
                Structure:
              </Typography>
              <Box display="flex" justifyContent="center" alignItems="center" minHeight={200}>
                <MoleculeStructure 
                  molecule={molecule}
                  width={300}
                  height={200}
                />
              </Box>
            </Box>
          </Card>
        </Grid>

        {/* Properties section */}
        <Grid item xs={12} md={6}>
          <Card>
            <Box p={2}>
              <Typography variant="h6" gutterBottom>
                Properties:
              </Typography>
              <Box mb={2}>
                <Typography variant="body1" gutterBottom>
                  <strong>SMILES:</strong> {molecule.smiles}
                </Typography>
                <Typography variant="body1" gutterBottom>
                  <strong>Name:</strong> {moleculeName}
                </Typography>
                <PropertyTable
                  molecule={molecule}
                  showCalculatedIndicator
                />
              </Box>
            </Box>
          </Card>
        </Grid>
      </Grid>

      {/* Libraries section */}
      <Box mt={3}>
        <Typography variant="h6" gutterBottom>
          Libraries:
        </Typography>
        <Box display="flex" gap={1} flexWrap="wrap" alignItems="center">
          {libraries.length > 0 ? (
            libraries.map((library) => (
              <Chip 
                key={library.id}
                label={library.name}
                color="primary"
                variant="outlined"
              />
            ))
          ) : (
            <Typography variant="body2" color="text.secondary">
              This molecule is not part of any library
            </Typography>
          )}
          <Button
            startIcon={<LibraryBooks />}
            onClick={() => onAddToLibrary(molecule)}
            variant="outlined"
            size="small"
          >
            Add to Library
          </Button>
        </Box>
      </Box>

      {/* Experiments section */}
      <Box mt={3}>
        <Typography variant="h6" gutterBottom>
          Experiments:
        </Typography>
        <TableContainer component={Paper}>
          <Table size="small" aria-label="experiments table">
            <TableHead>
              <TableRow>
                <TableCell>Experiment</TableCell>
                <TableCell>Status</TableCell>
                <TableCell>Date</TableCell>
                <TableCell>Results</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {experiments.length > 0 ? (
                experiments.map((experiment) => (
                  <TableRow key={experiment.id}>
                    <TableCell>{experiment.name}</TableCell>
                    <TableCell>
                      <Badge 
                        label={experiment.status} 
                        color={getStatusColor(experiment.status)}
                        size="small"
                      />
                    </TableCell>
                    <TableCell>{formatDate(experiment.created_at)}</TableCell>
                    <TableCell>
                      {experiment.status === ExperimentStatus.RESULTS_AVAILABLE || 
                       experiment.status === ExperimentStatus.COMPLETED ? (
                        <Button size="small" variant="outlined">View</Button>
                      ) : (
                        'Pending'
                      )}
                    </TableCell>
                  </TableRow>
                ))
              ) : (
                <TableRow>
                  <TableCell colSpan={4} align="center">
                    <Typography variant="body2" color="text.secondary">
                      No experiments for this molecule
                    </Typography>
                  </TableCell>
                </TableRow>
              )}
            </TableBody>
          </Table>
        </TableContainer>
      </Box>

      {/* Action buttons */}
      <Box display="flex" gap={2} mt={3}>
        <Button
          startIcon={molecule.flag_status === FlagStatus.IMPORTANT ? <Flag /> : <FlagOutlined />}
          onClick={handleFlagChange}
          variant="outlined"
          color={molecule.flag_status === FlagStatus.IMPORTANT ? 'primary' : 'default'}
        >
          {molecule.flag_status === FlagStatus.IMPORTANT ? 'Flagged Important' : 'Flag Important'}
        </Button>
        <Button
          startIcon={<Science />}
          onClick={() => onAddToExperiment(molecule)}
          variant="contained"
          color="primary"
        >
          Add to Experiment
        </Button>
        <Button
          startIcon={<ArrowBack />}
          onClick={onBack}
          variant="outlined"
          color="default"
        >
          Back to List
        </Button>
      </Box>
    </DetailContainer>
  );
};

export default MoleculeDetail;