import React, { useMemo, useCallback } from 'react'; // React v18.2+
import { Box, Typography, IconButton, Tooltip } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { Visibility, Add, Flag, Star, StarBorder, Queue } from '@mui/icons-material'; // v5.13+
import Table from '../../../components/common/Table';
import MoleculeStructure from './MoleculeStructure';
import PropertyBadge from '../../../components/molecular/PropertyBadge';
import Button from '../../../components/common/Button';
import useMolecules from '../hooks/useMolecules';
import { Molecule, FlagStatus } from '../../../types/molecule';
import { getMoleculePropertyByName, getPropertyDisplayName } from '../../../utils/molecularUtils';
import theme from '../../../theme';

/**
 * Props interface for the MoleculeTable component
 */
interface MoleculeTableProps {
  /** Array of molecule objects to display */
  molecules: Molecule[];
  /** Whether the table is in a loading state */
  loading: boolean;
  /** Error message to display if data loading failed */
  error: string | null;
  /** Array of selected molecule IDs */
  selectedMolecules: string[];
  /** Callback when selection changes */
  onSelectionChange: (selectedIds: string[]) => void;
  /** Callback when view details button is clicked */
  onViewDetails: (molecule: Molecule) => void;
  /** Callback when add to queue button is clicked */
  onAddToQueue: (molecule: Molecule) => void;
  /** Callback when flag status is changed */
  onFlagChange: (molecule: Molecule, flagStatus: FlagStatus | null) => void;
  /** Pagination configuration */
  pagination: {
    page: number;
    pageSize: number;
    totalItems: number;
    onPageChange: (page: number) => void;
    onPageSizeChange: (pageSize: number) => void;
  };
  /** Array of property names to display as columns */
  displayProperties: string[];
  /** Whether rows can be selected with checkboxes */
  selectable: boolean;
  /** Whether to show action buttons */
  showActions: boolean;
  /** Additional CSS class for the table */
  className?: string;
  /** Additional inline styles for the table */
  style?: React.CSSProperties;
}

/**
 * A specialized table component for displaying molecular data with advanced features
 */
const MoleculeTable: React.FC<MoleculeTableProps> = ({
  molecules,
  loading,
  error,
  selectedMolecules,
  onSelectionChange,
  onViewDetails,
  onAddToQueue,
  onFlagChange,
  pagination,
  displayProperties,
  selectable,
  showActions,
  className,
  style,
}) => {
  /**
   * Creates a render function for the structure cell
   * @returns A React component that renders the molecule structure
   */
  const renderStructureCell = useCallback((row: Molecule) => (
    <MoleculeStructure molecule={row} width={80} height={60} interactive={false} />
  ), []);

  /**
   * Creates a render function for the SMILES cell
   * @returns A React component that renders the SMILES string
   */
  const renderSmilesCell = useCallback((row: Molecule) => (
    <Typography variant="body2" style={{ maxWidth: '200px', overflow: 'hidden', textOverflow: 'ellipsis', whiteSpace: 'nowrap' }}>
      {row.smiles}
    </Typography>
  ), []);

  /**
   * Creates a render function for the property cell
   * @param propertyName The name of the property to display
   * @returns A React component that renders the property value with a badge
   */
  const renderPropertyCell = useCallback((propertyName: string) => (row: Molecule) => {
    const propertyValue = getMoleculePropertyByName(row, propertyName);
    return propertyValue !== null && propertyValue !== undefined ? (
      <PropertyBadge propertyName={propertyName} value={propertyValue} showTooltip />
    ) : (
      <Typography variant="body2">N/A</Typography>
    );
  }, []);

  /**
   * Creates a render function for the actions cell
   * @param molecule The molecule object
   * @returns A React component that renders the action buttons
   */
  const renderActionsCell = useCallback((molecule: Molecule) => (
    <Box display="flex" gap={1}>
      <Tooltip title="View Details">
        <IconButton onClick={() => onViewDetails(molecule)} aria-label="View Details">
          <Visibility />
        </IconButton>
      </Tooltip>
      <Tooltip title="Add to Experiment Queue">
        <IconButton onClick={() => onAddToQueue(molecule)} aria-label="Add to Experiment Queue">
          <Queue />
        </IconButton>
      </Tooltip>
      <Tooltip title="Flag as Important">
        <IconButton onClick={() => onFlagChange(molecule, molecule.flag_status === FlagStatus.IMPORTANT ? null : FlagStatus.IMPORTANT)} aria-label="Flag as Important">
          {molecule.flag_status === FlagStatus.IMPORTANT ? <Star /> : <StarBorder />}
        </IconButton>
      </Tooltip>
    </Box>
  ), [onViewDetails, onAddToQueue, onFlagChange]);

  /**
   * Defines the table columns
   */
  const columns = useMemo(() => {
    const baseColumns = [
      {
        field: 'structure',
        headerName: 'Structure',
        width: 100,
        renderCell: renderStructureCell,
      },
      {
        field: 'smiles',
        headerName: 'SMILES',
        width: 250,
        renderCell: renderSmilesCell,
      },
      ...displayProperties.map(propertyName => ({
        field: propertyName,
        headerName: getPropertyDisplayName(propertyName),
        width: 150,
        renderCell: renderPropertyCell(propertyName),
      })),
    ];

    if (showActions) {
      baseColumns.push({
        field: 'actions',
        headerName: 'Actions',
        width: 150,
        renderCell: renderActionsCell,
        sortable: false,
      });
    }

    return baseColumns;
  }, [renderStructureCell, renderSmilesCell, renderPropertyCell, displayProperties, showActions]);

  return (
    <Table
      data={molecules}
      columns={columns}
      pagination={pagination}
      loading={loading}
      error={error}
      selectable={selectable}
      selectedRows={selectedMolecules}
      onSelectionChange={onSelectionChange}
      className={className}
      style={style}
    />
  );
};

export default MoleculeTable;