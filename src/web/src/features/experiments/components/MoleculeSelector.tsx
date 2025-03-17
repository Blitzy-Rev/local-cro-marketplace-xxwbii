import React, { useState, useEffect, useCallback, useMemo } from 'react'; // React v18.2+
import { Box, Grid, Typography, Divider, Paper, Chip, IconButton, Tooltip } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { Add, Remove, FilterList, CheckCircle, Cancel } from '@mui/icons-material'; // v5.13+

import useMolecules from '../../molecules/hooks/useMolecules';
import MoleculeTable from '../../molecules/components/MoleculeTable';
import MoleculeFilter from '../../molecules/components/MoleculeFilter';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import Button from '../../../components/common/Button';
import { Molecule, MoleculeFilter as MoleculeFilterType } from '../../../types/molecule';

/**
 * Props interface for the MoleculeSelector component
 */
interface MoleculeSelectorProps {
  /** Array of currently selected molecule IDs */
  selectedMoleculeIds: string[];
  /** Callback when selection changes */
  onSelectionChange: (moleculeIds: string[]) => void;
  /** Initial filter criteria for molecules */
  initialFilter?: MoleculeFilterType;
  /** Maximum number of molecules that can be selected */
  maxSelectionCount?: number;
  /** Whether to show only selected molecules */
  showSelectedOnly?: boolean;
  /** Additional CSS class name */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
}

// Define styled components for the MoleculeSelector UI
const SelectorContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  gap: '16px',
  width: '100%',
});

const FilterSection = styled(Paper)({
  padding: '16px',
  marginBottom: '16px',
});

const SelectionSummary = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  padding: '8px 16px',
  backgroundColor: theme.palette.background.paper,
  borderRadius: '4px',
  marginBottom: '16px',
}));

const SelectionActions = styled(Box)({
  display: 'flex',
  gap: '8px',
});

const SelectionCount = styled(Chip)({
  fontWeight: 'bold',
});

const WarningMessage = styled(Typography)(({ theme }) => ({
  color: theme.palette.warning.main,
  fontSize: '0.875rem',
  marginTop: '8px',
}));

/**
 * A component that provides an interface for selecting molecules to include in experiments
 */
const MoleculeSelector: React.FC<MoleculeSelectorProps> = ({
  selectedMoleculeIds = [],
  onSelectionChange,
  initialFilter = {},
  maxSelectionCount = 100,
  showSelectedOnly = false,
  className = '',
  style,
}) => {
  // LD1: Initialize the useMolecules hook with initialFilter prop
  const {
    molecules,
    loading,
    error,
    filter,
    selectedMolecules,
    totalMolecules,
    pagination,
    setFilter,
    selectMolecule,
    deselectMolecule,
    selectAllMolecules,
    deselectAllMolecules,
  } = useMolecules({ initialFilter, enableAutoFetch: true });

  // LD1: Set up state for showing selected molecules only
  const [viewSelectedOnly, setViewSelectedOnly] = useState(showSelectedOnly);

  // LD1: Handle selection changes by calling onSelectionChange prop
  const handleSelectionChange = useCallback(
    (selectedIds: string[]) => {
      onSelectionChange(selectedIds);
    },
    [onSelectionChange]
  );

  // LD1: Handle filter changes through the MoleculeFilter component
  const handleFilterChange = useCallback(
    (newFilter: MoleculeFilterType) => {
      setFilter(newFilter);
    },
    [setFilter]
  );

  // LD1: Toggle between showing all molecules and selected molecules only
  const handleToggleViewSelected = () => {
    setViewSelectedOnly((prev) => !prev);
  };

  // LD1: Select all molecules up to maxSelectionCount limit
  const handleSelectAll = () => {
    const availableMoleculeIds = molecules.map((molecule) => molecule.id);
    const selectableMoleculeIds = availableMoleculeIds.slice(0, maxSelectionCount);
    selectAllMolecules(selectableMoleculeIds);
  };

  // LD1: Clear all selected molecules
  const handleClearSelection = () => {
    deselectAllMolecules();
  };

  // LD1: Filter molecules based on viewSelectedOnly state
  const displayedMolecules = useMemo(() => {
    if (viewSelectedOnly) {
      return molecules.filter((molecule) => selectedMoleculeIds.includes(molecule.id));
    }
    return molecules;
  }, [molecules, selectedMoleculeIds, viewSelectedOnly]);

  // LD1: Calculate percentage of selection limit used
  const selectionPercentage = useMemo(() => {
    return (selectedMolecules.length / maxSelectionCount) * 100;
  }, [selectedMolecules.length, maxSelectionCount]);

  // LD1: Check if selection is approaching the maximum limit
  const isApproachingLimit = useMemo(() => {
    return selectionPercentage > 80 && selectionPercentage < 100;
  }, [selectionPercentage]);

  // LD1: Check if selection has reached the maximum limit
  const isAtLimit = useMemo(() => {
    return selectionPercentage >= 100;
  }, [selectionPercentage]);

  // LD1: Render the component with filter, table, and selection summary
  return (
    <SelectorContainer className={className} style={style}>
      {/* LD1: Render MoleculeFilter for filtering available molecules */}
      <FilterSection>
        <MoleculeFilter initialFilter={initialFilter} onFilterChange={handleFilterChange} />
      </FilterSection>

      {/* LD1: Render SelectionSummary with selection count and actions */}
      <SelectionSummary>
        <Box>
          {/* LD1: Display selected molecule count and selection limit information */}
          <SelectionCount label={`${selectedMolecules.length} / ${maxSelectionCount} Molecules Selected`} />
          {/* LD1: Show warning when approaching or reaching selection limit */}
          {isApproachingLimit && (
            <WarningMessage>
              Approaching selection limit ({maxSelectionCount} molecules)
            </WarningMessage>
          )}
          {isAtLimit && (
            <WarningMessage>
              Maximum selection limit reached ({maxSelectionCount} molecules)
            </WarningMessage>
          )}
        </Box>
        {/* LD1: Provide buttons for selecting all and clearing selection */}
        <SelectionActions>
          <Button variant="outlined" size="small" onClick={handleToggleViewSelected}>
            {viewSelectedOnly ? 'Show All Molecules' : 'Show Selected Only'}
          </Button>
          <Button
            variant="contained"
            size="small"
            onClick={handleSelectAll}
            disabled={isAtLimit}
          >
            Select All
          </Button>
          <Button variant="outlined" size="small" onClick={handleClearSelection}>
            Clear Selection
          </Button>
        </SelectionActions>
      </SelectionSummary>

      {/* LD1: Render MoleculeTable with appropriate props for selection */}
      <MoleculeTable
        molecules={displayedMolecules}
        loading={loading}
        error={error}
        selectable={true}
        selectedMolecules={selectedMolecules}
        onSelectionChange={handleSelectionChange}
        pagination={pagination}
        displayProperties={['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors']}
        showActions={false}
      />
    </SelectorContainer>
  );
};

export default MoleculeSelector;