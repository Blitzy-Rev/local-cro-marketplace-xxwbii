import React, { useMemo } from 'react';
import { Box, Typography, Chip } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { CalculateOutlined } from '@mui/icons-material'; // v5.13+
import Table, { TableProps } from '../../../components/common/Table';
import { Molecule, MoleculeProperty } from '../../../types/molecule';
import { 
  getPropertyDisplayName, 
  formatPropertyValue, 
  DEFAULT_MOLECULE_PROPERTIES 
} from '../../../utils/molecularUtils';

// Styled component for the calculated property indicator
const StyledChip = styled(Chip)({
  height: '20px',
  fontSize: '0.7rem',
  marginLeft: '8px',
});

/**
 * Props interface for the PropertyTable component
 */
interface PropertyTableProps {
  /** Molecule object containing properties to display */
  molecule: Molecule;
  /** Optional list of property names to display, if undefined shows all properties */
  propertiesToShow?: string[];
  /** Whether to show an indicator for calculated properties */
  showCalculatedIndicator?: boolean;
  /** Optional CSS class name */
  className?: string;
  /** Optional inline styles */
  style?: React.CSSProperties;
}

/**
 * Interface for property data in table format
 */
interface PropertyTableData {
  /** Unique identifier for the property row */
  id: string;
  /** Display name of the property */
  name: string;
  /** Formatted property value with units */
  value: string;
  /** Whether the property is calculated or imported */
  isCalculated: boolean;
  /** Raw property value for sorting */
  rawValue: any;
}

/**
 * A component that displays molecular properties in a structured table format
 * 
 * @param props Component props
 * @returns The rendered property table component
 */
const PropertyTable: React.FC<PropertyTableProps> = ({
  molecule,
  propertiesToShow,
  showCalculatedIndicator = true,
  className,
  style,
}) => {
  // Transform molecule properties into table data format
  const propertyData = useMemo(() => {
    if (!molecule || !molecule.properties || molecule.properties.length === 0) {
      return [];
    }

    // Filter properties if propertiesToShow is provided
    let filteredProperties = molecule.properties;
    if (propertiesToShow && propertiesToShow.length > 0) {
      filteredProperties = molecule.properties.filter(
        prop => propertiesToShow.includes(prop.property_name)
      );
    }

    // Transform properties to table data format
    return filteredProperties.map(prop => ({
      id: `${prop.property_name}`,
      name: getPropertyDisplayName(prop.property_name),
      value: formatPropertyValue(prop.property_name, prop.property_value, true),
      isCalculated: prop.is_calculated,
      rawValue: prop.property_value,
    }));
  }, [molecule, propertiesToShow]);

  // Define table columns
  const columns = useMemo(() => [
    {
      field: 'name',
      headerName: 'Property',
      width: '40%',
      sortable: true,
    },
    {
      field: 'value',
      headerName: 'Value',
      width: '60%',
      sortable: true,
      renderCell: (row: PropertyTableData) => (
        <Box display="flex" alignItems="center">
          {row.value}
          {showCalculatedIndicator && row.isCalculated && (
            <StyledChip
              size="small"
              color="info"
              label="Calculated"
              icon={<CalculateOutlined fontSize="small" />}
            />
          )}
        </Box>
      ),
      // Add valueGetter for proper sorting by raw value
      valueGetter: (row: PropertyTableData) => row.rawValue,
    },
  ], [showCalculatedIndicator]);

  // Handle empty state
  if (!molecule || !molecule.properties || propertyData.length === 0) {
    return (
      <Box className={className} style={style} p={2} textAlign="center">
        <Typography variant="body2" color="textSecondary">
          No properties available for this molecule
        </Typography>
      </Box>
    );
  }

  // Render table with property data
  return (
    <Table
      data={propertyData}
      columns={columns}
      className={className}
      style={style}
      sortable={true}
      defaultSortField="name"
      defaultSortDirection="asc"
      emptyMessage="No properties available"
      dense={true}
      idField="id"
      stickyHeader={false}
    />
  );
};

export default PropertyTable;