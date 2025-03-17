import React, { useMemo } from 'react';
import { Box, Typography, Paper, Alert } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+

import Table, { TableColumn } from '../../../components/common/Table';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import { isValidSMILES } from '../../../utils/molecularUtils';
import { CUSTOM_PROPERTY_PREFIX } from '../../../utils/csvUtils';

// Define the props interface for the component
interface PreviewTableProps {
  /** Preview data from CSV file */
  preview: any[];
  /** Original CSV headers */
  headers: string[];
  /** Current mapping configuration of CSV columns to system properties */
  currentMapping: Array<{csvColumn: string, systemProperty: string, isCustom: boolean}>;
  /** Custom names for columns mapped to custom properties */
  customPropertyNames?: {[key: string]: string};
  /** Maximum number of rows to display in preview */
  maxRows?: number;
}

// Styled components for consistent appearance
const PreviewContainer = styled(Paper)(({ theme }) => ({
  padding: '24px',
  marginBottom: '24px',
  borderRadius: '8px'
}));

const PreviewTitle = styled(Typography)(({ theme }) => ({
  marginBottom: '16px',
  fontWeight: '500'
}));

const PreviewDescription = styled(Typography)(({ theme }) => ({
  marginBottom: '24px',
  color: theme.palette.text.secondary
}));

/**
 * Component for displaying a preview of CSV data with applied column mappings
 * during the CSV import process. It shows a sample of rows from the uploaded CSV file
 * with headers transformed according to the current mapping configuration.
 */
const PreviewTable: React.FC<PreviewTableProps> = ({
  preview,
  headers,
  currentMapping,
  customPropertyNames = {},
  maxRows = 5,
}) => {
  // Memoize table columns based on current mapping
  const columns = useMemo((): TableColumn[] => {
    if (!currentMapping || currentMapping.length === 0) return [];
    
    return currentMapping.map(mapping => {
      const isSmilesColumn = mapping.systemProperty.toLowerCase() === 'smiles';
      
      return {
        field: mapping.systemProperty,
        headerName: getDisplayNameForProperty(mapping.systemProperty, mapping.isCustom, mapping.csvColumn),
        width: isSmilesColumn ? 150 : undefined,
        align: 'left',
        sortable: true,
        renderCell: isSmilesColumn 
          ? (row) => renderSMILESCell(row)
          : undefined,
      };
    });
  }, [currentMapping, customPropertyNames]);
  
  // Memoize table data with transformed preview data
  const data = useMemo(() => {
    if (!preview || !currentMapping || currentMapping.length === 0) return [];
    
    // Transform the preview data according to the mapping
    return preview.slice(0, maxRows).map((row, index) => {
      const transformedRow: any = { id: index };
      
      currentMapping.forEach(mapping => {
        transformedRow[mapping.systemProperty] = row[mapping.csvColumn];
      });
      
      return transformedRow;
    });
  }, [preview, currentMapping, maxRows]);
  
  /**
   * Custom renderer for SMILES column cells
   * Displays molecular structure visualization for valid SMILES strings
   */
  const renderSMILESCell = (row: any) => {
    const smiles = row.smiles;
    if (!smiles) return '-';
    
    const isValid = isValidSMILES(smiles);
    
    if (isValid) {
      return (
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <MoleculeViewer 
            smiles={smiles} 
            width={80} 
            height={60} 
            showControls={false}
            showProperties={false}
          />
          <Typography variant="caption" sx={{ wordBreak: 'break-all' }}>
            {smiles}
          </Typography>
        </Box>
      );
    }
    
    return (
      <Typography variant="caption" color="error">
        {smiles} (Invalid SMILES)
      </Typography>
    );
  };
  
  /**
   * Gets display name for a system property
   * 
   * @param systemProperty The property name in the system
   * @param isCustom Whether the property is a custom property
   * @param csvColumn The original CSV column name
   * @returns Formatted display name for the property
   */
  const getDisplayNameForProperty = (systemProperty: string, isCustom: boolean, csvColumn: string) => {
    if (isCustom) {
      // Extract the key without the prefix for custom properties
      const customPropertyKey = systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX) 
        ? systemProperty.substring(CUSTOM_PROPERTY_PREFIX.length) 
        : systemProperty;
      
      return customPropertyNames[customPropertyKey] || csvColumn;
    }
    
    if (systemProperty.toLowerCase() === 'smiles') {
      return 'Structure';
    }
    
    // Format system property name for display (capitalize, replace underscores with spaces)
    return systemProperty
      .split('_')
      .map(word => word.charAt(0).toUpperCase() + word.slice(1))
      .join(' ');
  };
  
  return (
    <PreviewContainer>
      <PreviewTitle variant="h6">Data Preview</PreviewTitle>
      <PreviewDescription variant="body2">
        Preview of your data with the current column mapping applied. Showing up to {maxRows} rows.
      </PreviewDescription>
      
      {(!currentMapping || currentMapping.length === 0) ? (
        <Alert severity="info">
          Please create a mapping to see a preview of your data.
        </Alert>
      ) : (
        <Table 
          data={data}
          columns={columns}
          emptyMessage="No preview data available"
          stickyHeader
          dense
        />
      )}
    </PreviewContainer>
  );
};

export default PreviewTable;