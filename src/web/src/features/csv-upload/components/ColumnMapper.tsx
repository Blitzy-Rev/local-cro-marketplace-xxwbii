import React, { useState, useEffect, useMemo, useCallback } from 'react'; // react ^18.2.0
import {
  Box,
  Typography,
  Paper,
  Alert,
  Grid,
  Divider,
  TextField,
  InputAdornment,
  CircularProgress
} from '@mui/material'; // v5.13.0
import {
  AutoAwesome,
  Error,
  Warning
} from '@mui/icons-material'; // v5.13.0
import { styled } from '@mui/material/styles'; // v5.13.0

import { useCSVImport } from '../hooks/useCSVImport';
import PreviewTable from './PreviewTable';
import Dropdown from '../../../components/common/Dropdown';
import Table from '../../../components/common/Table';
import Button from '../../../components/common/Button';
import { validateMappingConfiguration, CUSTOM_PROPERTY_PREFIX } from '../../../utils/csvUtils';
import { useToast } from '../../../hooks/useToast';

/**
 * Interface for the props of the ColumnMapper component.
 * Defines the callback functions that are used to communicate
 * with the parent component about the mapping status.
 */
interface ColumnMapperProps {
  /** Callback function called when mapping is complete and submitted */
  onMappingComplete: () => void;
  /** Callback function called when mapping validity changes */
  onMappingChange: (isValid: boolean) => void;
}

/**
 * Interface for mapping validation errors.
 * Defines the structure of the error objects that are used
 * to display validation errors to the user.
 */
interface ValidationError {
  /** CSV column with error */
  column: string;
  /** Error message */
  message: string;
}

/**
 * Interface for tracking custom property names.
 * Defines the structure for storing custom property names
 * associated with CSV columns.
 */
interface CustomPropertyNames {
  [key: string]: string;
}

// Styled component for the mapping container
const MappingContainer = styled(Paper)(({ theme }) => ({
  padding: '24px',
  marginBottom: '24px',
  borderRadius: '8px'
}));

// Styled component for the mapping title
const MappingTitle = styled(Typography)(({ theme }) => ({
  marginBottom: '16px',
  fontWeight: '500'
}));

// Styled component for the mapping description
const MappingDescription = styled(Typography)(({ theme }) => ({
  marginBottom: '24px',
  color: theme.palette.text.secondary
}));

// Styled component for the mapping table
const MappingTable = styled(Box)(({ theme }) => ({
  marginBottom: '24px'
}));

// Styled component for the action container
const ActionContainer = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  marginTop: '24px',
  marginBottom: '24px'
}));

/**
 * Component for mapping CSV columns to system properties.
 * This component provides an interactive interface for users to select
 * appropriate property mappings, validates the mapping configuration,
 * and displays a preview of the data with applied mappings.
 *
 * @param {ColumnMapperProps} props - The props for the component, including callback functions.
 * @returns {JSX.Element} - The rendered ColumnMapper component.
 */
const ColumnMapper: React.FC<ColumnMapperProps> = ({ onMappingComplete, onMappingChange }) => {
  // Extract necessary state and functions from useCSVImport hook
  const {
    headers,
    preview,
    rowCount,
    currentMapping,
    mappingOptions,
    availableProperties,
    handleUpdateMapping,
    handleGenerateMapping,
    handleSubmitMapping,
    isLoading,
    errors
  } = useCSVImport();

  // Initialize state for mapping validation errors
  const [validationErrors, setValidationErrors] = useState<ValidationError[]>([]);

  // Initialize state for custom property names
  const [customPropertyNames, setCustomPropertyNames] = useState<CustomPropertyNames>({});

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  /**
   * Updates the mapping when a property selection changes.
   * This function is called when the user selects a new property
   * from the dropdown for a specific CSV column.
   *
   * @param {string} csvColumn - The name of the CSV column being mapped.
   * @param {string} systemProperty - The name of the system property to map to.
   */
  const handlePropertyChange = (csvColumn: string, systemProperty: string) => {
    // Determine if the selected property is a custom property
    const isCustom = systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX);

    // Call handleUpdateMapping from useCSVImport with appropriate parameters
    handleUpdateMapping(csvColumn, systemProperty, isCustom);

    // If custom property, initialize custom property name if not already set
    if (isCustom && !customPropertyNames[csvColumn]) {
      setCustomPropertyNames(prevNames => ({ ...prevNames, [csvColumn]: csvColumn }));
    }
  };

  /**
   * Updates the custom property name for a column.
   * This function is called when the user edits the custom property name
   * in the TextField for a specific CSV column.
   *
   * @param {string} csvColumn - The name of the CSV column being mapped.
   * @param {string} customName - The new custom property name.
   */
  const handleCustomPropertyChange = (csvColumn: string, customName: string) => {
    // Update customPropertyNames state with new name for the column
    setCustomPropertyNames(prevNames => ({ ...prevNames, [csvColumn]: customName }));
  };

  /**
   * Triggers automatic mapping generation.
   * This function is called when the user clicks the "Auto Map" button.
   */
  const handleAutoMap = () => {
    // Call handleGenerateMapping from useCSVImport
    handleGenerateMapping();

    // Show success toast notification
    showToast({
      type: 'success',
      message: 'Mapping automatically generated',
    });
  };

  /**
   * Validates and submits the current mapping.
   * This function is called when the user clicks the "Import Molecules" button.
   */
  const handleSubmit = async () => {
    // Validate the current mapping
    const validationResult = validateMappingConfiguration(currentMapping);

    // If invalid, show error toast with validation errors
    if (!validationResult.isValid) {
      setValidationErrors(validationResult.errors);
      showToast({
        type: 'error',
        message: 'Invalid mapping configuration. Please correct the errors.',
      });
      return;
    }

    // If valid, call handleSubmitMapping from useCSVImport
    try {
      await handleSubmitMapping();

      // Call onMappingComplete callback when submission is successful
      onMappingComplete();
    } catch (error: any) {
      // Handle any errors during submission
      showToast({
        type: 'error',
        message: error.message || 'Failed to submit mapping',
      });
    }
  };

  // Use useEffect to validate mapping when it changes
  useEffect(() => {
    const validationResult = validateMappingConfiguration(currentMapping);
    setValidationErrors(validationResult.errors);
  }, [currentMapping]);

  // Use useEffect to notify parent component about mapping validity changes
  useEffect(() => {
    onMappingChange(validationErrors.length === 0);
  }, [validationErrors, onMappingChange]);

  return (
    <MappingContainer>
      <MappingTitle variant="h6">Map CSV Columns</MappingTitle>
      <MappingDescription variant="body2">
        Map the columns from your CSV file to the appropriate system properties.
      </MappingDescription>

      {errors.length > 0 && (
        <Alert severity="error">
          {errors.map((error, index) => (
            <div key={index}>{error}</div>
          ))}
        </Alert>
      )}

      <MappingTable>
        <Table
          data={headers.map((header, index) => ({
            id: index,
            csvColumn: header,
            systemProperty: currentMapping.find(item => item.csvColumn === header)?.systemProperty || '',
          }))}
          columns={[
            {
              field: 'csvColumn',
              headerName: 'CSV Header',
              width: 200,
            },
            {
              field: 'systemProperty',
              headerName: 'System Property',
              width: 300,
              renderCell: (row) => {
                const csvColumn = row.csvColumn;
                const systemProperty = currentMapping.find(item => item.csvColumn === csvColumn)?.systemProperty || '';
                const mappingOption = mappingOptions.find(item => item.csvColumn === csvColumn);
                const isCustom = systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX);

                return (
                  <Grid container spacing={2} alignItems="center">
                    <Grid item xs={6}>
                      {mappingOption ? (
                        <Dropdown
                          id={`property-select-${csvColumn}`}
                          name={`property-select-${csvColumn}`}
                          label="Select Property"
                          value={systemProperty}
                          onChange={(e) => handlePropertyChange(csvColumn, e.target.value)}
                          options={mappingOption.mappingOptions.map(option => ({
                            value: option.value,
                            label: option.label,
                          }))}
                          fullWidth
                        />
                      ) : (
                        <div>No mapping options available</div>
                      )}
                    </Grid>
                    {isCustom && (
                      <Grid item xs={6}>
                        <TextField
                          label="Custom Property Name"
                          value={customPropertyNames[csvColumn] || ''}
                          onChange={(e) => handleCustomPropertyChange(csvColumn, e.target.value)}
                          fullWidth
                          size="small"
                          InputProps={{
                            startAdornment: <InputAdornment position="start">Custom:</InputAdornment>,
                          }}
                        />
                      </Grid>
                    )}
                  </Grid>
                );
              },
            },
          ]}
        />
      </MappingTable>

      <ActionContainer>
        <Button
          variant="outlined"
          startIcon={<AutoAwesome />}
          onClick={handleAutoMap}
        >
          Auto Map
        </Button>
        <Button
          variant="contained"
          color="primary"
          onClick={handleSubmit}
          disabled={isLoading || validationErrors.length > 0}
        >
          {isLoading ? <CircularProgress size={24} color="inherit" /> : 'Import Molecules'}
        </Button>
      </ActionContainer>

      <PreviewTable
        preview={preview}
        headers={headers}
        currentMapping={currentMapping}
        customPropertyNames={customPropertyNames}
      />
    </MappingContainer>
  );
};

export default ColumnMapper;