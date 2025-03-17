import React, { useState, useEffect, useMemo, useCallback } from 'react';
import { Box, Typography, Paper, Grid, Button, IconButton, Tooltip } from '@mui/material';
import { Add, Delete } from '@mui/icons-material';

import Input from '../../../components/common/Input';
import Table from '../../../components/common/Table';
import { ResultDataCreate } from '../../../types/result';
import { Molecule } from '../../../types/molecule';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';

// Interface for component props
interface StructuredDataEntryProps {
  molecules: Molecule[];
  value: ResultDataCreate[];
  onChange: (data: ResultDataCreate[]) => void;
  experimentType: string;
  error?: string;
}

// Interface for data fields based on experiment type
interface DataField {
  name: string;
  label: string;
  unit: string;
  required: boolean;
}

// Interface for molecule data entry row
interface MoleculeDataEntry {
  molecule: Molecule;
  dataEntries: ResultDataCreate[];
}

/**
 * A component for entering structured experimental result data for molecules.
 * Provides a tabular interface for entering property values, units, and measurements
 * for each molecule in an experiment, supporting the result upload workflow.
 */
const StructuredDataEntry: React.FC<StructuredDataEntryProps> = ({
  molecules,
  value,
  onChange,
  experimentType,
  error
}) => {
  // State for structured data entries
  const [dataEntries, setDataEntries] = useState<ResultDataCreate[]>(value || []);
  
  // Map to organize data entries by molecule for easier rendering
  const [moleculeDataMap, setMoleculeDataMap] = useState<Map<string, ResultDataCreate[]>>(new Map());
  
  // Map to track validation errors
  const [validationErrors, setValidationErrors] = useState<Map<string, string>>(new Map());

  // Effect to update the state when value prop changes
  useEffect(() => {
    setDataEntries(value || []);
  }, [value]);

  // Effect to organize data by molecule whenever dataEntries change
  useEffect(() => {
    const newMap = new Map<string, ResultDataCreate[]>();
    
    // Initialize map with empty arrays for each molecule
    molecules.forEach(molecule => {
      newMap.set(molecule.id, []);
    });
    
    // Populate map with data entries for each molecule
    dataEntries.forEach(entry => {
      const moleculeEntries = newMap.get(entry.molecule_id) || [];
      newMap.set(entry.molecule_id, [...moleculeEntries, entry]);
    });
    
    setMoleculeDataMap(newMap);
  }, [dataEntries, molecules]);

  /**
   * Returns the appropriate data fields based on experiment type
   */
  const getDataFieldsForExperimentType = useCallback((experimentType: string): DataField[] => {
    // Common data fields for different experiment types
    switch (experimentType.toLowerCase()) {
      case 'adme panel':
        return [
          { name: 'solubility', label: 'Solubility', unit: 'mg/mL', required: true },
          { name: 'permeability', label: 'Permeability', unit: '10^-6 cm/s', required: true },
          { name: 'stability', label: 'Metabolic Stability', unit: 't1/2 (h)', required: true },
          { name: 'toxicity', label: 'Toxicity', unit: 'IC50 (μM)', required: true },
        ];
      case 'binding assay':
        return [
          { name: 'binding_affinity', label: 'Binding Affinity', unit: '%', required: true },
          { name: 'ic50', label: 'IC50', unit: 'nM', required: true },
          { name: 'ki', label: 'Ki', unit: 'nM', required: false },
        ];
      case 'toxicity assay':
        return [
          { name: 'toxicity_score', label: 'Toxicity Score', unit: '%', required: true },
          { name: 'cell_viability', label: 'Cell Viability', unit: '%', required: true },
          { name: 'ld50', label: 'LD50', unit: 'mg/kg', required: false },
        ];
      case 'solubility test':
        return [
          { name: 'aqueous_solubility', label: 'Aqueous Solubility', unit: 'mg/mL', required: true },
          { name: 'ph_dependent_solubility', label: 'pH-Dependent Solubility', unit: 'mg/mL', required: false },
          { name: 'solubility_classification', label: 'Solubility Classification', unit: '', required: true },
        ];
      default:
        // Default fields if experiment type is not recognized
        return [
          { name: 'result_value', label: 'Result Value', unit: '', required: true },
          { name: 'measurement', label: 'Measurement', unit: '', required: false },
        ];
    }
  }, []);

  // Get data fields based on experiment type
  const dataFields = useMemo(() => {
    return getDataFieldsForExperimentType(experimentType);
  }, [experimentType, getDataFieldsForExperimentType]);

  /**
   * Handles changes to data value inputs
   */
  const handleDataValueChange = useCallback((moleculeId: string, dataName: string, event: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = event.target.value;
    
    setDataEntries(prevEntries => {
      // Find the existing entry or create a new one
      const existingEntryIndex = prevEntries.findIndex(
        entry => entry.molecule_id === moleculeId && entry.data_name === dataName
      );
      
      const updatedEntries = [...prevEntries];
      
      if (existingEntryIndex >= 0) {
        // Update existing entry
        updatedEntries[existingEntryIndex] = {
          ...updatedEntries[existingEntryIndex],
          data_value: newValue
        };
      } else {
        // Create new entry
        updatedEntries.push({
          molecule_id: moleculeId,
          data_name: dataName,
          data_value: newValue,
          data_unit: dataFields.find(field => field.name === dataName)?.unit || ''
        });
      }
      
      // Clear validation error if value is now valid
      if (newValue) {
        const errorKey = `${moleculeId}-${dataName}`;
        setValidationErrors(prev => {
          const newErrors = new Map(prev);
          newErrors.delete(errorKey);
          return newErrors;
        });
      }
      
      return updatedEntries;
    });
    
    // Validate data and call onChange with updated entries
    setTimeout(() => {
      const updatedEntries = dataEntries.map(entry => 
        entry.molecule_id === moleculeId && entry.data_name === dataName
          ? { ...entry, data_value: newValue }
          : entry
      );
      
      if (!dataEntries.some(entry => 
        entry.molecule_id === moleculeId && entry.data_name === dataName
      )) {
        updatedEntries.push({
          molecule_id: moleculeId,
          data_name: dataName,
          data_value: newValue,
          data_unit: dataFields.find(field => field.name === dataName)?.unit || ''
        });
      }
      
      onChange(updatedEntries);
    }, 0);
  }, [dataEntries, dataFields, onChange]);

  /**
   * Handles changes to data unit inputs
   */
  const handleDataUnitChange = useCallback((moleculeId: string, dataName: string, event: React.ChangeEvent<HTMLInputElement>) => {
    const newUnit = event.target.value;
    
    setDataEntries(prevEntries => {
      // Find the existing entry or create a new one
      const existingEntryIndex = prevEntries.findIndex(
        entry => entry.molecule_id === moleculeId && entry.data_name === dataName
      );
      
      const updatedEntries = [...prevEntries];
      
      if (existingEntryIndex >= 0) {
        // Update existing entry
        updatedEntries[existingEntryIndex] = {
          ...updatedEntries[existingEntryIndex],
          data_unit: newUnit
        };
      } else {
        // Create new entry with empty value
        updatedEntries.push({
          molecule_id: moleculeId,
          data_name: dataName,
          data_value: '',
          data_unit: newUnit
        });
      }
      
      return updatedEntries;
    });
    
    // Call onChange with updated entries
    setTimeout(() => {
      const updatedEntries = dataEntries.map(entry => 
        entry.molecule_id === moleculeId && entry.data_name === dataName
          ? { ...entry, data_unit: newUnit }
          : entry
      );
      
      if (!dataEntries.some(entry => 
        entry.molecule_id === moleculeId && entry.data_name === dataName
      )) {
        updatedEntries.push({
          molecule_id: moleculeId,
          data_name: dataName,
          data_value: '',
          data_unit: newUnit
        });
      }
      
      onChange(updatedEntries);
    }, 0);
  }, [dataEntries, onChange]);

  /**
   * Adds a new custom data field for a molecule
   */
  const handleAddDataField = useCallback((moleculeId: string) => {
    // Create a new data entry with custom field
    const newEntry: ResultDataCreate = {
      molecule_id: moleculeId,
      data_name: `custom_field_${Date.now()}`,
      data_value: '',
      data_unit: ''
    };
    
    setDataEntries(prevEntries => [...prevEntries, newEntry]);
    onChange([...dataEntries, newEntry]);
  }, [dataEntries, onChange]);

  /**
   * Removes a data field for a molecule
   */
  const handleRemoveDataField = useCallback((moleculeId: string, dataName: string) => {
    setDataEntries(prevEntries => 
      prevEntries.filter(entry => 
        !(entry.molecule_id === moleculeId && entry.data_name === dataName)
      )
    );
    
    // Call onChange with updated entries
    const updatedEntries = dataEntries.filter(entry => 
      !(entry.molecule_id === moleculeId && entry.data_name === dataName)
    );
    onChange(updatedEntries);
    
    // Remove any validation errors for this field
    const errorKey = `${moleculeId}-${dataName}`;
    setValidationErrors(prev => {
      const newErrors = new Map(prev);
      newErrors.delete(errorKey);
      return newErrors;
    });
  }, [dataEntries, onChange]);

  /**
   * Validates all data entries
   */
  const validateData = useCallback((): boolean => {
    const newErrors = new Map<string, string>();
    
    // Check for required fields
    molecules.forEach(molecule => {
      const moleculeEntries = moleculeDataMap.get(molecule.id) || [];
      
      // Check each required field
      dataFields.filter(field => field.required).forEach(field => {
        const entry = moleculeEntries.find(e => e.data_name === field.name);
        
        if (!entry || !entry.data_value) {
          const errorKey = `${molecule.id}-${field.name}`;
          newErrors.set(errorKey, `${field.label} is required`);
        }
      });
      
      // Validate numeric values
      moleculeEntries.forEach(entry => {
        if (entry.data_value && isNaN(Number(entry.data_value))) {
          const field = dataFields.find(f => f.name === entry.data_name);
          const fieldName = field ? field.label : entry.data_name;
          const errorKey = `${molecule.id}-${entry.data_name}`;
          newErrors.set(errorKey, `${fieldName} must be a number`);
        }
      });
    });
    
    setValidationErrors(newErrors);
    return newErrors.size === 0;
  }, [molecules, moleculeDataMap, dataFields]);

  // Call validateData whenever data entries change
  useEffect(() => {
    validateData();
  }, [dataEntries, validateData]);

  /**
   * Organizes data entries by molecule for display
   */
  const organizeDataByMolecule = useCallback((): MoleculeDataEntry[] => {
    return molecules.map(molecule => ({
      molecule,
      dataEntries: moleculeDataMap.get(molecule.id) || []
    }));
  }, [molecules, moleculeDataMap]);

  // Create table columns configuration
  const columns = useMemo(() => [
    {
      field: 'molecule',
      headerName: 'Molecule',
      width: 180,
      renderCell: (row: MoleculeDataEntry) => (
        <Box sx={{ width: 150, height: 100 }}>
          <MoleculeViewer 
            molecule={row.molecule} 
            width={150} 
            height={100} 
            showControls={false}
            showProperties={false}
          />
        </Box>
      )
    },
    {
      field: 'smiles',
      headerName: 'SMILES',
      width: 250,
      renderCell: (row: MoleculeDataEntry) => (
        <Typography variant="body2" noWrap sx={{ maxWidth: 220 }} title={row.molecule.smiles}>
          {row.molecule.smiles}
        </Typography>
      )
    },
    {
      field: 'data',
      headerName: 'Result Data',
      width: 500,
      renderCell: (row: MoleculeDataEntry) => (
        <Box sx={{ width: '100%' }}>
          <Grid container spacing={2}>
            {/* Render input fields for each data field based on experiment type */}
            {dataFields.map((field) => {
              const entry = row.dataEntries.find(e => e.data_name === field.name);
              const errorKey = `${row.molecule.id}-${field.name}`;
              const hasError = validationErrors.has(errorKey);
              const errorMessage = validationErrors.get(errorKey);
              
              return (
                <Grid item xs={12} key={field.name}>
                  <Grid container spacing={1} alignItems="center">
                    <Grid item xs={4}>
                      <Typography variant="body2">{field.label}:</Typography>
                    </Grid>
                    <Grid item xs={5}>
                      <Input
                        type="text"
                        value={entry?.data_value || ''}
                        onChange={(e) => handleDataValueChange(row.molecule.id, field.name, e)}
                        error={hasError ? errorMessage : false}
                        placeholder={`Enter ${field.label}`}
                        required={field.required}
                        size="small"
                        margin="none"
                      />
                    </Grid>
                    <Grid item xs={3}>
                      <Input
                        type="text"
                        value={entry?.data_unit || field.unit}
                        onChange={(e) => handleDataUnitChange(row.molecule.id, field.name, e)}
                        placeholder="Unit"
                        size="small"
                        margin="none"
                      />
                    </Grid>
                  </Grid>
                </Grid>
              );
            })}
            
            {/* Render custom data fields */}
            {row.dataEntries
              .filter(entry => !dataFields.some(field => field.name === entry.data_name))
              .map((entry) => (
                <Grid item xs={12} key={entry.data_name}>
                  <Grid container spacing={1} alignItems="center">
                    <Grid item xs={4}>
                      <Input
                        type="text"
                        value={entry.data_name.replace('custom_field_', 'Custom Field ')}
                        onChange={(e) => {
                          const newDataName = e.target.value;
                          handleDataValueChange(
                            row.molecule.id,
                            entry.data_name,
                            { target: { value: entry.data_value } } as React.ChangeEvent<HTMLInputElement>
                          );
                        }}
                        placeholder="Custom Field"
                        size="small"
                        margin="none"
                      />
                    </Grid>
                    <Grid item xs={4}>
                      <Input
                        type="text"
                        value={entry.data_value}
                        onChange={(e) => handleDataValueChange(row.molecule.id, entry.data_name, e)}
                        placeholder="Value"
                        size="small"
                        margin="none"
                      />
                    </Grid>
                    <Grid item xs={3}>
                      <Input
                        type="text"
                        value={entry.data_unit || ''}
                        onChange={(e) => handleDataUnitChange(row.molecule.id, entry.data_name, e)}
                        placeholder="Unit"
                        size="small"
                        margin="none"
                      />
                    </Grid>
                    <Grid item xs={1}>
                      <IconButton
                        size="small"
                        onClick={() => handleRemoveDataField(row.molecule.id, entry.data_name)}
                        aria-label="Remove field"
                      >
                        <Delete fontSize="small" />
                      </IconButton>
                    </Grid>
                  </Grid>
                </Grid>
              ))}
            
            {/* Add custom field button */}
            <Grid item xs={12}>
              <Button
                startIcon={<Add />}
                size="small"
                onClick={() => handleAddDataField(row.molecule.id)}
                sx={{ mt: 1 }}
              >
                Add Custom Field
              </Button>
            </Grid>
          </Grid>
        </Box>
      )
    }
  ], [dataFields, validationErrors, handleDataValueChange, handleDataUnitChange, handleAddDataField, handleRemoveDataField]);

  // Organize data for rendering
  const tableData = useMemo(() => organizeDataByMolecule(), [organizeDataByMolecule]);

  return (
    <Paper elevation={1} sx={{ p: 2, mb: 3 }}>
      <Typography variant="h6" gutterBottom>
        Structured Data Entry
      </Typography>
      
      <Box sx={{ mb: 2 }}>
        <Typography variant="body2" color="text.secondary">
          Enter result data for each molecule. Fields marked with * are required.
        </Typography>
      </Box>
      
      <Table 
        data={tableData}
        columns={columns}
        error={error}
        stickyHeader
        maxHeight={600}
      />
      
      {error && (
        <Typography variant="body2" color="error" sx={{ mt: 2 }}>
          {error}
        </Typography>
      )}
    </Paper>
  );
};

export default StructuredDataEntry;