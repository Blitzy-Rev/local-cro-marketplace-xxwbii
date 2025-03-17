import React, { useState, useEffect, useCallback } from 'react';
import { Box, Typography, Paper, Grid, Divider } from '@mui/material'; // @mui/material v5.13.0
import Input from '../../../components/common/Input';
import { ExperimentType } from '../../../types/experiment';

/**
 * Props interface for the ExperimentParameters component
 */
interface ExperimentParametersProps {
  experimentType: ExperimentType | null;
  parameters: Record<string, string>;
  onChange: (parameters: Record<string, string>) => void;
  errors: Record<string, string>;
  disabled: boolean;
}

/**
 * Interface defining default parameters for different experiment types
 */
interface DefaultParameters {
  bindingAssay: Record<string, string>;
  admePanel: Record<string, string>;
  toxicityAssay: Record<string, string>;
}

/**
 * A component that renders and manages experiment parameters based on the selected experiment type
 */
const ExperimentParameters: React.FC<ExperimentParametersProps> = ({
  experimentType,
  parameters = {},
  onChange,
  errors = {},
  disabled = false,
}) => {
  // Define default parameters for different experiment types
  const DEFAULT_PARAMETERS: DefaultParameters = {
    bindingAssay: {
      concentration: '10',
      temperature: '25',
    },
    admePanel: {
      concentration: '5',
      time: '24',
    },
    toxicityAssay: {
      concentration: '1',
      exposureTime: '48',
      cellType: 'HepG2',
    },
  };

  // Local state for parameters
  const [localParameters, setLocalParameters] = useState<Record<string, string>>(parameters || {});

  // Update local parameters when props change
  useEffect(() => {
    setLocalParameters(parameters || {});
  }, [parameters]);

  // Set default parameters when experiment type changes
  useEffect(() => {
    if (!experimentType) return;

    let defaultParams: Record<string, string> = {};
    const name = experimentType.name.toLowerCase();
    const category = experimentType.category.toLowerCase();

    if (name.includes('binding') || category.includes('binding')) {
      defaultParams = DEFAULT_PARAMETERS.bindingAssay;
    } else if (name.includes('adme') || category.includes('adme')) {
      defaultParams = DEFAULT_PARAMETERS.admePanel;
    } else if (name.includes('tox') || category.includes('tox')) {
      defaultParams = DEFAULT_PARAMETERS.toxicityAssay;
    }

    // Merge existing parameters with defaults for missing values
    const mergedParams = { ...defaultParams };
    
    // Only add values from current parameters that aren't already set
    Object.keys(mergedParams).forEach(key => {
      if (parameters && parameters[key]) {
        mergedParams[key] = parameters[key];
      }
    });
    
    // Add any additional parameters that might exist
    if (parameters) {
      Object.keys(parameters).forEach(key => {
        if (!mergedParams[key]) {
          mergedParams[key] = parameters[key];
        }
      });
    }
    
    // Only update if there are changes
    if (JSON.stringify(mergedParams) !== JSON.stringify(parameters)) {
      onChange(mergedParams);
    }
  }, [experimentType, onChange, parameters]);

  // Handle parameter changes
  const handleParameterChange = useCallback(
    (e: React.ChangeEvent<HTMLInputElement>) => {
      const { name, value } = e.target;
      const updatedParameters = { ...localParameters, [name]: value };
      setLocalParameters(updatedParameters);
      onChange(updatedParameters);
    },
    [localParameters, onChange]
  );

  // If no experiment type is selected, don't render anything
  if (!experimentType) {
    return null;
  }

  return (
    <Paper sx={{ p: 3, mt: 2 }}>
      <Typography variant="h6" gutterBottom>
        Experiment Parameters
      </Typography>
      <Divider sx={{ mb: 2 }} />
      
      <Grid container spacing={2}>
        {/* Binding Assay Parameters */}
        {(experimentType.name.toLowerCase().includes('binding') || 
          experimentType.category.toLowerCase().includes('binding')) && (
          <>
            <Grid item xs={12} sm={6}>
              <Input
                label="Concentration"
                name="concentration"
                value={localParameters.concentration || ''}
                onChange={handleParameterChange}
                error={errors.concentration}
                disabled={disabled}
                type="number"
                endAdornment="μM"
                fullWidth
              />
            </Grid>
            <Grid item xs={12} sm={6}>
              <Input
                label="Temperature"
                name="temperature"
                value={localParameters.temperature || ''}
                onChange={handleParameterChange}
                error={errors.temperature}
                disabled={disabled}
                type="number"
                endAdornment="°C"
                fullWidth
              />
            </Grid>
          </>
        )}

        {/* ADME Panel Parameters */}
        {(experimentType.name.toLowerCase().includes('adme') || 
          experimentType.category.toLowerCase().includes('adme')) && (
          <>
            <Grid item xs={12} sm={6}>
              <Input
                label="Concentration"
                name="concentration"
                value={localParameters.concentration || ''}
                onChange={handleParameterChange}
                error={errors.concentration}
                disabled={disabled}
                type="number"
                endAdornment="μM"
                fullWidth
              />
            </Grid>
            <Grid item xs={12} sm={6}>
              <Input
                label="Time"
                name="time"
                value={localParameters.time || ''}
                onChange={handleParameterChange}
                error={errors.time}
                disabled={disabled}
                type="number"
                endAdornment="hours"
                fullWidth
              />
            </Grid>
          </>
        )}

        {/* Toxicity Assay Parameters */}
        {(experimentType.name.toLowerCase().includes('tox') || 
          experimentType.category.toLowerCase().includes('tox')) && (
          <>
            <Grid item xs={12} sm={6}>
              <Input
                label="Concentration"
                name="concentration"
                value={localParameters.concentration || ''}
                onChange={handleParameterChange}
                error={errors.concentration}
                disabled={disabled}
                type="number"
                endAdornment="μM"
                fullWidth
              />
            </Grid>
            <Grid item xs={12} sm={6}>
              <Input
                label="Exposure Time"
                name="exposureTime"
                value={localParameters.exposureTime || ''}
                onChange={handleParameterChange}
                error={errors.exposureTime}
                disabled={disabled}
                type="number"
                endAdornment="hours"
                fullWidth
              />
            </Grid>
            <Grid item xs={12}>
              <Input
                label="Cell Type"
                name="cellType"
                value={localParameters.cellType || ''}
                onChange={handleParameterChange}
                error={errors.cellType}
                disabled={disabled}
                fullWidth
              />
            </Grid>
          </>
        )}

        {/* Generic Parameters for other experiment types */}
        {!experimentType.name.toLowerCase().includes('binding') && 
         !experimentType.category.toLowerCase().includes('binding') && 
         !experimentType.name.toLowerCase().includes('adme') && 
         !experimentType.category.toLowerCase().includes('adme') && 
         !experimentType.name.toLowerCase().includes('tox') && 
         !experimentType.category.toLowerCase().includes('tox') && (
          <Grid item xs={12}>
            <Typography variant="body2" color="textSecondary" sx={{ mb: 2 }}>
              Configure parameters for {experimentType.name}
            </Typography>
            
            <Grid container spacing={2}>
              <Grid item xs={12} sm={6}>
                <Input
                  label="Concentration"
                  name="concentration"
                  value={localParameters.concentration || ''}
                  onChange={handleParameterChange}
                  error={errors.concentration}
                  disabled={disabled}
                  type="number"
                  endAdornment="μM"
                  fullWidth
                />
              </Grid>
              <Grid item xs={12} sm={6}>
                <Input
                  label="Temperature"
                  name="temperature"
                  value={localParameters.temperature || ''}
                  onChange={handleParameterChange}
                  error={errors.temperature}
                  disabled={disabled}
                  type="number"
                  endAdornment="°C"
                  fullWidth
                />
              </Grid>
              <Grid item xs={12}>
                <Input
                  label="Notes"
                  name="notes"
                  value={localParameters.notes || ''}
                  onChange={handleParameterChange}
                  error={errors.notes}
                  disabled={disabled}
                  multiline
                  rows={3}
                  fullWidth
                />
              </Grid>
            </Grid>
          </Grid>
        )}
      </Grid>
    </Paper>
  );
};

export default ExperimentParameters;