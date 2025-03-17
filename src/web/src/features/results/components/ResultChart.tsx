import React, { useMemo, useState, useEffect } from 'react'; // react version ^18.2.0
import { Box, Typography, FormControl, InputLabel, Select, MenuItem, FormControlLabel, Switch } from '@mui/material'; // @mui/material version ^5.13.0

// Internal imports
import { ResultData } from '../../../types/result';
import PropertyChart from '../../../components/data-visualization/PropertyChart';
import Card from '../../../components/common/Card';
import { getPropertyDisplayName, getPropertyUnit } from '../../../utils/molecularUtils';

/**
 * Props for the ResultChart component
 */
interface ResultChartProps {
  /** Array of structured result data for molecules */
  data: ResultData[];
  /** Optional title for the chart section */
  title?: string;
  /** Optional height for the chart */
  height?: number;
  /** Whether to allow users to select which property to visualize */
  allowPropertySelection?: boolean;
  /** Whether to allow users to select the chart type */
  allowChartTypeSelection?: boolean;
  /** Default property to visualize */
  defaultProperty?: string;
  /** Default chart type to use */
  defaultChartType?: string;
  /** Optional callback when a molecule data point is selected */
  onMoleculeSelect?: (moleculeId: string) => void;
}

/**
 * A component that visualizes experimental result data in chart form, allowing researchers to analyze
 * trends and patterns in molecular property data from CRO experiment results.
 */
const ResultChart: React.FC<ResultChartProps> = ({
  data,
  title = 'Result Visualization',
  height = 400,
  allowPropertySelection = true,
  allowChartTypeSelection = true,
  defaultProperty,
  defaultChartType = 'bar',
  onMoleculeSelect
}) => {
  // Extract available properties from data
  const availableProperties = useMemo(() => getAvailableProperties(data), [data]);
  
  // Set state for chart configuration
  const [selectedProperty, setSelectedProperty] = useState<string>(defaultProperty || '');
  const [chartType, setChartType] = useState<string>(defaultChartType);
  const [showLegend, setShowLegend] = useState<boolean>(true);
  const [normalizeValues, setNormalizeValues] = useState<boolean>(false);
  
  // Initialize selected property when availableProperties changes
  useEffect(() => {
    if (availableProperties.length > 0 && (!selectedProperty || !availableProperties.includes(selectedProperty))) {
      setSelectedProperty(defaultProperty || availableProperties[0]);
    }
  }, [availableProperties, defaultProperty, selectedProperty]);
  
  // Process data for the chart
  const chartData = useMemo(() => {
    if (!selectedProperty || !data || data.length === 0) return [];
    return processChartData(data, selectedProperty, normalizeValues);
  }, [data, selectedProperty, normalizeValues]);
  
  // Event handlers
  const handlePropertyChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setSelectedProperty(event.target.value);
  };
  
  const handleChartTypeChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    setChartType(event.target.value);
  };
  
  const handleLegendToggle = (event: React.ChangeEvent<HTMLInputElement>) => {
    setShowLegend(event.target.checked);
  };
  
  const handleNormalizeToggle = (event: React.ChangeEvent<HTMLInputElement>) => {
    setNormalizeValues(event.target.checked);
  };
  
  const handlePointClick = (point: any) => {
    if (onMoleculeSelect && point.originalData && point.originalData.molecule_id) {
      onMoleculeSelect(point.originalData.molecule_id);
    }
  };
  
  // If no data or properties available, show empty state
  if (!data || data.length === 0 || availableProperties.length === 0) {
    return (
      <Card>
        <Box p={2}>
          <Typography variant="h6">{title}</Typography>
          <Box display="flex" justifyContent="center" alignItems="center" height={height} p={2}>
            <Typography color="textSecondary">
              No data available for visualization
            </Typography>
          </Box>
        </Box>
      </Card>
    );
  }
  
  return (
    <Card>
      <Box p={2}>
        <Typography variant="h6">{title}</Typography>
        
        {/* Controls section */}
        <Box mt={2} display="flex" flexWrap="wrap" gap={2} alignItems="center">
          {/* Property selection */}
          {allowPropertySelection && (
            <FormControl size="small" style={{ minWidth: 200 }}>
              <InputLabel id="property-select-label">Property</InputLabel>
              <Select
                labelId="property-select-label"
                id="property-select"
                value={selectedProperty}
                label="Property"
                onChange={handlePropertyChange}
              >
                {availableProperties.map((property) => (
                  <MenuItem key={property} value={property}>
                    {getPropertyDisplayName(property)}
                  </MenuItem>
                ))}
              </Select>
            </FormControl>
          )}
          
          {/* Chart type selection */}
          {allowChartTypeSelection && (
            <FormControl size="small" style={{ minWidth: 150 }}>
              <InputLabel id="chart-type-select-label">Chart Type</InputLabel>
              <Select
                labelId="chart-type-select-label"
                id="chart-type-select"
                value={chartType}
                label="Chart Type"
                onChange={handleChartTypeChange}
              >
                <MenuItem value="bar">Bar Chart</MenuItem>
                <MenuItem value="line">Line Chart</MenuItem>
                <MenuItem value="area">Area Chart</MenuItem>
              </Select>
            </FormControl>
          )}
          
          {/* Toggle switches */}
          <FormControlLabel
            control={<Switch checked={showLegend} onChange={handleLegendToggle} />}
            label="Show Legend"
          />
          <FormControlLabel
            control={<Switch checked={normalizeValues} onChange={handleNormalizeToggle} />}
            label="Normalize Values"
          />
        </Box>
        
        {/* Chart component */}
        <Box mt={2}>
          <PropertyChart
            data={chartData}
            chartType={chartType}
            propertyName={selectedProperty}
            title={getPropertyDisplayName(selectedProperty)}
            xAxisLabel="Molecules"
            yAxisLabel={`${getPropertyDisplayName(selectedProperty)} ${getPropertyUnit(selectedProperty) || ''}`}
            height={height}
            showLegend={showLegend}
            onPointClick={handlePointClick}
          />
        </Box>
      </Box>
    </Card>
  );
};

/**
 * Extracts unique property names from result data
 * @param data Array of result data
 * @returns Array of unique property names
 */
function getAvailableProperties(data: ResultData[]): string[] {
  if (!data || data.length === 0) return [];
  
  // Extract all data_name values and filter to unique values
  const properties = data.map(item => item.data_name);
  return Array.from(new Set(properties)).sort();
}

/**
 * Processes raw result data into a format suitable for the PropertyChart component
 * @param data Array of result data
 * @param selectedProperty Property name to visualize
 * @param normalize Whether to normalize values for comparison
 * @returns Processed data for chart visualization
 */
function processChartData(data: ResultData[], selectedProperty: string, normalize: boolean = false): any[] {
  // Filter data to only include entries with the selected property
  const filteredData = data.filter(item => item.data_name === selectedProperty);
  if (filteredData.length === 0) return [];
  
  // Group data by molecule_id to organize by molecule
  const groupedByMolecule: Record<string, ResultData[]> = {};
  filteredData.forEach(item => {
    if (!groupedByMolecule[item.molecule_id]) {
      groupedByMolecule[item.molecule_id] = [];
    }
    groupedByMolecule[item.molecule_id].push(item);
  });
  
  // Calculate normalization factor if needed
  let maxValue = 1;
  if (normalize) {
    maxValue = Math.max(...filteredData.map(item => {
      const value = typeof item.data_value === 'number' 
        ? item.data_value 
        : parseFloat(item.data_value.toString()) || 0;
      return value;
    }));
    maxValue = maxValue || 1; // Ensure we don't divide by zero
  }
  
  // Transform data for PropertyChart
  const chartData = Object.entries(groupedByMolecule).map(([moleculeId, items], index) => {
    const item = items[0]; // Take the first item for this molecule
    
    // Convert data_value to number
    let value = typeof item.data_value === 'number' 
      ? item.data_value 
      : parseFloat(item.data_value.toString()) || 0;
    
    // Normalize if required
    if (normalize && maxValue > 0) {
      value = value / maxValue;
    }
    
    // Use molecule SMILES or ID as label
    const label = item.molecule?.smiles || `Molecule ${index + 1}`;
    
    return {
      molecule_id: moleculeId,
      [selectedProperty]: value, // Use the property name as the key
      label,
      data_unit: item.data_unit,
      originalData: item
    };
  });
  
  // Sort data by value for proper rendering
  return chartData.sort((a, b) => a[selectedProperty] - b[selectedProperty]);
}

export default ResultChart;