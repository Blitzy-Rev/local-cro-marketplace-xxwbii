import React, { useState, useEffect, useCallback, useMemo } from 'react';
import {
  Box,
  Grid,
  Typography,
  IconButton,
  Chip,
  Collapse,
  Divider,
  FormControlLabel,
  Switch,
  SelectChangeEvent
} from '@mui/material';
import { styled } from '@mui/material/styles';
import {
  Add,
  Remove,
  FilterList,
  Clear,
  ExpandMore,
  ExpandLess
} from '@mui/icons-material';
import { DatePicker } from '@mui/x-date-pickers'; // @mui/x-date-pickers v6.0+

import Input from '../../../components/common/Input';
import Dropdown from '../../../components/common/Dropdown';
import Button from '../../../components/common/Button';
import useDebounce from '../../../hooks/useDebounce';

import {
  MoleculeFilter as MoleculeFilterType,
  PropertyFilter,
  FilterOperator,
  FlagStatus
} from '../../../types/molecule';

import {
  DEFAULT_MOLECULE_PROPERTIES,
  PROPERTY_DISPLAY_NAMES,
  getPropertyDisplayName
} from '../../../utils/molecularUtils';

import {
  validateSMILES,
  validateMolecularProperty,
  PROPERTY_RANGES
} from '../../../utils/validation';

// Interface for the MoleculeFilter component props
interface MoleculeFilterProps {
  initialFilter?: MoleculeFilterType;
  onFilterChange: (filter: MoleculeFilterType) => void;
  className?: string;
  style?: React.CSSProperties;
}

// Define styled components for the filter UI
const FilterContainer = styled(Box)(({ theme }) => ({
  padding: '16px',
  borderRadius: '8px',
  border: '1px solid',
  borderColor: theme.palette.divider,
  backgroundColor: theme.palette.background.paper
}));

const FilterHeader = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  marginBottom: '16px'
});

const PropertyFilterItem = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  gap: '8px',
  marginBottom: '8px'
});

const FilterActions = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  marginTop: '16px'
});

/**
 * A component that provides an interactive interface for filtering molecules by various criteria
 * including SMILES pattern, properties, flag status, and date ranges
 */
const MoleculeFilter: React.FC<MoleculeFilterProps> = ({
  initialFilter,
  onFilterChange,
  className,
  style
}) => {
  // Initialize filter state with provided initial filter or defaults
  const [filter, setFilter] = useState<MoleculeFilterType>(
    initialFilter || {
      smiles_pattern: '',
      property_filters: [],
      flag_status: null,
      created_after: null,
      created_before: null,
      sort_by: 'created_at',
      sort_desc: true
    }
  );
  
  // State for toggling advanced filter visibility
  const [showAdvancedFilters, setShowAdvancedFilters] = useState(false);
  
  // Create debounced filter to prevent excessive API calls during user input
  const debouncedFilter = useDebounce(filter, 500);
  
  // When debounced filter changes, notify parent component
  useEffect(() => {
    onFilterChange(debouncedFilter);
  }, [debouncedFilter, onFilterChange]);
  
  // Generate property options for dropdown from DEFAULT_MOLECULE_PROPERTIES
  const propertyOptions = useMemo(() => 
    DEFAULT_MOLECULE_PROPERTIES.map(prop => ({
      value: prop,
      label: getPropertyDisplayName(prop)
    })),
    []
  );
  
  // Generate options for filter operators
  const operatorOptions = [
    { value: FilterOperator.EQUALS, label: 'Equals' },
    { value: FilterOperator.NOT_EQUALS, label: 'Not Equals' },
    { value: FilterOperator.GREATER_THAN, label: 'Greater Than' },
    { value: FilterOperator.LESS_THAN, label: 'Less Than' },
    { value: FilterOperator.BETWEEN, label: 'Between' },
    { value: FilterOperator.CONTAINS, label: 'Contains' }
  ];
  
  // Generate options for flag status
  const flagStatusOptions = [
    { value: '', label: 'Any' },
    { value: FlagStatus.IMPORTANT, label: 'Important' },
    { value: FlagStatus.FAVORITE, label: 'Favorite' },
    { value: FlagStatus.REVIEW, label: 'Review' },
    { value: FlagStatus.ARCHIVED, label: 'Archived' }
  ];
  
  // Generate options for sorting
  const sortOptions = [
    { value: 'created_at', label: 'Creation Date' },
    ...propertyOptions
  ];
  
  // Handler for updating SMILES pattern filter
  const handleSmilesChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    setFilter(prev => ({
      ...prev,
      smiles_pattern: value
    }));
  };
  
  // Handler for updating flag status filter
  const handleFlagStatusChange = (event: SelectChangeEvent<any>) => {
    const value = event.target.value;
    setFilter(prev => ({
      ...prev,
      flag_status: value === '' ? null : value
    }));
  };
  
  // Handler for updating sort field
  const handleSortChange = (event: SelectChangeEvent<any>) => {
    const value = event.target.value;
    setFilter(prev => ({
      ...prev,
      sort_by: value
    }));
  };
  
  // Handler for updating sort direction
  const handleSortDirectionChange = (event: React.ChangeEvent<HTMLInputElement>) => {
    const checked = event.target.checked;
    setFilter(prev => ({
      ...prev,
      sort_desc: checked
    }));
  };
  
  // Handler for updating date range filters
  const handleDateChange = (type: 'created_after' | 'created_before', date: Date | null) => {
    setFilter(prev => ({
      ...prev,
      [type]: date ? date.toISOString() : null
    }));
  };
  
  // Handler for adding a new property filter
  const handleAddPropertyFilter = () => {
    // Default to first property in the list and equals operator
    const newPropertyFilter: PropertyFilter = {
      name: DEFAULT_MOLECULE_PROPERTIES[0],
      operator: FilterOperator.EQUALS,
      value: ''
    };
    
    setFilter(prev => ({
      ...prev,
      property_filters: [...(prev.property_filters || []), newPropertyFilter]
    }));
  };
  
  // Handler for removing a property filter
  const handleRemovePropertyFilter = (index: number) => {
    setFilter(prev => ({
      ...prev,
      property_filters: (prev.property_filters || []).filter((_, i) => i !== index)
    }));
  };
  
  // Handler for updating a property filter
  const handlePropertyFilterChange = (index: number, field: keyof PropertyFilter, value: any) => {
    setFilter(prev => {
      const updatedFilters = [...(prev.property_filters || [])];
      updatedFilters[index] = {
        ...updatedFilters[index],
        [field]: value
      };
      
      // If the property name or operator changes, reset the value
      if (field === 'name' || field === 'operator') {
        updatedFilters[index].value = '';
        // Clear min and max if operator is no longer "between"
        if (field === 'operator' && value !== FilterOperator.BETWEEN) {
          updatedFilters[index].min_value = undefined;
          updatedFilters[index].max_value = undefined;
        }
        // Set up min/max values for "between" operator
        if (field === 'operator' && value === FilterOperator.BETWEEN) {
          const propertyName = updatedFilters[index].name;
          const range = PROPERTY_RANGES[propertyName as keyof typeof PROPERTY_RANGES];
          if (range) {
            updatedFilters[index].min_value = range.min;
            updatedFilters[index].max_value = range.max;
          }
        }
      }
      
      return {
        ...prev,
        property_filters: updatedFilters
      };
    });
  };
  
  // Handler for clearing all filters
  const handleClearFilters = () => {
    setFilter({
      smiles_pattern: '',
      property_filters: [],
      flag_status: null,
      created_after: null,
      created_before: null,
      sort_by: 'created_at',
      sort_desc: true
    });
  };
  
  // Toggle advanced filter visibility
  const toggleAdvancedFilters = () => {
    setShowAdvancedFilters(prev => !prev);
  };
  
  return (
    <FilterContainer className={className} style={style}>
      <FilterHeader>
        <Typography variant="h6" component="div" display="flex" alignItems="center">
          <FilterList sx={{ mr: 1 }} />
          Filter Molecules
        </Typography>
        <IconButton
          onClick={toggleAdvancedFilters}
          aria-label={showAdvancedFilters ? "Hide advanced filters" : "Show advanced filters"}
        >
          {showAdvancedFilters ? <ExpandLess /> : <ExpandMore />}
        </IconButton>
      </FilterHeader>
      
      {/* Basic SMILES pattern filter */}
      <Box mb={2}>
        <Input
          label="SMILES Pattern"
          placeholder="Enter SMILES or substructure pattern"
          value={filter.smiles_pattern || ''}
          onChange={handleSmilesChange}
          fullWidth
          helperText="Filter molecules containing this SMILES pattern"
        />
      </Box>
      
      {/* Advanced filters section */}
      <Collapse in={showAdvancedFilters}>
        <Divider sx={{ my: 2 }} />
        
        {/* Property filters */}
        <Box mb={2}>
          <Typography variant="subtitle2" gutterBottom>
            Property Filters
          </Typography>
          
          {/* Render each property filter */}
          {(filter.property_filters || []).map((propertyFilter, index) => (
            <PropertyFilterItem key={index}>
              {/* Property name dropdown */}
              <Dropdown
                id={`property-name-${index}`}
                name={`property-name-${index}`}
                label="Property"
                value={propertyFilter.name}
                options={propertyOptions}
                onChange={(e) => handlePropertyFilterChange(index, 'name', e.target.value)}
                size="small"
                fullWidth={false}
                style={{ minWidth: '150px' }}
              />
              
              {/* Operator dropdown */}
              <Dropdown
                id={`property-operator-${index}`}
                name={`property-operator-${index}`}
                label="Operator"
                value={propertyFilter.operator}
                options={operatorOptions}
                onChange={(e) => handlePropertyFilterChange(index, 'operator', e.target.value)}
                size="small"
                fullWidth={false}
                style={{ minWidth: '150px' }}
              />
              
              {/* Value input (or min/max for "between" operator) */}
              {propertyFilter.operator === FilterOperator.BETWEEN ? (
                <>
                  <Input
                    label="Min"
                    type="number"
                    value={propertyFilter.min_value !== undefined ? propertyFilter.min_value : ''}
                    onChange={(e) => {
                      const value = e.target.value === '' ? undefined : Number(e.target.value);
                      handlePropertyFilterChange(index, 'min_value', value);
                    }}
                    size="small"
                    style={{ width: '120px' }}
                  />
                  <Input
                    label="Max"
                    type="number"
                    value={propertyFilter.max_value !== undefined ? propertyFilter.max_value : ''}
                    onChange={(e) => {
                      const value = e.target.value === '' ? undefined : Number(e.target.value);
                      handlePropertyFilterChange(index, 'max_value', value);
                    }}
                    size="small"
                    style={{ width: '120px' }}
                  />
                </>
              ) : (
                <Input
                  label="Value"
                  type={
                    propertyFilter.operator === FilterOperator.CONTAINS 
                      ? 'text' 
                      : 'number'
                  }
                  value={propertyFilter.value !== undefined ? propertyFilter.value : ''}
                  onChange={(e) => {
                    const value = e.target.value;
                    handlePropertyFilterChange(
                      index, 
                      'value', 
                      propertyFilter.operator === FilterOperator.CONTAINS ? value : 
                        value === '' ? undefined : Number(value)
                    );
                  }}
                  size="small"
                  style={{ width: '150px' }}
                />
              )}
              
              {/* Remove button */}
              <IconButton 
                size="small" 
                aria-label="Remove filter"
                onClick={() => handleRemovePropertyFilter(index)}
              >
                <Remove />
              </IconButton>
            </PropertyFilterItem>
          ))}
          
          {/* Add property filter button */}
          <Button
            variant="outlined"
            size="small"
            startIcon={<Add />}
            onClick={handleAddPropertyFilter}
            sx={{ mt: 1 }}
          >
            Add Property Filter
          </Button>
        </Box>
        
        <Grid container spacing={2}>
          {/* Flag status filter */}
          <Grid item xs={12} sm={6}>
            <Dropdown
              id="flag-status"
              name="flag-status"
              label="Flag Status"
              value={filter.flag_status || ''}
              options={flagStatusOptions}
              onChange={handleFlagStatusChange}
              fullWidth
            />
          </Grid>
          
          {/* Sort by filter */}
          <Grid item xs={12} sm={6}>
            <Box display="flex" alignItems="center" gap={1}>
              <Dropdown
                id="sort-by"
                name="sort-by"
                label="Sort By"
                value={filter.sort_by || 'created_at'}
                options={sortOptions}
                onChange={handleSortChange}
                fullWidth
              />
              <FormControlLabel
                control={
                  <Switch
                    checked={filter.sort_desc !== false}
                    onChange={handleSortDirectionChange}
                  />
                }
                label="Descending"
              />
            </Box>
          </Grid>
          
          {/* Date range filters */}
          <Grid item xs={12} sm={6}>
            <DatePicker
              label="Created After"
              value={filter.created_after ? new Date(filter.created_after) : null}
              onChange={(date) => handleDateChange('created_after', date)}
              slotProps={{ textField: { fullWidth: true } }}
            />
          </Grid>
          <Grid item xs={12} sm={6}>
            <DatePicker
              label="Created Before"
              value={filter.created_before ? new Date(filter.created_before) : null}
              onChange={(date) => handleDateChange('created_before', date)}
              slotProps={{ textField: { fullWidth: true } }}
            />
          </Grid>
        </Grid>
      </Collapse>
      
      {/* Filter actions */}
      <FilterActions>
        <Box>
          {(!!filter.smiles_pattern || 
            (filter.property_filters && filter.property_filters.length > 0) ||
            filter.flag_status ||
            filter.created_after ||
            filter.created_before) && (
            <Button
              variant="outlined"
              size="small"
              color="secondary"
              startIcon={<Clear />}
              onClick={handleClearFilters}
            >
              Clear Filters
            </Button>
          )}
        </Box>
        
        <Box>
          {/* Additional action buttons would go here if needed */}
        </Box>
      </FilterActions>
    </FilterContainer>
  );
};

export default MoleculeFilter;