import React from 'react';
import {
  Select,
  MenuItem,
  FormControl,
  InputLabel,
  FormHelperText,
  Box,
  Chip,
  OutlinedInput,
  SelectChangeEvent
} from '@mui/material';
import { styled } from '@mui/material/styles';
import theme from '../../theme';

/**
 * Interface for dropdown option items
 */
interface DropdownOption {
  value: string | number;
  label: string;
  disabled?: boolean;
  icon?: React.ReactNode;
}

/**
 * Props interface for the Dropdown component extending Material-UI SelectProps
 */
interface DropdownProps {
  id: string;
  name: string;
  label: string;
  value: any;
  onChange: (event: SelectChangeEvent<any>, child: React.ReactNode) => void;
  required?: boolean;
  disabled?: boolean;
  fullWidth?: boolean;
  error?: string;
  helperText?: string;
  variant?: 'outlined' | 'filled' | 'standard';
  size?: 'small' | 'medium';
  className?: string;
  style?: React.CSSProperties;
  multiple?: boolean;
  options: DropdownOption[];
  placeholder?: string;
  clearable?: boolean;
  onClear?: () => void;
}

// Styled components
const StyledFormControl = styled(FormControl)(({ fullWidth }) => ({
  marginBottom: '16px',
  width: fullWidth ? '100%' : 'auto',
}));

const StyledSelect = styled(Select)(() => ({
  '& .MuiOutlinedInput-notchedOutline': {
    borderRadius: '8px',
  },
  '&.Mui-focused .MuiOutlinedInput-notchedOutline': {
    borderWidth: '1px',
    borderColor: theme.palette.primary.main,
  },
}));

const StyledMenuItem = styled(MenuItem)(({ size }) => ({
  display: 'flex',
  alignItems: 'center',
  gap: '8px',
  minHeight: size === 'small' ? '32px' : '40px',
}));

/**
 * A customizable dropdown component extending Material-UI Select with consistent styling
 * and additional features. Used throughout the application for dropdown selections,
 * filters, and form fields.
 */
const Dropdown: React.FC<DropdownProps> = ({
  id,
  name,
  label,
  value,
  onChange,
  required = false,
  disabled = false,
  fullWidth = true,
  error,
  helperText,
  variant = 'outlined',
  size = 'medium',
  className,
  style,
  multiple = false,
  options = [],
  placeholder,
  clearable = false,
  onClear,
}) => {
  // Handle empty value display with placeholder
  const displayEmpty = Boolean(placeholder);
  
  // Check if value is empty (considering both single and multiple cases)
  const isValueEmpty = multiple 
    ? !value || (Array.isArray(value) && value.length === 0) 
    : (value === undefined || value === null || value === '');
  
  // Render the value (for multiple select with chips)
  const renderValue = (selected: any) => {
    // Show placeholder for empty values
    if (isValueEmpty && displayEmpty) {
      return <em style={{ opacity: 0.6 }}>{placeholder}</em>;
    }
    
    // For multiple selection, render chips
    if (multiple && Array.isArray(selected)) {
      return (
        <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 0.5 }}>
          {selected.map((val) => {
            const option = options.find((opt) => opt.value === val);
            return (
              <Chip 
                key={val} 
                label={option?.label || val} 
                size={size === 'small' ? 'small' : 'medium'}
                onDelete={clearable ? () => {
                  const newValue = selected.filter(v => v !== val);
                  onChange(
                    { target: { value: newValue } } as SelectChangeEvent<any>,
                    null
                  );
                } : undefined}
              />
            );
          })}
        </Box>
      );
    }
    
    // For single selection, show the label
    const option = options.find((opt) => opt.value === selected);
    return option?.label || selected;
  };

  // Handle clearing the selection
  const handleClear = (e: React.MouseEvent) => {
    e.stopPropagation();
    if (multiple) {
      onChange({ target: { value: [] } } as SelectChangeEvent<any>, null);
    } else {
      onChange({ target: { value: '' } } as SelectChangeEvent<any>, null);
    }
    onClear?.();
  };

  // Handle keyboard interaction for the clear button
  const handleClearKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' || e.key === ' ') {
      e.preventDefault();
      if (multiple) {
        onChange({ target: { value: [] } } as SelectChangeEvent<any>, null);
      } else {
        onChange({ target: { value: '' } } as SelectChangeEvent<any>, null);
      }
      onClear?.();
    }
  };

  // Prepare the value for the Select component
  const selectValue = isValueEmpty 
    ? (multiple ? [] : '') 
    : value;

  return (
    <StyledFormControl 
      variant={variant}
      size={size}
      fullWidth={fullWidth}
      error={Boolean(error)}
      disabled={disabled}
      className={className}
      style={style}
    >
      {label && (
        <InputLabel id={`${id}-label`} required={required}>
          {label}
        </InputLabel>
      )}
      
      <StyledSelect
        labelId={`${id}-label`}
        id={id}
        name={name}
        value={selectValue}
        onChange={onChange}
        input={<OutlinedInput label={label} />}
        multiple={multiple}
        displayEmpty={displayEmpty}
        renderValue={renderValue}
        size={size}
        MenuProps={{
          PaperProps: {
            style: {
              maxHeight: 300,
            },
          },
          // Improved menu positioning
          anchorOrigin: {
            vertical: 'bottom',
            horizontal: 'left',
          },
          transformOrigin: {
            vertical: 'top',
            horizontal: 'left',
          },
          // Better keyboard navigation
          disableRestoreFocus: true,
        }}
        // Add clear button to the end if clearable
        endAdornment={
          clearable && !isValueEmpty ? (
            <Box 
              component="div" 
              sx={{ 
                position: 'absolute',
                right: 32, // Position before the dropdown arrow
                top: '50%',
                transform: 'translateY(-50%)',
                display: 'flex',
                alignItems: 'center',
                justifyContent: 'center',
                width: '16px',
                height: '16px',
                borderRadius: '50%',
                cursor: 'pointer',
                color: theme.palette.text.secondary,
                '&:hover': {
                  color: theme.palette.text.primary,
                  backgroundColor: theme.palette.action.hover,
                },
                '&:focus': {
                  outline: '2px solid',
                  outlineColor: theme.palette.primary.main,
                  outlineOffset: '2px',
                },
                zIndex: 1,
                fontSize: '14px',
              }}
              onClick={handleClear}
              onKeyDown={handleClearKeyDown}
              role="button"
              aria-label="Clear selection"
              tabIndex={0}
            >
              ✕
            </Box>
          ) : null
        }
        // Accessibility support
        aria-describedby={error ? `${id}-error` : helperText ? `${id}-helper` : undefined}
      >
        {options.map((option) => (
          <StyledMenuItem 
            key={option.value} 
            value={option.value} 
            disabled={option.disabled}
            size={size}
          >
            {option.icon && (
              <Box component="span" sx={{ display: 'flex', alignItems: 'center', marginRight: '8px' }}>
                {option.icon}
              </Box>
            )}
            {option.label}
          </StyledMenuItem>
        ))}
      </StyledSelect>
      
      {(error || helperText) && (
        <FormHelperText id={error ? `${id}-error` : `${id}-helper`}>
          {error || helperText}
        </FormHelperText>
      )}
    </StyledFormControl>
  );
};

export default Dropdown;