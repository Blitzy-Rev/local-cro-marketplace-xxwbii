import React, { useState } from 'react';
import { TextField, TextFieldProps, InputAdornment, IconButton, InputProps } from '@mui/material';
import { styled } from '@mui/material/styles';
import Visibility from '@mui/icons-material/Visibility'; // @mui/icons-material v5.11.16
import VisibilityOff from '@mui/icons-material/VisibilityOff'; // @mui/icons-material v5.11.16
import theme from '../../theme';

// Props interface for the Input component extending Material-UI TextFieldProps
interface CustomInputProps {
  label?: string;
  name?: string;
  value?: string;
  placeholder?: string;
  type?: 'text' | 'password' | 'email' | 'number' | 'tel' | 'url' | 'search';
  helperText?: string;
  error?: string | boolean;
  fullWidth?: boolean;
  required?: boolean;
  disabled?: boolean;
  multiline?: boolean;
  rows?: number;
  maxRows?: number;
  variant?: 'outlined' | 'filled' | 'standard';
  size?: 'small' | 'medium';
  margin?: 'none' | 'dense' | 'normal';
  startAdornment?: React.ReactNode;
  endAdornment?: React.ReactNode;
  onChange?: (e: React.ChangeEvent<HTMLInputElement>) => void;
  onBlur?: () => void;
  onFocus?: () => void;
  inputRef?: React.RefObject<HTMLInputElement>;
  className?: string;
  style?: React.CSSProperties;
  InputProps?: InputProps;
}

// Styled TextField component with custom styling based on theme
const StyledTextField = styled(TextField)({
  '& .MuiOutlinedInput-root': {
    borderRadius: '8px',
    '&:hover .MuiOutlinedInput-notchedOutline': {
      borderColor: theme.palette.primary.main,
    },
    '&.Mui-focused .MuiOutlinedInput-notchedOutline': {
      borderWidth: '2px',
    },
  },
  '& .MuiInputLabel-root': {
    color: theme.palette.text.secondary,
    '&.Mui-focused': {
      color: theme.palette.primary.main,
    },
  },
  '& .MuiFormHelperText-root': {
    marginLeft: '0',
    '&.Mui-error': {
      color: theme.palette.error.main,
    },
  },
});

/**
 * A customizable input component extending Material-UI TextField with additional styling options and features
 */
const Input: React.FC<CustomInputProps> = ({
  label,
  name,
  value,
  placeholder,
  type = 'text',
  helperText,
  error,
  fullWidth = true,
  required = false,
  disabled = false,
  multiline = false,
  rows,
  maxRows,
  variant = 'outlined',
  size = 'medium',
  margin = 'normal',
  startAdornment,
  endAdornment,
  onChange,
  onBlur,
  onFocus,
  inputRef,
  className,
  style,
  InputProps = {},
  ...rest
}) => {
  // State for password visibility toggle
  const [showPassword, setShowPassword] = useState(false);

  // Toggles password visibility for password inputs
  const handleTogglePasswordVisibility = () => {
    setShowPassword(!showPassword);
  };

  // Clone InputProps to avoid mutating props
  const mergedInputProps = { ...InputProps };

  // Add start adornment if provided
  if (startAdornment) {
    mergedInputProps.startAdornment = (
      <InputAdornment position="start">{startAdornment}</InputAdornment>
    );
  }

  // Add end adornment with conditional password toggle
  if (type === 'password' || endAdornment) {
    mergedInputProps.endAdornment = (
      <InputAdornment position="end">
        {endAdornment}
        {type === 'password' && (
          <IconButton
            aria-label="toggle password visibility"
            onClick={handleTogglePasswordVisibility}
            onMouseDown={(e) => e.preventDefault()} // Prevent blur
            edge="end"
          >
            {showPassword ? <VisibilityOff /> : <Visibility />}
          </IconButton>
        )}
      </InputAdornment>
    );
  }

  // Determine error state and message
  const hasError = Boolean(error);
  const errorMessage = typeof error === 'string' ? error : helperText;

  // Handle onBlur and onFocus to match our interface
  const handleBlur = (e: React.FocusEvent<HTMLInputElement>) => {
    if (onBlur) onBlur();
  };

  const handleFocus = (e: React.FocusEvent<HTMLInputElement>) => {
    if (onFocus) onFocus();
  };

  return (
    <StyledTextField
      label={label}
      name={name}
      value={value}
      placeholder={placeholder}
      type={type === 'password' ? (showPassword ? 'text' : 'password') : type}
      helperText={errorMessage}
      error={hasError}
      fullWidth={fullWidth}
      required={required}
      disabled={disabled}
      multiline={multiline}
      rows={rows}
      maxRows={maxRows}
      variant={variant}
      size={size}
      margin={margin}
      onChange={onChange}
      onBlur={handleBlur}
      onFocus={handleFocus}
      InputProps={mergedInputProps}
      inputRef={inputRef}
      className={className}
      style={style}
      {...rest}
    />
  );
};

export default Input;