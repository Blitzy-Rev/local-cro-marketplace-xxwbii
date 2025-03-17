import { UserRole } from '../types/user';
import { ExperimentStatus } from '../types/experiment';

// Validation patterns
export const EMAIL_REGEX = /^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$/;
export const PASSWORD_REGEX = /^(?=.*[a-z])(?=.*[A-Z])(?=.*\d)(?=.*[@$!%*?&])[A-Za-z\d@$!%*?&]{10,}$/;
export const SMILES_REGEX = /^[A-Za-z0-9@+\-\[\]\(\)\\/#$.%=~&:,*]+$/;

// Configuration constants
export const MAX_CSV_SIZE_MB = 10;
export const REQUIRED_CSV_COLUMNS = ['SMILES'];
export const PROPERTY_RANGES = {
  molecular_weight: { min: 0, max: 2000 },
  logp: { min: -10, max: 10 },
  solubility: { min: 0, max: 1000 },
  activity: { min: 0, max: 100 }
};

/**
 * Custom error class for validation errors
 */
export class ValidationError extends Error {
  public errors: Record<string, string>;

  /**
   * Creates a new ValidationError instance
   * @param message - Error message
   * @param errors - Object containing field-specific error messages
   */
  constructor(message: string, errors: Record<string, string> = {}) {
    super(message);
    this.name = 'ValidationError';
    this.errors = errors;
    
    // This is needed to properly capture the stack trace in TypeScript
    if (Error.captureStackTrace) {
      Error.captureStackTrace(this, ValidationError);
    }
  }
}

/**
 * Validates a SMILES string using regex pattern and basic structure checks
 * @param smiles - SMILES string to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validateSMILES(smiles: string): { isValid: boolean; error: string | null } {
  if (!smiles || smiles.trim() === '') {
    return { isValid: false, error: 'SMILES string is required' };
  }

  if (!SMILES_REGEX.test(smiles)) {
    return { isValid: false, error: 'Invalid SMILES format' };
  }

  // Additional structure validation
  // Check for balanced parentheses and brackets
  const openParenCount = (smiles.match(/\(/g) || []).length;
  const closeParenCount = (smiles.match(/\)/g) || []).length;
  const openBracketCount = (smiles.match(/\[/g) || []).length;
  const closeBracketCount = (smiles.match(/\]/g) || []).length;

  if (openParenCount !== closeParenCount) {
    return { isValid: false, error: 'Unbalanced parentheses in SMILES string' };
  }

  if (openBracketCount !== closeBracketCount) {
    return { isValid: false, error: 'Unbalanced brackets in SMILES string' };
  }

  return { isValid: true, error: null };
}

/**
 * Validates an email address format
 * @param email - Email address to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validateEmail(email: string): { isValid: boolean; error: string | null } {
  if (!email || email.trim() === '') {
    return { isValid: false, error: 'Email is required' };
  }

  if (!EMAIL_REGEX.test(email)) {
    return { isValid: false, error: 'Invalid email format' };
  }

  return { isValid: true, error: null };
}

/**
 * Validates a password against security requirements
 * @param password - Password to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validatePassword(password: string): { isValid: boolean; error: string | null } {
  if (!password) {
    return { isValid: false, error: 'Password is required' };
  }

  if (!PASSWORD_REGEX.test(password)) {
    return { 
      isValid: false, 
      error: 'Password must be at least 10 characters and include uppercase, lowercase, number, and special character' 
    };
  }

  return { isValid: true, error: null };
}

/**
 * Validates that password and confirm password match
 * @param password - Password value
 * @param confirmPassword - Confirm password value
 * @returns Validation result with isValid flag and optional error message
 */
export function validatePasswordMatch(
  password: string, 
  confirmPassword: string
): { isValid: boolean; error: string | null } {
  if (!password || !confirmPassword) {
    return { isValid: false, error: 'Both password and confirm password are required' };
  }

  if (password !== confirmPassword) {
    return { isValid: false, error: 'Passwords do not match' };
  }

  return { isValid: true, error: null };
}

/**
 * Validates that a user role is valid
 * @param role - Role value to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validateUserRole(role: string): { isValid: boolean; error: string | null } {
  if (!role) {
    return { isValid: false, error: 'Role is required' };
  }

  const validRoles = Object.values(UserRole);
  if (!validRoles.includes(role as UserRole)) {
    return { isValid: false, error: 'Invalid user role' };
  }

  return { isValid: true, error: null };
}

/**
 * Validates user registration data
 * @param data - Registration form data
 * @returns Validation result with isValid flag and field-specific errors
 */
export function validateRegistrationData(data: {
  email: string;
  password: string;
  confirm_password: string;
  role: string;
}): { isValid: boolean; errors: Record<string, string> } {
  const errors: Record<string, string> = {};
  let isValid = true;

  // Validate email
  const emailValidation = validateEmail(data.email);
  if (!emailValidation.isValid) {
    errors.email = emailValidation.error || 'Invalid email';
    isValid = false;
  }

  // Validate password
  const passwordValidation = validatePassword(data.password);
  if (!passwordValidation.isValid) {
    errors.password = passwordValidation.error || 'Invalid password';
    isValid = false;
  }

  // Validate password match
  const passwordMatchValidation = validatePasswordMatch(data.password, data.confirm_password);
  if (!passwordMatchValidation.isValid) {
    errors.confirm_password = passwordMatchValidation.error || 'Passwords do not match';
    isValid = false;
  }

  // Validate role
  const roleValidation = validateUserRole(data.role);
  if (!roleValidation.isValid) {
    errors.role = roleValidation.error || 'Invalid role';
    isValid = false;
  }

  return { isValid, errors };
}

/**
 * Validates user login credentials
 * @param credentials - Login credentials
 * @returns Validation result with isValid flag and field-specific errors
 */
export function validateLoginCredentials(credentials: {
  email: string;
  password: string;
}): { isValid: boolean; errors: Record<string, string> } {
  const errors: Record<string, string> = {};
  let isValid = true;

  // Validate email
  const emailValidation = validateEmail(credentials.email);
  if (!emailValidation.isValid) {
    errors.email = emailValidation.error || 'Invalid email';
    isValid = false;
  }

  // Check if password is provided
  if (!credentials.password) {
    errors.password = 'Password is required';
    isValid = false;
  }

  return { isValid, errors };
}

/**
 * Validates a molecular property value against expected ranges
 * @param propertyName - Name of the property
 * @param propertyValue - Value to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validateMolecularProperty(
  propertyName: string,
  propertyValue: any
): { isValid: boolean; error: string | null } {
  // Get property range definition if it exists
  const propertyRange = (PROPERTY_RANGES as Record<string, { min: number; max: number }>)[propertyName.toLowerCase()];
  
  // If we don't have a defined range for this property, consider it valid
  if (!propertyRange) {
    return { isValid: true, error: null };
  }

  // Ensure value is a number
  const numValue = Number(propertyValue);
  if (isNaN(numValue)) {
    return { isValid: false, error: 'Property value must be a number' };
  }

  // Check if value is within the valid range
  const { min, max } = propertyRange;
  if (numValue < min || numValue > max) {
    return { isValid: false, error: `Value must be between ${min} and ${max}` };
  }

  return { isValid: true, error: null };
}

/**
 * Validates that an experiment status is valid
 * @param status - Status value to validate
 * @returns Validation result with isValid flag and optional error message
 */
export function validateExperimentStatus(status: string): { isValid: boolean; error: string | null } {
  if (!status) {
    return { isValid: false, error: 'Status is required' };
  }

  const validStatuses = Object.values(ExperimentStatus);
  if (!validStatuses.includes(status as ExperimentStatus)) {
    return { isValid: false, error: 'Invalid experiment status' };
  }

  return { isValid: true, error: null };
}

/**
 * Validates experiment creation data
 * @param data - Experiment creation form data
 * @returns Validation result with isValid flag and field-specific errors
 */
export function validateExperimentData(data: {
  name?: string;
  type_id?: string;
  molecule_ids?: string[];
  parameters?: Array<{ parameter_name: string; parameter_value: any }>;
}): { isValid: boolean; errors: Record<string, string> } {
  const errors: Record<string, string> = {};
  let isValid = true;

  // Validate experiment name
  if (!data.name || data.name.trim() === '') {
    errors.name = 'Experiment name is required';
    isValid = false;
  }

  // Validate experiment type
  if (!data.type_id) {
    errors.type_id = 'Experiment type is required';
    isValid = false;
  }

  // Validate molecule selection
  if (!data.molecule_ids || data.molecule_ids.length === 0) {
    errors.molecule_ids = 'At least one molecule must be selected';
    isValid = false;
  }

  // Validate parameters if provided
  if (data.parameters && data.parameters.length > 0) {
    const paramErrors: string[] = [];
    
    data.parameters.forEach((param, index) => {
      if (!param.parameter_name || param.parameter_name.trim() === '') {
        paramErrors.push(`Parameter #${index + 1}: Name is required`);
      }
      
      if (param.parameter_value === undefined || param.parameter_value === null) {
        paramErrors.push(`Parameter #${index + 1}: Value is required`);
      }
    });
    
    if (paramErrors.length > 0) {
      errors.parameters = paramErrors.join('; ');
      isValid = false;
    }
  }

  return { isValid, errors };
}

/**
 * Validates that CSV headers include required columns
 * @param headers - Array of CSV header names
 * @param requiredColumns - Array of required column names
 * @returns Validation result with isValid flag and array of missing columns
 */
export function validateCSVHeaders(
  headers: string[],
  requiredColumns: string[] = REQUIRED_CSV_COLUMNS
): { isValid: boolean; missingColumns: string[] } {
  const missingColumns: string[] = [];
  const headerUppercase = headers.map(h => h.toUpperCase());
  
  // Check if each required column exists in headers
  requiredColumns.forEach(col => {
    if (!headerUppercase.includes(col.toUpperCase())) {
      missingColumns.push(col);
    }
  });
  
  return { isValid: missingColumns.length === 0, missingColumns };
}

/**
 * Validates that a CSV file size is within limits
 * @param sizeInBytes - File size in bytes
 * @param maxSizeMB - Maximum allowed size in MB
 * @returns Validation result with isValid flag and optional error message
 */
export function validateCSVFileSize(
  sizeInBytes: number,
  maxSizeMB: number = MAX_CSV_SIZE_MB
): { isValid: boolean; error: string | null } {
  const maxSizeBytes = maxSizeMB * 1024 * 1024; // Convert MB to bytes
  
  if (sizeInBytes > maxSizeBytes) {
    return { 
      isValid: false, 
      error: `File size exceeds maximum allowed size of ${maxSizeMB}MB` 
    };
  }
  
  return { isValid: true, error: null };
}

/**
 * Validates that all required fields are present in an object
 * @param data - Object to validate
 * @param requiredFields - Array of required field names
 * @returns Validation result with isValid flag and array of missing fields
 */
export function validateRequiredFields(
  data: Record<string, any>,
  requiredFields: string[]
): { isValid: boolean; missingFields: string[] } {
  const missingFields: string[] = [];
  
  requiredFields.forEach(field => {
    if (data[field] === undefined || data[field] === null) {
      missingFields.push(field);
    }
  });
  
  return { isValid: missingFields.length === 0, missingFields };
}

/**
 * Validates that a string length is within specified limits
 * @param value - String to validate
 * @param minLength - Minimum allowed length
 * @param maxLength - Maximum allowed length
 * @returns Validation result with isValid flag and optional error message
 */
export function validateStringLength(
  value: string,
  minLength: number,
  maxLength: number
): { isValid: boolean; error: string | null } {
  if (typeof value !== 'string') {
    return { isValid: false, error: 'Value must be a string' };
  }
  
  const length = value.length;
  
  if (length < minLength) {
    return { isValid: false, error: `Value must be at least ${minLength} characters` };
  }
  
  if (length > maxLength) {
    return { isValid: false, error: `Value must be no more than ${maxLength} characters` };
  }
  
  return { isValid: true, error: null };
}

/**
 * Validates that a numeric value is within specified range
 * @param value - Number to validate
 * @param minValue - Minimum allowed value
 * @param maxValue - Maximum allowed value
 * @returns Validation result with isValid flag and optional error message
 */
export function validateNumericRange(
  value: any,
  minValue: number,
  maxValue: number
): { isValid: boolean; error: string | null } {
  const numValue = Number(value);
  
  if (isNaN(numValue)) {
    return { isValid: false, error: 'Value must be a number' };
  }
  
  if (numValue < minValue || numValue > maxValue) {
    return { isValid: false, error: `Value must be between ${minValue} and ${maxValue}` };
  }
  
  return { isValid: true, error: null };
}