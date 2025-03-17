/**
 * Utility functions for formatting data in the Molecular Data Management and CRO Integration Platform frontend.
 * Provides reusable formatting logic for dates, currency, file sizes, SMILES strings,
 * experiment statuses, and other data types displayed throughout the application.
 */

import { format, formatDistance } from 'date-fns'; // date-fns v2.30.0
import { ExperimentStatus } from '../types/experiment';
import { SubmissionStatus } from '../types/submission';
import { ResultStatus } from '../types/result';
import { FlagStatus } from '../types/molecule';

/**
 * Formats a date string or timestamp into a human-readable date format
 * @param date - The date to format (string, number, or Date object)
 * @param formatString - Optional custom format string (defaults to 'MMM d, yyyy')
 * @returns Formatted date string
 */
export const formatDate = (
  date: string | number | Date,
  formatString: string = 'MMM d, yyyy'
): string => {
  if (!date) return '';
  
  try {
    const dateObj = typeof date === 'string' || typeof date === 'number' 
      ? new Date(date) 
      : date;
    
    if (isNaN(dateObj.getTime())) {
      return '';
    }
    
    return format(dateObj, formatString);
  } catch (error) {
    console.error('Error formatting date:', error);
    return '';
  }
};

/**
 * Formats a date string or timestamp into a human-readable date and time format
 * @param date - The date to format (string, number, or Date object)
 * @returns Formatted date and time string
 */
export const formatDateTime = (date: string | number | Date): string => {
  return formatDate(date, 'MMM d, yyyy h:mm a');
};

/**
 * Formats a date as a relative time string (e.g., '2 days ago')
 * @param date - The date to format (string, number, or Date object)
 * @returns Relative time string
 */
export const formatRelativeTime = (date: string | number | Date): string => {
  if (!date) return '';
  
  try {
    const dateObj = typeof date === 'string' || typeof date === 'number' 
      ? new Date(date) 
      : date;
    
    if (isNaN(dateObj.getTime())) {
      return '';
    }
    
    return formatDistance(dateObj, new Date(), { addSuffix: true });
  } catch (error) {
    console.error('Error formatting relative time:', error);
    return '';
  }
};

/**
 * Formats a number as a currency string with dollar sign and decimal places
 * @param amount - The amount to format
 * @param currencyCode - The currency code (defaults to USD)
 * @returns Formatted currency string
 */
export const formatCurrency = (amount: number, currencyCode: string = 'USD'): string => {
  if (typeof amount !== 'number' || isNaN(amount)) {
    return '';
  }
  
  try {
    return new Intl.NumberFormat('en-US', {
      style: 'currency',
      currency: currencyCode,
    }).format(amount);
  } catch (error) {
    console.error('Error formatting currency:', error);
    return `$${amount.toFixed(2)}`;
  }
};

/**
 * Formats a file size in bytes to a human-readable format (KB, MB, GB)
 * @param bytes - The file size in bytes
 * @returns Formatted file size string
 */
export const formatFileSize = (bytes: number): string => {
  if (typeof bytes !== 'number' || isNaN(bytes) || bytes < 0) {
    return '';
  }
  
  if (bytes === 0) return '0 Bytes';
  
  const k = 1024;
  const sizes = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
  const i = Math.floor(Math.log(bytes) / Math.log(k));
  
  // Format with appropriate precision (more precision for smaller values)
  const value = bytes / Math.pow(k, i);
  const precision = i > 0 ? (value < 10 ? 2 : 1) : 0;
  
  return `${value.toFixed(precision)} ${sizes[i]}`;
};

/**
 * Formats a decimal number as a percentage string
 * @param value - The decimal value to format as percentage (e.g., 0.75 for 75%)
 * @param decimalPlaces - Number of decimal places (defaults to 1)
 * @returns Formatted percentage string
 */
export const formatPercentage = (value: number, decimalPlaces: number = 1): string => {
  if (typeof value !== 'number' || isNaN(value)) {
    return '';
  }
  
  // Convert to percentage and round to specified decimal places
  const percentage = value * 100;
  return `${percentage.toFixed(decimalPlaces)}%`;
};

/**
 * Formats an experiment status enum value into a human-readable label
 * @param status - The experiment status enum value
 * @returns Human-readable status label
 */
export const formatExperimentStatus = (status: ExperimentStatus): string => {
  const statusLabels: Record<ExperimentStatus, string> = {
    [ExperimentStatus.DRAFT]: 'Draft',
    [ExperimentStatus.QUEUED]: 'Queued',
    [ExperimentStatus.SUBMITTED]: 'Submitted',
    [ExperimentStatus.REJECTED]: 'Rejected',
    [ExperimentStatus.QUOTE_PENDING]: 'Quote Pending',
    [ExperimentStatus.QUOTE_REJECTED]: 'Quote Rejected',
    [ExperimentStatus.IN_PROGRESS]: 'In Progress',
    [ExperimentStatus.RESULTS_PENDING]: 'Results Pending',
    [ExperimentStatus.RESULTS_AVAILABLE]: 'Results Available',
    [ExperimentStatus.RESULTS_REJECTED]: 'Results Rejected',
    [ExperimentStatus.COMPLETED]: 'Completed',
    [ExperimentStatus.CANCELLED]: 'Cancelled'
  };
  
  return statusLabels[status] || 'Unknown';
};

/**
 * Formats a submission status enum value into a human-readable label
 * @param status - The submission status enum value
 * @returns Human-readable status label
 */
export const formatSubmissionStatus = (status: SubmissionStatus): string => {
  const statusLabels: Record<SubmissionStatus, string> = {
    [SubmissionStatus.PENDING]: 'Pending',
    [SubmissionStatus.REJECTED]: 'Rejected',
    [SubmissionStatus.QUOTE_PROVIDED]: 'Quote Provided',
    [SubmissionStatus.QUOTE_REJECTED]: 'Quote Rejected',
    [SubmissionStatus.APPROVED]: 'Approved',
    [SubmissionStatus.IN_PROGRESS]: 'In Progress',
    [SubmissionStatus.COMPLETED]: 'Completed',
    [SubmissionStatus.CANCELLED]: 'Cancelled'
  };
  
  return statusLabels[status] || 'Unknown';
};

/**
 * Formats a result status enum value into a human-readable label
 * @param status - The result status enum value
 * @returns Human-readable status label
 */
export const formatResultStatus = (status: ResultStatus): string => {
  const statusLabels: Record<ResultStatus, string> = {
    [ResultStatus.PENDING]: 'Pending',
    [ResultStatus.UPLOADED]: 'Uploaded',
    [ResultStatus.APPROVED]: 'Approved',
    [ResultStatus.REJECTED]: 'Rejected'
  };
  
  return statusLabels[status] || 'Unknown';
};

/**
 * Formats a molecule flag status enum value into a human-readable label
 * @param status - The flag status enum value
 * @returns Human-readable flag label
 */
export const formatFlagStatus = (status: FlagStatus | null): string => {
  if (!status) return 'None';
  
  const statusLabels: Record<FlagStatus, string> = {
    [FlagStatus.IMPORTANT]: 'Important',
    [FlagStatus.FAVORITE]: 'Favorite',
    [FlagStatus.REVIEW]: 'Review',
    [FlagStatus.ARCHIVED]: 'Archived'
  };
  
  return statusLabels[status] || 'None';
};

/**
 * Truncates a SMILES string to a specified length with ellipsis if needed
 * @param smiles - The SMILES string to truncate
 * @param maxLength - Maximum length before truncation (defaults to 30)
 * @returns Truncated SMILES string
 */
export const truncateSMILES = (smiles: string, maxLength: number = 30): string => {
  if (!smiles) return '';
  
  if (smiles.length <= maxLength) {
    return smiles;
  }
  
  return `${smiles.substring(0, maxLength - 3)}...`;
};

/**
 * Formats a molecular property value with appropriate precision and unit
 * @param value - The property value (number or numeric string)
 * @param propertyName - The name of the property (determines formatting precision)
 * @param unit - Optional unit to append to the formatted value
 * @returns Formatted property value with unit
 */
export const formatPropertyValue = (
  value: number | string,
  propertyName: string,
  unit?: string
): string => {
  if (value === null || value === undefined || value === '') {
    return '';
  }
  
  // Convert string to number if it's a numeric string
  const numValue = typeof value === 'string' ? parseFloat(value) : value;
  
  if (typeof numValue !== 'number' || isNaN(numValue)) {
    return String(value);
  }
  
  // Determine appropriate precision based on property name
  let precision = 2; // Default precision
  
  // Common property name variations (case-insensitive)
  const propertyLower = propertyName.toLowerCase();
  
  if (propertyLower.includes('mw') || propertyLower.includes('molecular weight')) {
    precision = 2;
  } else if (propertyLower.includes('logp') || propertyLower.includes('clogp')) {
    precision = 2;
  } else if (propertyLower.includes('psa') || propertyLower.includes('polar surface area')) {
    precision = 1;
  } else if (propertyLower.includes('hba') || propertyLower.includes('hbd') || 
             propertyLower.includes('rotatable bonds') || propertyLower.includes('rotbonds')) {
    precision = 0; // Integer values
  } else if (propertyLower.includes('solubility')) {
    precision = 2;
  } else if (propertyLower.includes('ic50') || propertyLower.includes('ec50') || propertyLower.includes('ki')) {
    precision = 3; // More precision for binding constants
  }
  
  // Format the number with appropriate precision
  const formattedValue = numValue.toFixed(precision);
  
  // Remove trailing zeros after decimal point if they exist
  const cleanValue = formattedValue.replace(/\.0+$/, '').replace(/(\.\d*[1-9])0+$/, '$1');
  
  // Append unit if provided
  return unit ? `${cleanValue} ${unit}` : cleanValue;
};

/**
 * Formats a duration in days into a human-readable format
 * @param days - Number of days
 * @returns Formatted duration string
 */
export const formatDuration = (days: number): string => {
  if (typeof days !== 'number' || isNaN(days)) {
    return '';
  }
  
  if (days === 1) {
    return '1 day';
  }
  
  return `${days} days`;
};