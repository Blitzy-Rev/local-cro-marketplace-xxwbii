import React, { useState, useCallback } from 'react'; // React 18.2+
import { Box, Typography, Paper, Alert, CircularProgress } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { CloudUpload } from '@mui/icons-material'; // v5.13+
import FileUpload from '../../../components/common/FileUpload';
import ProgressBar from '../../../components/common/ProgressBar';
import { useCSVImport } from '../hooks/useCSVImport';
import { validateCSVFile } from '../../../utils/csvUtils';
import { useToast } from '../../../hooks/useToast';

/**
 * Interface defining the props for the FileUploader component
 */
interface FileUploaderProps {
  /** Callback function called when file upload is complete */
  onUploadComplete: (headers: string[], rowCount: number) => void;
  /** Additional CSS class for styling */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
}

/**
 * A styled Paper component for the upload container
 */
const UploadContainer = styled(Paper)({
  padding: '24px',
  borderRadius: '8px',
  marginBottom: '24px',
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
});

/**
 * A styled Typography component for the upload title
 */
const UploadTitle = styled(Typography)(({ theme }) => ({
  marginBottom: '16px',
  fontWeight: 500,
  color: theme.palette.text.primary,
}));

/**
 * A styled Box component for the progress container
 */
const ProgressContainer = styled(Box)({
  width: '100%',
  marginTop: '16px',
});

/**
 * A specialized file upload component for CSV files containing molecular data.
 * It handles file selection, drag-and-drop functionality, file validation, and provides feedback on the upload process.
 * This component is a key part of the CSV import workflow for molecular data ingestion.
 */
const FileUploader: React.FC<FileUploaderProps> = ({ onUploadComplete, className, style }) => {
  // Destructure state and functions from the useCSVImport hook
  const {
    isLoading,
    errors,
    handleFileUpload,
    progress,
    isUploading
  } = useCSVImport();

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  /**
   * Handles file selection from the FileUpload component
   * @param files - Array of selected files
   */
  const handleFilesSelected = useCallback(async (files: File[]) => {
    // Check if files array is empty
    if (!files || files.length === 0) {
      showToast({
        type: 'warning',
        message: 'No file selected',
      });
      return;
    }

    // Get the first file from the array
    const file = files[0];

    // Validate the CSV file
    const validationResult = await validateCSVFile(file);
    if (!validationResult.isValid) {
      showToast({
        type: 'error',
        message: validationResult.error || 'Invalid CSV file',
      });
      return;
    }

    // Call the handleFileUpload function from the useCSVImport hook
    await handleFileUpload(file);
  }, [handleFileUpload, showToast]);

  /**
   * Handles validation failures from the FileUpload component
   * @param error - Error message from file validation
   */
  const handleFileValidationFail = useCallback((error: string) => {
    showToast({
      type: 'error',
      message: error || 'Invalid CSV file',
    });
  }, [showToast]);

  /**
   * Custom validation function for CSV files
   * @param file - The CSV file to validate
   * @returns Promise<{ isValid: boolean; error?: string }> - Validation result
   */
  const customValidation = useCallback(async (file: File) => {
    // Add custom validation logic here if needed
    return { isValid: true };
  }, []);

  return (
    <UploadContainer className={className} style={style}>
      <UploadTitle>Upload CSV File</UploadTitle>
      <FileUpload
        accept={['.csv', 'text/csv']}
        onFilesSelected={handleFilesSelected}
        onFileValidationFail={handleFileValidationFail}
        customValidation={customValidation}
        multiple={false}
        disabled={isLoading || isUploading}
        helperText="Upload a CSV file containing molecular data. Ensure the file contains a 'SMILES' column."
        label="Select CSV File"
      />
      {(isLoading || isUploading) && (
        <CircularProgress
          size={24}
          sx={{
            position: 'absolute',
            top: '50%',
            left: '50%',
            marginTop: '-12px',
            marginLeft: '-12px',
          }}
        />
      )}
      {progress > 0 && (
        <ProgressContainer>
          <ProgressBar value={progress} showPercentage />
        </ProgressContainer>
      )}
      {errors.length > 0 && (
        <Alert severity="error">
          {errors.join(', ')}
        </Alert>
      )}
    </UploadContainer>
  );
};

export default FileUploader;