import React, { useState, useCallback } from 'react'; // React 18.2+
import { Box, Typography, Paper, Grid, CircularProgress } from '@mui/material'; // v5.13+
import FileUpload from '../../../components/common/FileUpload';
import Button from '../../../components/common/Button';
import Input from '../../../components/common/Input';
import StructuredDataEntry from './StructuredDataEntry';
import { useCROInterface } from '../hooks/useCROInterface';
import { useToast } from '../../../hooks/useToast';
import { ResultCreate, ResultDataCreate } from '../../../types/result';
import { Submission } from '../../../types/submission';

/**
 * Interface defining the props for the ResultUploadForm component.
 */
interface ResultUploadFormProps {
  /** The submission to upload results for. */
  submission: Submission;
  /** Callback function called when results are successfully uploaded. */
  onSuccess: () => void;
  /** Callback function called when form submission is cancelled. */
  onCancel: () => void;
}

/**
 * A form component for CRO users to upload experimental results for a completed submission.
 * It provides file upload functionality, structured data entry for molecules, and submission of results to the backend.
 * @param {ResultUploadFormProps} { submission, onSuccess, onCancel } - The props for the component.
 * @returns {JSX.Element} - The rendered component.
 */
const ResultUploadForm: React.FC<ResultUploadFormProps> = ({ submission, onSuccess, onCancel }) => {
  // State variables for form fields
  const [files, setFiles] = useState<File[]>([]);
  const [notes, setNotes] = useState('');
  const [structuredData, setStructuredData] = useState<ResultDataCreate[]>([]);

  // State variables for loading status and validation errors
  const [isLoading, setIsLoading] = useState(false);
  const [errors, setErrors] = useState<{ files?: string; notes?: string; structuredData?: string }>({});

  // Extract molecules from submission experiment data
  const molecules = submission.experiment?.molecules?.map(m => m.molecule) || [];

  // Custom hook for CRO interface functionality
  const { handleUploadResult, handleAddResultData } = useCROInterface();

  // Custom hook for displaying toast notifications
  const { showToast } = useToast();

  /**
   * Custom validation for uploaded files
   * @param {File} file - The file to validate.
   * @returns {Promise<{ isValid: boolean; error?: string }>} - A promise that resolves to an object with a boolean indicating whether the file is valid and an optional error message.
   */
  const handleFileValidation = async (file: File): Promise<{ isValid: boolean; error?: string }> => {
    // Allowed result file types
    const allowedTypes = ['application/pdf', 'application/vnd.ms-excel', 'text/csv', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'];
    // Maximum file size (10MB)
    const maxSize = 10 * 1024 * 1024;

    if (!allowedTypes.includes(file.type)) {
      return { isValid: false, error: 'Invalid file type. Allowed types are PDF, Excel, and CSV.' };
    }

    if (file.size > maxSize) {
      return { isValid: false, error: 'File size exceeds the maximum allowed size of 10MB.' };
    }

    return { isValid: true };
  };

  /**
   * Handler for when files are selected
   * @param {File[]} selectedFiles - The selected files.
   * @returns {void}
   */
  const handleFilesSelected = (selectedFiles: File[]): void => {
    setFiles(selectedFiles);
    setErrors(prev => ({ ...prev, files: undefined }));
  };

  /**
   * Handler for when structured data changes
   * @param {ResultDataCreate[]} data - The structured data.
   * @returns {void}
   */
  const handleStructuredDataChange = (data: ResultDataCreate[]): void => {
    setStructuredData(data);
    setErrors(prev => ({ ...prev, structuredData: undefined }));
  };

  /**
   * Handler for when notes field changes
   * @param {React.ChangeEvent<HTMLInputElement>} e - The change event.
   * @returns {void}
   */
  const handleNotesChange = (e: React.ChangeEvent<HTMLInputElement>): void => {
    setNotes(e.target.value);
    setErrors(prev => ({ ...prev, notes: undefined }));
  };

  /**
   * Validates the form fields before submission
   * @returns {boolean} - Whether the form is valid.
   */
  const validateForm = (): boolean => {
    const newErrors: { files?: string; notes?: string; structuredData?: string } = {};

    if (files.length === 0) {
      newErrors.files = 'At least one result file is required';
    }

    if (!notes) {
      newErrors.notes = 'Notes are required';
    }

    if (structuredData.length === 0 && molecules.length > 0) {
      newErrors.structuredData = 'Structured data is required for at least one molecule';
    }

    setErrors(newErrors);
    return Object.keys(newErrors).length === 0;
  };

  /**
   * Handler for form submission
   * @param {React.FormEvent} e - The form event.
   * @returns {Promise<void>}
   */
  const handleSubmit = async (e: React.FormEvent): Promise<void> => {
    e.preventDefault();

    if (!validateForm()) {
      return;
    }

    setIsLoading(true);

    try {
      // Create result data object
      const resultData: ResultCreate = {
        submission_id: submission.id,
        notes: notes,
        files: files.map(file => ({
          file: file,
          file_name: file.name,
          file_type: file.type
        }))
      };

      // Call handleUploadResult from useCROInterface hook
      const result = await handleUploadResult(resultData);

      // Add structured data for each entry
      if (structuredData.length > 0) {
        for (const data of structuredData) {
          await handleAddResultData(result.id, data);
        }
      }

      // Show success toast notification
      showToast({ type: 'success', message: 'Results uploaded successfully!' });

      // Call onSuccess callback
      onSuccess();
    } catch (error: any) {
      // Show error toast notification
      showToast({ type: 'error', message: error.message || 'Failed to upload results' });
      console.error("Result upload failed:", error);
    } finally {
      setIsLoading(false);
    }
  };

  return (
    <Paper elevation={1} sx={{ p: 2, mb: 3 }}>
      <Typography variant="h6" gutterBottom>
        Upload Results
      </Typography>
      <form onSubmit={handleSubmit}>
        <FileUpload
          label="Result Files"
          buttonText="Select Files"
          accept={['application/pdf', 'application/vnd.ms-excel', 'text/csv', 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet']}
          multiple
          onFilesSelected={handleFilesSelected}
          onFileValidationFail={(error) => console.error("File validation failed:", error)}
          customValidation={handleFileValidation}
          disabled={isLoading}
        />
        {errors.files && (
          <Typography variant="body2" color="error" sx={{ mt: 1 }}>
            {errors.files}
          </Typography>
        )}
        <Input
          label="Notes"
          multiline
          rows={3}
          fullWidth
          placeholder="Additional notes about the results"
          value={notes}
          onChange={handleNotesChange}
          error={errors.notes}
          disabled={isLoading}
          sx={{ mt: 2 }}
        />
        <StructuredDataEntry
          molecules={molecules}
          value={structuredData}
          onChange={handleStructuredDataChange}
          experimentType={submission.experiment?.type?.name || 'Unknown'}
          error={errors.structuredData}
        />
        <Box sx={{ display: 'flex', justifyContent: 'flex-end', mt: 3 }}>
          <Button onClick={onCancel} disabled={isLoading} sx={{ mr: 2 }}>
            Cancel
          </Button>
          <Button type="submit" variant="contained" color="primary" disabled={isLoading} loading={isLoading}>
            Submit Results
          </Button>
        </Box>
      </form>
    </Paper>
  );
};

export default ResultUploadForm;