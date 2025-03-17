import { useState, useEffect, useCallback, useRef } from 'react'; // react ^18.2.0
import {
  uploadCSV,
  mapCSVColumns,
  getCSVProcessingStatus,
  cancelCSVProcessing,
  getAvailableProperties
} from '../../../api/csv';
import {
  getCSVPreview,
  createColumnMappingOptions,
  generateMappingFromHeaders,
  validateMappingConfiguration,
  CUSTOM_PROPERTY_PREFIX
} from '../../../utils/csvUtils';
import { useToast } from '../../../hooks/useToast';

/**
 * Interface for CSV import state
 */
interface CSVImportState {
  /** The uploaded CSV file */
  file: File | null;
  /** CSV file headers */
  headers: string[];
  /** Preview data from CSV file */
  preview: any[];
  /** Total number of rows in CSV file */
  rowCount: number;
  /** ID of uploaded file on server */
  fileId: string;
  /** ID of processing job on server */
  jobId: string;
  /** Current status of processing (pending, processing, completed, failed) */
  status: string;
  /** Processing progress percentage (0-100) */
  progress: number;
  /** Summary of processing results */
  summary: any;
  /** Current mapping of CSV columns to system properties */
  currentMapping: Array<{ csvColumn: string; systemProperty: string; isCustom: boolean }>;
  /** Available mapping options for each CSV column */
  mappingOptions: Array<{ csvColumn: string; mappingOptions: Array<{ value: string; label: string; isCustom: boolean }> }>;
  /** Available system properties for mapping */
  availableProperties: Array<{ name: string; display_name: string; data_type: string; unit?: string }>;
  /** Whether any operation is in progress */
  isLoading: boolean;
  /** Error messages from operations */
  errors: string[];
}

/**
 * Return type of useCSVImport hook
 */
interface CSVImportHook extends CSVImportState {
  /** All state properties */
  // ...CSVImportState; // This is already included by extending CSVImportState
  /** Function to handle file upload */
  handleFileUpload: (file: File) => Promise<void>;
  /** Function to update column mapping */
  handleUpdateMapping: (csvColumn: string, systemProperty: string, isCustom: boolean) => void;
  /** Function to auto-generate mapping from headers */
  handleGenerateMapping: () => void;
  /** Function to submit mapping and start processing */
  handleSubmitMapping: () => Promise<void>;
  /** Function to cancel ongoing processing */
  handleCancelProcessing: () => Promise<void>;
  /** Function to reset all state */
  handleReset: () => void;
}

/**
 * Custom hook for managing CSV import workflow
 * @returns State and functions for CSV import workflow
 */
export const useCSVImport = (): CSVImportHook => {
  // Initialize state for file, headers, preview data, and row count
  const [file, setFile] = useState<File | null>(null);
  const [headers, setHeaders] = useState<string[]>([]);
  const [preview, setPreview] = useState<any[]>([]);
  const [rowCount, setRowCount] = useState<number>(0);

  // Initialize state for fileId, jobId, status, progress, and summary
  const [fileId, setFileId] = useState<string>('');
  const [jobId, setJobId] = useState<string>('');
  const [status, setStatus] = useState<string>('');
  const [progress, setProgress] = useState<number>(0);
  const [summary, setSummary] = useState<any>(null);

  // Initialize state for mapping-related data (currentMapping, mappingOptions, availableProperties)
  const [currentMapping, setCurrentMapping] = useState<Array<{ csvColumn: string; systemProperty: string; isCustom: boolean }>>([]);
  const [mappingOptions, setMappingOptions] = useState<Array<{ csvColumn: string; mappingOptions: Array<{ value: string; label: string; isCustom: boolean }> }>>([]);
  const [availableProperties, setAvailableProperties] = useState<Array<{ name: string; display_name: string; data_type: string; unit?: string }>>([]);

  // Initialize state for loading status and errors
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [errors, setErrors] = useState<string[]>([]);

  // Create polling interval reference with useRef
  const pollingInterval = useRef<number | null>(null);

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  // Fetch available system properties on mount
  useEffect(() => {
    const fetchProperties = async () => {
      try {
        setIsLoading(true);
        const response = await getAvailableProperties();
        setAvailableProperties(response.properties);
      } catch (error: any) {
        setErrors(prevErrors => [...prevErrors, error.message || 'Failed to fetch available properties']);
        showToast({
          type: 'error',
          message: error.message || 'Failed to fetch available properties',
        });
      } finally {
        setIsLoading(false);
      }
    };

    fetchProperties();
  }, [showToast]);

  // Create handleFileUpload function to upload CSV file and get preview data
  const handleFileUpload = useCallback(async (file: File) => {
    try {
      setIsLoading(true);
      setErrors([]);

      const response = await uploadCSV(file);
      setFileId(response.fileId);
      setHeaders(response.headers);
      setRowCount(response.rowCount);

      const previewData = await getCSVPreview(file);
      setPreview(previewData.data);

      const generatedMapping = generateMappingFromHeaders(response.headers, availableProperties);
      setCurrentMapping(generatedMapping.map(item => ({ ...item, isCustom: item.systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX) })));

      const columnMappingOptions = createColumnMappingOptions(response.headers, availableProperties);
      setMappingOptions(columnMappingOptions);

      showToast({
        type: 'success',
        message: 'File uploaded and preview generated successfully',
      });
    } catch (error: any) {
      setErrors(prevErrors => [...prevErrors, error.message || 'File upload failed']);
      showToast({
        type: 'error',
        message: error.message || 'File upload failed',
      });
    } finally {
      setIsLoading(false);
      setFile(file);
    }
  }, [availableProperties, showToast]);

  // Create handleUpdateMapping function to update column mapping
  const handleUpdateMapping = useCallback((csvColumn: string, systemProperty: string, isCustom: boolean) => {
    setCurrentMapping(prevMapping =>
      prevMapping.map(item =>
        item.csvColumn === csvColumn ? { ...item, systemProperty, isCustom } : item
      )
    );
  }, []);

  // Create handleGenerateMapping function to auto-generate mapping from headers
  const handleGenerateMapping = useCallback(() => {
    const generatedMapping = generateMappingFromHeaders(headers, availableProperties);
    setCurrentMapping(generatedMapping.map(item => ({ ...item, isCustom: item.systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX) })));
  }, [headers, availableProperties]);

  // Create handleSubmitMapping function to submit mapping and start processing
  const handleSubmitMapping = useCallback(async () => {
    try {
      setIsLoading(true);
      setErrors([]);

      const validationResult = validateMappingConfiguration(currentMapping);
      if (!validationResult.isValid) {
        setErrors(validationResult.errors.map(err => err.message));
        showToast({
          type: 'error',
          message: 'Invalid mapping configuration. Please correct the errors.',
        });
        return;
      }

      const response = await mapCSVColumns({
        fileId: fileId,
        mapping: currentMapping.map(item => ({ csvColumn: item.csvColumn, systemProperty: item.systemProperty }))
      });

      setJobId(response.jobId);
      setStatus(response.status);

      showToast({
        type: 'success',
        message: 'Mapping submitted and processing started',
      });
    } catch (error: any) {
      setErrors(prevErrors => [...prevErrors, error.message || 'Failed to submit mapping']);
      showToast({
        type: 'error',
        message: error.message || 'Failed to submit mapping',
      });
    } finally {
      setIsLoading(false);
    }
  }, [fileId, currentMapping, showToast]);

  // Create handleCheckStatus function to poll for processing status
  const handleCheckStatus = useCallback(async () => {
    if (!jobId) return;

    try {
      const response = await getCSVProcessingStatus(jobId);
      setStatus(response.status);
      setProgress(response.progress);
      setSummary(response.summary);
    } catch (error: any) {
      setErrors(prevErrors => [...prevErrors, error.message || 'Failed to check processing status']);
      showToast({
        type: 'error',
        message: error.message || 'Failed to check processing status',
      });
    }
  }, [jobId, showToast]);

  // Create handleCancelProcessing function to cancel ongoing processing
  const handleCancelProcessing = useCallback(async () => {
    if (!jobId) return;

    try {
      setIsLoading(true);
      setErrors([]);

      await cancelCSVProcessing(jobId);
      setStatus('cancelled');
      showToast({
        type: 'info',
        message: 'Processing cancelled',
      });
    } catch (error: any) {
      setErrors(prevErrors => [...prevErrors, error.message || 'Failed to cancel processing']);
      showToast({
        type: 'error',
        message: error.message || 'Failed to cancel processing',
      });
    } finally {
      setIsLoading(false);
    }
  }, [jobId, showToast]);

  // Create handleReset function to reset all state
  const handleReset = useCallback(() => {
    setFile(null);
    setHeaders([]);
    setPreview([]);
    setRowCount(0);
    setFileId('');
    setJobId('');
    setStatus('');
    setProgress(0);
    setSummary(null);
    setCurrentMapping([]);
    setMappingOptions([]);
    setErrors([]);
  }, []);

  // Set up polling effect when jobId changes
  useEffect(() => {
    if (jobId) {
      pollingInterval.current = window.setInterval(handleCheckStatus, 2000); // Poll every 2 seconds
    } else {
      if (pollingInterval.current) {
        clearInterval(pollingInterval.current);
        pollingInterval.current = null;
      }
    }

    // Clean up polling interval on unmount
    return () => {
      if (pollingInterval.current) {
        clearInterval(pollingInterval.current);
      }
    };
  }, [jobId, handleCheckStatus]);

  // Return all state variables and functions
  return {
    file,
    headers,
    preview,
    rowCount,
    fileId,
    jobId,
    status,
    progress,
    summary,
    currentMapping,
    mappingOptions,
    availableProperties,
    isLoading,
    errors,
    handleFileUpload,
    handleUpdateMapping,
    handleGenerateMapping,
    handleSubmitMapping,
    handleCancelProcessing,
    handleReset,
  };
};