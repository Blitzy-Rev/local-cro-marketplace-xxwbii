import React, { useState, useRef, useCallback, useEffect } from 'react'; // React 18.2+
import { Box, Typography, Paper, CircularProgress, Alert } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { CloudUpload, InsertDriveFile } from '@mui/icons-material'; // v5.13+
import Button from './Button';
import { useToast } from '../../hooks/useToast';

/**
 * Props interface for the FileUpload component
 */
interface FileUploadProps {
  /** Array of accepted file types (e.g., ['.csv', 'text/csv']) */
  accept?: string[];
  /** Maximum file size in bytes (default: 10MB) */
  maxSize?: number;
  /** Whether multiple files can be selected (default: false) */
  multiple?: boolean;
  /** Whether the upload functionality is disabled (default: false) */
  disabled?: boolean;
  /** Whether to show a preview of selected files (default: true) */
  showPreview?: boolean;
  /** Callback function called when files are selected and validated */
  onFilesSelected: (files: File[]) => void;
  /** Callback function called when file validation fails */
  onFileValidationFail?: (error: string) => void;
  /** Custom validation function for additional file validation */
  customValidation?: (file: File) => Promise<{ isValid: boolean; error?: string }>;
  /** Helper text displayed in the upload area (default: 'Drag and drop files here, or click to select files') */
  helperText?: string;
  /** Label for the upload component (default: 'Upload Files') */
  label?: string;
  /** Text for the upload button (default: 'Select Files') */
  buttonText?: string;
  /** Whether the component is in a loading state (default: false) */
  isLoading?: boolean;
  /** Additional CSS class for styling */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
}

/**
 * A styled Paper component for the upload container
 */
const UploadContainer = styled(Paper, {
  shouldForwardProp: (prop) => prop !== 'dragActive' && prop !== 'disabled',
})<{ dragActive?: boolean; disabled?: boolean }>(({ theme, dragActive, disabled }) => ({
  padding: theme.spacing(2),
  borderRadius: theme.spacing(1),
  border: `2px dashed`,
  borderColor: dragActive ? theme.palette.primary.main : theme.palette.divider,
  backgroundColor: dragActive ? 'rgba(0, 0, 0, 0.04)' : 'transparent',
  transition: 'all 0.3s ease',
  cursor: 'pointer',
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  minHeight: '200px',
  position: 'relative',
  '&:hover': {
    borderColor: !disabled && theme.palette.primary.main,
    backgroundColor: !disabled && 'rgba(0, 0, 0, 0.04)',
  },
  opacity: disabled ? 0.6 : 1,
  pointerEvents: disabled ? 'none' : 'auto',
}));

/**
 * A styled input component for the hidden file input
 */
const HiddenInput = styled('input')({
  display: 'none',
});

/**
 * A styled Box component for the file preview
 */
const FilePreview = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  marginTop: '16px',
  padding: '8px',
  borderRadius: '4px',
  backgroundColor: 'rgba(0, 0, 0, 0.04)',
  width: '100%',
});

/**
 * A reusable file upload component with drag-and-drop functionality
 */
const FileUpload: React.FC<FileUploadProps> = ({
  accept = [],
  maxSize = 10485760, // 10MB
  multiple = false,
  disabled = false,
  showPreview = true,
  onFilesSelected,
  onFileValidationFail = () => {},
  customValidation = null,
  helperText = 'Drag and drop files here, or click to select files',
  label = 'Upload Files',
  buttonText = 'Select Files',
  isLoading = false,
  className,
  style,
}) => {
  // State for drag and drop active
  const [dragActive, setDragActive] = useState<boolean>(false);
  // State for selected files
  const [selectedFiles, setSelectedFiles] = useState<File[]>([]);
  // State for validation error
  const [error, setError] = useState<string>('');
  // Ref for the hidden file input
  const inputRef = useRef<HTMLInputElement>(null);
  // Toast hook for displaying notifications
  const { showToast } = useToast();

  /**
   * Handles drag events to update the dragActive state
   * @param e React.DragEvent<HTMLDivElement>
   */
  const handleDrag = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    if (e.type === 'dragenter' || e.type === 'dragover') {
      setDragActive(true);
    } else if (e.type === 'dragleave' || e.type === 'drop') {
      setDragActive(false);
    }
  };

  /**
   * Handles file drop events
   * @param e React.DragEvent<HTMLDivElement>
   */
  const handleDrop = (e: React.DragEvent<HTMLDivElement>) => {
    e.preventDefault();
    e.stopPropagation();
    setDragActive(false);
    const files = [...e.dataTransfer.files];
    handleFiles(files);
  };

  /**
   * Handles file selection from the file input
   * @param e React.ChangeEvent<HTMLInputElement>
   */
  const handleChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    if (e.target.files) {
      const files = [...e.target.files];
      handleFiles(files);
      // Reset the input value to allow selecting the same file again
      e.target.value = '';
    }
  };

  /**
   * Handles click on the upload button
   */
  const handleButtonClick = () => {
    if (!disabled && inputRef.current) {
      inputRef.current.click();
    }
  };

  /**
   * Validates files against accept types, size limit, and custom validation
   * @param files File[]
   * @returns Promise<{isValid: boolean, files?: File[], error?: string}>
   */
  const validateFiles = useCallback(async (files: File[]) => {
    if (!files || files.length === 0) {
      return { isValid: false, error: 'No files selected' };
    }

    for (const file of files) {
      if (accept.length > 0 && !accept.includes(file.type)) {
        return { isValid: false, error: `File type not accepted: ${file.name}` };
      }

      if (file.size > maxSize) {
        return { isValid: false, error: `File is too large: ${file.name}` };
      }

      if (customValidation) {
        const validationResult = await customValidation(file);
        if (!validationResult.isValid) {
          return { isValid: false, error: validationResult.error || `Validation failed for: ${file.name}` };
        }
      }
    }

    return { isValid: true, files };
  }, [accept, maxSize, customValidation]);

  /**
   * Processes files after selection or drop
   * @param fileList FileList | File[]
   */
  const handleFiles = useCallback(async (fileList: FileList | File[]) => {
    let files: File[] = [];
    if (fileList instanceof FileList) {
      files = Array.from(fileList);
    } else {
      files = fileList;
    }

    if (!multiple && files.length > 1) {
      files = [files[0]];
    }

    const validationResult = await validateFiles(files);

    if (!validationResult.isValid) {
      setError(validationResult.error);
      onFileValidationFail(validationResult.error);
      showToast({ type: 'error', message: validationResult.error });
    } else {
      setError('');
      setSelectedFiles(validationResult.files);
      onFilesSelected(validationResult.files);
    }
  }, [multiple, validateFiles, onFilesSelected, onFileValidationFail, showToast]);

  /**
   * Clears selected files and error state
   */
  const clearFiles = useCallback(() => {
    setSelectedFiles([]);
    setError('');
    if (inputRef.current) {
      inputRef.current.value = '';
    }
  }, []);

  useEffect(() => {
    // Clear files when disabled state changes
    if (disabled) {
      clearFiles();
    }
  }, [disabled, clearFiles]);

  return (
    <Box className={className} style={style}>
      {label && <Typography variant="subtitle1" gutterBottom>{label}</Typography>}
      <UploadContainer
        dragActive={dragActive}
        disabled={disabled}
        onClick={handleButtonClick}
        onDragEnter={handleDrag}
        onDragOver={handleDrag}
        onDragLeave={handleDrag}
        onDrop={handleDrop}
      >
        <HiddenInput
          type="file"
          multiple={multiple}
          accept={accept.join(',')}
          onChange={handleChange}
          ref={inputRef}
          disabled={disabled}
        />
        <CloudUpload color="primary" sx={{ fontSize: 60 }} />
        <Typography variant="body2" color="textSecondary" align="center">
          {helperText}
        </Typography>
        <Button onClick={handleButtonClick} disabled={disabled} loading={isLoading}>
          {buttonText}
        </Button>
        {isLoading && (
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
      </UploadContainer>
      {error && (
        <Alert severity="error" sx={{ mt: 2 }}>
          {error}
        </Alert>
      )}
      {showPreview && selectedFiles.length > 0 && (
        <Box mt={2}>
          {selectedFiles.map((file) => (
            <FilePreview key={file.name}>
              <InsertDriveFile color="action" sx={{ mr: 1 }} />
              <Typography variant="body2">{file.name}</Typography>
            </FilePreview>
          ))}
        </Box>
      )}
    </Box>
  );
};

export default FileUpload;