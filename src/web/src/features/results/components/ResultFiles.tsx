# src/web/src/features/results/components/ResultFiles.tsx
```typescript
import React, { useState, useCallback } from 'react'; // react ^18.2.0
import {
  Box,
  Typography,
  Button,
  IconButton,
  List,
  ListItem,
  ListItemText,
  ListItemIcon,
  ListItemSecondaryAction,
  Tooltip,
  CircularProgress,
  Divider
} from '@mui/material'; // @mui/material v5.13.0
import {
  Download,
  Delete,
  InsertDriveFile,
  PictureAsPdf,
  Description,
  Image
} from '@mui/icons-material'; // @mui/icons-material v5.11.16
import Card from '../../../components/common/Card';
import useResults from '../hooks/useResults';
import { useToast } from '../../../hooks/useToast';
import { ResultFile } from '../../../types/result';

/**
 * Interface defining the props for the ResultFiles component.
 */
interface ResultFilesProps {
  /** Array of result files to display */
  files: ResultFile[];
  /** ID of the result these files belong to */
  resultId: string;
  /** Whether to show delete buttons for files */
  allowDelete?: boolean;
  /** Callback function when a file is deleted */
  onFileDeleted?: () => void;
  /** Optional custom title for the component */
  title?: string;
}

/**
 * A component that displays and manages files associated with experimental results.
 * It provides functionality for viewing, downloading, and optionally deleting result files uploaded by CRO users.
 *
 * @param {ResultFilesProps} props - The props for the component.
 * @returns {JSX.Element} The rendered component.
 */
const ResultFiles: React.FC<ResultFilesProps> = ({
  files,
  resultId,
  allowDelete = false,
  onFileDeleted = () => { },
  title = 'Result Files'
}) => {
  // State to track loading state for each file during download or delete operations
  const [loading, setLoading] = useState<Record<string, boolean>>({});

  // Custom hook for result file operations
  const { downloadResultFile, deleteResultFile } = useResults();

  // Hook for displaying toast notifications
  const { showToast } = useToast();

  /**
   * Handles downloading a file.
   * Sets loading state for the file to true, attempts to download the file,
   * shows success or error toast notifications, and sets loading state to false.
   *
   * @param {string} fileId - The ID of the file to download.
   * @param {string} fileName - The name of the file to download.
   */
  const handleDownload = useCallback(async (fileId: string, fileName: string) => {
    setLoading(prev => ({ ...prev, [fileId]: true }));
    try {
      const success = await downloadResultFile(fileId, fileName);
      if (success) {
        showToast({
          type: 'success',
          message: `File "${fileName}" downloaded successfully.`
        });
      } else {
        showToast({
          type: 'error',
          message: `Failed to download file "${fileName}".`
        });
      }
    } catch (error) {
      console.error('Download error:', error);
      showToast({
        type: 'error',
        message: `Failed to download file "${fileName}".`
      });
    } finally {
      setLoading(prev => ({ ...prev, [fileId]: false }));
    }
  }, [downloadResultFile, showToast]);

  /**
   * Handles deleting a file.
   * Sets loading state for the file to true, attempts to delete the file,
   * shows success or error toast notifications, and calls the onFileDeleted callback.
   *
   * @param {string} fileId - The ID of the file to delete.
   */
  const handleDelete = useCallback(async (fileId: string) => {
    setLoading(prev => ({ ...prev, [fileId]: true }));
    try {
      const success = await deleteResultFile(fileId);
      if (success) {
        showToast({
          type: 'success',
          message: 'File deleted successfully.'
        });
        onFileDeleted();
      } else {
        showToast({
          type: 'error',
          message: 'Failed to delete file.'
        });
      }
    } catch (error) {
      console.error('Delete error:', error);
      showToast({
        type: 'error',
        message: 'Failed to delete file.'
      });
    } finally {
      setLoading(prev => ({ ...prev, [fileId]: false }));
    }
  }, [deleteResultFile, showToast, onFileDeleted]);

  /**
   * Returns the appropriate icon based on file type.
   *
   * @param {string} fileType - The type of the file.
   * @returns {JSX.Element} The Material-UI icon component.
   */
  const getFileIcon = (fileType: string): JSX.Element => {
    if (fileType?.includes('pdf')) {
      return <PictureAsPdf />;
    }
    if (fileType?.includes('image')) {
      return <Image />;
    }
    if (fileType?.includes('text') || fileType?.includes('csv')) {
      return <Description />;
    }
    return <InsertDriveFile />;
  };

  /**
   * Formats file size in bytes to human-readable format.
   *
   * @param {number} bytes - The file size in bytes.
   * @returns {string} The formatted file size with appropriate unit.
   */
  const formatFileSize = (bytes: number): string => {
    if (bytes === 0) return '0 Bytes';
    const units = ['Bytes', 'KB', 'MB', 'GB', 'TB'];
    const k = 1024;
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + units[i];
  };

  return (
    <Card title={title}>
      {files?.length === 0 ? (
        <Typography variant="body2" color="textSecondary" align="center">
          No files available
        </Typography>
      ) : (
        <List>
          {files?.map((file, index) => (
            <React.Fragment key={file.id}>
              <ListItem>
                <ListItemIcon>
                  {getFileIcon(file.file_type)}
                </ListItemIcon>
                <ListItemText
                  primary={file.file_name}
                  secondary={`${formatFileSize(file.file_size)} - Uploaded: ${file.uploaded_at}`}
                />
                <ListItemSecondaryAction>
                  <Tooltip title="Download">
                    <IconButton
                      edge="end"
                      aria-label="download"
                      onClick={() => handleDownload(file.id, file.file_name)}
                      disabled={loading[file.id]}
                    >
                      {loading[file.id] ? (
                        <CircularProgress size={24} />
                      ) : (
                        <Download />
                      )}
                    </IconButton>
                  </Tooltip>
                  {allowDelete && (
                    <Tooltip title="Delete">
                      <IconButton
                        edge="end"
                        aria-label="delete"
                        onClick={() => handleDelete(file.id)}
                        disabled={loading[file.id]}
                      >
                        {loading[file.id] ? (
                          <CircularProgress size={24} />
                        ) : (
                          <Delete />
                        )}
                      </IconButton>
                    </Tooltip>
                  )}
                </ListItemSecondaryAction>
              </ListItem>
              {index < files.length - 1 && <Divider />}
            </React.Fragment>
          ))}
        </List>
      )}
    </Card>
  );
};

export default ResultFiles;