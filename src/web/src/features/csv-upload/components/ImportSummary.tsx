import React, { useMemo } from 'react'; // React 18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom ^6.10.0
import {
  Box,
  Typography,
  Paper,
  Alert,
  Grid,
  Divider,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Chip,
} from '@mui/material'; // @mui/material ^5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles ^5.13.0
import {
  CheckCircle,
  Error,
  Warning,
  Info,
  ArrowForward,
  Refresh,
} from '@mui/icons-material'; // @mui/icons-material ^5.13.0
import { useCSVImport } from '../hooks/useCSVImport';
import Button from '../../../components/common/Button';
import { useToast } from '../../../hooks/useToast';
import { ImportSummaryProps, ImportSummaryData } from '../../../types/csv';

/**
 * Styled component for the summary container
 */
const SummaryContainer = styled(Paper)(({ theme }) => ({
  padding: theme.spacing(3), // 24px
  marginBottom: theme.spacing(3), // 24px
  borderRadius: theme.spacing(1), // 8px
}));

/**
 * Styled component for the summary title
 */
const SummaryTitle = styled(Typography)(({ theme }) => ({
  marginBottom: theme.spacing(2), // 16px
  fontWeight: 500,
}));

/**
 * Styled component for the summary description
 */
const SummaryDescription = styled(Typography)(({ theme }) => ({
  marginBottom: theme.spacing(3), // 24px
  color: theme.palette.text.secondary,
}));

/**
 * Styled component for the statistics grid
 */
const StatGrid = styled(Grid)(({ theme }) => ({
  marginBottom: theme.spacing(3), // 24px
}));

/**
 * Styled component for the statistics item
 */
const StatItem = styled(Box)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  padding: theme.spacing(2), // 16px
}));

/**
 * Styled component for the statistics value
 */
const StatValue = styled(Typography)(({ theme }) => ({
  fontSize: '24px',
  fontWeight: 'bold',
  marginBottom: theme.spacing(1), // 8px
}));

/**
 * Styled component for the statistics label
 */
const StatLabel = styled(Typography)(({ theme }) => ({
  color: theme.palette.text.secondary,
}));

/**
 * Styled component for the detail section
 */
const DetailSection = styled(Box)(({ theme }) => ({
  marginBottom: theme.spacing(3), // 24px
}));

/**
 * Styled component for the detail title
 */
const DetailTitle = styled(Typography)(({ theme }) => ({
  fontWeight: 500,
  marginBottom: theme.spacing(1), // 8px
}));

/**
 * Styled component for the action container
 */
const ActionContainer = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  marginTop: theme.spacing(3), // 24px
}));

/**
 * Component for displaying CSV import results summary
 * @param props - Props for the component
 * @returns JSX.Element
 */
const ImportSummary: React.FC<ImportSummaryProps> = ({ onStartNew }) => {
  // Extract summary, file, and handleReset from useCSVImport hook
  const { summary, file, handleReset } = useCSVImport();

  // Initialize navigate function from useNavigate hook
  const navigate = useNavigate();

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  /**
   * Navigates to the molecules page filtered by the imported library
   */
  const handleViewMolecules = () => {
    // Navigate to '/molecules' with query parameter for the imported library
    navigate(`/molecules?library_id=${summary?.importedLibraryId}`);

    // Show success toast notification
    showToast({
      type: 'success',
      message: 'Navigating to imported molecules',
    });
  };

  /**
   * Resets the import process and calls the onStartNew callback
   */
  const handleStartNew = () => {
    // Call handleReset from useCSVImport hook
    handleReset();

    // Call onStartNew callback from props
    onStartNew();
  };

  /**
   * Formats processing time in a readable format
   * @param milliseconds - Processing time in milliseconds
   * @returns Formatted duration string
   */
  const formatDuration = (milliseconds: number): string => {
    // Convert milliseconds to seconds
    const seconds = milliseconds / 1000;

    // If less than 60 seconds, return seconds with appropriate unit
    if (seconds < 60) {
      return `${seconds.toFixed(2)} seconds`;
    }

    // If 60 seconds or more, convert to minutes and seconds format
    const minutes = Math.floor(seconds / 60);
    const remainingSeconds = seconds % 60;

    // Return formatted string
    return `${minutes} minutes and ${remainingSeconds.toFixed(2)} seconds`;
  };

  // Use useMemo to calculate success rate and other statistics
  const { successRate, hasErrors, hasWarnings, formattedProcessingTime } = useMemo(() => {
    // Calculate percentage of successful imports (successfulImports / totalMolecules * 100)
    const successRate = summary ? (summary.successfulImports / summary.totalMolecules) * 100 : 0;

    // Check if there are any errors (summary?.errors?.length > 0)
    const hasErrors = summary?.errors?.length > 0;

    // Check if there are any warnings (summary?.warnings?.length > 0)
    const hasWarnings = summary?.warnings?.length > 0;

    // Format processing time using formatDuration function
    const formattedProcessingTime = summary ? formatDuration(summary.processingTime) : 'N/A';

    return { successRate, hasErrors, hasWarnings, formattedProcessingTime };
  }, [summary]);

  return (
    <>
      {/* Check if summary is available */}
      {!summary ? (
        // If no summary, render Alert component with message to complete import process
        <Alert severity="info">Please complete the CSV import process to view the summary.</Alert>
      ) : (
        // Otherwise, render SummaryContainer with title 'Import Summary'
        <SummaryContainer>
          {/* Render title 'Import Summary' */}
          <SummaryTitle variant="h5">Import Summary</SummaryTitle>

          {/* Render success message with CheckCircle icon if successRate is 100% */}
          {successRate === 100 && (
            <Alert icon={<CheckCircle />} severity="success">
              All molecules imported successfully!
            </Alert>
          )}

          {/* Render partial success message with Warning icon if successRate is between 0 and 100% */}
          {successRate > 0 && successRate < 100 && (
            <Alert icon={<Warning />} severity="warning">
              Some molecules were not imported due to errors.
            </Alert>
          )}

          {/* Render failure message with Error icon if successRate is 0% */}
          {successRate === 0 && (
            <Alert icon={<Error />} severity="error">
              No molecules were imported due to errors.
            </Alert>
          )}

          {/* Render StatGrid with statistics for total molecules, successful imports, errors, and warnings */}
          <StatGrid container spacing={3}>
            <Grid item xs={12} sm={6} md={3}>
              <StatItem>
                <StatValue>{summary.totalMolecules}</StatValue>
                <StatLabel>Total Molecules</StatLabel>
              </StatItem>
            </Grid>
            <Grid item xs={12} sm={6} md={3}>
              <StatItem>
                <StatValue>{summary.successfulImports}</StatValue>
                <StatLabel>Successful Imports</StatLabel>
              </StatItem>
            </Grid>
            <Grid item xs={12} sm={6} md={3}>
              <StatItem>
                <StatValue>{summary.errors.length}</StatValue>
                <StatLabel>Errors</StatLabel>
              </StatItem>
            </Grid>
            <Grid item xs={12} sm={6} md={3}>
              <StatItem>
                <StatValue>{summary.warnings.length}</StatValue>
                <StatLabel>Warnings</StatLabel>
              </StatItem>
            </Grid>
          </StatGrid>

          <Divider />

          {/* Render DetailSection with processing details (file name, size, processing time) */}
          <DetailSection>
            <DetailTitle variant="subtitle1">Processing Details</DetailTitle>
            <Typography variant="body2">
              File Name: {file?.name}
            </Typography>
            <Typography variant="body2">
              File Size: {(file?.size / 1024).toFixed(2)} KB
            </Typography>
            <Typography variant="body2">
              Processing Time: {formattedProcessingTime}
            </Typography>
          </DetailSection>

          {/* If hasErrors, render error details section with list of errors */}
          {hasErrors && (
            <DetailSection>
              <DetailTitle variant="subtitle1">Error Details</DetailTitle>
              <List>
                {summary.errors.map((error: { row: number; message: string }, index: number) => (
                  <ListItem key={index} alignItems="flex-start">
                    <ListItemIcon>
                      <Error color="error" />
                    </ListItemIcon>
                    <ListItemText
                      primary={`Row ${error.row}`}
                      secondary={error.message}
                    />
                  </ListItem>
                ))}
              </List>
            </DetailSection>
          )}

          {/* If hasWarnings, render warning details section with list of warnings */}
          {hasWarnings && (
            <DetailSection>
              <DetailTitle variant="subtitle1">Warning Details</DetailTitle>
              <List>
                {summary.warnings.map((warning: { row: number; message: string }, index: number) => (
                  <ListItem key={index} alignItems="flex-start">
                    <ListItemIcon>
                      <Warning color="warning" />
                    </ListItemIcon>
                    <ListItemText
                      primary={`Row ${warning.row}`}
                      secondary={warning.message}
                    />
                  </ListItem>
                ))}
              </List>
            </DetailSection>
          )}

          <Divider />

          {/* Render ActionContainer with buttons for viewing molecules and starting new import */}
          <ActionContainer>
            <Button
              variant="outlined"
              startIcon={<ArrowForward />}
              onClick={handleViewMolecules}
            >
              View Imported Molecules
            </Button>
            <Button
              variant="contained"
              color="primary"
              startIcon={<Refresh />}
              onClick={handleStartNew}
            >
              Start New Import
            </Button>
          </ActionContainer>
        </SummaryContainer>
      )}
    </>
  );
};

export default ImportSummary;