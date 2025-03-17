import React from 'react';
import { Box, Grid, Typography, Chip, Divider } from '@mui/material'; // @mui/material v5.13.0
import { CheckCircle, Cancel, Schedule, Info } from '@mui/icons-material'; // @mui/icons-material v5.11.16

import Card from '../../../components/common/Card';
import { ResultDetailed, ResultStatus } from '../../../types/result';
import { formatDate } from '../../../utils/formatters';

/**
 * Props for the ResultSummary component
 */
interface ResultSummaryProps {
  /**
   * The detailed result data to display
   */
  result: ResultDetailed;
  /**
   * Optional custom title for the component
   */
  title?: string;
}

/**
 * Returns the appropriate icon component based on result status
 */
const getStatusIcon = (status: ResultStatus): JSX.Element => {
  switch (status) {
    case ResultStatus.APPROVED:
      return <CheckCircle fontSize="small" />;
    case ResultStatus.REJECTED:
      return <Cancel fontSize="small" />;
    case ResultStatus.UPLOADED:
      return <Schedule fontSize="small" />;
    case ResultStatus.PENDING:
      return <Info fontSize="small" />;
    default:
      return <Info fontSize="small" />;
  }
};

/**
 * Returns the appropriate color for status chip based on result status
 */
const getStatusColor = (status: ResultStatus): 'success' | 'error' | 'warning' | 'info' | 'default' => {
  switch (status) {
    case ResultStatus.APPROVED:
      return 'success';
    case ResultStatus.REJECTED:
      return 'error';
    case ResultStatus.UPLOADED:
      return 'warning';
    case ResultStatus.PENDING:
      return 'info';
    default:
      return 'default';
  }
};

/**
 * Returns a human-readable label for the result status
 */
const getStatusLabel = (status: ResultStatus): string => {
  switch (status) {
    case ResultStatus.APPROVED:
      return 'Approved';
    case ResultStatus.REJECTED:
      return 'Rejected';
    case ResultStatus.UPLOADED:
      return 'Uploaded';
    case ResultStatus.PENDING:
      return 'Pending';
    default:
      return 'Unknown';
  }
};

/**
 * A component that displays summary information about an experimental result
 */
const ResultSummary: React.FC<ResultSummaryProps> = ({ result, title = 'Result Summary' }) => {
  // Get molecule count (unique molecule IDs in data points)
  const moleculeCount = result.data_points 
    ? new Set(result.data_points.map(dp => dp.molecule_id)).size 
    : 0;
    
  // Get file count
  const fileCount = result.files?.length || 0;
  
  return (
    <Card header={<Typography variant="h6">{title}</Typography>}>
      <Grid container spacing={3}>
        {/* Left column: Submission information */}
        <Grid item xs={12} md={6}>
          <Typography variant="subtitle1" gutterBottom fontWeight="bold">
            Submission Information
          </Typography>
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Submission ID
            </Typography>
            <Typography variant="body1">
              #{result.submission_id}
            </Typography>
          </Box>
          
          {result.submission?.experiment_id && (
            <Box sx={{ mb: 2 }}>
              <Typography variant="body2" color="text.secondary">
                Experiment ID
              </Typography>
              <Typography variant="body1">
                #{result.submission.experiment_id}
              </Typography>
            </Box>
          )}
          
          {result.submission?.status && (
            <Box sx={{ mb: 2 }}>
              <Typography variant="body2" color="text.secondary">
                Submission Status
              </Typography>
              <Typography variant="body1">
                {result.submission.status}
              </Typography>
            </Box>
          )}
          
          <Divider sx={{ my: 2 }} />
          
          <Typography variant="subtitle1" gutterBottom fontWeight="bold">
            Result Details
          </Typography>
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Result ID
            </Typography>
            <Typography variant="body1">
              #{result.id}
            </Typography>
          </Box>
          
          {result.notes && (
            <Box sx={{ mb: 2 }}>
              <Typography variant="body2" color="text.secondary">
                Notes
              </Typography>
              <Typography variant="body1">
                {result.notes}
              </Typography>
            </Box>
          )}
        </Grid>
        
        {/* Right column: Status and timeline */}
        <Grid item xs={12} md={6}>
          <Box sx={{ mb: 2, display: 'flex', justifyContent: 'flex-start' }}>
            <Chip 
              icon={getStatusIcon(result.status)}
              label={getStatusLabel(result.status)}
              color={getStatusColor(result.status)}
              sx={{ fontWeight: 'bold' }}
            />
          </Box>
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Uploaded
            </Typography>
            <Typography variant="body1">
              {formatDate(result.uploaded_at)}
            </Typography>
          </Box>
          
          {result.approved_at && (
            <Box sx={{ mb: 2 }}>
              <Typography variant="body2" color="text.secondary">
                {result.status === ResultStatus.APPROVED ? 'Approved' : 'Rejected'}
              </Typography>
              <Typography variant="body1">
                {formatDate(result.approved_at)}
              </Typography>
            </Box>
          )}
          
          <Divider sx={{ my: 2 }} />
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Molecules
            </Typography>
            <Typography variant="body1">
              {moleculeCount > 0 ? `${moleculeCount} molecules` : 'No molecules'}
            </Typography>
          </Box>
          
          <Box sx={{ mb: 2 }}>
            <Typography variant="body2" color="text.secondary">
              Files
            </Typography>
            <Typography variant="body1">
              {fileCount > 0 ? `${fileCount} files` : 'No files'}
            </Typography>
          </Box>
        </Grid>
      </Grid>
    </Card>
  );
};

export default ResultSummary;