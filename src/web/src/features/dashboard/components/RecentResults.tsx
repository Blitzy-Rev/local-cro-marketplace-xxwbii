import React, { useEffect, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import { 
  Box, 
  Typography, 
  Chip, 
  Button, 
  List, 
  ListItem, 
  ListItemText, 
  ListItemSecondaryAction 
} from '@mui/material'; // @mui/material v5.13+
import { LabResults as LabResultsIcon } from '@mui/icons-material'; // @mui/icons-material v5.13+

import Card from '../../../components/common/Card';
import useResults from '../../results/hooks/useResults';
import { Result, ResultStatus } from '../../../types/result';
import { formatDate } from '../../../utils/formatters';

/**
 * Props interface for the RecentResults component
 */
interface RecentResultsProps {
  /** Callback function when the View All button is clicked */
  onViewAll: () => void;
}

/**
 * Component that displays a summary of recent experimental results on the dashboard
 * Shows a limited number of the most recent results with their status and
 * provides navigation to the full results page.
 */
const RecentResults: React.FC<RecentResultsProps> = ({ onViewAll }) => {
  const navigate = useNavigate();
  
  // Set up useResults hook with filter for recent results
  const { results, loading, error } = useResults({
    initialFilter: {
      page: 1,
      page_size: 5,
      sort_by: 'uploaded_at',
      sort_desc: true
    }
  });
  
  // Helper function to determine color for status chips
  const getStatusColor = (status: string) => {
    switch (status) {
      case ResultStatus.UPLOADED:
        return { color: 'primary' as const, variant: 'filled' as const };
      case ResultStatus.APPROVED:
        return { color: 'success' as const, variant: 'outlined' as const };
      case ResultStatus.REJECTED:
        return { color: 'error' as const, variant: 'outlined' as const };
      case ResultStatus.PENDING:
        return { color: 'warning' as const, variant: 'outlined' as const };
      default:
        return { color: 'default' as const, variant: 'outlined' as const };
    }
  };
  
  // Handler for when a result item is clicked
  const handleResultClick = (resultId: string) => {
    navigate(`/results/${resultId}`);
  };
  
  // Memoize the number of new results
  const newResultsCount = useMemo(() => {
    return results.filter(r => r.status === ResultStatus.UPLOADED).length;
  }, [results]);
  
  // Memoize card header
  const header = useMemo(() => (
    <Box sx={{ p: 2, display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
      <Box sx={{ display: 'flex', alignItems: 'center' }}>
        <LabResultsIcon color="primary" sx={{ mr: 1 }} />
        <Typography variant="h6">Recent Results</Typography>
      </Box>
      {newResultsCount > 0 && (
        <Chip 
          label={`${newResultsCount} new`} 
          color="primary" 
          size="small" 
        />
      )}
    </Box>
  ), [newResultsCount]);
  
  // Memoize card footer
  const footer = useMemo(() => (
    <Box sx={{ p: 1, display: 'flex', justifyContent: 'flex-end' }}>
      <Button onClick={onViewAll} color="primary">
        View All Results
      </Button>
    </Box>
  ), [onViewAll]);
  
  return (
    <Card 
      header={header}
      footer={footer}
      elevation={1}
      sx={{ height: '100%' }}
    >
      <Box sx={{ px: 2, pb: 1 }}>
        {loading ? (
          <Box sx={{ py: 2, display: 'flex', justifyContent: 'center' }}>
            <Typography>Loading recent results...</Typography>
          </Box>
        ) : error ? (
          <Box sx={{ py: 2, display: 'flex', justifyContent: 'center' }}>
            <Typography color="error">Error loading results: {error}</Typography>
          </Box>
        ) : results.length === 0 ? (
          <Box sx={{ py: 2, display: 'flex', justifyContent: 'center' }}>
            <Typography>No results available</Typography>
          </Box>
        ) : (
          <List disablePadding>
            {results.map((result: Result) => {
              const statusColor = getStatusColor(result.status);
              const isNew = result.status === ResultStatus.UPLOADED;
              
              // Format the result identifier for display
              let displayName = "Result";
              if (result.submission && result.submission.experiment_id) {
                displayName = `Experiment ${result.submission.experiment_id.substring(0, 8)}...`;
              }
              
              // Format the status for display (capitalize first letter)
              const statusDisplay = result.status.charAt(0).toUpperCase() + result.status.slice(1);
              
              return (
                <ListItem
                  key={result.id}
                  divider
                  button
                  onClick={() => handleResultClick(result.id)}
                  sx={{ 
                    py: 1,
                    backgroundColor: isNew ? 'rgba(25, 118, 210, 0.08)' : 'transparent',
                    '&:hover': {
                      backgroundColor: isNew ? 'rgba(25, 118, 210, 0.12)' : 'rgba(0, 0, 0, 0.04)'
                    }
                  }}
                >
                  <ListItemText
                    primary={displayName}
                    secondary={formatDate(result.uploaded_at)}
                  />
                  <ListItemSecondaryAction>
                    <Chip
                      label={statusDisplay}
                      size="small"
                      color={statusColor.color}
                      variant={statusColor.variant}
                    />
                  </ListItemSecondaryAction>
                </ListItem>
              );
            })}
          </List>
        )}
      </Box>
    </Card>
  );
};

export default RecentResults;