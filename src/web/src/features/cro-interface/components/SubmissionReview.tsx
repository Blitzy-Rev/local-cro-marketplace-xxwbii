# src/web/src/features/cro-interface/components/SubmissionReview.tsx
```typescript
import React, { useState, useEffect } from 'react'; // React 18.2+
import { Box, Typography, Grid, Divider, Paper, Chip, List, ListItem, ListItemText } from '@mui/material'; // @mui/material v5.13.0
import { ArrowBack, Science, Assignment, CheckCircle, Cancel } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Loading from '../../../components/common/Loading';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import PropertyBadge from '../../../components/molecular/PropertyBadge';
import useToast from '../../../hooks/useToast';
import { useCROInterface } from '../hooks/useCROInterface';
import QuoteForm from './QuoteForm';
import { getSubmission, SubmissionDetailed, SubmissionStatus } from '../../../types/submission';

/**
 * Interface defining the props for the SubmissionReview component
 */
interface SubmissionReviewProps {
  /** ID of the submission to review */
  submissionId: string;
  /** Callback function called after successful quote submission */
  onQuoteSubmitted: () => void;
  /** Callback function called when submission is declined */
  onDecline: () => void;
  /** Callback function called when navigating back to submission list */
  onBack: () => void;
}

/**
 * A component that displays detailed information about a submission and allows CRO users to provide a quote or decline
 * @param {SubmissionReviewProps} props - The props for the component
 * @returns {JSX.Element} Rendered component
 */
const SubmissionReview: React.FC<SubmissionReviewProps> = ({ submissionId, onQuoteSubmitted, onDecline, onBack }) => {
  // LD1: Extract submissionId, onQuoteSubmitted, onDecline, and onBack from props
  
  // LD1: Initialize state for submission data, loading status, and error handling
  const [submission, setSubmission] = useState<SubmissionDetailed | null>(null);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string | null>(null);
  
  // LD1: Initialize state for showing quote form
  const [showQuoteForm, setShowQuoteForm] = useState<boolean>(false);
  
  // LD1: Get showToast function from useToast hook
  const { showToast } = useToast();
  
  // LD1: Define handleProvideQuote function from useCROInterface hook
  const { handleProvideQuote } = useCROInterface();

  // LD1: Define fetchSubmissionData function to load submission details
  const fetchSubmissionData = async () => {
    setLoading(true);
    setError(null);
    try {
      const data = await getSubmission(submissionId);
      setSubmission(data);
    } catch (e: any) {
      setError(e?.message || 'Failed to load submission details.');
    } finally {
      setLoading(false);
    }
  };

  // LD1: Set up effect to fetch submission data when component mounts or submissionId changes
  useEffect(() => {
    fetchSubmissionData();
  }, [submissionId]);

  // LD1: Define handleShowQuoteForm function to display the quote form
  const handleShowQuoteForm = () => {
    setShowQuoteForm(true);
  };

  // LD1: Define handleQuoteSubmitted function to handle successful quote submission
  const handleQuoteSubmitted = () => {
    setShowQuoteForm(false);
    onQuoteSubmitted();
  };

  // LD1: Define handleDeclineSubmission function to decline the submission
  const handleDeclineSubmission = () => {
    onDecline();
  };

  // LD1: Define handleBack function to navigate back to the submission list
  const handleBack = () => {
    onBack();
  };

  // LD1: Render loading state while data is being fetched
  if (loading) {
    return <Loading fullScreen message="Loading submission details..." />;
  }

  // LD1: Render error message if data fetching fails
  if (error) {
    return (
      <Box textAlign="center" mt={4}>
        <Typography color="error">{error}</Typography>
        <Button onClick={fetchSubmissionData}>Retry</Button>
      </Box>
    );
  }

  // LD1: Render submission details including:
  //   - Header with back button and submission ID
  //   - Experiment details section with name, type, and parameters
  //   - Molecules section with structure visualization and properties
  //   - Action buttons for providing quote or declining submission
  return (
    <Paper elevation={3} style={{ padding: '20px' }}>
      {/* Header */}
      <Box display="flex" alignItems="center" mb={2}>
        <Button onClick={handleBack} startIcon={<ArrowBack />}>
          Back
        </Button>
        <Typography variant="h6" style={{ marginLeft: '20px', flexGrow: 1 }}>
          Submission ID: {submission?.id}
        </Typography>
        {submission?.status && (
          <Chip
            label={submission.status}
            color={
              submission.status === SubmissionStatus.PENDING
                ? 'primary'
                : submission.status === SubmissionStatus.APPROVED
                ? 'success'
                : submission.status === SubmissionStatus.REJECTED
                ? 'error'
                : 'info'
            }
          />
        )}
      </Box>

      <Divider style={{ margin: '10px 0' }} />

      {/* Experiment Details */}
      <Grid container spacing={2}>
        <Grid item xs={6}>
          <Card>
            <CardContent>
              <Typography variant="h6" gutterBottom>
                Experiment Details
              </Typography>
              <Typography variant="body1">
                Name: {submission?.experiment?.name}
              </Typography>
              <Typography variant="body1">
                Type: {submission?.experiment?.type?.name}
              </Typography>
              <Typography variant="body2" color="textSecondary">
                Client: {submission?.experiment?.creator?.email}
              </Typography>
              <Typography variant="body2" color="textSecondary">
                Submission Date: {submission?.submitted_at}
              </Typography>
              <Typography variant="subtitle2" style={{ marginTop: '10px' }}>
                Parameters:
              </Typography>
              <List>
                {submission?.experiment?.parameters?.map((param) => (
                  <ListItem key={param.parameter_name}>
                    <ListItemText
                      primary={param.parameter_name}
                      secondary={param.parameter_value}
                    />
                  </ListItem>
                ))}
              </List>
            </CardContent>
          </Card>
        </Grid>

        {/* Molecules */}
        <Grid item xs={6}>
          <Card>
            <CardContent>
              <Typography variant="h6" gutterBottom>
                Molecules
              </Typography>
              <Typography variant="body1">
                {submission?.experiment?.molecule_count} molecules
              </Typography>
              <List>
                {submission?.experiment?.molecules?.map((molecule) => (
                  <ListItem key={molecule.molecule_id} alignItems="flex-start">
                    <MoleculeViewer
                      smiles={molecule.molecule?.smiles}
                      width={100}
                      height={100}
                      interactive={false}
                    />
                    <ListItemText
                      primary={molecule.molecule?.smiles}
                      secondary={
                        <Box>
                          {molecule.molecule?.properties?.map((property) => (
                            <PropertyBadge
                              key={property.id}
                              propertyName={property.property_name}
                              value={property.property_value}
                            />
                          ))}
                        </Box>
                      }
                    />
                  </ListItem>
                ))}
              </List>
            </CardContent>
          </Card>
        </Grid>
      </Grid>

      <Divider style={{ margin: '20px 0' }} />

      {/* Action Buttons */}
      <Box display="flex" justifyContent="flex-end">
        <Button variant="contained" color="primary" onClick={handleShowQuoteForm}>
          Provide Quote
        </Button>
        <Button onClick={handleDeclineSubmission} style={{ marginLeft: '10px' }}>
          Decline Submission
        </Button>
      </Box>

      {/* Conditionally render QuoteForm component when showQuoteForm is true */}
      {showQuoteForm && submission && (
        <QuoteForm
          submission={submission}
          onSubmitted={handleQuoteSubmitted}
          onCancel={() => setShowQuoteForm(false)}
        />
      )}
    </Paper>
  );
};

export default SubmissionReview;