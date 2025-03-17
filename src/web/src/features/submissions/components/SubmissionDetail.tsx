import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import { useParams, useNavigate, Link } from 'react-router-dom'; // react-router-dom v6.11.0
import { Box, Typography, Divider, Grid, Paper, Chip, Stack, Alert } from '@mui/material'; // @mui/material v5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13.0
import { Science, Assignment, Cancel, ArrowBack, Visibility } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Dialog from '../../../components/common/Dialog';
import Badge from '../../../components/common/Badge';
import SubmissionTimeline from './SubmissionTimeline';
import QuoteApproval from './QuoteApproval';
import useSubmissions from '../hooks/useSubmissions';
import useToast from '../../../hooks/useToast';
import { Submission, SubmissionDetailed, SubmissionStatus } from '../../../types/submission';
import useAuth from '../../auth/hooks/useAuth';
import { formatDate, formatCurrency } from '../../../utils/formatters';

/**
 * Interface defining the props for the SubmissionDetail component.
 */
interface SubmissionDetailProps {
  /**
   * ID of the submission to display, can be provided directly or via URL params.
   */
  submissionId?: string | undefined;
  /**
   * Submission object if already fetched by parent component.
   */
  submission?: SubmissionDetailed | undefined;
  /**
   * Callback function when back button is clicked.
   */
  onBack?: () => void;
  /**
   * Callback function when submission status changes.
   */
  onStatusChange?: () => void;
}

/**
 * Styled component for the detail container.
 */
const DetailContainer = styled(Box)(({ theme }) => ({
  marginBottom: theme.spacing(3),
}));

/**
 * Styled component for the section title.
 */
const SectionTitle = styled(Typography)(({ theme }) => ({
  fontWeight: '600',
  marginBottom: theme.spacing(1),
}));

/**
 * Styled component for individual detail items.
 */
const DetailItem = styled(Box)(({ theme }) => ({
  display: 'flex',
  marginBottom: theme.spacing(1),
}));

/**
 * Styled component for detail labels.
 */
const DetailLabel = styled(Typography)(({ theme }) => ({
  fontWeight: '500',
  width: '180px',
  color: theme.palette.text.secondary,
}));

/**
 * Styled component for detail values.
 */
const DetailValue = styled(Typography)(({ theme }) => ({
  flex: '1',
}));

/**
 * Styled component for the action container.
 */
const ActionContainer = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  marginTop: theme.spacing(3),
}));

/**
 * Component that displays detailed information about a submission to a CRO.
 */
const SubmissionDetail: React.FC<SubmissionDetailProps> = ({
  submissionId,
  submission,
  onBack = () => {},
  onStatusChange = () => {},
}) => {
  // LD1: Initialize state for confirmation dialogs
  const [cancelDialogOpen, setCancelDialogOpen] = useState(false);
  const [cancelReason, setCancelReason] = useState('');

  // LD1: Initialize state for submission data, loading, and error
  const [submissionData, setSubmissionData] = useState<SubmissionDetailed | undefined>(submission);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // LD1: Get submission ID from props or URL params
  const { id } = useParams<{ id: string }>();
  const submissionIdToUse = submissionId || id;

  // LD1: Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Get user authentication state and role from useAuth hook
  const { user, isAuthenticated } = useAuth();

  // LD1: Get submission operations from useSubmissions hook
  const { fetchSubmissionDetails, cancelSubmission, updateSubmissionStatus } = useSubmissions();

  // LD1: Get toast notification function from useToast hook
  const { showToast } = useToast();

  /**
   * LD1: Handles fetching submission details.
   */
  const handleFetchSubmission = useCallback(() => {
    if (!submissionIdToUse) return;

    setLoading(true);
    setError(null);

    fetchSubmissionDetails(submissionIdToUse)
      .then((data) => {
        setSubmissionData(data);
      })
      .catch((err) => {
        setError(err?.message || 'Failed to fetch submission details.');
      })
      .finally(() => {
        setLoading(false);
      });
  }, [fetchSubmissionDetails, submissionIdToUse]);

  /**
   * LD1: Handles cancelling the submission.
   */
  const handleCancelSubmission = useCallback(() => {
    if (!submissionIdToUse) return;

    cancelSubmission(submissionIdToUse, cancelReason)
      .then(() => {
        setCancelDialogOpen(false);
        showToast({
          type: 'success',
          message: 'Submission cancelled successfully!',
        });
        onStatusChange();
      })
      .catch((err) => {
        showToast({
          type: 'error',
          message: err?.message || 'Failed to cancel submission.',
        });
      });
  }, [cancelReason, cancelSubmission, onStatusChange, showToast, submissionIdToUse]);

  /**
   * LD1: Handles opening the cancel dialog.
   */
  const handleOpenCancelDialog = useCallback(() => {
    setCancelDialogOpen(true);
  }, []);

  /**
   * LD1: Handles changes to the cancel reason.
   * @param event - The change event.
   */
  const handleCancelReasonChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setCancelReason(event.target.value);
  }, []);

  /**
   * LD1: Handles quote responded.
   */
  const handleQuoteResponded = useCallback(() => {
    handleFetchSubmission();
    onStatusChange();
  }, [handleFetchSubmission, onStatusChange]);

  /**
   * LD1: Handles view results.
   */
  const handleViewResults = useCallback(() => {
    if (!submissionIdToUse) return;
    navigate(`/results?submissionId=${submissionIdToUse}`);
  }, [navigate, submissionIdToUse]);

  /**
   * LD1: Handles back.
   */
  const handleBack = useCallback(() => {
    onBack();
    navigate('/submissions');
  }, [navigate, onBack]);

  // LD1: Fetch submission details when component mounts or ID changes
  useEffect(() => {
    if (submissionIdToUse && !submission) {
      handleFetchSubmission();
    }
  }, [handleFetchSubmission, submission, submissionIdToUse]);

  // LD1: Determine available actions based on submission status and user role
  const isPharma = user?.role === 'pharma';
  const isCRO = user?.role === 'cro';

  const showApproveRejectQuote = submissionData?.status === SubmissionStatus.QUOTE_PROVIDED && isPharma;
  const showCancelButton = submissionData?.status === SubmissionStatus.PENDING || submissionData?.status === SubmissionStatus.APPROVED;
  const showViewResultsButton = submissionData?.status === SubmissionStatus.COMPLETED;

  // LD1: Render loading state if submission is being fetched
  if (loading) {
    return <Typography>Loading submission details...</Typography>;
  }

  // LD1: Render error state if submission fetch failed
  if (error) {
    return <Alert severity="error">{error}</Alert>;
  }

  // LD1: Render submission details including experiment info, CRO info, status, and timeline
  return (
    <Card>
      <DetailContainer>
        <Button onClick={handleBack} startIcon={<ArrowBack />}>
          Back to Submissions
        </Button>
        {submissionData && (
          <>
            <SectionTitle variant="h5" component="h2">
              {submissionData.experiment?.name} - {submissionData.experiment?.type?.name}
              <Badge label={submissionData.status} color="primary" style={{ marginLeft: '8px' }} />
            </SectionTitle>

            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <SectionTitle variant="h6" component="h3">
                  Experiment Details
                </SectionTitle>
                <DetailItem>
                  <DetailLabel>Experiment Type:</DetailLabel>
                  <DetailValue>{submissionData.experiment?.type?.name}</DetailValue>
                </DetailItem>
                <DetailItem>
                  <DetailLabel>Created By:</DetailLabel>
                  <DetailValue>{submissionData.experiment?.creator?.email}</DetailValue>
                </DetailItem>
              </Grid>

              <Grid item xs={12} md={6}>
                <SectionTitle variant="h6" component="h3">
                  CRO Details
                </SectionTitle>
                <DetailItem>
                  <DetailLabel>CRO:</DetailLabel>
                  <DetailValue>{submissionData.cro?.email}</DetailValue>
                </DetailItem>
              </Grid>
            </Grid>

            <SectionTitle variant="h6" component="h3">
              Submission Details
            </SectionTitle>
            <DetailItem>
              <DetailLabel>Submitted At:</DetailLabel>
              <DetailValue>{formatDate(submissionData.submitted_at)}</DetailValue>
            </DetailItem>
            {submissionData.price && (
              <DetailItem>
                <DetailLabel>Price:</DetailLabel>
                <DetailValue>{formatCurrency(submissionData.price)}</DetailValue>
              </DetailItem>
            )}
            {submissionData.notes && (
              <DetailItem>
                <DetailLabel>Notes:</DetailLabel>
                <DetailValue>{submissionData.notes}</DetailValue>
              </DetailItem>
            )}

            <SectionTitle variant="h6" component="h3">
              Timeline
            </SectionTitle>
            <SubmissionTimeline submission={submissionData} />

            {showApproveRejectQuote && (
              <QuoteApproval submission={submissionData} onQuoteResponded={handleQuoteResponded} />
            )}

            <ActionContainer>
              {showCancelButton && (
                <Button
                  variant="outlined"
                  color="error"
                  onClick={handleOpenCancelDialog}
                  startIcon={<Cancel />}
                >
                  Cancel Submission
                </Button>
              )}
              {showViewResultsButton && (
                <Button
                  variant="contained"
                  color="primary"
                  onClick={handleViewResults}
                  startIcon={<Visibility />}
                >
                  View Results
                </Button>
              )}
            </ActionContainer>
          </>
        )}
      </DetailContainer>

      <Dialog
        open={cancelDialogOpen}
        onClose={() => setCancelDialogOpen(false)}
        title="Cancel Submission"
        contentText="Please provide a reason for cancelling the submission:"
        confirmButtonText="Cancel Submission"
        cancelButtonText="Keep Submission"
        onConfirm={handleCancelSubmission}
        onCancel={() => setCancelDialogOpen(false)}
      >
        <TextField
          label="Cancellation Reason"
          multiline
          rows={4}
          fullWidth
          variant="outlined"
          value={cancelReason}
          onChange={handleCancelReasonChange}
        />
      </Dialog>
    </Card>
  );
};

export default SubmissionDetail;