import React, { useState, useCallback } from 'react'; // React 18.2+
import {
  Box,
  Typography,
  Paper,
  Divider,
  TextField,
  Stack,
  Alert,
} from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import {
  CheckCircle,
  Cancel,
  AttachMoney,
  Schedule,
} from '@mui/icons-material'; // v5.13+
import Button from '../../../components/common/Button';
import Dialog from '../../../components/common/Dialog';
import useSubmissions from '../hooks/useSubmissions';
import useToast from '../../../hooks/useToast';
import { Submission, QuoteResponse } from '../../../types/submission';
import { formatCurrency } from '../../../utils/formatters';

/**
 * Interface defining the props for the QuoteApproval component.
 */
interface QuoteApprovalProps {
  /**
   * The submission object containing quote details.
   */
  submission: Submission;
  /**
   * Callback function to be executed when the quote is approved or rejected.
   */
  onQuoteResponded: () => void;
}

/**
 * Styled Paper component for the quote container.
 */
const QuoteContainer = styled(Paper)(({ theme }) => ({
  padding: theme.spacing(3),
  marginBottom: theme.spacing(3),
  border: `1px solid ${theme.palette.divider}`,
  borderRadius: '8px',
}));

/**
 * Styled Box component for the quote header.
 */
const QuoteHeader = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  marginBottom: theme.spacing(2),
}));

/**
 * Styled Box component for the quote details.
 */
const QuoteDetails = styled(Box)(({ theme }) => ({
  marginBottom: theme.spacing(2),
}));

/**
 * Styled Box component for the quote actions.
 */
const QuoteActions = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'flex-end',
  gap: theme.spacing(2),
  marginTop: theme.spacing(2),
}));

/**
 * Styled Box component for individual detail items.
 */
const DetailItem = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  marginBottom: theme.spacing(1),
}));

/**
 * Styled Box component for detail icons.
 */
const DetailIcon = styled(Box)(({ theme }) => ({
  marginRight: theme.spacing(1),
  color: theme.palette.primary.main,
}));

/**
 * Component for displaying and responding to quotes from CROs.
 *
 * @param {QuoteApprovalProps} props - The props for the component.
 * @returns {JSX.Element} The rendered component.
 */
const QuoteApproval: React.FC<QuoteApprovalProps> = ({ submission, onQuoteResponded = () => {} }) => {
  // LD1: Extract submission data from props
  const { id, price, turnaround_days, notes } = submission;

  // LD1: Initialize state for rejection dialog and rejection notes
  const [rejectDialogOpen, setRejectDialogOpen] = useState(false);
  const [rejectionNotes, setRejectionNotes] = useState('');

  // LD1: Get respondToQuote function from useSubmissions hook
  const { respondToQuote } = useSubmissions();

  // LD1: Get showToast function from useToast hook
  const { showToast } = useToast();

  /**
   * LD1: Handles the approval of the quote.
   */
  const handleApproveQuote = useCallback(() => {
    if (!id) return;

    respondToQuote(id, { approved: true })
      .then(() => {
        showToast({
          type: 'success',
          message: 'Quote approved successfully!',
        });
        onQuoteResponded();
      })
      .catch((error: any) => {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to approve quote.',
        });
      });
  }, [id, respondToQuote, showToast, onQuoteResponded]);

  /**
   * LD1: Handles opening the rejection dialog.
   */
  const handleOpenRejectDialog = useCallback(() => {
    setRejectDialogOpen(true);
  }, []);

  /**
   * LD1: Handles the rejection of the quote.
   */
  const handleRejectQuote = useCallback(() => {
    if (!id) return;

    respondToQuote(id, { approved: false, notes: rejectionNotes })
      .then(() => {
        setRejectDialogOpen(false);
        showToast({
          type: 'success',
          message: 'Quote rejected successfully!',
        });
        onQuoteResponded();
      })
      .catch((error: any) => {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to reject quote.',
        });
      });
  }, [id, rejectionNotes, respondToQuote, showToast, onQuoteResponded]);

  /**
   * LD1: Handles changes to the rejection notes.
   *
   * @param {React.ChangeEvent<HTMLInputElement>} event - The change event.
   */
  const handleNotesChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setRejectionNotes(event.target.value);
  }, []);

  // LD1: Render quote details including price, turnaround time, and notes
  // LD1: Render approval and rejection buttons
  // LD1: Render rejection confirmation dialog when open
  return (
    <QuoteContainer>
      <QuoteHeader>
        <Typography variant="h6" component="div">
          Quote Details
        </Typography>
      </QuoteHeader>
      <QuoteDetails>
        <DetailItem>
          <DetailIcon>
            <AttachMoney />
          </DetailIcon>
          <Typography variant="body1">
            Price: {formatCurrency(price || 0)}
          </Typography>
        </DetailItem>
        <DetailItem>
          <DetailIcon>
            <Schedule />
          </DetailIcon>
          <Typography variant="body1">
            Turnaround Time: {turnaround_days} days
          </Typography>
        </DetailItem>
        {notes && (
          <DetailItem>
            <Typography variant="body2">
              CRO Notes: {notes}
            </Typography>
          </DetailItem>
        )}
      </QuoteDetails>
      <Divider />
      <QuoteActions>
        <Button
          variant="contained"
          color="success"
          startIcon={<CheckCircle />}
          onClick={handleApproveQuote}
        >
          Approve Quote
        </Button>
        <Button
          variant="outlined"
          color="error"
          startIcon={<Cancel />}
          onClick={handleOpenRejectDialog}
        >
          Reject Quote
        </Button>
      </QuoteActions>
      <Dialog
        open={rejectDialogOpen}
        onClose={() => setRejectDialogOpen(false)}
        title="Reject Quote"
        contentText="Please provide a reason for rejecting the quote:"
        confirmButtonText="Reject"
        cancelButtonText="Cancel"
        onConfirm={handleRejectQuote}
        onCancel={() => setRejectDialogOpen(false)}
      >
        <TextField
          label="Rejection Notes"
          multiline
          rows={4}
          fullWidth
          variant="outlined"
          value={rejectionNotes}
          onChange={handleNotesChange}
        />
      </Dialog>
    </QuoteContainer>
  );
};

export default QuoteApproval;