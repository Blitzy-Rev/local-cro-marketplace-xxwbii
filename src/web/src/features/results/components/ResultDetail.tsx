# src/web/src/features/results/components/ResultDetail.tsx
```typescript
import React, { useState, useEffect, useCallback } from 'react'; // react ^18.2.0
import { useParams, useNavigate } from 'react-router-dom'; // react-router-dom ^6.11.2
import {
  Box,
  Grid,
  Typography,
  Divider,
  Tabs,
  Tab,
  CircularProgress,
  TextField
} from '@mui/material'; // @mui/material v5.13.0
import {
  CheckCircle,
  Cancel,
  CloudDownload,
  Assessment,
  Message
} from '@mui/icons-material'; // @mui/icons-material v5.11.16

import { ResultDetailed, ResultStatus } from '../../../types/result';
import useResults from '../hooks/useResults';
import ResultSummary from './ResultSummary';
import ResultFiles from './ResultFiles';
import MoleculeResults from './MoleculeResults';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import { useToast } from '../../../hooks/useToast';
import { useDialog } from '../../../hooks/useDialog';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Interface defining the props for the ResultDetail component.
 */
interface ResultDetailProps {
  /** ID of the result to display, optional as it can be obtained from URL params */
  resultId?: string;
  /** Optional callback when a result is approved */
  onApprove?: (resultId: string) => void;
  /** Optional callback when a result is rejected */
  onReject?: (resultId: string) => void;
  /** Optional callback for navigating back to results list */
  onBack?: () => void;
}

/**
 * A component that displays detailed information about an experimental result.
 * It includes result summary, associated files, and molecule-specific data.
 * This component serves as the main view for reviewing results uploaded by CRO users.
 *
 * @param {ResultDetailProps} props - The props for the component.
 * @returns {JSX.Element} The rendered component.
 */
const ResultDetail: React.FC<ResultDetailProps> = ({ resultId: propResultId, onApprove, onReject, onBack }) => {
  // Get result ID from URL params if not provided as prop
  const { id } = useParams<{ id: string }>();
  const resultId = propResultId || id;

  // React Router hook for navigation
  const navigate = useNavigate();

  // Custom hook for result operations
  const { fetchResultDetail, approveResult, rejectResult, exportResult, analyzeResult, requestAdditionalData, currentResult: result, loading, error, clearError } = useResults();

  // Hook for displaying toast notifications
  const { showToast } = useToast();

  // Hook for displaying dialog confirmations
  const { showDialog } = useDialog();

  // Hook for checking user permissions
  const { canApproveResults } = usePermissions();

  // State for managing the active tab
  const [activeTab, setActiveTab] = useState(0);

  // State for managing approval/rejection notes
  const [approvalNotes, setApprovalNotes] = useState('');

  // State for managing loading state during approval/rejection actions
  const [actionLoading, setActionLoading] = useState(false);

  /**
   * Handles changing the active tab.
   *
   * @param {React.SyntheticEvent} event - The event object.
   * @param {number} newValue - The index of the new tab.
   */
  const handleTabChange = (event: React.SyntheticEvent, newValue: number) => {
    setActiveTab(newValue);
  };

  /**
   * Handles approving the result.
   * Shows a confirmation dialog with notes input, then calls the approveResult function.
   */
  const handleApprove = useCallback(async () => {
    showDialog({
      title: 'Approve Result',
      message: 'Are you sure you want to approve this result? Please add any relevant notes.',
      confirmText: 'Approve',
      cancelText: 'Cancel',
      inputLabel: 'Approval Notes',
      onConfirm: async (notes: string) => {
        setActionLoading(true);
        try {
          await approveResult(resultId, notes);
          showToast({ type: 'success', message: 'Result approved successfully.' });
          if (onApprove) {
            onApprove(resultId);
          }
        } finally {
          setActionLoading(false);
          loadResultData(); // Refresh data after action
        }
      }
    });
  }, [resultId, approveResult, showToast, showDialog, onApprove, loadResultData]);

  /**
   * Handles rejecting the result.
   * Shows a confirmation dialog with notes input, then calls the rejectResult function.
   */
  const handleReject = useCallback(async () => {
    showDialog({
      title: 'Reject Result',
      message: 'Are you sure you want to reject this result? Please add a reason for rejection.',
      confirmText: 'Reject',
      cancelText: 'Cancel',
      inputLabel: 'Rejection Reason',
      onConfirm: async (notes: string) => {
        setActionLoading(true);
        try {
          await rejectResult(resultId, notes);
          showToast({ type: 'success', message: 'Result rejected successfully.' });
          if (onReject) {
            onReject(resultId);
          }
        } finally {
          setActionLoading(false);
          loadResultData(); // Refresh data after action
        }
      }
    });
  }, [resultId, rejectResult, showToast, showDialog, onReject, loadResultData]);

  /**
   * Handles exporting the result data.
   * Calls the exportResult function and shows a success or error toast notification.
   * @param {string} format - The format to export the result data in.
   */
  const handleExport = useCallback(async (format: string) => {
    setActionLoading(true);
    try {
      await exportResult(resultId, format as 'csv' | 'xlsx' | 'json');
      showToast({ type: 'success', message: `Result data exported successfully as ${format}.` });
    } catch (error) {
      console.error('Export error:', error);
      showToast({ type: 'error', message: `Failed to export result data as ${format}.` });
    } finally {
      setActionLoading(false);
    }
  }, [resultId, exportResult, showToast]);

  /**
   * Handles analyzing the result data.
   * Calls the analyzeResult function and shows the analysis results in a dialog.
   */
  const handleAnalyze = useCallback(async () => {
    setActionLoading(true);
    try {
      await analyzeResult(resultId);
      showToast({ type: 'success', message: 'Result data analyzed successfully.' });
    } catch (error) {
      console.error('Analyze error:', error);
      showToast({ type: 'error', message: 'Failed to analyze result data.' });
    } finally {
      setActionLoading(false);
    }
  }, [resultId, analyzeResult, showToast]);

  /**
   * Handles requesting additional data from the CRO.
   * Shows a dialog with notes input, then calls the requestAdditionalData function.
   */
  const handleRequestAdditionalData = useCallback(async () => {
    showDialog({
      title: 'Request Additional Data',
      message: 'Please specify what additional data you require from the CRO.',
      confirmText: 'Send Request',
      cancelText: 'Cancel',
      inputLabel: 'Request Details',
      onConfirm: async (notes: string) => {
        setActionLoading(true);
        try {
          await requestAdditionalData(resultId, notes);
          showToast({ type: 'success', message: 'Request for additional data sent successfully.' });
        } finally {
          setActionLoading(false);
        }
      }
    });
  }, [resultId, requestAdditionalData, showToast, showDialog]);

  /**
   * Handles navigating back to the results list.
   * Uses the onBack prop if provided, otherwise uses the navigate function.
   */
  const handleBack = useCallback(() => {
    if (onBack) {
      onBack();
    } else {
      navigate(-1); // Go back to the previous page
    }
  }, [navigate, onBack]);

  /**
   * Loads the result data.
   * Calls the fetchResultDetail function to load the result data.
   */
  const loadResultData = useCallback(async () => {
    setLoading(true);
    setError(null);
    try {
      if (resultId) {
        await fetchResultDetail(resultId);
      }
    } catch (e: any) {
      setError(e?.message || 'Failed to load result data.');
      showToast({
        type: 'error',
        message: e?.message || 'Failed to load result data.'
      });
    } finally {
      setLoading(false);
    }
  }, [fetchResultDetail, resultId, showToast]);

  // Load result data on component mount and when resultId changes
  useEffect(() => {
    loadResultData();
  }, [resultId, loadResultData]);

  return (
    <Box>
      {loading && (
        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', minHeight: 200 }}>
          <CircularProgress />
        </Box>
      )}

      {error && (
        <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', flexDirection: 'column', minHeight: 200 }}>
          <Typography color="error">{error}</Typography>
          <Button onClick={loadResultData}>Retry</Button>
        </Box>
      )}

      {result && (
        <>
          {/* Header section with back button, title, and action buttons */}
          <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
            <Button variant="outlined" onClick={handleBack}>
              Back to Results
            </Button>
            <Typography variant="h5">Result Details</Typography>
            <Box>
              {canApproveResults() && result.status !== ResultStatus.APPROVED && (
                <Button
                  variant="contained"
                  color="success"
                  onClick={handleApprove}
                  disabled={actionLoading}
                  startIcon={<CheckCircle />}
                  sx={{ mr: 1 }}
                >
                  Approve
                </Button>
              )}
              {canApproveResults() && result.status !== ResultStatus.REJECTED && (
                <Button
                  variant="contained"
                  color="error"
                  onClick={handleReject}
                  disabled={actionLoading}
                  startIcon={<Cancel />}
                  sx={{ mr: 1 }}
                >
                  Reject
                </Button>
              )}
              <Button
                variant="outlined"
                onClick={() => handleExport('csv')}
                disabled={actionLoading}
                startIcon={<CloudDownload />}
                sx={{ mr: 1 }}
              >
                Export CSV
              </Button>
              <Button
                variant="outlined"
                onClick={handleAnalyze}
                disabled={actionLoading}
                startIcon={<Assessment />}
              >
                Analyze
              </Button>
              <Button
                variant="outlined"
                onClick={handleRequestAdditionalData}
                disabled={actionLoading}
                startIcon={<Message />}
              >
                Request Data
              </Button>
            </Box>
          </Box>

          {/* ResultSummary component with result data */}
          <ResultSummary result={result} />

          {/* Tabs component with 'Files', 'Molecule Results', and 'Notes' tabs */}
          <Tabs value={activeTab} onChange={handleTabChange} aria-label="result details tabs">
            <Tab label="Files" />
            <Tab label="Molecule Results" />
            <Tab label="Notes" />
          </Tabs>

          {/* Tab panels that display different aspects of the result */}
          {activeTab === 0 && (
            <Box sx={{ mt: 3 }}>
              <ResultFiles files={result.files || []} resultId={result.id} />
            </Box>
          )}
          {activeTab === 1 && (
            <Box sx={{ mt: 3 }}>
              <MoleculeResults result={result} />
            </Box>
          )}
          {activeTab === 2 && (
            <Box sx={{ mt: 3 }}>
              <Card>
                <Box sx={{ p: 2 }}>
                  <Typography variant="subtitle1">Notes</Typography>
                  <Typography variant="body2">{result.notes || 'No notes available'}</Typography>
                </Box>
              </Card>
            </Box>
          )}
        </>
      )}
    </Box>
  );
};

export default ResultDetail;