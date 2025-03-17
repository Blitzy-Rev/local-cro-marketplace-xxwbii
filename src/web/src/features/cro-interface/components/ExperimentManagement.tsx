import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import { Box, Typography, Paper, Grid, Divider, Chip, Slider, TextField, Dialog, DialogTitle, DialogContent, DialogActions } from '@mui/material'; // v5.13+
import { Science, Assignment, CheckCircle, Update, CloudUpload } from '@mui/icons-material'; // v5.13+

import useCROInterface from '../hooks/useCROInterface';
import ResultUploadForm from './ResultUploadForm';
import Table from '../../../components/common/Table';
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Input from '../../../components/common/Input';
import Loading from '../../../components/common/Loading';
import useToast from '../../../hooks/useToast';
import useDialog from '../../../hooks/useDialog';
import { Submission, SubmissionStatus } from '../../../types/submission';
import { ExperimentStatus } from '../../../types/experiment';

/**
 * Interface for the status update form data
 */
interface StatusUpdateFormData {
  completion: number;
  notes: string;
}

/**
 * A component for CRO users to manage in-progress experiments and update their status
 * @returns Rendered component
 */
const ExperimentManagement: React.FC = () => {
  // LD1: Initialize useCROInterface hook to get submissions and management functions
  const { submissions, loading, error, filters, setFilters, refreshSubmissions, handleUpdateStatus, handleUploadResult } = useCROInterface();

  // LD1: Initialize useToast hook for displaying notifications
  const { showToast } = useToast();

  // LD1: Initialize state for selected experiment, status update dialog, and result upload dialog
  const [selectedExperiment, setSelectedExperiment] = useState<Submission | null>(null);
  const [statusDialogOpen, setStatusDialogOpen] = useState<boolean>(false);
  const [resultDialogOpen, setResultDialogOpen] = useState<boolean>(false);

  // LD1: Initialize state for status update form data (completion percentage and notes)
  const [statusFormData, setStatusFormData] = useState<StatusUpdateFormData>({ completion: 0, notes: '' });

  // LD1: Set up filter to show only IN_PROGRESS submissions by default
  useEffect(() => {
    setFilters({ status: ExperimentStatus.IN_PROGRESS });
  }, [setFilters]);

  // LD1: Define table columns for experiment display with appropriate fields and formatting
  const columns = useMemo(() => [
    {
      field: 'experiment.name',
      headerName: 'Experiment Name',
      width: 200,
      valueGetter: (params: { row: Submission }) => params.row.experiment?.name,
    },
    {
      field: 'experiment.type.name',
      headerName: 'Experiment Type',
      width: 150,
      valueGetter: (params: { row: Submission }) => params.row.experiment?.type?.name,
    },
    {
      field: 'submitted_at',
      headerName: 'Submitted Date',
      width: 150,
      valueGetter: (params: { row: Submission }) => new Date(params.row.submitted_at).toLocaleDateString(),
    },
    {
      field: 'status',
      headerName: 'Status',
      width: 150,
      renderCell: (params: { row: Submission }) => (
        <Chip label={params.row.status} color="primary" />
      ),
    },
    {
      field: 'actions',
      headerName: 'Actions',
      width: 200,
      renderCell: (params: { row: Submission }) => (
        <>
          <Button
            startIcon={<Update />}
            onClick={() => handleStatusUpdate(params.row)}
            size="small"
          >
            Update Status
          </Button>
          <Button
            startIcon={<CloudUpload />}
            onClick={() => handleUploadResults(params.row)}
            size="small"
          >
            Upload Results
          </Button>
        </>
      ),
    },
  ], [handleStatusUpdate, handleUploadResults]);

  // LD1: Define handleStatusUpdate function to open status update dialog for an experiment
  const handleStatusUpdate = useCallback((experiment: Submission) => {
    setSelectedExperiment(experiment);
    const completion = experiment.experiment?.parameters?.find(p => p.parameter_name === 'completion')?.parameter_value || 0;
    setStatusFormData({ completion: Number(completion), notes: '' });
    setStatusDialogOpen(true);
  }, []);

  // LD1: Define handleStatusUpdateSubmit function to submit status updates
  const handleStatusUpdateSubmit = useCallback(async () => {
    if (!selectedExperiment) return;

    const statusUpdate = {
      notes: statusFormData.notes,
      status: statusFormData.completion === 100 ? ExperimentStatus.RESULTS_PENDING : ExperimentStatus.IN_PROGRESS,
    };

    try {
      await handleUpdateStatus(selectedExperiment.id, statusUpdate);
      showToast({ type: 'success', message: 'Status updated successfully!' });
      refreshSubmissions();
    } catch (error: any) {
      showToast({ type: 'error', message: error.message || 'Failed to update status' });
    } finally {
      closeDialogs();
    }
  }, [selectedExperiment, statusFormData, handleUpdateStatus, showToast, refreshSubmissions]);

  // LD1: Define handleUploadResults function to open result upload dialog for an experiment
  const handleUploadResults = useCallback((experiment: Submission) => {
    setSelectedExperiment(experiment);
    setResultDialogOpen(true);
  }, []);

  // LD1: Define handleResultUploadSuccess function to handle successful result uploads
  const handleResultUploadSuccess = useCallback(() => {
    showToast({ type: 'success', message: 'Results uploaded successfully!' });
    refreshSubmissions();
    closeDialogs();
  }, [showToast, refreshSubmissions]);

  // LD1: Define closeDialogs function to close all open dialogs
  const closeDialogs = useCallback(() => {
    setStatusDialogOpen(false);
    setResultDialogOpen(false);
    setSelectedExperiment(null);
  }, []);

  // Handler for completion slider change
  const handleCompletionChange = useCallback((event: Event, newValue: number | number[]) => {
    setStatusFormData(prev => ({ ...prev, completion: Number(newValue) }));
  }, []);

  // Handler for notes text field change
  const handleNotesChange = useCallback((event: React.ChangeEvent<HTMLInputElement>) => {
    setStatusFormData(prev => ({ ...prev, notes: event.target.value }));
  }, []);

  // LD1: Render loading state while data is being fetched
  if (loading) {
    return <Loading message="Loading in-progress experiments..." />;
  }

  // LD1: Render error message if data fetching fails
  if (error) {
    return (
      <Card>
        <Typography color="error">Error: {error}</Typography>
        <Button onClick={refreshSubmissions}>Retry</Button>
      </Card>
    );
  }

  // LD1: Render table of in-progress experiments with status indicators
  return (
    <Paper sx={{ p: 2, height: '100%' }}>
      <Typography variant="h6" gutterBottom>
        In-Progress Experiments
      </Typography>
      {submissions && submissions.length > 0 ? (
        <Table
          data={submissions.filter(sub => sub.status === ExperimentStatus.IN_PROGRESS)}
          columns={columns}
          loading={loading}
          error={error}
        />
      ) : (
        <Typography>No in-progress experiments</Typography>
      )}

      {/* LD1: Render status update dialog when open */}
      <Dialog open={statusDialogOpen} onClose={closeDialogs}>
        <DialogTitle>Update Experiment Status</DialogTitle>
        <DialogContent>
          <Typography variant="subtitle1">
            Experiment: {selectedExperiment?.experiment?.name}
          </Typography>
          <Typography variant="body2">
            Type: {selectedExperiment?.experiment?.type?.name}
          </Typography>
          <Divider sx={{ my: 2 }} />
          <Typography variant="body1">Completion:</Typography>
          <Slider
            value={statusFormData.completion}
            onChange={handleCompletionChange}
            valueLabelDisplay="auto"
            aria-labelledby="completion-slider"
          />
          <TextField
            label="Notes"
            multiline
            rows={4}
            fullWidth
            value={statusFormData.notes}
            onChange={handleNotesChange}
            sx={{ mt: 2 }}
          />
        </DialogContent>
        <DialogActions>
          <Button onClick={closeDialogs}>Cancel</Button>
          <Button onClick={handleStatusUpdateSubmit}>Update</Button>
        </DialogActions>
      </Dialog>

      {/* LD1: Render result upload dialog when open */}
      <Dialog open={resultDialogOpen} onClose={closeDialogs} fullWidth maxWidth="md">
        <DialogTitle>Upload Results</DialogTitle>
        <DialogContent>
          {selectedExperiment && (
            <ResultUploadForm
              submission={selectedExperiment}
              onSuccess={handleResultUploadSuccess}
              onCancel={closeDialogs}
            />
          )}
        </DialogContent>
      </Dialog>
    </Paper>
  );
};

export default ExperimentManagement;