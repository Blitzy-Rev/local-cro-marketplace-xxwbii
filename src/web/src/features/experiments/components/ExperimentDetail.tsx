import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import { useParams, useNavigate, Link } from 'react-router-dom'; // react-router-dom v6.10.0
import { Box, Typography, Divider, Grid, Chip, Paper, Stack, Alert, Tabs, Tab } from '@mui/material'; // @mui/material v5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13.0
import { ArrowBack, Edit, Delete, PlayArrow, Send, Cancel, Science, Assignment } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Loading from '../../../components/common/Loading';
import Dialog from '../../../components/common/Dialog';
import ExperimentStatus from './ExperimentStatus';
import ExperimentParameters from './ExperimentParameters';
import MoleculeSelector from './MoleculeSelector';
import SubmissionDetail from '../../submissions/components/SubmissionDetail';
import useExperiments from '../hooks/useExperiments';
import usePermissions from '../../../hooks/usePermissions';
import useToast from '../../../hooks/useToast';
import { Experiment, ExperimentDetail as ExperimentDetailType, ExperimentStatus as ExperimentStatusEnum } from '../../../types/experiment';
import { formatDate } from '../../../utils/formatters';

// Interface defining the props for the ExperimentDetail component
interface ExperimentDetailProps {
  experimentId: string;
  onBack?: () => void;
  onStatusChange?: () => void;
  className?: string;
}

// Styled components for layout and visual enhancements
const DetailContainer = styled(Box)({
  width: '100%',
  marginBottom: theme.spacing(2),
});

const DetailHeader = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  marginBottom: theme.spacing(2),
});

const DetailSection = styled(Box)({
  marginBottom: theme.spacing(3),
});

const ActionButtons = styled(Box)({
  display: 'flex',
  justifyContent: 'flex-end',
  gap: theme.spacing(2),
  marginTop: theme.spacing(2),
});

const TabPanel = styled(Box)({
  padding: theme.spacing(3),
  backgroundColor: theme.palette.background.paper,
  borderRadius: '0 0 8px 8px',
  borderTop: 'none',
});

/**
 * Component that displays detailed information about an experiment
 */
const ExperimentDetail: React.FC<ExperimentDetailProps> = ({
  experimentId,
  onBack = () => {},
  onStatusChange = () => {},
  className = '',
}) => {
  // LD1: Get navigation function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Get URL parameters using useParams hook
  const { id } = useParams<{ id: string }>();

  // LD1: Determine the experiment ID to use (from props or URL params)
  const experimentIdToUse = experimentId || id;

  // LD1: Get experiment data and functions from useExperiments hook
  const {
    currentExperiment,
    loading,
    error,
    fetchExperimentById,
    updateExperiment,
    deleteExperiment,
    addMoleculesToExperiment,
    removeMoleculesFromExperiment,
    queueExperiment,
    submitExperiment,
    cancelExperiment,
  } = useExperiments();

  // LD1: Get permission checking functions from usePermissions hook
  const { isPharma, isCRO, can } = usePermissions();

  // LD1: Get toast notification function from useToast hook
  const { showToast } = useToast();

  // LD1: Set up state for active tab, edit mode, and confirmation dialogs
  const [activeTab, setActiveTab] = useState('0');
  const [editMode, setEditMode] = useState(false);
  const [deleteDialogOpen, setDeleteDialogOpen] = useState(false);
  const [cancelDialogOpen, setCancelDialogOpen] = useState(false);
  const [queueDialogOpen, setQueueDialogOpen] = useState(false);
  const [submitDialogOpen, setSubmitDialogOpen] = useState(false);
  const [selectedMoleculeIds, setSelectedMoleculeIds] = useState<string[]>([]);
  const [editData, setEditData] = useState({ name: '', parameters: {} });

  // LD1: Set up effect to fetch experiment data on component mount
  useEffect(() => {
    if (experimentIdToUse) {
      fetchExperimentById(experimentIdToUse);
    }
  }, [experimentIdToUse, fetchExperimentById]);

  // LD1: Set up effect to update edit data and selected molecule IDs when experiment data changes
  useEffect(() => {
    if (currentExperiment) {
      setEditData({
        name: currentExperiment.name,
        parameters: currentExperiment.parameters?.reduce((obj, item) => {
          obj[item.parameter_name] = item.parameter_value;
          return obj;
        }, {}) || {},
      });
      setSelectedMoleculeIds(currentExperiment.molecules?.map(m => m.molecule_id) || []);
    }
  }, [currentExperiment]);

  // LD1: Create handler for back button click
  const handleBack = useCallback(() => {
    onBack();
    navigate('/experiments');
  }, [navigate, onBack]);

  // LD1: Create handler for tab change
  const handleTabChange = useCallback((event: React.SyntheticEvent, newValue: string) => {
    setActiveTab(newValue);
  }, []);

  // LD1: Create handler for edit button click
  const handleEdit = useCallback(() => {
    setEditMode(true);
  }, []);

  // LD1: Create handler for cancel edit button click
  const handleCancelEdit = useCallback(() => {
    setEditMode(false);
    setEditData({
      name: currentExperiment?.name || '',
      parameters: currentExperiment?.parameters?.reduce((obj, item) => {
        obj[item.parameter_name] = item.parameter_value;
        return obj;
      }, {}) || {},
    });
  }, [currentExperiment]);

  // LD1: Create handler for save button click
  const handleSave = useCallback(() => {
    if (currentExperiment) {
      updateExperiment(currentExperiment.id, {
        name: editData.name,
        parameters: Object.entries(editData.parameters).map(([key, value]) => ({
          parameter_name: key,
          parameter_value: value,
        })),
      })
        .then(() => {
          setEditMode(false);
          showToast({ type: 'success', message: 'Experiment updated successfully!' });
          onStatusChange();
        })
        .catch((err) => {
          showToast({ type: 'error', message: err?.message || 'Failed to update experiment.' });
        });
    }
  }, [currentExperiment, editData, updateExperiment, showToast, onStatusChange]);

  // LD1: Create handler for delete button click
  const handleDeleteClick = useCallback(() => {
    setDeleteDialogOpen(true);
  }, []);

  // LD1: Create handler for confirm delete button click
  const handleConfirmDelete = useCallback(() => {
    if (currentExperiment) {
      deleteExperiment(currentExperiment.id)
        .then(() => {
          setDeleteDialogOpen(false);
          showToast({ type: 'success', message: 'Experiment deleted successfully!' });
          navigate('/experiments');
        })
        .catch((err) => {
          showToast({ type: 'error', message: err?.message || 'Failed to delete experiment.' });
        });
    }
  }, [currentExperiment, deleteExperiment, navigate, showToast]);

  // LD1: Create handler for queue button click
  const handleQueueClick = useCallback(() => {
    setQueueDialogOpen(true);
  }, []);

  // LD1: Create handler for confirm queue button click
  const handleConfirmQueue = useCallback(() => {
    if (currentExperiment) {
      queueExperiment(currentExperiment.id)
        .then(() => {
          setQueueDialogOpen(false);
          showToast({ type: 'success', message: 'Experiment queued successfully!' });
          onStatusChange();
        })
        .catch((err) => {
          showToast({ type: 'error', message: err?.message || 'Failed to queue experiment.' });
        });
    }
  }, [currentExperiment, queueExperiment, showToast, onStatusChange]);

  // LD1: Create handler for submit button click
  const handleSubmitClick = useCallback(() => {
    setSubmitDialogOpen(true);
  }, []);

  // LD1: Create handler for confirm submit button click
  const handleConfirmSubmit = useCallback(() => {
    if (currentExperiment) {
      submitExperiment(currentExperiment.id)
        .then(() => {
          setSubmitDialogOpen(false);
          showToast({ type: 'success', message: 'Experiment submitted successfully!' });
          onStatusChange();
        })
        .catch((err) => {
          showToast({ type: 'error', message: err?.message || 'Failed to submit experiment.' });
        });
    }
  }, [currentExperiment, submitExperiment, showToast, onStatusChange]);

  // LD1: Create handler for cancel button click
  const handleCancelClick = useCallback(() => {
    setCancelDialogOpen(true);
  }, []);

  // LD1: Create handler for confirm cancel button click
  const handleConfirmCancel = useCallback(() => {
    if (currentExperiment) {
      cancelExperiment(currentExperiment.id)
        .then(() => {
          setCancelDialogOpen(false);
          showToast({ type: 'success', message: 'Experiment cancelled successfully!' });
          onStatusChange();
        })
        .catch((err) => {
          showToast({ type: 'error', message: err?.message || 'Failed to cancel experiment.' });
        });
    }
  }, [currentExperiment, cancelExperiment, showToast, onStatusChange]);

  // LD1: Create handler for molecule selection change
  const handleMoleculeSelectionChange = useCallback(
    (moleculeIds: string[]) => {
      if (currentExperiment) {
        const moleculesToAdd = moleculeIds.filter(id => !currentExperiment.molecules.find(m => m.molecule_id === id));
        const moleculesToRemove = currentExperiment.molecules.filter(m => !moleculeIds.includes(m.molecule_id)).map(m => m.molecule_id);

        if (moleculesToAdd.length > 0) {
          addMoleculesToExperiment(currentExperiment.id, moleculesToAdd)
            .then(() => {
              showToast({ type: 'success', message: 'Molecules added to experiment successfully!' });
              onStatusChange();
            })
            .catch((err) => {
              showToast({ type: 'error', message: err?.message || 'Failed to add molecules to experiment.' });
            });
        }

        if (moleculesToRemove.length > 0) {
          removeMoleculesFromExperiment(currentExperiment.id, moleculesToRemove)
            .then(() => {
              showToast({ type: 'success', message: 'Molecules removed from experiment successfully!' });
              onStatusChange();
            })
            .catch((err) => {
              showToast({ type: 'error', message: err?.message || 'Failed to remove molecules from experiment.' });
            });
        }
      }
    },
    [currentExperiment, addMoleculesToExperiment, removeMoleculesFromExperiment, showToast, onStatusChange]
  );

  // LD1: Create handler for edit data change
  const handleEditDataChange = useCallback((field: string, value: any) => {
    setEditData(prev => ({ ...prev, [field]: value }));
  }, []);

  // LD1: Create handler for submission status change
  const handleSubmissionStatusChange = useCallback(() => {
    fetchExperimentById(experimentIdToUse);
    onStatusChange();
  }, [fetchExperimentById, onStatusChange, experimentIdToUse]);

  // LD1: Render loading state if data is being fetched
  if (loading) {
    return <Loading message="Loading experiment details..." />;
  }

  // LD1: Render error state if there was an error fetching data
  if (error) {
    return (
      <Alert severity="error">
        {error}
        <Button onClick={() => fetchExperimentById(experimentIdToUse)}>Retry</Button>
      </Alert>
    );
  }

  // LD1: Render experiment not found message if no experiment data is available
  if (!currentExperiment) {
    return <Typography>Experiment not found.</Typography>;
  }

  // LD1: Determine if the current user can edit the experiment
  const canEditExperiment = (): boolean => {
    return isPharma() && currentExperiment.status === ExperimentStatusEnum.DRAFT;
  };

  // LD1: Determine if the current user can queue the experiment
  const canQueueExperiment = (): boolean => {
    return isPharma() && currentExperiment.status === ExperimentStatusEnum.DRAFT && currentExperiment.molecules.length > 0;
  };

  // LD1: Determine if the current user can submit the experiment to a CRO
  const canSubmitExperiment = (): boolean => {
    return isPharma() && currentExperiment.status === ExperimentStatusEnum.QUEUED;
  };

  // LD1: Determine if the current user can cancel the experiment
  const canCancelExperiment = (): boolean => {
    return isPharma() && can('cancel_experiment') && currentExperiment.status !== ExperimentStatusEnum.COMPLETED && currentExperiment.status !== ExperimentStatusEnum.CANCELLED && currentExperiment.status !== ExperimentStatusEnum.RESULTS_AVAILABLE;
  };

  // LD1: Determine if the current user can delete the experiment
  const canDeleteExperiment = (): boolean => {
    return isPharma() && can('delete_experiment') && currentExperiment.status === ExperimentStatusEnum.DRAFT;
  };

  // LD1: Render experiment details with header, status, and action buttons
  return (
    <DetailContainer className={className}>
      <DetailHeader>
        <Button onClick={handleBack} startIcon={<ArrowBack />}>
          Back to Experiments
        </Button>
        <Box display="flex" alignItems="center">
          <Typography variant="h5" component="h2" sx={{ mr: 2 }}>
            {currentExperiment.name}
          </Typography>
          <ExperimentStatus status={currentExperiment.status} />
        </Box>
      </DetailHeader>

      <ActionButtons>
        {canEditExperiment() && (
          <Button variant="contained" startIcon={<Edit />} onClick={handleEdit}>
            Edit
          </Button>
        )}
        {canDeleteExperiment() && (
          <Button variant="outlined" color="error" startIcon={<Delete />} onClick={handleDeleteClick}>
            Delete
          </Button>
        )}
        {canQueueExperiment() && (
          <Button variant="contained" color="info" startIcon={<PlayArrow />} onClick={handleQueueClick}>
            Queue
          </Button>
        )}
        {canSubmitExperiment() && (
          <Button variant="contained" color="primary" startIcon={<Send />} onClick={handleSubmitClick}>
            Submit to CRO
          </Button>
        )}
        {canCancelExperiment() && (
          <Button variant="outlined" color="error" startIcon={<Cancel />} onClick={handleCancelClick}>
            Cancel Experiment
          </Button>
        )}
      </ActionButtons>

      <Tabs value={activeTab} onChange={handleTabChange} aria-label="experiment details tabs">
        <Tab label="Details" value="0" />
        <Tab label="Molecules" value="1" />
        <Tab label="Submissions" value="2" />
      </Tabs>

      {activeTab === '0' && (
        <TabPanel>
          <DetailSection>
            <Typography variant="h6" component="h3">
              Experiment Information
            </Typography>
            <Divider sx={{ mb: 2 }} />
            <Grid container spacing={3}>
              <Grid item xs={12} md={6}>
                <Typography variant="subtitle1">Name:</Typography>
                <Typography variant="body1">{currentExperiment.name}</Typography>
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="subtitle1">Type:</Typography>
                <Typography variant="body1">{currentExperiment.type?.name}</Typography>
              </Grid>
              <Grid item xs={12}>
                <Typography variant="subtitle1">Parameters:</Typography>
                <ExperimentParameters
                  experimentType={currentExperiment.type}
                  parameters={editData.parameters}
                  onChange={(newParams) => handleEditDataChange('parameters', newParams)}
                  disabled={!editMode}
                />
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="subtitle1">Created At:</Typography>
                <Typography variant="body1">{formatDate(currentExperiment.created_at)}</Typography>
              </Grid>
              <Grid item xs={12} md={6}>
                <Typography variant="subtitle1">Status:</Typography>
                <Typography variant="body1">{currentExperiment.status}</Typography>
              </Grid>
            </Grid>
          </DetailSection>
        </TabPanel>
      )}

      {activeTab === '1' && (
        <TabPanel>
          {editMode ? (
            <MoleculeSelector
              selectedMoleculeIds={selectedMoleculeIds}
              onSelectionChange={handleMoleculeSelectionChange}
            />
          ) : (
            <Stack direction="column" spacing={1}>
              {currentExperiment.molecules.length > 0 ? (
                currentExperiment.molecules.map((molecule) => (
                  <Chip key={molecule.molecule_id} label={molecule.molecule_id} />
                ))
              ) : (
                <Typography>No molecules added to this experiment.</Typography>
              )}
            </Stack>
          )}
        </TabPanel>
      )}

      {activeTab === '2' && (
        <TabPanel>
          {currentExperiment.submissions && currentExperiment.submissions.length > 0 ? (
            <SubmissionDetail submissionId={currentExperiment.submissions[0].id} onStatusChange={handleSubmissionStatusChange} />
          ) : (
            <Typography>No submissions yet.</Typography>
          )}
        </TabPanel>
      )}

      <Dialog
        open={deleteDialogOpen}
        onClose={() => setDeleteDialogOpen(false)}
        title="Delete Experiment"
        contentText="Are you sure you want to delete this experiment? This action cannot be undone."
        confirmButtonText="Delete"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmDelete}
        onCancel={() => setDeleteDialogOpen(false)}
      />

      <Dialog
        open={cancelDialogOpen}
        onClose={() => setCancelDialogOpen(false)}
        title="Cancel Experiment"
        contentText="Are you sure you want to cancel this experiment?"
        confirmButtonText="Cancel"
        cancelButtonText="Keep Experiment"
        onConfirm={handleConfirmCancel}
        onCancel={() => setCancelDialogOpen(false)}
      />

      <Dialog
        open={queueDialogOpen}
        onClose={() => setQueueDialogOpen(false)}
        title="Queue Experiment"
        contentText="Are you sure you want to queue this experiment?"
        confirmButtonText="Queue"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmQueue}
        onCancel={() => setQueueDialogOpen(false)}
      />

      <Dialog
        open={submitDialogOpen}
        onClose={() => setSubmitDialogOpen(false)}
        title="Submit Experiment"
        contentText="Are you sure you want to submit this experiment to a CRO?"
        confirmButtonText="Submit"
        cancelButtonText="Cancel"
        onConfirm={handleConfirmSubmit}
        onCancel={() => setSubmitDialogOpen(false)}
      />
    </DetailContainer>
  );
};

export default ExperimentDetail;