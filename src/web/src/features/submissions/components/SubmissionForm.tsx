import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import Input from '../../../components/common/Input'; // Reusable input component for form fields
import Dropdown from '../../../components/common/Dropdown'; // Reusable dropdown component for selecting experiment and CRO
import Button from '../../../components/common/Button'; // Reusable button component for form submission
import useSubmissions from '../hooks/useSubmissions'; // Custom hook for submission management
import useExperiments from '../../experiments/hooks/useExperiments'; // Custom hook for fetching experiment data
import useToast from '../../../hooks/useToast'; // Custom hook for displaying toast notifications
import {
  SubmissionCreate,
  SubmissionUpdate,
  Submission,
  SubmissionDetailCreate,
} from '../../../types/submission'; // Type definition for submission creation data
import { UserRole } from '../../../types/user'; // Enum for user roles to filter CRO users
import { getUsers } from '../../../api/admin'; // API function to fetch users for CRO selection
import {
  Box,
  Grid,
  Typography,
  Paper,
  Divider,
} from '@mui/material'; // Material-UI components v5.13+

/**
 * Interface for the props of the SubmissionForm component.
 * @property initialData - The initial submission data, if editing an existing submission.
 * @property experimentId - The ID of the experiment to create a submission for.
 * @property onSubmit - A callback function to be called when the form is submitted.
 * @property onCancel - A callback function to be called when the form is cancelled.
 */
interface SubmissionFormProps {
  initialData?: Submission | null;
  experimentId?: number | null;
  onSubmit: (submission: Submission) => void;
  onCancel: () => void;
}

/**
 * Interface for the form state.
 * @property experiment_id - The ID of the experiment.
 * @property cro_id - The ID of the CRO.
 * @property notes - Additional notes for the submission.
 * @property details - Submission details.
 */
interface FormState {
  experiment_id: number | null;
  cro_id: number | null;
  notes: string;
  details: SubmissionDetailCreate[];
}

/**
 * A form component for creating and submitting experiments to Contract Research Organizations (CROs).
 * This component allows users to select an experiment, choose a CRO, add submission details, and submit the form.
 * It supports both creation of new submissions and editing of existing ones.
 * @param props - The props for the component.
 * @returns A React component.
 */
const SubmissionForm: React.FC<SubmissionFormProps> = (props) => {
  // LD1: Initialize form state with initial data or default values
  const [formState, setFormState] = useState<FormState>({
    experiment_id: props.experimentId || props.initialData?.experiment_id || null,
    cro_id: props.initialData?.cro_id || null,
    notes: props.initialData?.notes || '',
    details: props.initialData?.details || [],
  });

  // LD1: Initialize state for validation errors
  const [errors, setErrors] = useState<Record<string, string>>({});

  // LD1: Initialize state for submission status
  const [isSubmitting, setIsSubmitting] = useState<boolean>(false);

  // LD1: Initialize state for CRO users
  const [croUsers, setCroUsers] = useState<Array<{ value: number; label: string }>>([]);

  // LD1: Initialize state for loading CRO users
  const [isLoadingCroUsers, setIsLoadingCroUsers] = useState<boolean>(true);

  // LD1: Access submission-related functionality using the useSubmissions hook
  const { createSubmission, updateSubmission } = useSubmissions();

  // LD1: Access experiment data using the useExperiments hook
  const { experiments, fetchExperiments } = useExperiments();

  // LD1: Access toast notification functionality using the useToast hook
  const { showToast } = useToast();

  // LD1: Initialize form state when props change
  useEffect(() => {
    setFormState({
      experiment_id: props.experimentId || props.initialData?.experiment_id || null,
      cro_id: props.initialData?.cro_id || null,
      notes: props.initialData?.notes || '',
      details: props.initialData?.details || [],
    });
  }, [props.experimentId, props.initialData]);

  // LD1: Fetch experiments and CRO users on component mount
  useEffect(() => {
    fetchExperiments();
    fetchCROUsers();
  }, [fetchExperiments]);

  // LD1: Fetch CRO users for the dropdown
  const fetchCROUsers = useCallback(async () => {
    setIsLoadingCroUsers(true);
    try {
      const response = await getUsers({ role: UserRole.CRO });
      setCroUsers(
        response.data.items.map((user) => ({
          value: user.id,
          label: user.email,
        }))
      );
    } catch (error: any) {
      showToast({
        type: 'error',
        message: error?.message || 'Failed to load CRO users',
      });
    } finally {
      setIsLoadingCroUsers(false);
    }
  }, [showToast]);

  // LD1: Handle changes to form input fields
  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const { name, value } = e.target;
    setFormState((prevState) => ({
      ...prevState,
      [name]: value,
    }));
    setErrors((prevErrors) => ({ ...prevErrors, [name]: '' }));
  };

  // LD1: Handle changes to dropdown selections
  const handleDropdownChange = (name: string, value: any) => {
    setFormState((prevState) => ({
      ...prevState,
      [name]: value,
    }));
    setErrors((prevErrors) => ({ ...prevErrors, [name]: '' }));
  };

  // LD1: Handle changes to submission details
  const handleDetailChange = (index: number, field: string, value: string) => {
    const updatedDetails = [...formState.details];
    updatedDetails[index] = {
      ...updatedDetails[index],
      [field]: value,
    };
    setFormState((prevState) => ({
      ...prevState,
      details: updatedDetails,
    }));
  };

  // LD1: Add a new empty detail to the form
  const addDetail = () => {
    setFormState((prevState) => ({
      ...prevState,
      details: [...prevState.details, { detail_name: '', detail_value: '' }],
    }));
  };

  // LD1: Remove a detail from the form
  const removeDetail = (index: number) => {
    const updatedDetails = [...formState.details];
    updatedDetails.splice(index, 1);
    setFormState((prevState) => ({
      ...prevState,
      details: updatedDetails,
    }));
  };

  // LD1: Validate the form data
  const validateForm = () => {
    let isValid = true;
    const newErrors: Record<string, string> = {};

    if (!formState.experiment_id) {
      newErrors.experiment_id = 'Experiment is required';
      isValid = false;
    }

    if (!formState.cro_id) {
      newErrors.cro_id = 'CRO is required';
      isValid = false;
    }

    formState.details.forEach((detail, index) => {
      if (!detail.detail_name) {
        newErrors[`detail_name_${index}`] = 'Detail name is required';
        isValid = false;
      }
      if (!detail.detail_value) {
        newErrors[`detail_value_${index}`] = 'Detail value is required';
        isValid = false;
      }
    });

    setErrors(newErrors);
    return isValid;
  };

  // LD1: Handle form submission
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    if (!validateForm()) {
      return;
    }

    setIsSubmitting(true);

    const submissionData: SubmissionCreate | SubmissionUpdate = {
      experiment_id: formState.experiment_id!.toString(),
      cro_id: formState.cro_id!,
      notes: formState.notes,
      details: formState.details,
    };

    try {
      if (props.initialData) {
        await updateSubmission(props.initialData.id, submissionData as SubmissionUpdate);
        showToast({ type: 'success', message: 'Submission updated successfully!' });
        props.onSubmit({ ...props.initialData, ...submissionData } as Submission);
      } else {
        const newSubmission = await createSubmission(submissionData as SubmissionCreate);
        showToast({ type: 'success', message: 'Submission created successfully!' });
        props.onSubmit(newSubmission);
      }
    } catch (error: any) {
      showToast({ type: 'error', message: error?.message || 'An error occurred' });
    } finally {
      setIsSubmitting(false);
    }
  };

  // LD1: Render the component
  return (
    <Paper elevation={3} style={{ padding: '20px', width: '100%' }}>
      <Typography variant="h6" gutterBottom>
        {props.initialData ? 'Edit Submission' : 'Create Submission'}
      </Typography>
      <Divider style={{ marginBottom: '20px' }} />
      <form onSubmit={handleSubmit}>
        <Grid container spacing={2}>
          <Grid item xs={12} md={6}>
            <Dropdown
              id="experiment_id"
              name="experiment_id"
              label="Experiment"
              value={formState.experiment_id || ''}
              onChange={(e) => handleDropdownChange('experiment_id', Number(e.target.value))}
              options={experiments.map((exp) => ({
                value: exp.id,
                label: exp.name,
              }))}
              error={errors.experiment_id}
              helperText={errors.experiment_id}
              required
              disabled={props.initialData !== undefined}
            />
          </Grid>
          <Grid item xs={12} md={6}>
            <Dropdown
              id="cro_id"
              name="cro_id"
              label="CRO"
              value={formState.cro_id || ''}
              onChange={(e) => handleDropdownChange('cro_id', Number(e.target.value))}
              options={croUsers.map((user) => ({
                value: user.value,
                label: user.label,
              }))}
              error={errors.cro_id}
              helperText={errors.cro_id}
              required
              disabled={isLoadingCroUsers}
            />
          </Grid>
          <Grid item xs={12}>
            <Input
              label="Notes"
              name="notes"
              value={formState.notes}
              onChange={handleInputChange}
              multiline
              rows={4}
              fullWidth
            />
          </Grid>
          <Grid item xs={12}>
            <Typography variant="subtitle1">Submission Details</Typography>
            {formState.details.map((detail, index) => (
              <Grid container spacing={2} key={index} alignItems="center">
                <Grid item xs={5}>
                  <Input
                    label={`Detail Name ${index + 1}`}
                    name={`detail_name_${index}`}
                    value={detail.detail_name}
                    onChange={(e) => handleDetailChange(index, 'detail_name', e.target.value)}
                    error={errors[`detail_name_${index}`]}
                    helperText={errors[`detail_name_${index}`]}
                    fullWidth
                  />
                </Grid>
                <Grid item xs={5}>
                  <Input
                    label={`Detail Value ${index + 1}`}
                    name={`detail_value_${index}`}
                    value={detail.detail_value}
                    onChange={(e) => handleDetailChange(index, 'detail_value', e.target.value)}
                    error={errors[`detail_value_${index}`]}
                    helperText={errors[`detail_value_${index}`]}
                    fullWidth
                  />
                </Grid>
                <Grid item xs={2} style={{ textAlign: 'center' }}>
                  <Button
                    type="button"
                    onClick={() => removeDetail(index)}
                    color="error"
                    size="small"
                  >
                    Remove
                  </Button>
                </Grid>
              </Grid>
            ))}
            <Button type="button" onClick={addDetail}>
              Add Detail
            </Button>
          </Grid>
          <Grid item xs={12} style={{ textAlign: 'right' }}>
            <Button onClick={props.onCancel} disabled={isSubmitting}>
              Cancel
            </Button>
            <Button type="submit" variant="contained" color="primary" disabled={isSubmitting}>
              {isSubmitting ? 'Submitting...' : 'Submit'}
            </Button>
          </Grid>
        </Grid>
      </form>
    </Paper>
  );
};

SubmissionForm.defaultProps = {
  initialData: null,
  experimentId: null,
};

export default SubmissionForm;