import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom ^6.11.2
import {
  Box,
  Paper,
  Typography,
  Grid,
  Divider,
  MenuItem,
  FormControl,
  FormHelperText,
  CircularProgress,
} from '@mui/material'; // @mui/material ^5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles ^5.13.0

import { useExperiments } from '../hooks/useExperiments';
import MoleculeSelector from './MoleculeSelector';
import ExperimentParameters from './ExperimentParameters';
import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import { ExperimentCreate, ExperimentType } from '../../../types/experiment';
import useToast from '../../../hooks/useToast';

/**
 * Props interface for the CreateExperimentForm component
 */
interface CreateExperimentFormProps {
  /** Callback function called when experiment is successfully created */
  onSuccess: (experimentId: string) => void;
  /** Initial data for the form, if any */
  initialData?: Partial<ExperimentCreate>;
  /** Additional CSS class name */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
}

/**
 * Interface for the form's internal state
 */
interface FormData {
  /** Name of the experiment */
  name: string;
  /** ID of the selected experiment type */
  type_id: string;
  /** Experiment parameters as key-value pairs */
  parameters: Record<string, string>;
  /** IDs of selected molecules */
  molecule_ids: string[];
}

/**
 * Interface for form validation errors
 */
interface FormErrors {
  /** Error message for experiment name */
  name: string;
  /** Error message for experiment type */
  type_id: string;
  /** Error messages for parameters */
  parameters: Record<string, string>;
  /** Error message for molecule selection */
  molecule_ids: string;
}

// Define styled components for consistent styling
const FormContainer = styled(Paper)({
  padding: '24px',
  marginBottom: '24px',
  width: '100%',
});

const FormSection = styled(Box)({
  marginBottom: '24px',
});

const FormActions = styled(Box)({
  display: 'flex',
  justifyContent: 'flex-end',
  gap: '16px',
  marginTop: '32px',
});

/**
 * A form component for creating new experiments
 */
const CreateExperimentForm: React.FC<CreateExperimentFormProps> = ({
  onSuccess,
  initialData,
  className,
  style,
}) => {
  // LD1: Initialize state for form data
  const [formData, setFormData] = useState<FormData>({
    name: initialData?.name || '',
    type_id: initialData?.type_id || '',
    parameters: initialData?.parameters || {},
    molecule_ids: initialData?.molecule_ids || [],
  });

  // LD1: Initialize state for form validation errors
  const [errors, setErrors] = useState<FormErrors>({
    name: '',
    type_id: '',
    parameters: {},
    molecule_ids: '',
  });

  // LD1: Initialize state for form submission status
  const [submitting, setSubmitting] = useState(false);

  // LD1: Get experiment types, loading state, and createExperiment function from useExperiments hook
  const { experimentTypes, loading, fetchExperimentTypes, createExperiment } =
    useExperiments();

  // LD1: Get showToast function from useToast hook
  const { showToast } = useToast();

  // LD1: Get navigate function from useNavigate hook
  const navigate = useNavigate();

  // LD1: Initialize state for selected experiment type
  const [selectedExperimentType, setSelectedExperimentType] = useState<ExperimentType | null>(null);

  // LD1: Fetch experiment types on component mount
  useEffect(() => {
    fetchExperimentTypes();
  }, [fetchExperimentTypes]);

  // LD1: Update selectedExperimentType when type_id or experimentTypes change
  useEffect(() => {
    if (experimentTypes && formData.type_id) {
      const type = experimentTypes.find((t) => t.id === formData.type_id) || null;
      setSelectedExperimentType(type);
    }
  }, [experimentTypes, formData.type_id]);

  // LD1: Handle experiment type selection change
  const handleTypeChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const typeId = e.target.value;
    setFormData((prev) => ({
      ...prev,
      type_id: typeId,
      parameters: {}, // Reset parameters when type changes
    }));
    setErrors((prev) => ({
      ...prev,
      type_id: '',
      parameters: {},
    }));
  };

  // LD1: Handle experiment name change
  const handleNameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const name = e.target.value;
    setFormData((prev) => ({
      ...prev,
      name: name,
    }));
    setErrors((prev) => ({
      ...prev,
      name: '',
    }));
  };

  // LD1: Handle experiment parameters change
  const handleParametersChange = useCallback((parameters: Record<string, string>) => {
    setFormData((prev) => ({
      ...prev,
      parameters: parameters,
    }));
    setErrors((prev) => ({
      ...prev,
      parameters: {},
    }));
  }, []);

  // LD1: Handle molecule selection change
  const handleMoleculeSelectionChange = (moleculeIds: string[]) => {
    setFormData((prev) => ({
      ...prev,
      molecule_ids: moleculeIds,
    }));
    setErrors((prev) => ({
      ...prev,
      molecule_ids: '',
    }));
  };

  // LD1: Validate form data before submission
  const validateForm = (): boolean => {
    let isValid = true;
    const newErrors: FormErrors = {
      name: '',
      type_id: '',
      parameters: {},
      molecule_ids: '',
    };

    if (!formData.name) {
      newErrors.name = 'Experiment name is required';
      isValid = false;
    }

    if (!formData.type_id) {
      newErrors.type_id = 'Experiment type is required';
      isValid = false;
    }

    if (formData.molecule_ids.length === 0) {
      newErrors.molecule_ids = 'At least one molecule must be selected';
      isValid = false;
    }

    setErrors(newErrors);
    return isValid;
  };

  // LD1: Handle form submission
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    if (!validateForm()) {
      return;
    }

    setSubmitting(true);
    try {
      const experimentData: ExperimentCreate = {
        name: formData.name,
        type_id: formData.type_id,
        parameters: Object.entries(formData.parameters).map(([key, value]) => ({
          parameter_name: key,
          parameter_value: value,
        })),
        molecule_ids: formData.molecule_ids,
      };

      const newExperiment = await createExperiment(experimentData);
      if (newExperiment && newExperiment.id) {
        showToast({
          type: 'success',
          message: 'Experiment created successfully!',
        });
        onSuccess(newExperiment.id);
      } else {
        showToast({
          type: 'error',
          message: 'Failed to create experiment.',
        });
      }
    } catch (error: any) {
      showToast({
        type: 'error',
        message: error.message || 'Failed to create experiment.',
      });
    } finally {
      setSubmitting(false);
    }
  };

  // LD1: Handle cancel action
  const handleCancel = () => {
    navigate('/experiments');
  };

  // LD1: Render form with sections for experiment details, type selection, parameters, and molecule selection
  return (
    <FormContainer className={className} style={style}>
      <Typography variant="h5" gutterBottom>
        Create New Experiment
      </Typography>
      <Divider sx={{ mb: 2 }} />
      <form onSubmit={handleSubmit}>
        {/* LD1: Experiment Details Section */}
        <FormSection>
          <Typography variant="h6">Experiment Details</Typography>
          <Input
            label="Experiment Name"
            value={formData.name}
            onChange={handleNameChange}
            error={errors.name}
            fullWidth
            required
          />
        </FormSection>

        {/* LD1: Experiment Type Selection */}
        <FormSection>
          <Typography variant="h6">Experiment Type</Typography>
          {loading ? (
            <CircularProgress />
          ) : (
            <FormControl fullWidth error={!!errors.type_id}>
              <Input
                select
                label="Select Experiment Type"
                value={formData.type_id}
                onChange={handleTypeChange}
                required
                error={!!errors.type_id}
                helperText={errors.type_id}
              >
                {experimentTypes &&
                  experimentTypes.map((type) => (
                    <MenuItem key={type.id} value={type.id}>
                      {type.name}
                    </MenuItem>
                  ))}
              </Input>
              {errors.type_id && (
                <FormHelperText error>{errors.type_id}</FormHelperText>
              )}
            </FormControl>
          )}
        </FormSection>

        {/* LD1: Experiment Parameters Section */}
        {selectedExperimentType && (
          <FormSection>
            <ExperimentParameters
              experimentType={selectedExperimentType}
              parameters={formData.parameters}
              onChange={handleParametersChange}
              errors={errors.parameters}
              disabled={submitting}
            />
          </FormSection>
        )}

        {/* LD1: Molecule Selection Section */}
        <FormSection>
          <Typography variant="h6">Select Molecules</Typography>
          <MoleculeSelector
            selectedMoleculeIds={formData.molecule_ids}
            onSelectionChange={handleMoleculeSelectionChange}
            maxSelectionCount={10}
          />
          {errors.molecule_ids && (
            <FormHelperText error>{errors.molecule_ids}</FormHelperText>
          )}
        </FormSection>

        {/* LD1: Form Actions */}
        <FormActions>
          <Button onClick={handleCancel} disabled={submitting}>
            Cancel
          </Button>
          <Button type="submit" variant="contained" disabled={submitting}>
            {submitting ? <CircularProgress size={24} /> : 'Save'}
          </Button>
        </FormActions>
      </form>
    </FormContainer>
  );
};

export default CreateExperimentForm;