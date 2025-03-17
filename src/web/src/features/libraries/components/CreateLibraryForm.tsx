import React, { useState, useCallback } from 'react'; // React v18.2+
import { Box, Typography, Divider, Grid } from '@mui/material'; // @mui/material v5.13+
import { useForm, Controller } from 'react-hook-form'; // react-hook-form v7.43+
import { yupResolver } from '@hookform/resolvers/yup'; // @hookform/resolvers/yup v3.1+
import * as yup from 'yup'; // yup v1.1+

import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import Card from '../../../components/common/Card';
import MoleculeSelector from '../../experiments/components/MoleculeSelector';
import useLibraries from '../hooks/useLibraries';
import { LibraryCreate } from '../../../types/library';
import { Molecule } from '../../../types/molecule';
import useToast from '../../../hooks/useToast';

/**
 * Interface defining the properties for the CreateLibraryForm component
 */
interface CreateLibraryFormProps {
  /** Callback function called when library is successfully created */
  onSuccess: (libraryId: string) => void;
  /** Callback function called when form is cancelled */
  onCancel: () => void;
  /** Optional initial molecules to include in the library */
  initialMolecules?: Molecule[];
}

/**
 * Validation schema for the CreateLibraryForm using Yup
 */
const validationSchema = yup.object({
  name: yup.string().required('Library name is required').max(100, 'Name must be less than 100 characters'),
  description: yup.string().max(500, 'Description must be less than 500 characters')
});

/**
 * A form component for creating new molecule libraries
 */
const CreateLibraryForm: React.FC<CreateLibraryFormProps> = ({
  onSuccess,
  onCancel,
  initialMolecules = []
}) => {
  // LD1: Initialize state variables for managing selected molecules, molecule selector visibility, and submission status
  const [selectedMolecules, setSelectedMolecules] = useState<Molecule[]>(initialMolecules || []);
  const [isSelectorOpen, setIsSelectorOpen] = useState<boolean>(false);
  const [isSubmitting, setIsSubmitting] = useState<boolean>(false);

  // LD1: Use the useLibraries hook to access library management functions
  const { createLibrary } = useLibraries();

  // LD1: Use the useToast hook to display toast notifications
  const { showToast } = useToast();

  // LD1: Use the useForm hook to manage form state and validation
  const { control, handleSubmit, formState, reset } = useForm<LibraryCreate>({
    resolver: yupResolver(validationSchema),
    defaultValues: {
      name: '',
      description: ''
    }
  });

  // LD1: Function to open the molecule selector dialog
  const handleOpenSelector = () => {
    setIsSelectorOpen(true);
  };

  // LD1: Function to close the molecule selector dialog
  const handleCloseSelector = () => {
    setIsSelectorOpen(false);
  };

  // LD1: Function to handle when molecules are selected from the selector
  const handleMoleculesSelected = (molecules: Molecule[]) => {
    setSelectedMolecules(molecules);
    handleCloseSelector();
  };

  // LD1: Function to handle form submission
  const onSubmit = async (formData: LibraryCreate) => {
    setIsSubmitting(true);
    try {
      // IE1: Create library data object with form values and selected molecule IDs
      const libraryData: LibraryCreate = {
        name: formData.name,
        description: formData.description,
        molecule_ids: selectedMolecules.map(molecule => molecule.id)
      };

      // IE1: Call createLibrary function from useLibraries hook
      const newLibrary = await createLibrary(libraryData);

      // IE1: Show success toast notification
      showToast({
        type: 'success',
        message: 'Library created successfully!'
      });

      // IE1: Call onSuccess callback with created library ID
      onSuccess(newLibrary.id);

      // IE1: Reset form and selected molecules
      reset();
      setSelectedMolecules([]);
    } catch (error: any) {
      // IE1: Show error toast notification
      showToast({
        type: 'error',
        message: error.message || 'Failed to create library.'
      });
    } finally {
      // IE1: Set isSubmitting state to false regardless of outcome
      setIsSubmitting(false);
    }
  };

  return (
    <Card>
      <Typography variant="h6" component="div" sx={{ padding: '16px' }}>
        Create New Library
      </Typography>
      <Divider />
      <Box sx={{ padding: '16px' }}>
        <Controller
          name="name"
          control={control}
          render={({ field, fieldState }) => (
            <Input
              {...field}
              label="Library Name"
              error={fieldState.error?.message}
              helperText={fieldState.error?.message}
              required
              fullWidth
            />
          )}
        />
        <Controller
          name="description"
          control={control}
          render={({ field }) => (
            <Input
              {...field}
              label="Description"
              multiline
              rows={3}
              fullWidth
            />
          )}
        />
        <Typography variant="subtitle1" sx={{ mt: 2 }}>
          Molecules: {selectedMolecules.length} selected
        </Typography>
        <Button variant="outlined" onClick={handleOpenSelector}>
          Select Molecules
        </Button>
      </Box>
      <Divider />
      <Box sx={{ display: 'flex', justifyContent: 'flex-end', padding: '16px' }}>
        <Button onClick={onCancel} sx={{ mr: 1 }}>
          Cancel
        </Button>
        <Button
          variant="contained"
          onClick={handleSubmit(onSubmit)}
          disabled={isSubmitting}
        >
          {isSubmitting ? 'Creating...' : 'Create Library'}
        </Button>
      </Box>

      {/* Molecule Selector Dialog */}
      {isSelectorOpen && (
        <MoleculeSelector
          selectedMoleculeIds={selectedMolecules.map(molecule => molecule.id)}
          onSelectionChange={(selectedIds) => {
            const selected = selectedMolecules.filter(molecule => selectedIds.includes(molecule.id));
            handleMoleculesSelected(selected);
          }}
          onCancel={handleCloseSelector}
        />
      )}
    </Card>
  );
};

export default CreateLibraryForm;