import { useState, useEffect, useCallback } from 'react'; // react 18.2+
import { useAppDispatch, useAppSelector } from '../../store';
import {
  Experiment,
  ExperimentDetail,
  ExperimentCreate,
  ExperimentFilter,
  ExperimentMoleculeOperation,
  ExperimentStatusUpdate,
  ExperimentType
} from '../../types/experiment';
import {
  fetchExperiments,
  fetchMyExperiments,
  fetchExperimentById,
  fetchExperimentTypes,
  createNewExperiment,
  updateExistingExperiment,
  removeExperiment,
  addMoleculesToExperiment,
  removeMoleculesFromExperiment,
  updateStatus,
  queueExistingExperiment,
  submitExistingExperiment,
  cancelExistingExperiment,
  setExperimentFilter,
  setCurrentPage,
  setPageSize,
  clearExperimentError
} from '../../store/experiments/experimentsSlice';
import { useToast } from '../../hooks/useToast';

/**
 * A hook that provides functionality for managing experiments
 * @returns Experiment management functions and state: { experiments, experimentTypes, currentExperiment, loading, error, totalExperiments, currentPage, pageSize, filter, fetchExperiments, fetchMyExperiments, fetchExperimentById, fetchExperimentTypes, createExperiment, updateExperiment, deleteExperiment, addMoleculesToExperiment, removeMoleculesFromExperiment, updateExperimentStatus, queueExperiment, submitExperiment, cancelExperiment, setFilter, setPage, setSize, clearError }
 */
export const useExperiments = () => {
  // Get the Redux dispatch function using useAppDispatch
  const dispatch = useAppDispatch();

  // Select experiment-related state from the Redux store using useAppSelector
  const experiments = useAppSelector((state) => state.experiments.experiments);
  const experimentTypes = useAppSelector((state) => state.experiments.experimentTypes);
  const currentExperiment = useAppSelector((state) => state.experiments.currentExperiment);
  const loading = useAppSelector((state) => state.experiments.loading);
  const error = useAppSelector((state) => state.experiments.error);
  const totalExperiments = useAppSelector((state) => state.experiments.totalExperiments);
  const currentPage = useAppSelector((state) => state.experiments.currentPage);
  const pageSize = useAppSelector((state) => state.experiments.pageSize);
  const filter = useAppSelector((state) => state.experiments.filter);

  // Get the showToast function from useToast hook
  const { showToast } = useToast();

  // Create a memoized setFilter function that dispatches setExperimentFilter action
  const setFilter = useCallback((filter: ExperimentFilter) => {
    dispatch(setExperimentFilter(filter));
  }, [dispatch]);

  // Create a memoized setPage function that dispatches setCurrentPage action
  const setPage = useCallback((page: number) => {
    dispatch(setCurrentPage(page));
  }, [dispatch]);

  // Create a memoized setSize function that dispatches setPageSize action
  const setSize = useCallback((pageSize: number) => {
    dispatch(setPageSize(pageSize));
  }, [dispatch]);

    // Create a memoized clearError function that dispatches clearExperimentError action
    const clearError = useCallback(() => {
      dispatch(clearExperimentError());
    }, [dispatch]);

  // Create a memoized fetchExperiments function that dispatches the fetchExperiments thunk
  const fetchExperimentsData = useCallback(() => {
    dispatch(fetchExperiments(filter));
  }, [dispatch, filter]);

  // Create a memoized fetchMyExperiments function that dispatches the fetchMyExperiments thunk
  const fetchMyExperimentsData = useCallback(() => {
    dispatch(fetchMyExperiments(filter));
  }, [dispatch, filter]);

  // Create a memoized fetchExperimentById function that dispatches the fetchExperimentById thunk
  const fetchExperimentByIdData = useCallback((id: string) => {
    dispatch(fetchExperimentById(id));
  }, [dispatch]);

  // Create a memoized fetchExperimentTypes function that dispatches the fetchExperimentTypes thunk
  const fetchExperimentTypesData = useCallback(() => {
    dispatch(fetchExperimentTypes());
  }, [dispatch]);

  // Create a memoized createExperiment function that dispatches the createNewExperiment thunk and shows success/error toast
  const createExperiment = useCallback((experimentData: ExperimentCreate) => {
    dispatch(createNewExperiment(experimentData))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment created successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to create experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized updateExperiment function that dispatches the updateExistingExperiment thunk and shows success/error toast
  const updateExperiment = useCallback((id: string, data: ExperimentUpdate) => {
    dispatch(updateExistingExperiment({ id, data }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment updated successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to update experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized deleteExperiment function that dispatches the removeExperiment thunk and shows success/error toast
  const deleteExperiment = useCallback((id: string) => {
    dispatch(removeExperiment(id))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment deleted successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to delete experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized addMoleculesToExperiment function that dispatches the addMolecules thunk and shows success/error toast
  const addMoleculesToExperimentData = useCallback((id: string, moleculeIds: string[]) => {
    dispatch(addMolecules({ id, moleculeIds }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Molecules added to experiment successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to add molecules to experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized removeMoleculesFromExperiment function that dispatches the removeMolecules thunk and shows success/error toast
  const removeMoleculesFromExperimentData = useCallback((id: string, moleculeIds: string[]) => {
    dispatch(removeMolecules({ id, moleculeIds }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Molecules removed from experiment successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to remove molecules from experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized updateExperimentStatus function that dispatches the updateStatus thunk and shows success/error toast
  const updateExperimentStatusData = useCallback((id: string, status: string, notes?: string) => {
    dispatch(updateStatus({ id, status, notes }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment status updated successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to update experiment status: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized queueExperiment function that dispatches the queueExistingExperiment thunk and shows success/error toast
  const queueExperimentData = useCallback((id: string, notes?: string) => {
    dispatch(queueExistingExperiment({ id, notes }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment queued successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to queue experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized submitExperiment function that dispatches the submitExistingExperiment thunk and shows success/error toast
  const submitExperimentData = useCallback((id: string, notes?: string) => {
    dispatch(submitExistingExperiment({ id, notes }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment submitted successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to submit experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Create a memoized cancelExperiment function that dispatches the cancelExistingExperiment thunk and shows success/error toast
  const cancelExperimentData = useCallback((id: string, notes?: string) => {
    dispatch(cancelExistingExperiment({ id, notes }))
      .unwrap()
      .then(() => {
        showToast({
          type: 'success',
          message: 'Experiment cancelled successfully!',
        });
      })
      .catch((error) => {
        showToast({
          type: 'error',
          message: `Failed to cancel experiment: ${error}`,
        });
      });
  }, [dispatch, showToast]);

  // Return an object containing all the state and functions
  return {
    experiments,
    experimentTypes,
    currentExperiment,
    loading,
    error,
    totalExperiments,
    currentPage,
    pageSize,
    filter,
    fetchExperiments: fetchExperimentsData,
    fetchMyExperiments: fetchMyExperimentsData,
    fetchExperimentById: fetchExperimentByIdData,
    fetchExperimentTypes: fetchExperimentTypesData,
    createExperiment: createExperiment,
    updateExperiment: updateExperiment,
    deleteExperiment: deleteExperiment,
    addMoleculesToExperiment: addMoleculesToExperimentData,
    removeMoleculesFromExperiment: removeMoleculesFromExperimentData,
    updateExperimentStatus: updateExperimentStatusData,
    queueExperiment: queueExperimentData,
    submitExperiment: submitExperimentData,
    cancelExperiment: cancelExperimentData,
    setFilter: setFilter,
    setPage: setPage,
    setSize: setSize,
    clearError: clearError
  };
};