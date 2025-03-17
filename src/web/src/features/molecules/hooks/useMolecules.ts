import { useState, useEffect, useCallback, useMemo } from 'react'; // react v18.2+
import { useAppDispatch, useAppSelector } from '../../../store';
import {
  fetchMolecules,
  updateMoleculeFlagStatus,
  bulkUpdateMoleculeFlags,
  addMoleculesToLibrary,
  addMoleculesToExperiment,
  selectMolecules,
  selectMoleculeFilter,
  selectSelectedMolecules,
  selectMoleculeLoading,
  selectMoleculeError,
  selectTotalMolecules,
  selectCurrentPage,
  selectPageSize,
  setFilter,
  setSelectedMolecules,
  addSelectedMolecule,
  removeSelectedMolecule,
  clearSelectedMolecules,
  setCurrentPage,
  setPageSize
} from '../../../store/molecules/moleculesSlice';
import useDebounce from '../../../hooks/useDebounce';
import usePagination from '../../../hooks/usePagination';
import useToast from '../../../hooks/useToast';
import { Molecule, MoleculeFilter, PropertyFilter } from '../../../types/molecule';
import { getMoleculePropertyByName, sortMoleculesByProperty, filterMoleculesByProperty } from '../../../utils/molecularUtils';

/**
 * A hook that provides comprehensive functionality for managing molecule data.
 * @param options - Configuration options for the hook.
 * @returns An object containing molecule management state and functions.
 */
const useMolecules = (options: { initialFilter?: MoleculeFilter; initialPageSize?: number; enableAutoFetch?: boolean } = {}) => {
  // LD1: Initialize Redux dispatch and selectors for molecule state
  const dispatch = useAppDispatch();
  const molecules = useAppSelector(selectMolecules);
  const loading = useAppSelector(selectMoleculeLoading);
  const error = useAppSelector(selectMoleculeError);
  const selectedMolecules = useAppSelector(selectSelectedMolecules);
  const totalMolecules = useAppSelector(selectTotalMolecules);
  const currentPage = useAppSelector(selectCurrentPage);
  const pageSize = useAppSelector(selectPageSize);
  const initialFilter = options.initialFilter || {};
  const initialPageSize = options.initialPageSize || 10;
  const enableAutoFetch = options.enableAutoFetch !== false;

  // LD1: Set up local state for filter with default values
  const [filter, setLocalFilter] = useState<MoleculeFilter>(initialFilter);

  // LD1: Create debounced filter to prevent excessive API calls
  const debouncedFilter = useDebounce(filter, 300);

  // LD1: Initialize pagination hook with total molecules count
  const pagination = usePagination({
    initialPage: currentPage,
    initialPageSize: initialPageSize,
    totalItems: totalMolecules,
  });

  // LD1: Initialize toast notification hook for user feedback
  const { showToast } = useToast();

  // LD1: Set up effect to fetch molecules when filter or pagination changes
  useEffect(() => {
    if (enableAutoFetch) {
      dispatch(
        fetchMolecules({
          filter: debouncedFilter,
          page: pagination.page,
          pageSize: pagination.pageSize,
        })
      );
    }
  }, [debouncedFilter, pagination.page, pagination.pageSize, dispatch, enableAutoFetch]);

  // LD1: Implement refreshMolecules function to manually trigger data fetch
  const refreshMolecules = useCallback(() => {
    dispatch(
      fetchMolecules({
        filter: debouncedFilter,
        page: pagination.page,
        pageSize: pagination.pageSize,
      })
    );
  }, [debouncedFilter, pagination.page, pagination.pageSize, dispatch]);

  // LD1: Implement molecule selection functions (select, deselect, selectAll, deselectAll)
  const selectMolecule = useCallback(
    (id: string) => {
      dispatch(addSelectedMolecule(id));
    },
    [dispatch]
  );

  const deselectMolecule = useCallback(
    (id: string) => {
      dispatch(removeSelectedMolecule(id));
    },
    [dispatch]
  );

  const selectAllMolecules = useCallback(
    (ids: string[]) => {
      dispatch(setSelectedMolecules(ids));
    },
    [dispatch]
  );

  const deselectAllMolecules = useCallback(() => {
    dispatch(clearSelectedMolecules());
  }, [dispatch]);

  // LD1: Implement updateMoleculeFlag function to update flag status
  const updateMoleculeFlag = useCallback(
    async (id: string, flagStatus: string) => {
      try {
        await dispatch(updateMoleculeFlagStatus({ id, flagStatus }));
        showToast({
          type: 'success',
          message: 'Molecule flag updated successfully',
        });
        refreshMolecules();
      } catch (e: any) {
        showToast({
          type: 'error',
          message: e.message || 'Failed to update molecule flag',
        });
      }
    },
    [dispatch, showToast, refreshMolecules]
  );

  // LD1: Implement addToLibrary function to add molecules to a library
  const addToLibrary = useCallback(
    async (libraryId: string) => {
      try {
        await dispatch(
          addMoleculesToLibrary({
            moleculeIds: selectedMolecules,
            libraryId: libraryId,
          })
        );
        showToast({
          type: 'success',
          message: 'Molecules added to library successfully',
        });
        dispatch(clearSelectedMolecules());
      } catch (e: any) {
        showToast({
          type: 'error',
          message: e.message || 'Failed to add molecules to library',
        });
      }
    },
    [dispatch, selectedMolecules, showToast]
  );

  // LD1: Implement addToExperiment function to add molecules to an experiment
  const addToExperiment = useCallback(
    async (experimentId: string) => {
      try {
        await dispatch(
          addMoleculesToExperiment({
            moleculeIds: selectedMolecules,
            experimentId: experimentId,
          })
        );
        showToast({
          type: 'success',
          message: 'Molecules added to experiment successfully',
        });
        dispatch(clearSelectedMolecules());
      } catch (e: any) {
        showToast({
          type: 'error',
          message: e.message || 'Failed to add molecules to experiment',
        });
      }
    },
    [dispatch, selectedMolecules, showToast]
  );

  // LD1: Return comprehensive object with all molecule management state and functions
  return {
    molecules,
    loading,
    error,
    filter: debouncedFilter,
    selectedMolecules,
    totalMolecules,
    pagination,
    setFilter: (newFilter: MoleculeFilter) => {
      dispatch(setCurrentPage(1));
      setLocalFilter(newFilter);
      dispatch(setFilter(newFilter));
    },
    refreshMolecules,
    selectMolecule,
    deselectMolecule,
    selectAllMolecules,
    deselectAllMolecules,
    updateMoleculeFlag,
    addToLibrary,
    addToExperiment,
  };
};

export default useMolecules;