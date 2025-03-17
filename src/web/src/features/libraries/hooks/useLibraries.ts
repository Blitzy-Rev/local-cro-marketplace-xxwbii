# src/web/src/features/libraries/hooks/useLibraries.ts
```typescript
import { useState, useEffect, useCallback, useMemo } from 'react'; // react 18.2+
import { useQuery, useMutation, useQueryClient } from 'react-query'; // react-query 4.28+
import { useAppDispatch, useAppSelector } from '../../../store';
import { 
  Library, 
  LibraryDetail, 
  LibraryCreate, 
  LibraryUpdate, 
  LibraryFilter,
  LibraryOperationType
} from '../../../types/library';
import { Molecule } from '../../../types/molecule';
import { 
  selectLibraries, 
  selectCurrentLibrary, 
  selectLibraryLoading, 
  selectLibraryError, 
  selectTotalLibraries, 
  selectCurrentPage, 
  selectPageSize 
} from '../../../store/libraries/librariesSlice';
import { 
  fetchLibraries, 
  fetchMyLibraries,
  fetchLibraryById,
  fetchLibraryDetailById,
  createNewLibrary, 
  updateExistingLibrary, 
  removeLibrary, 
  addMoleculesToExistingLibrary, 
  removeMoleculesFromExistingLibrary,
  updateLibraryMoleculesOperation,
  setCurrentLibrary, 
  clearCurrentLibrary,
  setCurrentPage,
  setPageSize
} from '../../../store/libraries/librariesSlice';
import useDragDrop from './useDragDrop';
import usePagination from '../../../hooks/usePagination';
import { useToast } from '../../../hooks/useToast';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * A hook that provides comprehensive functionality for managing molecule libraries
 * @param options - Configuration options for the hook
 * @returns An object containing all library management state and functions
 */
const useLibraries = (options: {
  filter?: LibraryFilter;
  initialPage?: number;
  initialPageSize?: number;
  fetchMyLibrariesOnly?: boolean;
} = {}) => {
  // LD1: Initialize Redux dispatch and selectors for library state
  const dispatch = useAppDispatch();
  const libraries = useAppSelector(selectLibraries);
  const currentLibrary = useAppSelector(selectCurrentLibrary);
  const loading = useAppSelector(selectLibraryLoading);
  const error = useAppSelector(selectLibraryError);
  const totalLibraries = useAppSelector(selectTotalLibraries);

  // LD1: Initialize React Query client for cache invalidation
  const queryClient = useQueryClient();

  // LD1: Initialize toast notifications with useToast
  const { showToast } = useToast();

  // LD1: Check user permissions with usePermissions
  const { canManageLibraries } = usePermissions();

  // LD1: Set up pagination with usePagination hook
  const pagination = usePagination({
    initialPage: options.initialPage,
    initialPageSize: options.initialPageSize,
    totalItems: totalLibraries,
  });

  // LD1: Set up drag and drop functionality with useDragDrop hook
  const dragDrop = useDragDrop({
    onDrop: async (molecule: Molecule, libraryId: string, sourceLibraryId: string | null) => {
      if (!canManageLibraries()) {
        showToast({
          type: 'error',
          message: 'You do not have permission to manage libraries.'
        });
        return;
      }

      try {
        await dispatch(updateLibraryMoleculesOperation({
          libraryId: libraryId,
          operation: LibraryOperationType.ADD,
          moleculeIds: [molecule.id]
        })).unwrap();

        // IE1: Invalidate queries related to the target library to reflect changes
        queryClient.invalidateQueries(['libraries', libraryId]);
        queryClient.invalidateQueries(['molecules']);

        showToast({
          type: 'success',
          message: `Molecule added to library successfully.`
        });
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error.message || 'Failed to add molecule to library.'
        });
      }
    }
  });

  // LD1: Fetch libraries based on filter and pagination using Redux thunks
  useEffect(() => {
    const fetchAction = options.fetchMyLibrariesOnly ? fetchMyLibraries : fetchLibraries;
    dispatch(fetchAction({ 
      page: pagination.page, 
      pageSize: pagination.pageSize,
      filter: options.filter
    }));
  }, [dispatch, pagination.page, pagination.pageSize, options.filter, options.fetchMyLibrariesOnly]);

  // LD1: Create mutation functions for library operations (create, update, delete)
  const createLibraryMutation = useMutation(
    (libraryData: LibraryCreate) => dispatch(createNewLibrary(libraryData)).unwrap(),
    {
      onSuccess: () => {
        // IE1: Invalidate queries related to libraries to reflect changes
        queryClient.invalidateQueries('libraries');
        showToast({ type: 'success', message: 'Library created successfully.' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: error.message || 'Failed to create library.' });
      },
    }
  );

  const updateLibraryMutation = useMutation(
    (params: { id: string; data: LibraryUpdate }) => dispatch(updateExistingLibrary(params)).unwrap(),
    {
      onSuccess: () => {
        // IE1: Invalidate queries related to libraries to reflect changes
        queryClient.invalidateQueries('libraries');
        showToast({ type: 'success', message: 'Library updated successfully.' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: error.message || 'Failed to update library.' });
      },
    }
  );

  const deleteLibraryMutation = useMutation(
    (id: string) => dispatch(removeLibrary(id)).unwrap(),
    {
      onSuccess: () => {
        // IE1: Invalidate queries related to libraries to reflect changes
        queryClient.invalidateQueries('libraries');
        showToast({ type: 'success', message: 'Library deleted successfully.' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: error.message || 'Failed to delete library.' });
      },
    }
  );

  // LD1: Create mutation functions for molecule operations (add, remove)
  const addMoleculesMutation = useMutation(
    (params: { libraryId: string; moleculeIds: string[] }) => dispatch(addMoleculesToExistingLibrary(params)).unwrap(),
    {
      onSuccess: () => {
        // IE1: Invalidate queries related to libraries to reflect changes
        queryClient.invalidateQueries('libraries');
        showToast({ type: 'success', message: 'Molecules added to library successfully.' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: error.message || 'Failed to add molecules to library.' });
      },
    }
  );

  const removeMoleculesMutation = useMutation(
    (params: { libraryId: string; moleculeIds: string[] }) => dispatch(removeMoleculesFromExistingLibrary(params)).unwrap(),
    {
      onSuccess: () => {
        // IE1: Invalidate queries related to libraries to reflect changes
        queryClient.invalidateQueries('libraries');
        showToast({ type: 'success', message: 'Molecules removed from library successfully.' });
      },
      onError: (error: any) => {
        showToast({ type: 'error', message: error.message || 'Failed to remove molecules from library.' });
      },
    }
  );

  // LD1: Set up handlers for drag and drop operations
  const fetchLibrary = useCallback((id: string) => {
    dispatch(fetchLibraryById(id));
  }, [dispatch]);

  // LD1: Implement library refresh functionality
  const refreshLibraries = useCallback(() => {
    queryClient.invalidateQueries('libraries');
    dispatch(fetchLibraries({ 
      page: pagination.page, 
      pageSize: pagination.pageSize,
      filter: options.filter
    }));
  }, [queryClient, dispatch, pagination.page, pagination.pageSize, options.filter]);

  // LD1: Return an object containing all library management state and functions
  return {
    libraries,
    currentLibrary,
    loading,
    error,
    totalLibraries,
    pagination,
    dragDrop,
    fetchLibrary,
    createLibrary: createLibraryMutation.mutate,
    updateLibrary: updateLibraryMutation.mutate,
    deleteLibrary: deleteLibraryMutation.mutate,
    addMolecules: addMoleculesMutation.mutate,
    removeMolecules: removeMoleculesMutation.mutate,
    setCurrentLibrary: (library: LibraryDetail) => dispatch(setCurrentLibrary(library)),
    clearCurrentLibrary: () => dispatch(clearCurrentLibrary()),
    refreshLibraries,
  };
};

export default useLibraries;