import { useState, useEffect, useCallback, useMemo } from 'react'; // react ^18.2.0
import { useDispatch, useSelector } from 'react-redux'; // react-redux ^8.0.5
import { saveAs } from 'file-saver'; // file-saver ^2.0.5

import useDebounce from '../../../hooks/useDebounce';
import { 
  Result, 
  ResultDetailed, 
  ResultFilter,
  ResultApproval,
  ResultListResponse,
  ResultDataCreate
} from '../../../types/result';
import {
  fetchResults,
  fetchResultsBySubmission,
  fetchResultDetail,
  approveResult,
  rejectResult,
  batchApproveResults,
  batchRejectResults,
  exportResultData,
  analyzeResultData
} from '../../../store/results/resultsSlice';
import {
  selectResults,
  selectCurrentResult,
  selectResultFilter,
  selectResultLoading,
  selectResultError,
  selectTotalResults,
  selectCurrentPage,
  selectPageSize,
  setFilter,
  setCurrentPage,
  setPageSize,
  clearResultError
} from '../../../store/results/resultsSlice';
import {
  getResultFile,
  deleteResultFile,
  addResultData,
  getResultData
} from '../../../api/results';

/**
 * Interface for the options that can be passed to useResults hook
 */
interface UseResultsOptions {
  /** Initial filter criteria for results */
  initialFilter?: ResultFilter;
  /** Debounce delay in milliseconds for filter changes */
  debounceMs?: number;
  /** Number of results per page */
  pageSize?: number;
}

/**
 * Interface for pagination information
 */
interface PaginationInfo {
  /** Total number of results matching the current filter */
  totalResults: number;
  /** Current page number (1-based) */
  currentPage: number;
  /** Number of results per page */
  pageSize: number;
  /** Total number of pages */
  totalPages: number;
}

/**
 * Custom hook that provides comprehensive functionality for managing experimental results
 * in the Molecular Data Management and CRO Integration Platform.
 * 
 * This hook handles:
 * - Fetching and filtering results
 * - Pagination
 * - Result approval/rejection workflows
 * - Batch operations
 * - File downloads and management
 * - Exporting results in different formats
 * - Result data analysis
 * - Error handling
 * 
 * @param options Configuration options for the hook
 * @returns Object containing result data, loading state, error state, pagination info, and result management functions
 */
const useResults = (options?: UseResultsOptions) => {
  const {
    initialFilter = {},
    debounceMs = 300,
    pageSize = 10
  } = options || {};

  // Initialize Redux dispatch and selectors
  const dispatch = useDispatch();
  const results = useSelector(selectResults);
  const currentResult = useSelector(selectCurrentResult);
  const storeFilter = useSelector(selectResultFilter);
  const loading = useSelector(selectResultLoading);
  const error = useSelector(selectResultError);
  const totalResults = useSelector(selectTotalResults);
  const currentPage = useSelector(selectCurrentPage);
  const storePageSize = useSelector(selectPageSize);

  // Local filter state with initial values merged from store and options
  const [filter, setLocalFilter] = useState<ResultFilter>({
    ...storeFilter,
    ...initialFilter,
    page: initialFilter.page || storeFilter.page || 1,
    page_size: initialFilter.page_size || pageSize
  });

  // Debounce filter changes to prevent excessive API calls
  const debouncedFilter = useDebounce(filter, debounceMs);

  // Calculate total pages based on total results and page size
  const totalPages = useMemo(() => {
    return Math.ceil(totalResults / (filter.page_size || pageSize));
  }, [totalResults, filter.page_size, pageSize]);

  // Update global state when debounced filter changes
  useEffect(() => {
    dispatch(setFilter(debouncedFilter));
    
    // Update page size in store if changed
    if (debouncedFilter.page_size && debouncedFilter.page_size !== storePageSize) {
      dispatch(setPageSize(debouncedFilter.page_size));
    }
    
    // Update current page in store if changed
    if (debouncedFilter.page && debouncedFilter.page !== currentPage) {
      dispatch(setCurrentPage(debouncedFilter.page));
    }
  }, [dispatch, debouncedFilter, storePageSize, currentPage]);

  // Fetch results when filter changes
  useEffect(() => {
    if (debouncedFilter.submission_id) {
      // If filtering by submission, use specialized endpoint
      dispatch(fetchResultsBySubmission({
        submissionId: debouncedFilter.submission_id,
        page: debouncedFilter.page,
        pageSize: debouncedFilter.page_size
      }));
    } else {
      // Otherwise fetch all results with general filter
      dispatch(fetchResults({
        filter: debouncedFilter,
        page: debouncedFilter.page,
        pageSize: debouncedFilter.page_size
      }));
    }
  }, [dispatch, debouncedFilter]);

  /**
   * Updates the filter with new criteria
   * @param partialFilter New filter criteria to apply (partial)
   */
  const updateFilter = useCallback((partialFilter: Partial<ResultFilter>) => {
    setLocalFilter(prevFilter => ({
      ...prevFilter,
      ...partialFilter,
      // Reset to page 1 when filter changes (except when explicitly changing page)
      page: 'page' in partialFilter ? partialFilter.page : 1
    }));
  }, []);

  /**
   * Changes the current page
   * @param page New page number (1-based)
   */
  const setPage = useCallback((page: number) => {
    updateFilter({ page });
  }, [updateFilter]);

  /**
   * Changes the page size
   * @param newPageSize New number of items per page
   */
  const updatePageSize = useCallback((newPageSize: number) => {
    updateFilter({ page: 1, page_size: newPageSize });
  }, [updateFilter]);

  /**
   * Fetches detailed information for a specific result
   * @param resultId ID of the result to fetch
   * @returns Promise resolving to the result details
   */
  const fetchResultDetails = useCallback((resultId: string) => {
    return dispatch(fetchResultDetail(resultId));
  }, [dispatch]);

  /**
   * Approves a result
   * @param resultId ID of the result to approve
   * @param notes Optional notes explaining the approval
   * @returns Promise resolving when approval is complete
   */
  const approveResultFn = useCallback((resultId: string, notes?: string) => {
    return dispatch(approveResult({ resultId, notes }));
  }, [dispatch]);

  /**
   * Rejects a result
   * @param resultId ID of the result to reject
   * @param notes Optional notes explaining the rejection
   * @returns Promise resolving when rejection is complete
   */
  const rejectResultFn = useCallback((resultId: string, notes?: string) => {
    return dispatch(rejectResult({ resultId, notes }));
  }, [dispatch]);

  /**
   * Approves multiple results in a batch operation
   * @param resultIds Array of result IDs to approve
   * @param notes Optional notes explaining the approval
   * @returns Promise resolving when approval is complete
   */
  const batchApproveResultsFn = useCallback((resultIds: string[], notes?: string) => {
    return dispatch(batchApproveResults({ resultIds, notes }));
  }, [dispatch]);

  /**
   * Rejects multiple results in a batch operation
   * @param resultIds Array of result IDs to reject
   * @param notes Optional notes explaining the rejection
   * @returns Promise resolving when rejection is complete
   */
  const batchRejectResultsFn = useCallback((resultIds: string[], notes?: string) => {
    return dispatch(batchRejectResults({ resultIds, notes }));
  }, [dispatch]);

  /**
   * Downloads a file associated with a result
   * @param fileId ID of the file to download
   * @param fileName Name to save the file as
   * @returns Promise resolving to true if download was successful, false otherwise
   */
  const downloadResultFile = useCallback(async (fileId: string, fileName: string): Promise<boolean> => {
    try {
      const blob = await getResultFile(fileId);
      saveAs(blob, fileName);
      return true;
    } catch (error) {
      console.error('Error downloading file:', error);
      return false;
    }
  }, []);

  /**
   * Deletes a file associated with a result
   * @param fileId ID of the file to delete
   * @returns Promise resolving to true if deletion was successful, false otherwise
   */
  const deleteResultFileFn = useCallback(async (fileId: string): Promise<boolean> => {
    try {
      await deleteResultFile(fileId);
      return true;
    } catch (error) {
      console.error('Error deleting file:', error);
      return false;
    }
  }, []);

  /**
   * Exports result data in various formats
   * @param resultId ID of the result to export
   * @param format Format to export (csv, xlsx, or json)
   * @returns Promise resolving to the exported blob
   */
  const exportResult = useCallback((resultId: string, format: 'csv' | 'xlsx' | 'json') => {
    return dispatch(exportResultData({ resultId, format }))
      .unwrap()
      .then((blob) => {
        const fileName = `result_${resultId}.${format}`;
        saveAs(blob, fileName);
        return blob;
      });
  }, [dispatch]);

  /**
   * Generates analysis of result data
   * @param resultId ID of the result to analyze
   * @returns Promise resolving to the analysis data
   */
  const analyzeResult = useCallback((resultId: string) => {
    return dispatch(analyzeResultData(resultId));
  }, [dispatch]);

  /**
   * Adds structured data to a result
   * @param resultId ID of the result to add data to
   * @param dataPoint Data point to add
   * @returns Promise resolving to true if data was added successfully, false otherwise
   */
  const addResultDataFn = useCallback(async (resultId: string, dataPoint: ResultDataCreate): Promise<boolean> => {
    try {
      await addResultData(resultId, dataPoint);
      // Refresh result details after adding data
      dispatch(fetchResultDetail(resultId));
      return true;
    } catch (error) {
      console.error('Error adding result data:', error);
      return false;
    }
  }, [dispatch]);

  /**
   * Gets structured data for a result
   * @param resultId ID of the result to get data for
   * @param options Optional filter options
   * @returns Promise resolving to the result data
   */
  const getResultDataFn = useCallback(async (resultId: string, options?: { molecule_id?: string }) => {
    try {
      const response = await getResultData(resultId, options);
      return response.data;
    } catch (error) {
      console.error('Error getting result data:', error);
      return null;
    }
  }, []);

  /**
   * Clears any result error message
   */
  const clearError = useCallback(() => {
    dispatch(clearResultError());
  }, [dispatch]);

  /**
   * Requests additional data from the CRO for a result
   * @param resultId ID of the result to request additional data for
   * @param notes Notes explaining what additional data is needed
   * @returns Promise resolving to true if request was sent successfully
   */
  const requestAdditionalData = useCallback(async (resultId: string, notes: string): Promise<boolean> => {
    try {
      // This would be implemented with a real API call in production
      // For now, we'll log the request and return success
      console.log(`Requesting additional data for result ${resultId}: ${notes}`);
      
      // In a real implementation, this would dispatch an action to create a notification
      // or message for the CRO, for example:
      // await dispatch(createResultRequest({ resultId, notes }));
      
      return true;
    } catch (error) {
      console.error('Error requesting additional data:', error);
      return false;
    }
  }, []);

  // Return all necessary data and functions
  return {
    // Data
    results,
    currentResult,
    loading,
    error,
    
    // Pagination
    pagination: {
      totalResults,
      currentPage,
      pageSize: filter.page_size || pageSize,
      totalPages
    },
    
    // Filter
    filter,
    setFilter: updateFilter,
    
    // Pagination functions
    setPage,
    setPageSize: updatePageSize,
    
    // Result operations
    fetchResultDetail: fetchResultDetails,
    approveResult: approveResultFn,
    rejectResult: rejectResultFn,
    batchApproveResults: batchApproveResultsFn,
    batchRejectResults: batchRejectResultsFn,
    downloadResultFile,
    deleteResultFile: deleteResultFileFn,
    exportResult,
    analyzeResult,
    addResultData: addResultDataFn,
    getResultData: getResultDataFn,
    clearError,
    requestAdditionalData
  };
};

export default useResults;