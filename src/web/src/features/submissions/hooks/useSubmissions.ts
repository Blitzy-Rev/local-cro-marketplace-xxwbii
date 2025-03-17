import { useState, useEffect, useCallback, useMemo } from 'react'; // react v18.2.0
import { useAppDispatch, useAppSelector } from '../../../store';
import {
  fetchSubmissions,
  fetchSubmissionById,
  fetchSubmissionsByExperiment,
  fetchSubmissionsByCRO,
  fetchSubmissionsByStatus,
  createNewSubmission,
  updateExistingSubmission,
  updateSubmissionStatusThunk,
  provideQuoteForSubmission,
  respondToSubmissionQuote,
  cancelExistingSubmission,
  setSubmissionFilter,
  setCurrentPage,
  setPageSize,
  clearSubmissionError
} from '../../../store/submissions/submissionsSlice';
import {
  selectSubmissions,
  selectCurrentSubmission,
  selectSubmissionLoading,
  selectSubmissionError,
  selectTotalSubmissions,
  selectCurrentPage,
  selectPageSize,
  selectSubmissionFilter
} from '../../../store/submissions/submissionsSlice';
import useDebounce from '../../../hooks/useDebounce';
import usePagination from '../../../hooks/usePagination';
import useToast from '../../../hooks/useToast';
import {
  Submission,
  SubmissionDetailed,
  SubmissionFilter,
  SubmissionCreate,
  SubmissionUpdate,
  SubmissionStatus,
  QuoteProvide,
  QuoteResponse
} from '../../../types/submission';

/**
 * A hook that provides comprehensive functionality for managing submission data.
 *
 * @param options - Configuration options for the hook.
 * @returns An object containing submission management state and functions.
 */
const useSubmissions = (options: {
  initialFilter?: SubmissionFilter;
  initialPageSize?: number;
  enableAutoFetch?: boolean;
} = {}) => {
  // LD1: Destructure options with default values
  const { initialFilter = {}, initialPageSize = 10, enableAutoFetch = true } = options;

  // LD1: Initialize Redux dispatch and selectors for submission state
  const dispatch = useAppDispatch();
  const submissions = useAppSelector(selectSubmissions);
  const currentSubmission = useAppSelector(selectCurrentSubmission);
  const loading = useAppSelector(selectSubmissionLoading);
  const error = useAppSelector(selectSubmissionError);
  const totalSubmissions = useAppSelector(selectTotalSubmissions);
  const currentPage = useAppSelector(selectCurrentPage);
  const pageSize = useAppSelector(selectPageSize);
  const currentFilter = useAppSelector(selectSubmissionFilter);

  // LD1: Set up local state for filter with default values from options
  const [filter, setLocalFilter] = useState<SubmissionFilter>(initialFilter);

  // LD1: Create debounced filter to prevent excessive API calls
  const debouncedFilter = useDebounce(filter, 300);

  // LD1: Initialize pagination hook with total submissions count
  const pagination = usePagination({
    initialPageSize: initialPageSize,
    totalItems: totalSubmissions,
    initialPage: currentPage,
  });

  // LD1: Initialize toast notification hook for user feedback
  const { showToast } = useToast();

  // LD1: Set up effect to fetch submissions when filter or pagination changes
  useEffect(() => {
    // IE1: Check if auto-fetching is enabled or not
    if (enableAutoFetch) {
      dispatch(
        fetchSubmissions({
          ...debouncedFilter,
          page: pagination.page,
          pageSize: pagination.pageSize,
        })
      );
    }
  }, [debouncedFilter, pagination.page, pagination.pageSize, dispatch, enableAutoFetch]);

  // LD1: Implement refreshSubmissions function to manually trigger data fetch
  const refreshSubmissions = useCallback(() => {
    dispatch(
      fetchSubmissions({
        ...debouncedFilter,
        page: pagination.page,
        pageSize: pagination.pageSize,
      })
    );
  }, [debouncedFilter, pagination.page, pagination.pageSize, dispatch]);

  // LD1: Implement fetchSubmissionDetails function to fetch detailed submission data
  const fetchSubmissionDetails = useCallback(
    (id: string) => {
      dispatch(fetchSubmissionById(id));
    },
    [dispatch]
  );

  // LD1: Implement createSubmission function to create a new submission
  const createSubmission = useCallback(
    async (submissionData: SubmissionCreate) => {
      try {
        await dispatch(createNewSubmission(submissionData)).unwrap();
        showToast({
          type: 'success',
          message: 'Submission created successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to create submission',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement updateSubmission function to update an existing submission
  const updateSubmission = useCallback(
    async (id: string, submissionData: SubmissionUpdate) => {
      try {
        await dispatch(updateExistingSubmission({ id, data: submissionData })).unwrap();
        showToast({
          type: 'success',
          message: 'Submission updated successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to update submission',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement updateSubmissionStatus function to update submission status
  const updateSubmissionStatus = useCallback(
    async (id: string, status: SubmissionStatus, notes?: string) => {
      try {
        await dispatch(updateSubmissionStatusThunk({ id, status, notes })).unwrap();
        showToast({
          type: 'success',
          message: 'Submission status updated successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to update submission status',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement provideQuote function for CRO to provide a quote
  const provideQuote = useCallback(
    async (id: string, quoteData: QuoteProvide) => {
      try {
        await dispatch(provideQuoteForSubmission({ id, quoteData })).unwrap();
        showToast({
          type: 'success',
          message: 'Quote provided successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to provide quote',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement respondToQuote function for pharma user to respond to a quote
  const respondToQuote = useCallback(
    async (id: string, responseData: QuoteResponse) => {
      try {
        await dispatch(respondToSubmissionQuote({ id, responseData })).unwrap();
        showToast({
          type: 'success',
          message: 'Quote response submitted successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to respond to quote',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement cancelSubmission function to cancel a submission
  const cancelSubmission = useCallback(
    async (id: string, notes?: string) => {
      try {
        await dispatch(cancelExistingSubmission({ id, notes })).unwrap();
        showToast({
          type: 'success',
          message: 'Submission cancelled successfully!',
        });
        refreshSubmissions();
      } catch (error: any) {
        showToast({
          type: 'error',
          message: error?.message || 'Failed to cancel submission',
        });
      }
    },
    [dispatch, showToast, refreshSubmissions]
  );

  // LD1: Implement fetchByExperiment function to fetch submissions by experiment
  const fetchByExperiment = useCallback(
    (experimentId: string) => {
      dispatch(fetchSubmissionsByExperiment(experimentId));
    },
    [dispatch]
  );

  // LD1: Implement fetchByCRO function to fetch submissions assigned to CRO
  const fetchByCRO = useCallback(
    () => {
      dispatch(fetchSubmissionsByCRO());
    },
    [dispatch]
  );

  // LD1: Implement fetchByStatus function to fetch submissions by status
  const fetchByStatus = useCallback(
    (status: string) => {
      dispatch(fetchSubmissionsByStatus({ status }));
    },
    [dispatch]
  );

  // LD1: Implement setFilter function to update the filter in Redux store
  const setFilter = useCallback(
    (newFilter: SubmissionFilter) => {
      setLocalFilter(newFilter);
      dispatch(setSubmissionFilter(newFilter));
    },
    [dispatch]
  );

  // IE1: Implement setPage function to update the current page in Redux store
  const setPage = useCallback(
    (page: number) => {
      dispatch(setCurrentPage(page));
    },
    [dispatch]
  );

  // IE1: Implement setPageSize function to update the page size in Redux store
  const setPageSize = useCallback(
    (pageSize: number) => {
      dispatch(setPageSize(pageSize));
    },
    [dispatch]
  );

  // LD1: Clear Submission Error
  const clearError = useCallback(() => {
    dispatch(clearSubmissionError());
  }, [dispatch]);

  // LD1: Return comprehensive object with all submission management state and functions
  return {
    submissions,
    currentSubmission,
    loading,
    error,
    filter: currentFilter,
    totalSubmissions,
    pagination: {
      ...pagination,
      goToPage: setPage,
      setPageSize: setPageSize,
    },
    setFilter,
    refreshSubmissions,
    fetchSubmissionDetails,
    createSubmission,
    updateSubmission,
    updateSubmissionStatus,
    provideQuote,
    respondToQuote,
    cancelSubmission,
    fetchByExperiment,
    fetchByCRO,
    fetchByStatus,
    clearError
  };
};

export default useSubmissions;