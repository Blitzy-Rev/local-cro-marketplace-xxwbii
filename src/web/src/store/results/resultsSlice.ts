/**
 * Redux Toolkit slice for managing experimental results data in the Molecular Data Management and CRO Integration Platform.
 * This slice handles state related to results, including fetching, filtering, and managing result data from CRO submissions.
 * 
 * @module store/results/resultsSlice
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import { 
  ResultState, 
  Result, 
  ResultDetailed, 
  ResultFilter, 
  ResultApproval 
} from '../../types/result';
import {
  getResults,
  getResult,
  getResultsBySubmission,
  approveRejectResult,
  batchApproveRejectResults,
  exportResult,
  analyzeResult
} from '../../api/results';

/**
 * Initial state for the results slice
 */
const initialState: ResultState = {
  results: [],
  currentResult: null,
  loading: false,
  error: null,
  totalResults: 0,
  currentPage: 1,
  pageSize: 10,
  filter: {
    submission_id: null,
    status: null,
    uploaded_after: null,
    uploaded_before: null,
    sort_by: 'uploaded_at',
    sort_desc: true
  }
};

/**
 * Async thunk for fetching paginated results with optional filtering
 */
export const fetchResults = createAsyncThunk(
  'results/fetchResults',
  async (params: { filter?: ResultFilter, page?: number, pageSize?: number }, { rejectWithValue }) => {
    try {
      const response = await getResults({
        ...params.filter,
        page: params.page,
        page_size: params.pageSize
      });
      
      return {
        results: response.data.items,
        total: response.data.total
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch results');
    }
  }
);

/**
 * Async thunk for fetching results for a specific submission
 */
export const fetchResultsBySubmission = createAsyncThunk(
  'results/fetchResultsBySubmission',
  async (params: { submissionId: string, page?: number, pageSize?: number }, { rejectWithValue }) => {
    try {
      const response = await getResultsBySubmission(
        params.submissionId,
        { page: params.page, pageSize: params.pageSize }
      );
      
      return {
        results: response.data.items,
        total: response.data.total
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch results for submission');
    }
  }
);

/**
 * Async thunk for fetching detailed information about a single result
 */
export const fetchResultDetail = createAsyncThunk(
  'results/fetchResultDetail',
  async (resultId: string, { rejectWithValue }) => {
    try {
      const response = await getResult(resultId);
      return response.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch result details');
    }
  }
);

/**
 * Async thunk for approving a result
 */
export const approveResult = createAsyncThunk(
  'results/approveResult',
  async (params: { resultId: string, notes?: string }, { rejectWithValue }) => {
    try {
      const approvalData: ResultApproval = {
        approved: true,
        notes: params.notes
      };
      
      const response = await approveRejectResult(params.resultId, approvalData);
      return response.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to approve result');
    }
  }
);

/**
 * Async thunk for rejecting a result
 */
export const rejectResult = createAsyncThunk(
  'results/rejectResult',
  async (params: { resultId: string, notes?: string }, { rejectWithValue }) => {
    try {
      const approvalData: ResultApproval = {
        approved: false,
        notes: params.notes
      };
      
      const response = await approveRejectResult(params.resultId, approvalData);
      return response.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to reject result');
    }
  }
);

/**
 * Async thunk for batch approving multiple results
 */
export const batchApproveResults = createAsyncThunk(
  'results/batchApproveResults',
  async (params: { resultIds: string[], notes?: string }, { rejectWithValue }) => {
    try {
      const response = await batchApproveRejectResults({
        result_ids: params.resultIds,
        approved: true,
        notes: params.notes
      });
      
      return {
        successCount: params.resultIds.length,
        failures: []
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to batch approve results');
    }
  }
);

/**
 * Async thunk for batch rejecting multiple results
 */
export const batchRejectResults = createAsyncThunk(
  'results/batchRejectResults',
  async (params: { resultIds: string[], notes?: string }, { rejectWithValue }) => {
    try {
      const response = await batchApproveRejectResults({
        result_ids: params.resultIds,
        approved: false,
        notes: params.notes
      });
      
      return {
        successCount: params.resultIds.length,
        failures: []
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to batch reject results');
    }
  }
);

/**
 * Async thunk for exporting result data in various formats
 */
export const exportResultData = createAsyncThunk(
  'results/exportResultData',
  async (params: { resultId: string, format: 'csv' | 'xlsx' | 'json' }, { rejectWithValue }) => {
    try {
      const blob = await exportResult(params);
      return blob;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to export result data');
    }
  }
);

/**
 * Async thunk for analyzing result data
 */
export const analyzeResultData = createAsyncThunk(
  'results/analyzeResultData',
  async (resultId: string, { rejectWithValue }) => {
    try {
      const response = await analyzeResult(resultId);
      return response.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to analyze result data');
    }
  }
);

/**
 * Redux Toolkit slice for results state management
 */
export const resultsSlice = createSlice({
  name: 'results',
  initialState,
  reducers: {
    /**
     * Sets the filter for results
     */
    setFilter: (state, action: PayloadAction<Partial<ResultFilter>>) => {
      state.filter = {
        ...state.filter,
        ...action.payload
      };
      // Reset to first page when filter changes
      state.currentPage = 1;
    },
    
    /**
     * Sets the current page number
     */
    setCurrentPage: (state, action: PayloadAction<number>) => {
      state.currentPage = action.payload;
    },
    
    /**
     * Sets the page size
     */
    setPageSize: (state, action: PayloadAction<number>) => {
      state.pageSize = action.payload;
      // Reset to first page when page size changes
      state.currentPage = 1;
    },
    
    /**
     * Clears any error message
     */
    clearResultError: (state) => {
      state.error = null;
    }
  },
  extraReducers: (builder) => {
    // Fetch Results
    builder.addCase(fetchResults.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    
    builder.addCase(fetchResults.fulfilled, (state, action) => {
      state.loading = false;
      state.results = action.payload.results;
      state.totalResults = action.payload.total;
    });
    
    builder.addCase(fetchResults.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to fetch results';
    });
    
    // Fetch Results By Submission
    builder.addCase(fetchResultsBySubmission.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    
    builder.addCase(fetchResultsBySubmission.fulfilled, (state, action) => {
      state.loading = false;
      state.results = action.payload.results;
      state.totalResults = action.payload.total;
    });
    
    builder.addCase(fetchResultsBySubmission.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to fetch results for submission';
    });
    
    // Fetch Result Detail
    builder.addCase(fetchResultDetail.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    
    builder.addCase(fetchResultDetail.fulfilled, (state, action) => {
      state.loading = false;
      state.currentResult = action.payload;
    });
    
    builder.addCase(fetchResultDetail.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to fetch result details';
    });
    
    // Approve Result
    builder.addCase(approveResult.pending, (state) => {
      state.loading = true;
    });
    
    builder.addCase(approveResult.fulfilled, (state, action) => {
      state.loading = false;
      state.currentResult = action.payload;
      
      // Update the result in the results array if it exists
      const index = state.results.findIndex(result => result.id === action.payload.id);
      if (index !== -1) {
        state.results[index] = action.payload;
      }
    });
    
    builder.addCase(approveResult.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to approve result';
    });
    
    // Reject Result
    builder.addCase(rejectResult.pending, (state) => {
      state.loading = true;
    });
    
    builder.addCase(rejectResult.fulfilled, (state, action) => {
      state.loading = false;
      state.currentResult = action.payload;
      
      // Update the result in the results array if it exists
      const index = state.results.findIndex(result => result.id === action.payload.id);
      if (index !== -1) {
        state.results[index] = action.payload;
      }
    });
    
    builder.addCase(rejectResult.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to reject result';
    });
    
    // Batch Approve Results
    builder.addCase(batchApproveResults.pending, (state) => {
      state.loading = true;
    });
    
    builder.addCase(batchApproveResults.fulfilled, (state) => {
      state.loading = false;
      // Note: This doesn't update individual results in state,
      // typically would trigger a refetch via component after batch operation
    });
    
    builder.addCase(batchApproveResults.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to batch approve results';
    });
    
    // Batch Reject Results
    builder.addCase(batchRejectResults.pending, (state) => {
      state.loading = true;
    });
    
    builder.addCase(batchRejectResults.fulfilled, (state) => {
      state.loading = false;
      // Note: This doesn't update individual results in state,
      // typically would trigger a refetch via component after batch operation
    });
    
    builder.addCase(batchRejectResults.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to batch reject results';
    });
    
    // Analyze Result Data
    builder.addCase(analyzeResultData.pending, (state) => {
      state.loading = true;
    });
    
    builder.addCase(analyzeResultData.fulfilled, (state, action) => {
      state.loading = false;
      if (state.currentResult) {
        state.currentResult = {
          ...state.currentResult,
          analysis: action.payload
        };
      }
    });
    
    builder.addCase(analyzeResultData.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Failed to analyze result data';
    });
    
    // Note: exportResultData is handled differently since it downloads a file
    // and doesn't directly update state, so no extra reducers needed for it
  }
});

// Extract actions from the slice
export const { setFilter, setCurrentPage, setPageSize, clearResultError } = resultsSlice.actions;

// Define selectors
export const selectResults = (state: { results: ResultState }) => state.results.results;
export const selectCurrentResult = (state: { results: ResultState }) => state.results.currentResult;
export const selectResultFilter = (state: { results: ResultState }) => state.results.filter;
export const selectResultLoading = (state: { results: ResultState }) => state.results.loading;
export const selectResultError = (state: { results: ResultState }) => state.results.error;
export const selectTotalResults = (state: { results: ResultState }) => state.results.totalResults;
export const selectCurrentPage = (state: { results: ResultState }) => state.results.currentPage;
export const selectPageSize = (state: { results: ResultState }) => state.results.pageSize;

// Export the reducer
export default resultsSlice.reducer;