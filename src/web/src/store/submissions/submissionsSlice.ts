/**
 * Redux Toolkit slice for submission state management in the Molecular Data Management and CRO Integration Platform.
 * Handles fetching, creating, updating, and managing CRO submissions, including quote handling and status transitions.
 *
 * @packageDocumentation
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import {
  SubmissionState,
  Submission,
  SubmissionDetailed,
  SubmissionCreate,
  SubmissionFilter,
  SubmissionStatus,
  QuoteProvide,
  QuoteResponse,
  SubmissionStatusUpdate
} from '../../types/submission';
import {
  getSubmissions,
  getSubmission,
  createSubmission,
  updateSubmission,
  updateSubmissionStatus,
  provideQuote,
  respondToQuote,
  getSubmissionsByExperiment,
  getSubmissionsByCRO,
  getSubmissionsByStatus,
  cancelSubmission
} from '../../api/submissions';

/**
 * Initial state for the submissions slice
 */
const initialState: SubmissionState = {
  submissions: [],
  currentSubmission: null,
  loading: false,
  error: null,
  totalSubmissions: 0,
  currentPage: 1,
  pageSize: 10,
  filter: {}
};

/**
 * Async thunk that fetches submissions with optional filtering
 */
export const fetchSubmissions = createAsyncThunk(
  'submissions/fetchSubmissions',
  async (filter: SubmissionFilter = {}, thunkAPI) => {
    try {
      const state = thunkAPI.getState() as { submissions: SubmissionState };
      const { currentPage, pageSize } = state.submissions;
      
      const response = await getSubmissions({
        ...filter,
        page: currentPage,
        page_size: pageSize
      });
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to fetch submissions');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that fetches a single submission by ID
 */
export const fetchSubmissionById = createAsyncThunk(
  'submissions/fetchSubmissionById',
  async (id: string, thunkAPI) => {
    try {
      const response = await getSubmission(id);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to fetch submission');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that fetches submissions for a specific experiment
 */
export const fetchSubmissionsByExperiment = createAsyncThunk(
  'submissions/fetchSubmissionsByExperiment',
  async (experimentId: string, thunkAPI) => {
    try {
      const state = thunkAPI.getState() as { submissions: SubmissionState };
      const { filter, currentPage, pageSize } = state.submissions;
      
      const response = await getSubmissionsByExperiment(experimentId, {
        ...filter,
        page: currentPage,
        page_size: pageSize
      });
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to fetch submissions for experiment');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that fetches submissions assigned to the current CRO user
 */
export const fetchSubmissionsByCRO = createAsyncThunk(
  'submissions/fetchSubmissionsByCRO',
  async (filter: SubmissionFilter = {}, thunkAPI) => {
    try {
      const state = thunkAPI.getState() as { submissions: SubmissionState };
      const { currentPage, pageSize } = state.submissions;
      
      const response = await getSubmissionsByCRO({
        ...filter,
        page: currentPage,
        page_size: pageSize
      });
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to fetch CRO submissions');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that fetches submissions with a specific status
 */
export const fetchSubmissionsByStatus = createAsyncThunk(
  'submissions/fetchSubmissionsByStatus',
  async (payload: { status: string, filter?: SubmissionFilter }, thunkAPI) => {
    try {
      const { status, filter = {} } = payload;
      const state = thunkAPI.getState() as { submissions: SubmissionState };
      const { currentPage, pageSize } = state.submissions;
      
      const response = await getSubmissionsByStatus(status, {
        ...filter,
        page: currentPage,
        page_size: pageSize
      });
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || `Failed to fetch submissions with status ${status}`);
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that creates a new submission
 */
export const createNewSubmission = createAsyncThunk(
  'submissions/createNewSubmission',
  async (submissionData: SubmissionCreate, thunkAPI) => {
    try {
      const response = await createSubmission(submissionData);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to create submission');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that updates an existing submission
 */
export const updateExistingSubmission = createAsyncThunk(
  'submissions/updateExistingSubmission',
  async (updateData: { id: string, data: any }, thunkAPI) => {
    try {
      const { id, data } = updateData;
      const response = await updateSubmission(id, data);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to update submission');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that updates the status of a submission
 */
export const updateSubmissionStatusThunk = createAsyncThunk(
  'submissions/updateSubmissionStatus',
  async (payload: { id: string, status: string, notes?: string }, thunkAPI) => {
    try {
      const { id, status, notes } = payload;
      const statusData: SubmissionStatusUpdate = { status, notes };
      
      const response = await updateSubmissionStatus(id, statusData);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to update submission status');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that provides a quote for a submission
 */
export const provideQuoteForSubmission = createAsyncThunk(
  'submissions/provideQuoteForSubmission',
  async (payload: { id: string, quoteData: QuoteProvide }, thunkAPI) => {
    try {
      const { id, quoteData } = payload;
      const response = await provideQuote(id, quoteData);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to provide quote');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that responds to a quote for a submission
 */
export const respondToSubmissionQuote = createAsyncThunk(
  'submissions/respondToSubmissionQuote',
  async (payload: { id: string, responseData: QuoteResponse }, thunkAPI) => {
    try {
      const { id, responseData } = payload;
      const response = await respondToQuote(id, responseData);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to respond to quote');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Async thunk that cancels a submission
 */
export const cancelExistingSubmission = createAsyncThunk(
  'submissions/cancelExistingSubmission',
  async (payload: { id: string, notes?: string }, thunkAPI) => {
    try {
      const { id, notes } = payload;
      const response = await cancelSubmission(id, notes);
      
      if (!response.success) {
        return thunkAPI.rejectWithValue(response.error || 'Failed to cancel submission');
      }
      
      return response.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'An unknown error occurred');
    }
  }
);

/**
 * Submissions slice definition containing reducers and extra reducers for async thunks
 */
const submissionsSlice = createSlice({
  name: 'submissions',
  initialState,
  reducers: {
    setSubmissionFilter: (state, action: PayloadAction<SubmissionFilter>) => {
      state.filter = action.payload;
      // Reset to first page when filter changes
      state.currentPage = 1;
    },
    setCurrentPage: (state, action: PayloadAction<number>) => {
      state.currentPage = action.payload;
    },
    setPageSize: (state, action: PayloadAction<number>) => {
      state.pageSize = action.payload;
      // Reset to first page when page size changes
      state.currentPage = 1;
    },
    clearSubmissionError: (state) => {
      state.error = null;
    }
  },
  extraReducers: (builder) => {
    // fetchSubmissions
    builder.addCase(fetchSubmissions.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchSubmissions.fulfilled, (state, action) => {
      state.loading = false;
      state.submissions = action.payload.items;
      state.totalSubmissions = action.payload.total;
    });
    builder.addCase(fetchSubmissions.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // fetchSubmissionById
    builder.addCase(fetchSubmissionById.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchSubmissionById.fulfilled, (state, action) => {
      state.loading = false;
      state.currentSubmission = action.payload;
    });
    builder.addCase(fetchSubmissionById.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // fetchSubmissionsByExperiment
    builder.addCase(fetchSubmissionsByExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchSubmissionsByExperiment.fulfilled, (state, action) => {
      state.loading = false;
      state.submissions = action.payload.items;
      state.totalSubmissions = action.payload.total;
    });
    builder.addCase(fetchSubmissionsByExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // fetchSubmissionsByCRO
    builder.addCase(fetchSubmissionsByCRO.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchSubmissionsByCRO.fulfilled, (state, action) => {
      state.loading = false;
      state.submissions = action.payload.items;
      state.totalSubmissions = action.payload.total;
    });
    builder.addCase(fetchSubmissionsByCRO.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // fetchSubmissionsByStatus
    builder.addCase(fetchSubmissionsByStatus.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchSubmissionsByStatus.fulfilled, (state, action) => {
      state.loading = false;
      state.submissions = action.payload.items;
      state.totalSubmissions = action.payload.total;
    });
    builder.addCase(fetchSubmissionsByStatus.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // createNewSubmission
    builder.addCase(createNewSubmission.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(createNewSubmission.fulfilled, (state, action) => {
      state.loading = false;
      state.submissions.push(action.payload);
      state.totalSubmissions += 1;
    });
    builder.addCase(createNewSubmission.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // updateExistingSubmission
    builder.addCase(updateExistingSubmission.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateExistingSubmission.fulfilled, (state, action) => {
      state.loading = false;
      // Update in the list
      const index = state.submissions.findIndex(s => s.id === action.payload.id);
      if (index !== -1) {
        state.submissions[index] = action.payload;
      }
      // Update current submission if it's the same one
      if (state.currentSubmission && state.currentSubmission.id === action.payload.id) {
        state.currentSubmission = {
          ...state.currentSubmission,
          ...action.payload
        };
      }
    });
    builder.addCase(updateExistingSubmission.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // updateSubmissionStatusThunk
    builder.addCase(updateSubmissionStatusThunk.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateSubmissionStatusThunk.fulfilled, (state, action) => {
      state.loading = false;
      // Update in the list
      const index = state.submissions.findIndex(s => s.id === action.payload.id);
      if (index !== -1) {
        state.submissions[index] = action.payload;
      }
      // Update current submission if it's the same one
      if (state.currentSubmission && state.currentSubmission.id === action.payload.id) {
        state.currentSubmission = {
          ...state.currentSubmission,
          ...action.payload
        };
      }
    });
    builder.addCase(updateSubmissionStatusThunk.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // provideQuoteForSubmission
    builder.addCase(provideQuoteForSubmission.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(provideQuoteForSubmission.fulfilled, (state, action) => {
      state.loading = false;
      // Update in the list
      const index = state.submissions.findIndex(s => s.id === action.payload.id);
      if (index !== -1) {
        state.submissions[index] = action.payload;
      }
      // Update current submission if it's the same one
      if (state.currentSubmission && state.currentSubmission.id === action.payload.id) {
        state.currentSubmission = {
          ...state.currentSubmission,
          ...action.payload
        };
      }
    });
    builder.addCase(provideQuoteForSubmission.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // respondToSubmissionQuote
    builder.addCase(respondToSubmissionQuote.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(respondToSubmissionQuote.fulfilled, (state, action) => {
      state.loading = false;
      // Update in the list
      const index = state.submissions.findIndex(s => s.id === action.payload.id);
      if (index !== -1) {
        state.submissions[index] = action.payload;
      }
      // Update current submission if it's the same one
      if (state.currentSubmission && state.currentSubmission.id === action.payload.id) {
        state.currentSubmission = {
          ...state.currentSubmission,
          ...action.payload
        };
      }
    });
    builder.addCase(respondToSubmissionQuote.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // cancelExistingSubmission
    builder.addCase(cancelExistingSubmission.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(cancelExistingSubmission.fulfilled, (state, action) => {
      state.loading = false;
      // Update in the list
      const index = state.submissions.findIndex(s => s.id === action.payload.id);
      if (index !== -1) {
        state.submissions[index] = action.payload;
      }
      // Update current submission if it's the same one
      if (state.currentSubmission && state.currentSubmission.id === action.payload.id) {
        state.currentSubmission = {
          ...state.currentSubmission,
          ...action.payload
        };
      }
    });
    builder.addCase(cancelExistingSubmission.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });
  }
});

// Export actions and reducer
export const { setSubmissionFilter, setCurrentPage, setPageSize, clearSubmissionError } = submissionsSlice.actions;
export const submissionsActions = submissionsSlice.actions;

export default submissionsSlice.reducer;