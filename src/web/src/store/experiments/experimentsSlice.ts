/**
 * Redux Toolkit slice for experiment state management in the Molecular Data Management and CRO Integration Platform.
 * This slice handles fetching, creating, updating, and managing experiments, including adding molecules to experiments
 * and updating experiment statuses.
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import {
  ExperimentState,
  Experiment,
  ExperimentDetail,
  ExperimentCreate,
  ExperimentUpdate,
  ExperimentFilter,
  ExperimentMoleculeOperation,
  ExperimentStatusUpdate,
  ExperimentType
} from '../../types/experiment';
import {
  getExperiments,
  getMyExperiments,
  getExperimentById,
  getExperimentTypes,
  createExperiment,
  updateExperiment,
  deleteExperiment,
  addMoleculesToExperiment,
  removeMoleculesFromExperiment,
  updateExperimentStatus,
  queueExperiment,
  submitExperiment,
  cancelExperiment
} from '../../api/experiments';

// Fetch experiments with optional filtering
export const fetchExperiments = createAsyncThunk(
  'experiments/fetchExperiments',
  async (filter: ExperimentFilter, thunkAPI) => {
    try {
      const state = thunkAPI.getState() as { experiments: ExperimentState };
      const { currentPage, pageSize } = state.experiments;
      
      const response = await getExperiments(filter, currentPage, pageSize);
      return {
        experiments: response.data.data.items,
        total: response.data.data.total
      };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch experiments');
    }
  }
);

// Fetch experiments created by the current user
export const fetchMyExperiments = createAsyncThunk(
  'experiments/fetchMyExperiments',
  async (filter: ExperimentFilter, thunkAPI) => {
    try {
      const state = thunkAPI.getState() as { experiments: ExperimentState };
      const { currentPage, pageSize } = state.experiments;
      
      const response = await getMyExperiments(filter, currentPage, pageSize);
      return {
        experiments: response.data.data.items,
        total: response.data.data.total
      };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch user experiments');
    }
  }
);

// Fetch a single experiment by ID
export const fetchExperimentById = createAsyncThunk(
  'experiments/fetchExperimentById',
  async (id: string, thunkAPI) => {
    try {
      const response = await getExperimentById(id);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch experiment details');
    }
  }
);

// Fetch available experiment types
export const fetchExperimentTypes = createAsyncThunk(
  'experiments/fetchExperimentTypes',
  async (_, thunkAPI) => {
    try {
      const response = await getExperimentTypes();
      return response.data.data.items;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to fetch experiment types');
    }
  }
);

// Create a new experiment
export const createNewExperiment = createAsyncThunk(
  'experiments/createNewExperiment',
  async (experimentData: ExperimentCreate, thunkAPI) => {
    try {
      const response = await createExperiment(experimentData);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to create experiment');
    }
  }
);

// Update an existing experiment
export const updateExistingExperiment = createAsyncThunk(
  'experiments/updateExistingExperiment',
  async ({ id, data }: { id: string, data: ExperimentUpdate }, thunkAPI) => {
    try {
      const response = await updateExperiment(id, data);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to update experiment');
    }
  }
);

// Delete an experiment
export const removeExperiment = createAsyncThunk(
  'experiments/removeExperiment',
  async (id: string, thunkAPI) => {
    try {
      const response = await deleteExperiment(id);
      return { id, success: response.data.data.success };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to delete experiment');
    }
  }
);

// Add molecules to an experiment
export const addMolecules = createAsyncThunk(
  'experiments/addMolecules',
  async ({ id, moleculeIds }: { id: string, moleculeIds: string[] }, thunkAPI) => {
    try {
      const moleculeData: ExperimentMoleculeOperation = {
        molecule_ids: moleculeIds
      };
      const response = await addMoleculesToExperiment(id, moleculeData);
      return { id, ...response.data.data };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to add molecules to experiment');
    }
  }
);

// Remove molecules from an experiment
export const removeMolecules = createAsyncThunk(
  'experiments/removeMolecules',
  async ({ id, moleculeIds }: { id: string, moleculeIds: string[] }, thunkAPI) => {
    try {
      const moleculeData: ExperimentMoleculeOperation = {
        molecule_ids: moleculeIds
      };
      const response = await removeMoleculesFromExperiment(id, moleculeData);
      return { id, ...response.data.data };
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to remove molecules from experiment');
    }
  }
);

// Update the status of an experiment
export const updateStatus = createAsyncThunk(
  'experiments/updateStatus',
  async ({ id, status, notes }: { id: string, status: string, notes?: string }, thunkAPI) => {
    try {
      const statusData: ExperimentStatusUpdate = { status, notes };
      const response = await updateExperimentStatus(id, statusData);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to update experiment status');
    }
  }
);

// Queue an experiment
export const queueExistingExperiment = createAsyncThunk(
  'experiments/queueExistingExperiment',
  async ({ id, notes }: { id: string, notes?: string }, thunkAPI) => {
    try {
      const response = await queueExperiment(id, notes);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to queue experiment');
    }
  }
);

// Submit an experiment to CRO
export const submitExistingExperiment = createAsyncThunk(
  'experiments/submitExistingExperiment',
  async ({ id, notes }: { id: string, notes?: string }, thunkAPI) => {
    try {
      const response = await submitExperiment(id, notes);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to submit experiment');
    }
  }
);

// Cancel an experiment
export const cancelExistingExperiment = createAsyncThunk(
  'experiments/cancelExistingExperiment',
  async ({ id, notes }: { id: string, notes?: string }, thunkAPI) => {
    try {
      const response = await cancelExperiment(id, notes);
      return response.data.data;
    } catch (error) {
      return thunkAPI.rejectWithValue(error instanceof Error ? error.message : 'Failed to cancel experiment');
    }
  }
);

// Initial state for the experiments slice
const initialState: ExperimentState = {
  experiments: [],
  experimentTypes: [],
  currentExperiment: null,
  loading: false,
  error: null,
  totalExperiments: 0,
  currentPage: 1,
  pageSize: 10,
  filter: {}
};

// Create the experiments slice
const experimentsSlice = createSlice({
  name: 'experiments',
  initialState,
  reducers: {
    setExperimentFilter(state, action: PayloadAction<ExperimentFilter>) {
      state.filter = action.payload;
    },
    setCurrentPage(state, action: PayloadAction<number>) {
      state.currentPage = action.payload;
    },
    setPageSize(state, action: PayloadAction<number>) {
      state.pageSize = action.payload;
    },
    clearExperimentError(state) {
      state.error = null;
    }
  },
  extraReducers: (builder) => {
    // Fetch experiments
    builder.addCase(fetchExperiments.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchExperiments.fulfilled, (state, action) => {
      state.experiments = action.payload.experiments;
      state.totalExperiments = action.payload.total;
      state.loading = false;
    });
    builder.addCase(fetchExperiments.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Fetch my experiments
    builder.addCase(fetchMyExperiments.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchMyExperiments.fulfilled, (state, action) => {
      state.experiments = action.payload.experiments;
      state.totalExperiments = action.payload.total;
      state.loading = false;
    });
    builder.addCase(fetchMyExperiments.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Fetch experiment by ID
    builder.addCase(fetchExperimentById.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchExperimentById.fulfilled, (state, action) => {
      state.currentExperiment = action.payload;
      state.loading = false;
    });
    builder.addCase(fetchExperimentById.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Fetch experiment types
    builder.addCase(fetchExperimentTypes.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchExperimentTypes.fulfilled, (state, action) => {
      state.experimentTypes = action.payload;
      state.loading = false;
    });
    builder.addCase(fetchExperimentTypes.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Create new experiment
    builder.addCase(createNewExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(createNewExperiment.fulfilled, (state, action) => {
      state.experiments.push(action.payload);
      state.loading = false;
    });
    builder.addCase(createNewExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Update existing experiment
    builder.addCase(updateExistingExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateExistingExperiment.fulfilled, (state, action) => {
      const index = state.experiments.findIndex(exp => exp.id === action.payload.id);
      if (index !== -1) {
        state.experiments[index] = action.payload;
      }
      if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
        state.currentExperiment = { ...state.currentExperiment, ...action.payload };
      }
      state.loading = false;
    });
    builder.addCase(updateExistingExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Remove experiment
    builder.addCase(removeExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(removeExperiment.fulfilled, (state, action) => {
      if (action.payload.success) {
        state.experiments = state.experiments.filter(exp => exp.id !== action.payload.id);
        if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
          state.currentExperiment = null;
        }
      }
      state.loading = false;
    });
    builder.addCase(removeExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Add molecules to experiment
    builder.addCase(addMolecules.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(addMolecules.fulfilled, (state, action) => {
      // Note: We don't update the currentExperiment here because we don't have the full updated molecules list
      // Components should dispatch fetchExperimentById after this operation to get the updated data
      state.loading = false;
    });
    builder.addCase(addMolecules.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Remove molecules from experiment
    builder.addCase(removeMolecules.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(removeMolecules.fulfilled, (state, action) => {
      // Note: We don't update the currentExperiment here because we don't have the full updated molecules list
      // Components should dispatch fetchExperimentById after this operation to get the updated data
      state.loading = false;
    });
    builder.addCase(removeMolecules.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Update experiment status
    builder.addCase(updateStatus.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateStatus.fulfilled, (state, action) => {
      const index = state.experiments.findIndex(exp => exp.id === action.payload.id);
      if (index !== -1) {
        state.experiments[index] = action.payload;
      }
      if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
        state.currentExperiment = { ...state.currentExperiment, ...action.payload };
      }
      state.loading = false;
    });
    builder.addCase(updateStatus.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Queue experiment
    builder.addCase(queueExistingExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(queueExistingExperiment.fulfilled, (state, action) => {
      const index = state.experiments.findIndex(exp => exp.id === action.payload.id);
      if (index !== -1) {
        state.experiments[index] = action.payload;
      }
      if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
        state.currentExperiment = { ...state.currentExperiment, ...action.payload };
      }
      state.loading = false;
    });
    builder.addCase(queueExistingExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Submit experiment
    builder.addCase(submitExistingExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(submitExistingExperiment.fulfilled, (state, action) => {
      const index = state.experiments.findIndex(exp => exp.id === action.payload.id);
      if (index !== -1) {
        state.experiments[index] = action.payload;
      }
      if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
        state.currentExperiment = { ...state.currentExperiment, ...action.payload };
      }
      state.loading = false;
    });
    builder.addCase(submitExistingExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });

    // Cancel experiment
    builder.addCase(cancelExistingExperiment.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(cancelExistingExperiment.fulfilled, (state, action) => {
      const index = state.experiments.findIndex(exp => exp.id === action.payload.id);
      if (index !== -1) {
        state.experiments[index] = action.payload;
      }
      if (state.currentExperiment && state.currentExperiment.id === action.payload.id) {
        state.currentExperiment = { ...state.currentExperiment, ...action.payload };
      }
      state.loading = false;
    });
    builder.addCase(cancelExistingExperiment.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string;
    });
  }
});

// Export actions and reducer
export const { setExperimentFilter, setCurrentPage, setPageSize, clearExperimentError } = experimentsSlice.actions;
export const experimentsActions = experimentsSlice.actions;
export default experimentsSlice.reducer;