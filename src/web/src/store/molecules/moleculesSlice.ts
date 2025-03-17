/**
 * Redux slice for managing molecule data in the Molecular Data Management and CRO Integration Platform.
 * Handles state related to molecules, including fetching, filtering, selecting, and managing molecule data.
 * 
 * @packageDocumentation
 * @module store/molecules
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import { Molecule, MoleculeFilter, MoleculeState } from '../../types/molecule';
import { 
  getMolecules, 
  updateMoleculeFlag,
  bulkUpdateMolecules
} from '../../api/molecules';

/**
 * Initial state for the molecules slice
 */
const initialState: MoleculeState = {
  molecules: [],
  filter: {},
  selectedMolecules: [],
  loading: false,
  error: null,
  totalMolecules: 0,
  currentPage: 1,
  pageSize: 10
};

/**
 * Async thunk for fetching molecules with filtering and pagination
 */
export const fetchMolecules = createAsyncThunk(
  'molecules/fetchMolecules',
  async ({ filter, page, pageSize }: { filter: MoleculeFilter, page: number, pageSize: number }, { rejectWithValue }) => {
    try {
      const response = await getMolecules(filter, page, pageSize);
      return { 
        molecules: response.data.data.items, 
        total: response.data.data.total 
      };
    } catch (error: any) {
      return rejectWithValue(error.error || error.message || 'Failed to fetch molecules');
    }
  }
);

/**
 * Async thunk for updating a molecule's flag status
 */
export const updateMoleculeFlagStatus = createAsyncThunk(
  'molecules/updateFlagStatus',
  async ({ id, flagStatus }: { id: string, flagStatus: string }, { rejectWithValue }) => {
    try {
      const response = await updateMoleculeFlag(id, flagStatus);
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || error.message || 'Failed to update molecule flag');
    }
  }
);

/**
 * Async thunk for updating flag status of multiple molecules
 */
export const bulkUpdateMoleculeFlags = createAsyncThunk(
  'molecules/bulkUpdateFlags',
  async ({ moleculeIds, flagStatus }: { moleculeIds: string[], flagStatus: string }, { rejectWithValue }) => {
    try {
      const response = await bulkUpdateMolecules({
        operation_type: 'update_flag',
        molecule_ids: moleculeIds,
        flag_status: flagStatus
      });
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || error.message || 'Failed to update molecule flags');
    }
  }
);

/**
 * Async thunk for adding molecules to a library
 */
export const addMoleculesToLibrary = createAsyncThunk(
  'molecules/addToLibrary',
  async ({ moleculeIds, libraryId }: { moleculeIds: string[], libraryId: string }, { rejectWithValue }) => {
    try {
      const response = await bulkUpdateMolecules({
        operation_type: 'add_to_library',
        molecule_ids: moleculeIds,
        library_id: libraryId
      });
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || error.message || 'Failed to add molecules to library');
    }
  }
);

/**
 * Async thunk for adding molecules to an experiment
 */
export const addMoleculesToExperiment = createAsyncThunk(
  'molecules/addToExperiment',
  async ({ moleculeIds, experimentId }: { moleculeIds: string[], experimentId: string }, { rejectWithValue }) => {
    try {
      const response = await bulkUpdateMolecules({
        operation_type: 'add_to_experiment',
        molecule_ids: moleculeIds,
        experiment_id: experimentId
      });
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || error.message || 'Failed to add molecules to experiment');
    }
  }
);

/**
 * Redux slice for molecules
 * Contains reducers for molecule state management and handling async action results
 */
export const moleculesSlice = createSlice({
  name: 'molecules',
  initialState,
  reducers: {
    setFilter: (state, action: PayloadAction<MoleculeFilter>) => {
      state.filter = action.payload;
      state.currentPage = 1; // Reset to first page when filter changes
    },
    setSelectedMolecules: (state, action: PayloadAction<string[]>) => {
      state.selectedMolecules = action.payload;
    },
    addSelectedMolecule: (state, action: PayloadAction<string>) => {
      if (!state.selectedMolecules.includes(action.payload)) {
        state.selectedMolecules.push(action.payload);
      }
    },
    removeSelectedMolecule: (state, action: PayloadAction<string>) => {
      state.selectedMolecules = state.selectedMolecules.filter(id => id !== action.payload);
    },
    clearSelectedMolecules: (state) => {
      state.selectedMolecules = [];
    },
    setCurrentPage: (state, action: PayloadAction<number>) => {
      state.currentPage = action.payload;
    },
    setPageSize: (state, action: PayloadAction<number>) => {
      state.pageSize = action.payload;
      state.currentPage = 1; // Reset to first page when page size changes
    }
  },
  extraReducers: (builder) => {
    builder
      // fetchMolecules
      .addCase(fetchMolecules.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(fetchMolecules.fulfilled, (state, action) => {
        state.loading = false;
        state.molecules = action.payload.molecules;
        state.totalMolecules = action.payload.total;
      })
      .addCase(fetchMolecules.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to fetch molecules';
      })
      
      // updateMoleculeFlagStatus
      .addCase(updateMoleculeFlagStatus.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(updateMoleculeFlagStatus.fulfilled, (state, action) => {
        state.loading = false;
        // Update the molecule in the state
        const index = state.molecules.findIndex(molecule => molecule.id === action.payload.id);
        if (index !== -1) {
          state.molecules[index] = action.payload;
        }
      })
      .addCase(updateMoleculeFlagStatus.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to update molecule flag';
      })
      
      // bulkUpdateMoleculeFlags
      .addCase(bulkUpdateMoleculeFlags.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(bulkUpdateMoleculeFlags.fulfilled, (state) => {
        state.loading = false;
        // We'll need to refresh the molecules to see the updated flags
        // This will happen when the component re-fetches the data
      })
      .addCase(bulkUpdateMoleculeFlags.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to update molecule flags';
      })
      
      // addMoleculesToLibrary
      .addCase(addMoleculesToLibrary.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(addMoleculesToLibrary.fulfilled, (state) => {
        state.loading = false;
        // No need to update state as this operation doesn't change molecule data
        // Component should refetch data if needed
      })
      .addCase(addMoleculesToLibrary.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to add molecules to library';
      })
      
      // addMoleculesToExperiment
      .addCase(addMoleculesToExperiment.pending, (state) => {
        state.loading = true;
        state.error = null;
      })
      .addCase(addMoleculesToExperiment.fulfilled, (state) => {
        state.loading = false;
        // No need to update state as this operation doesn't change molecule data
        // Component should refetch data if needed
      })
      .addCase(addMoleculesToExperiment.rejected, (state, action) => {
        state.loading = false;
        state.error = action.payload as string || 'Failed to add molecules to experiment';
      })
  }
});

// Export selectors
export const selectMolecules = (state: { molecules: MoleculeState }) => state.molecules.molecules;
export const selectMoleculeFilter = (state: { molecules: MoleculeState }) => state.molecules.filter;
export const selectSelectedMolecules = (state: { molecules: MoleculeState }) => state.molecules.selectedMolecules;
export const selectMoleculeLoading = (state: { molecules: MoleculeState }) => state.molecules.loading;
export const selectMoleculeError = (state: { molecules: MoleculeState }) => state.molecules.error;
export const selectTotalMolecules = (state: { molecules: MoleculeState }) => state.molecules.totalMolecules;
export const selectCurrentPage = (state: { molecules: MoleculeState }) => state.molecules.currentPage;
export const selectPageSize = (state: { molecules: MoleculeState }) => state.molecules.pageSize;

// Export actions
export const { 
  setFilter, 
  setSelectedMolecules, 
  addSelectedMolecule, 
  removeSelectedMolecule, 
  clearSelectedMolecules,
  setCurrentPage,
  setPageSize
} = moleculesSlice.actions;

// Export reducer as default
export default moleculesSlice.reducer;