import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import { 
  Library, 
  LibraryDetail, 
  LibraryFilter, 
  LibraryState, 
  LibraryCreate, 
  LibraryUpdate,
  LibraryOperationType
} from '../../types/library';
import { 
  getLibraries,
  getMyLibraries,
  getLibraryById,
  getLibraryDetailById,
  createLibrary,
  updateLibrary,
  deleteLibrary,
  addMoleculesToLibrary,
  removeMoleculesFromLibrary,
  updateLibraryMolecules
} from '../../api/libraries';

// Initial state
const initialState: LibraryState = {
  libraries: [],
  currentLibrary: null,
  loading: false,
  error: null,
  totalLibraries: 0,
  currentPage: 1,
  pageSize: 10
};

/**
 * Async thunk for fetching libraries with filtering and pagination
 */
export const fetchLibraries = createAsyncThunk(
  'libraries/fetchLibraries',
  async (params: { filter?: LibraryFilter; page?: number; pageSize?: number }, { rejectWithValue }) => {
    try {
      const response = await getLibraries(params.filter, params.page, params.pageSize);
      return {
        libraries: response.data.data.items,
        total: response.data.data.total
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch libraries');
    }
  }
);

/**
 * Async thunk for fetching libraries created by the current user
 */
export const fetchMyLibraries = createAsyncThunk(
  'libraries/fetchMyLibraries',
  async (params: { page?: number; pageSize?: number }, { rejectWithValue }) => {
    try {
      const response = await getMyLibraries(params.page, params.pageSize);
      return {
        libraries: response.data.data.items,
        total: response.data.data.total
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch my libraries');
    }
  }
);

/**
 * Async thunk for fetching a single library by ID
 */
export const fetchLibraryById = createAsyncThunk(
  'libraries/fetchLibraryById',
  async (id: string, { rejectWithValue }) => {
    try {
      const response = await getLibraryById(id);
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch library');
    }
  }
);

/**
 * Async thunk for fetching detailed library information including contained molecules
 */
export const fetchLibraryDetailById = createAsyncThunk(
  'libraries/fetchLibraryDetailById',
  async (id: string, { rejectWithValue }) => {
    try {
      const response = await getLibraryDetailById(id);
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to fetch library details');
    }
  }
);

/**
 * Async thunk for creating a new library
 */
export const createNewLibrary = createAsyncThunk(
  'libraries/createNewLibrary',
  async (libraryData: LibraryCreate, { rejectWithValue }) => {
    try {
      const response = await createLibrary(libraryData);
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to create library');
    }
  }
);

/**
 * Async thunk for updating an existing library
 */
export const updateExistingLibrary = createAsyncThunk(
  'libraries/updateExistingLibrary',
  async (params: { id: string; data: LibraryUpdate }, { rejectWithValue }) => {
    try {
      const response = await updateLibrary(params.id, params.data);
      return response.data.data;
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to update library');
    }
  }
);

/**
 * Async thunk for deleting a library
 */
export const removeLibrary = createAsyncThunk(
  'libraries/removeLibrary',
  async (id: string, { rejectWithValue }) => {
    try {
      const response = await deleteLibrary(id);
      return { id, success: response.data.data.success };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to delete library');
    }
  }
);

/**
 * Async thunk for adding molecules to a library
 */
export const addMoleculesToExistingLibrary = createAsyncThunk(
  'libraries/addMoleculesToExistingLibrary',
  async (params: { libraryId: string; moleculeIds: string[] }, { rejectWithValue }) => {
    try {
      const response = await addMoleculesToLibrary(params.libraryId, params.moleculeIds);
      return { 
        libraryId: params.libraryId, 
        added_count: response.data.data.added_count 
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to add molecules to library');
    }
  }
);

/**
 * Async thunk for removing molecules from a library
 */
export const removeMoleculesFromExistingLibrary = createAsyncThunk(
  'libraries/removeMoleculesFromExistingLibrary',
  async (params: { libraryId: string; moleculeIds: string[] }, { rejectWithValue }) => {
    try {
      const response = await removeMoleculesFromLibrary(params.libraryId, params.moleculeIds);
      return { 
        libraryId: params.libraryId, 
        removed_count: response.data.data.removed_count 
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to remove molecules from library');
    }
  }
);

/**
 * Async thunk for performing add or remove operations on library molecules
 */
export const updateLibraryMoleculesOperation = createAsyncThunk(
  'libraries/updateLibraryMoleculesOperation',
  async (params: { 
    libraryId: string; 
    operation: LibraryOperationType;
    moleculeIds: string[] 
  }, { rejectWithValue }) => {
    try {
      const response = await updateLibraryMolecules(params.libraryId, {
        operation: params.operation,
        molecule_ids: params.moleculeIds
      });
      return { 
        libraryId: params.libraryId, 
        count: response.data.data.count,
        operation: params.operation
      };
    } catch (error: any) {
      return rejectWithValue(error.error || 'Failed to update library molecules');
    }
  }
);

// Create the slice
export const librariesSlice = createSlice({
  name: 'libraries',
  initialState,
  reducers: {
    setCurrentLibrary: (state, action: PayloadAction<LibraryDetail>) => {
      state.currentLibrary = action.payload;
    },
    clearCurrentLibrary: (state) => {
      state.currentLibrary = null;
    },
    setCurrentPage: (state, action: PayloadAction<number>) => {
      state.currentPage = action.payload;
    },
    setPageSize: (state, action: PayloadAction<number>) => {
      state.pageSize = action.payload;
    }
  },
  extraReducers: (builder) => {
    // fetchLibraries
    builder.addCase(fetchLibraries.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchLibraries.fulfilled, (state, action) => {
      state.loading = false;
      state.libraries = action.payload.libraries;
      state.totalLibraries = action.payload.total;
    });
    builder.addCase(fetchLibraries.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // fetchMyLibraries
    builder.addCase(fetchMyLibraries.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchMyLibraries.fulfilled, (state, action) => {
      state.loading = false;
      state.libraries = action.payload.libraries;
      state.totalLibraries = action.payload.total;
    });
    builder.addCase(fetchMyLibraries.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // fetchLibraryById
    builder.addCase(fetchLibraryById.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchLibraryById.fulfilled, (state, action) => {
      state.loading = false;
      
      // Update the library in the libraries array if it exists
      const index = state.libraries.findIndex(lib => lib.id === action.payload.id);
      if (index !== -1) {
        state.libraries[index] = action.payload;
      } else {
        // Add to the array if not found
        state.libraries.push(action.payload);
      }
    });
    builder.addCase(fetchLibraryById.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // fetchLibraryDetailById
    builder.addCase(fetchLibraryDetailById.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(fetchLibraryDetailById.fulfilled, (state, action) => {
      state.loading = false;
      state.currentLibrary = action.payload;
      
      // Also update the basic library information in the libraries array
      const basicLibrary: Library = {
        id: action.payload.id,
        name: action.payload.name,
        description: action.payload.description,
        created_by: action.payload.created_by,
        created_at: action.payload.created_at,
        updated_at: action.payload.updated_at,
        molecule_count: action.payload.molecules.length,
        creator: action.payload.creator
      };
      
      const index = state.libraries.findIndex(lib => lib.id === action.payload.id);
      if (index !== -1) {
        state.libraries[index] = basicLibrary;
      } else {
        // Add to the array if not found
        state.libraries.push(basicLibrary);
      }
    });
    builder.addCase(fetchLibraryDetailById.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // createNewLibrary
    builder.addCase(createNewLibrary.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(createNewLibrary.fulfilled, (state, action) => {
      state.loading = false;
      state.libraries.push(action.payload);
      state.totalLibraries += 1;
    });
    builder.addCase(createNewLibrary.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // updateExistingLibrary
    builder.addCase(updateExistingLibrary.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateExistingLibrary.fulfilled, (state, action) => {
      state.loading = false;
      
      // Update the library in the libraries array
      const index = state.libraries.findIndex(lib => lib.id === action.payload.id);
      if (index !== -1) {
        state.libraries[index] = {
          ...state.libraries[index],
          ...action.payload
        };
      }
      
      // Update the current library if it's the same one
      if (state.currentLibrary && state.currentLibrary.id === action.payload.id) {
        state.currentLibrary = {
          ...state.currentLibrary,
          ...action.payload
        };
      }
    });
    builder.addCase(updateExistingLibrary.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // removeLibrary
    builder.addCase(removeLibrary.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(removeLibrary.fulfilled, (state, action) => {
      state.loading = false;
      
      // Only remove if the operation was successful
      if (action.payload.success) {
        // Remove the library from the libraries array
        state.libraries = state.libraries.filter(lib => lib.id !== action.payload.id);
        state.totalLibraries -= 1;
        
        // Clear current library if it's the same one
        if (state.currentLibrary && state.currentLibrary.id === action.payload.id) {
          state.currentLibrary = null;
        }
      }
    });
    builder.addCase(removeLibrary.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // addMoleculesToExistingLibrary
    builder.addCase(addMoleculesToExistingLibrary.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(addMoleculesToExistingLibrary.fulfilled, (state, action) => {
      state.loading = false;
      
      // Update molecule count in the library
      const index = state.libraries.findIndex(lib => lib.id === action.payload.libraryId);
      if (index !== -1 && state.libraries[index].molecule_count !== undefined) {
        state.libraries[index].molecule_count += action.payload.added_count;
      }
      
      // If this is the current library, update its molecule count
      if (state.currentLibrary && state.currentLibrary.id === action.payload.libraryId) {
        state.currentLibrary.molecule_count = (state.currentLibrary.molecule_count || 0) + action.payload.added_count;
      }
    });
    builder.addCase(addMoleculesToExistingLibrary.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // removeMoleculesFromExistingLibrary
    builder.addCase(removeMoleculesFromExistingLibrary.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(removeMoleculesFromExistingLibrary.fulfilled, (state, action) => {
      state.loading = false;
      
      // Update molecule count in the library
      const index = state.libraries.findIndex(lib => lib.id === action.payload.libraryId);
      if (index !== -1 && state.libraries[index].molecule_count !== undefined) {
        state.libraries[index].molecule_count -= action.payload.removed_count;
      }
      
      // If this is the current library, update its molecule count
      if (state.currentLibrary && state.currentLibrary.id === action.payload.libraryId) {
        state.currentLibrary.molecule_count = (state.currentLibrary.molecule_count || 0) - action.payload.removed_count;
      }
    });
    builder.addCase(removeMoleculesFromExistingLibrary.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });

    // updateLibraryMoleculesOperation
    builder.addCase(updateLibraryMoleculesOperation.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(updateLibraryMoleculesOperation.fulfilled, (state, action) => {
      state.loading = false;
      
      // Update molecule count in the library based on operation type
      const index = state.libraries.findIndex(lib => lib.id === action.payload.libraryId);
      if (index !== -1 && state.libraries[index].molecule_count !== undefined) {
        if (action.payload.operation === LibraryOperationType.ADD) {
          state.libraries[index].molecule_count += action.payload.count;
        } else if (action.payload.operation === LibraryOperationType.REMOVE) {
          state.libraries[index].molecule_count -= action.payload.count;
        }
      }
      
      // Update current library count if it's the same one
      if (state.currentLibrary && state.currentLibrary.id === action.payload.libraryId) {
        if (action.payload.operation === LibraryOperationType.ADD) {
          state.currentLibrary.molecule_count = (state.currentLibrary.molecule_count || 0) + action.payload.count;
        } else if (action.payload.operation === LibraryOperationType.REMOVE) {
          state.currentLibrary.molecule_count = (state.currentLibrary.molecule_count || 0) - action.payload.count;
        }
      }
    });
    builder.addCase(updateLibraryMoleculesOperation.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Unknown error';
    });
  }
});

// Extract actions and reducer
export const { setCurrentLibrary, clearCurrentLibrary, setCurrentPage, setPageSize } = librariesSlice.actions;

// Selectors
export const selectLibraries = (state: { libraries: LibraryState }) => state.libraries.libraries;
export const selectCurrentLibrary = (state: { libraries: LibraryState }) => state.libraries.currentLibrary;
export const selectLibraryLoading = (state: { libraries: LibraryState }) => state.libraries.loading;
export const selectLibraryError = (state: { libraries: LibraryState }) => state.libraries.error;
export const selectTotalLibraries = (state: { libraries: LibraryState }) => state.libraries.totalLibraries;
export const selectCurrentPage = (state: { libraries: LibraryState }) => state.libraries.currentPage;
export const selectPageSize = (state: { libraries: LibraryState }) => state.libraries.pageSize;

// Export the slice and reducer
export default librariesSlice.reducer;