/**
 * Redux Toolkit slice for authentication state management in the Molecular Data Management and CRO Integration Platform.
 * This slice handles user authentication, registration, token refresh, and logout operations
 * without external dependencies.
 *
 * @module store/auth/authSlice
 */

import { createSlice, createAsyncThunk, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import { AuthState, User, LoginCredentials, RegisterData } from '../../types/auth';
import { 
  login, 
  register, 
  logout, 
  refreshToken 
} from '../../api/auth';
import { 
  setAuthToken, 
  setRefreshToken, 
  setUser, 
  clearAuth, 
  getInitialAuthState 
} from '../../utils/auth';

// Get initial state from localStorage via the utility function
const initialAuthState = getInitialAuthState();

// Initial state of the auth slice
const initialState: AuthState = {
  isAuthenticated: initialAuthState.isAuthenticated || false,
  user: initialAuthState.user || null,
  accessToken: initialAuthState.accessToken || null,
  refreshToken: initialAuthState.refreshToken || null,
  loading: false,
  error: null
};

/**
 * Async thunk for user login
 * Authenticates a user with provided credentials and stores tokens upon success
 * 
 * @param credentials - User login credentials (email/password)
 * @returns Promise resolving to user data and tokens on successful login
 */
export const loginUser = createAsyncThunk(
  'auth/login',
  async (credentials: LoginCredentials, { rejectWithValue }) => {
    try {
      const response = await login(credentials);
      
      // Store tokens and user data in localStorage
      setAuthToken(response.access_token);
      setRefreshToken(response.refresh_token);
      setUser(response.user);
      
      return { 
        user: response.user, 
        accessToken: response.access_token, 
        refreshToken: response.refresh_token 
      };
    } catch (error: any) {
      return rejectWithValue(
        error.error || 'Failed to login. Please check your credentials.'
      );
    }
  }
);

/**
 * Async thunk for user registration
 * Registers a new user with the provided data
 * 
 * @param registrationData - User registration information
 * @returns Promise resolving to registration result
 */
export const registerUser = createAsyncThunk(
  'auth/register',
  async (registrationData: RegisterData, { rejectWithValue }) => {
    try {
      const response = await register(registrationData);
      
      if (!response.success) {
        return rejectWithValue(response.message || 'Registration failed');
      }
      
      return { 
        success: response.success, 
        message: response.message 
      };
    } catch (error: any) {
      return rejectWithValue(
        error.error || 'Failed to register. Please try again.'
      );
    }
  }
);

/**
 * Async thunk for user logout
 * Logs out the user by invalidating their refresh token and clearing auth state
 * 
 * @returns Promise resolving when logout is complete
 */
export const logoutUser = createAsyncThunk(
  'auth/logout',
  async (_, { getState, rejectWithValue }) => {
    try {
      // For TypeScript to understand the state shape
      const state = getState() as { auth: AuthState };
      const refreshTokenValue = state.auth.refreshToken;
      
      if (refreshTokenValue) {
        await logout(refreshTokenValue);
      }
      
      // Clear auth data from localStorage
      clearAuth();
      
      return { success: true };
    } catch (error: any) {
      // Even if the logout API fails, we still clear the local auth state
      clearAuth();
      
      return rejectWithValue(
        error.error || 'Error during logout, but session was cleared locally.'
      );
    }
  }
);

/**
 * Async thunk for refreshing authentication token
 * Uses the refresh token to obtain a new access token
 * 
 * @returns Promise resolving to new tokens on successful refresh
 */
export const refreshUserToken = createAsyncThunk(
  'auth/refreshToken',
  async (_, { getState, rejectWithValue }) => {
    try {
      // For TypeScript to understand the state shape
      const state = getState() as { auth: AuthState };
      const currentRefreshToken = state.auth.refreshToken;
      
      if (!currentRefreshToken) {
        return rejectWithValue('No refresh token available.');
      }
      
      const response = await refreshToken(currentRefreshToken);
      
      // Store new tokens in localStorage
      setAuthToken(response.access_token);
      setRefreshToken(response.refresh_token);
      
      // Update user data if returned by the API
      if (response.user) {
        setUser(response.user);
      }
      
      return { 
        accessToken: response.access_token, 
        refreshToken: response.refresh_token,
        user: response.user
      };
    } catch (error: any) {
      // If refresh fails, clear auth state
      clearAuth();
      
      return rejectWithValue(
        error.error || 'Failed to refresh authentication token.'
      );
    }
  }
);

/**
 * Auth slice with reducers for managing authentication state
 */
const authSlice = createSlice({
  name: 'auth',
  initialState,
  reducers: {
    /**
     * Sets an error message in the auth state
     */
    setAuthError: (state, action: PayloadAction<string>) => {
      state.error = action.payload;
    },
    /**
     * Clears the error message in the auth state
     */
    clearAuthError: (state) => {
      state.error = null;
    }
  },
  extraReducers: (builder) => {
    // Login
    builder.addCase(loginUser.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(loginUser.fulfilled, (state, action) => {
      state.isAuthenticated = true;
      state.user = action.payload.user;
      state.accessToken = action.payload.accessToken;
      state.refreshToken = action.payload.refreshToken;
      state.loading = false;
      state.error = null;
    });
    builder.addCase(loginUser.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Login failed';
    });
    
    // Register
    builder.addCase(registerUser.pending, (state) => {
      state.loading = true;
      state.error = null;
    });
    builder.addCase(registerUser.fulfilled, (state) => {
      state.loading = false;
    });
    builder.addCase(registerUser.rejected, (state, action) => {
      state.loading = false;
      state.error = action.payload as string || 'Registration failed';
    });
    
    // Logout
    builder.addCase(logoutUser.pending, (state) => {
      state.loading = true;
    });
    builder.addCase(logoutUser.fulfilled, (state) => {
      state.isAuthenticated = false;
      state.user = null;
      state.accessToken = null;
      state.refreshToken = null;
      state.loading = false;
    });
    builder.addCase(logoutUser.rejected, (state, action) => {
      state.isAuthenticated = false;
      state.user = null;
      state.accessToken = null;
      state.refreshToken = null;
      state.loading = false;
      state.error = action.payload as string || 'Logout failed';
    });
    
    // Token Refresh
    builder.addCase(refreshUserToken.pending, (state) => {
      state.loading = true;
    });
    builder.addCase(refreshUserToken.fulfilled, (state, action) => {
      state.accessToken = action.payload.accessToken;
      state.refreshToken = action.payload.refreshToken;
      if (action.payload.user) {
        state.user = action.payload.user;
      }
      state.loading = false;
    });
    builder.addCase(refreshUserToken.rejected, (state, action) => {
      state.isAuthenticated = false;
      state.loading = false;
      state.error = action.payload as string || 'Token refresh failed';
    });
  }
});

// Export actions and reducer
export const { setAuthError, clearAuthError } = authSlice.actions;
export const authActions = authSlice.actions;

export default authSlice.reducer;