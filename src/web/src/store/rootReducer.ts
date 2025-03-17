/**
 * Root reducer file that combines all Redux slice reducers for the Molecular Data Management and CRO Integration Platform frontend.
 * This file creates the root state type and exports the combined reducer for use in the Redux store configuration.
 */

import { combineReducers } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+

// Import all slice reducers
import authReducer from './auth/authSlice';
import moleculesReducer from './molecules/moleculesSlice';
import librariesReducer from './libraries/librariesSlice';
import experimentsReducer from './experiments/experimentsSlice';
import submissionsReducer from './submissions/submissionsSlice';
import resultsReducer from './results/resultsSlice';
import uiReducer from './ui/uiSlice';
import notificationsReducer from './notifications/notificationsSlice';

/**
 * Root reducer that combines all slice reducers
 * Each slice reducer manages a specific part of the application state
 */
const rootReducer = combineReducers({
  auth: authReducer,             // Authentication state reducer
  molecules: moleculesReducer,   // Molecular data state reducer
  libraries: librariesReducer,   // Library management state reducer
  experiments: experimentsReducer, // Experiment management state reducer
  submissions: submissionsReducer, // CRO submission state reducer
  results: resultsReducer,       // Experimental results state reducer
  ui: uiReducer,                 // UI state reducer
  notifications: notificationsReducer // Notification state reducer
});

/**
 * RootState type derived from the rootReducer
 * This provides type safety for state selectors and components accessing the Redux store
 */
export type RootState = ReturnType<typeof rootReducer>;

export default rootReducer;