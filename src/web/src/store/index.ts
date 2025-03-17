/**
 * Redux store configuration file for the Molecular Data Management and CRO Integration Platform frontend.
 * This file creates and configures the Redux store with middleware, devtools, and typed hooks for use throughout the application.
 */

import { configureStore } from '@reduxjs/toolkit'; // @reduxjs/toolkit ^1.9+
import { TypedUseSelectorHook, useDispatch, useSelector } from 'react-redux'; // react-redux ^8.0+
import rootReducer, { RootState } from './rootReducer';

/**
 * Configure the Redux store with:
 * - The root reducer combining all slice reducers
 * - Default middleware with serializability check configuration
 * - Redux DevTools enabled in non-production environments for debugging
 */
const store = configureStore({
  reducer: rootReducer,
  middleware: (getDefaultMiddleware) => 
    getDefaultMiddleware({
      serializableCheck: {
        // Ignore non-serializable data in specific actions if needed
        ignoredActions: ['...']
      }
    }),
  devTools: process.env.NODE_ENV !== 'production'
});

/**
 * Type for the store's dispatch function
 * Provides type safety when dispatching actions
 */
export type AppDispatch = typeof store.dispatch;

/**
 * Custom typed version of useDispatch hook
 * Use this throughout the application instead of plain useDispatch
 * @returns Typed dispatch function
 */
export const useAppDispatch = () => useDispatch<AppDispatch>();

/**
 * Custom typed version of useSelector hook
 * Use this throughout the application instead of plain useSelector
 * Provides type safety when selecting from the Redux store
 */
export const useAppSelector: TypedUseSelectorHook<RootState> = useSelector;

export default store;