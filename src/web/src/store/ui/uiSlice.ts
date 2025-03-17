/**
 * Redux Toolkit slice for UI state management in the Molecular Data Management and CRO Integration Platform
 * 
 * This slice handles global UI state including:
 * - Sidebar visibility
 * - Theme mode (light/dark)
 * - Active view/page
 * - Toast notifications
 * - Modal dialogs
 */

import { createSlice, PayloadAction } from '@reduxjs/toolkit'; // @reduxjs/toolkit version 1.9+

/**
 * Interface for toast notification
 */
export interface Toast {
  id: string;
  type: 'success' | 'error' | 'warning' | 'info';
  message: string;
  duration: number;
  createdAt: number;
}

/**
 * Interface for the UI state
 */
export interface UIState {
  sidebarOpen: boolean;
  themeMode: 'light' | 'dark';
  activeView: string;
  toasts: Toast[];
  modals: Record<string, { open: boolean; props?: any }>;
}

/**
 * Initial state for the UI slice
 */
const initialState: UIState = {
  sidebarOpen: true,
  themeMode: 'light',
  activeView: 'dashboard',
  toasts: [],
  modals: {},
};

/**
 * UI slice definition with reducers
 */
const uiSlice = createSlice({
  name: 'ui',
  initialState,
  reducers: {
    /**
     * Sets the sidebar open/closed state
     */
    setSidebarOpen: (state, action: PayloadAction<boolean>) => {
      state.sidebarOpen = action.payload;
    },
    
    /**
     * Sets the theme mode (light/dark)
     */
    setThemeMode: (state, action: PayloadAction<'light' | 'dark'>) => {
      state.themeMode = action.payload;
    },
    
    /**
     * Sets the active view/page
     */
    setActiveView: (state, action: PayloadAction<string>) => {
      state.activeView = action.payload;
    },
    
    /**
     * Adds a toast notification to the state
     */
    addToast: (state, action: PayloadAction<Omit<Toast, 'createdAt'>>) => {
      const toast = {
        ...action.payload,
        createdAt: Date.now()
      };
      state.toasts.push(toast);
    },
    
    /**
     * Removes a toast notification from the state
     */
    removeToast: (state, action: PayloadAction<string>) => {
      state.toasts = state.toasts.filter(toast => toast.id !== action.payload);
    },
    
    /**
     * Clears all toast notifications
     */
    clearToasts: (state) => {
      state.toasts = [];
    },
    
    /**
     * Opens a modal dialog
     */
    openModal: (state, action: PayloadAction<{ id: string; props?: any }>) => {
      const { id, props } = action.payload;
      state.modals[id] = { open: true, props };
    },
    
    /**
     * Closes a modal dialog
     */
    closeModal: (state, action: PayloadAction<string>) => {
      const id = action.payload;
      if (state.modals[id]) {
        state.modals[id].open = false;
      }
    },
  },
});

// Extract the action creators and reducer
export const { 
  setSidebarOpen,
  setThemeMode,
  setActiveView,
  addToast,
  removeToast,
  clearToasts,
  openModal,
  closeModal
} = uiSlice.actions;

// Export all actions as a grouped object for convenience
export const uiActions = uiSlice.actions;

// Export the reducer as the default export
export default uiSlice.reducer;

// Selectors
export const selectSidebarOpen = (state: { ui: UIState }) => state.ui.sidebarOpen;
export const selectThemeMode = (state: { ui: UIState }) => state.ui.themeMode;
export const selectActiveView = (state: { ui: UIState }) => state.ui.activeView;
export const selectToasts = (state: { ui: UIState }) => state.ui.toasts;
export const selectModals = (state: { ui: UIState }) => state.ui.modals;