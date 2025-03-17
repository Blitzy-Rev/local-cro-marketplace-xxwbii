/**
 * A custom React hook that provides functionality for managing toast notifications throughout the application.
 * Leverages Redux for state management and provides a simple API for creating, displaying, and dismissing toast notifications.
 */

import { useCallback } from 'react'; // react 18.2+
import { v4 as uuidv4 } from 'uuid'; // uuid 9.0+
import { useAppDispatch, useAppSelector } from '../store';
import { addToast, removeToast, clearToasts } from '../store/ui/uiSlice';

/**
 * Interface defining the options for creating a toast notification
 */
interface ToastOptions {
  /** The type of toast notification */
  type: 'success' | 'error' | 'warning' | 'info';
  /** The content of the toast notification */
  message: string;
  /** How long the toast should be displayed in milliseconds (default: 5000) */
  duration?: number;
}

/**
 * A hook that provides functionality for managing toast notifications
 * @returns Toast management functions and state: { toasts, showToast, hideToast, clearAllToasts }
 */
export const useToast = () => {
  // Get the Redux dispatch function
  const dispatch = useAppDispatch();
  
  // Select the toasts array from the Redux store
  const toasts = useAppSelector((state) => state.ui.toasts);

  /**
   * Displays a toast notification
   * @param options Configuration options for the toast
   * @returns The unique ID of the created toast
   */
  const showToast = useCallback((options: ToastOptions) => {
    const { type, message, duration = 5000 } = options;
    const id = uuidv4(); // Generate unique ID for the toast
    
    dispatch(
      addToast({
        id,
        type,
        message,
        duration,
      })
    );
    
    return id; // Return ID for potential manual dismissal
  }, [dispatch]);

  /**
   * Hides a specific toast notification
   * @param id The unique ID of the toast to hide
   */
  const hideToast = useCallback((id: string) => {
    dispatch(removeToast(id));
  }, [dispatch]);

  /**
   * Clears all toast notifications
   */
  const clearAllToasts = useCallback(() => {
    dispatch(clearToasts());
  }, [dispatch]);

  // Return the toasts array and management functions
  return {
    toasts,
    showToast,
    hideToast,
    clearAllToasts,
  };
};