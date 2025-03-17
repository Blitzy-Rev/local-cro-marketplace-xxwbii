# src/web/src/App.tsx
```tsx
import React, { useState, useEffect, useMemo } from 'react'; // React v18.2+
import { useState, useEffect, useMemo } from 'react'; // react v18.2+
import { Provider } from 'react-redux'; // react-redux v8.0.5
import { RouterProvider } from 'react-router-dom'; // react-router-dom v6.11.1
import { ThemeProvider, CssBaseline } from '@mui/material'; // @mui/material v5.13.0
import { PaletteMode } from '@mui/material/styles'; // @mui/material/styles v5.13.0

// Internal imports
import router from './routes';
import store from './store';
import theme, { createCustomTheme, getStoredThemeMode, storeThemeMode } from './theme';
import ErrorBoundary from './components/common/ErrorBoundary';
import ToastContainer from './components/common/Toast';
import { ToastContext } from './contexts/ToastContext';

/**
 * Main application component that sets up providers and global components
 * @returns {JSX.Element} The rendered application
 */
const App: React.FC = () => {
  // LD1: Initialize theme mode state with useState, defaulting to stored preference or 'light'
  const [themeMode, setThemeMode] = useState<PaletteMode>(getStoredThemeMode() || 'light');

  // LD1: Create memoized custom theme based on current theme mode
  const customTheme = useMemo(() => createCustomTheme(themeMode), [themeMode]);

  // LD1: Set up effect to store theme mode preference when it changes
  useEffect(() => {
    storeThemeMode(themeMode);
  }, [themeMode]);

    // LD1: Set up toast notification state and handlers
    const [toasts, setToasts] = useState([]);

    const addToast = (message, type, duration) => {
        const id = String(Date.now());
        setToasts([...toasts, { id, type, message, duration }]);
    };

    const removeToast = (id) => {
        setToasts(toasts.filter((toast) => toast.id !== id));
    };

  // LD1: Render the component tree with all necessary providers
  return (
    <ErrorBoundary>
      {/* LD1: Provide Redux store with Provider */}
      <Provider store={store}>
        {/* LD1: Apply theme with ThemeProvider and CssBaseline */}
        <ThemeProvider theme={customTheme}>
          <CssBaseline />
            <ToastContext.Provider value={{ toasts, addToast, removeToast }}>
              {/* LD1: Render RouterProvider with the configured router */}
              <RouterProvider router={router} />

              {/* LD1: Include ToastContainer for application-wide notifications */}
              <ToastContainer toasts={toasts} onClose={removeToast} />
            </ToastContext.Provider>
        </ThemeProvider>
      </Provider>
    </ErrorBoundary>
  );
};

// IE3: Export the main App component as the default export
export default App;