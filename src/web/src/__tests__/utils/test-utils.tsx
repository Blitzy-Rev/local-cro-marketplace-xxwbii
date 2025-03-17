import React from 'react'; // react ^18.2.0
import { render, screen, waitFor, fireEvent, within } from '@testing-library/react'; // @testing-library/react ^14.0.0
import userEvent from '@testing-library/user-event'; // @testing-library/user-event ^14.4.3
import { Provider } from 'react-redux'; // react-redux ^8.0.5
import { RouterProvider, createMemoryRouter } from 'react-router-dom'; // react-router-dom ^6.10.0
import { ThemeProvider } from '@mui/material'; // @mui/material ^5.13.0
import { theme } from '@mui/material/styles'; // @mui/material/styles ^5.13.0
import { rest } from 'msw'; // msw ^1.2.1
import { store } from '../../store';
import { RootState } from '../../store/rootReducer';
import router from '../../routes';
import { server } from '../mocks/server';

/**
 * Custom render function that wraps components with necessary providers
 * @param {React.ReactElement} ui - The UI element to render
 * @param {object} options - Options for the render function
 * @returns {object} Rendered component with additional utilities
 */
const customRender = (ui: React.ReactElement, options = {}) => {
  // LD1: Create a wrapper component that provides Redux store, Router, and Theme providers
  const Wrapper = ({ children }: { children: React.ReactNode }) => {
    return (
      <Provider store={store}>
        <ThemeProvider theme={theme}>
          {children}
        </ThemeProvider>
      </Provider>
    );
  };

  // LD1: Render the component using RTL's render function with the wrapper
  return render(ui, { wrapper: Wrapper, ...options });
};

/**
 * Render a component with a memory router for testing routing
 * @param {React.ReactElement} ui - The UI element to render
 * @param {object} options - Options for the render function
 * @returns {object} Rendered component with router utilities
 */
const renderWithRouter = (ui: React.ReactElement, options = {}) => {
  // LD1: Create a memory router with the provided routes and initial entries
  const testRouter = createMemoryRouter([{
    path: '/',
    element: ui,
  }], {
    initialEntries: ['/'],
  });

  // LD1: Render the component with the memory router using customRender
  return customRender(<RouterProvider router={testRouter} />, options);
};

/**
 * Set up custom API mock handlers for specific test cases
 * @param {array} handlers - Array of MSW request handlers
 * @returns {void} No return value
 */
const setupApiMock = (handlers: any[]) => {
  // LD1: Use the MSW server to register the provided handlers
  server.use(...handlers);
};

/**
 * Create a mock Redux store with initial state for testing
 * @param {object} initialState - Initial state for the Redux store
 * @returns {object} Configured Redux store with initial state
 */
const createMockStore = (initialState: Partial<RootState>) => {
  // LD1: Create a new Redux store with the provided initial state
  const mockStore = {
    getState: () => ({
      ...initialState,
    }),
    dispatch: jest.fn(),
    subscribe: jest.fn(),
  };

  // LD1: Return the configured store for use in tests
  return mockStore;
};

// Re-export everything
export * from '@testing-library/react';
export { customRender as render, renderWithRouter, setupApiMock, createMockStore, screen, waitFor, fireEvent, userEvent, within, rest };