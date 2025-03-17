/**
 * Jest testing environment setup file
 * 
 * This file configures the global test settings for Jest, sets up the Mock Service Worker (MSW)
 * for API mocking, and extends Jest with additional matchers for improved testing capabilities.
 */

// Import @testing-library/jest-dom to extend Jest with custom DOM matchers
import '@testing-library/jest-dom'; // v5.16.5

// Import jest-fetch-mock for fetch mocking capabilities
import fetchMock from 'jest-fetch-mock'; // v3.0.3
fetchMock.enableMocks();

// Import MSW server for API mocking
import { server } from './mocks/server';

// Configure MSW to listen for API requests before all tests
// Using onUnhandledRequest: 'error' to identify API calls that aren't properly mocked
global.beforeAll(() => server.listen({ onUnhandledRequest: 'error' }));

// Reset handlers after each test for a clean state
global.afterEach(() => server.resetHandlers());

// Close the server after all tests are done
global.afterAll(() => server.close());

// Mock window.matchMedia for responsive design testing
// This is needed for testing components that use media queries
Object.defineProperty(window, 'matchMedia', {
  writable: true,
  value: jest.fn().mockImplementation(query => ({
    matches: false,
    media: query,
    onchange: null,
    addListener: jest.fn(),
    removeListener: jest.fn(),
    addEventListener: jest.fn(),
    removeEventListener: jest.fn(),
    dispatchEvent: jest.fn(),
  })),
});