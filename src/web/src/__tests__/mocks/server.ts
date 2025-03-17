import { setupServer } from 'msw/node'; // msw ^1.2.1
import handlers from './handlers';

/**
 * Creates a Mock Service Worker (MSW) server instance that intercepts network requests 
 * during testing and returns mock responses based on the defined handlers.
 * 
 * This server allows frontend tests to run without actual API dependencies, ensuring
 * consistent and predictable test data.
 * 
 * @see https://mswjs.io/docs/api/setup-server
 */
export const server = setupServer(...handlers);

// The server exposes the following methods:
// - server.listen() - Starts the server and begins intercepting requests
// - server.close() - Stops the server and cleans up any request interception
// - server.resetHandlers() - Resets any runtime request handlers to the initial set
// - server.use() - Registers runtime request handlers for specific test cases