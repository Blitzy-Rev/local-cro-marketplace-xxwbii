/**
 * Central export file for all TypeScript interfaces, types, and enums used in the
 * Molecular Data Management and CRO Integration Platform frontend application.
 * 
 * This barrel file consolidates type definitions to provide a single import point for components,
 * ensuring consistent type usage throughout the application and improving maintainability.
 * 
 * @packageDocumentation
 */

// Re-export all types from individual module type files
export * from './auth';
export * from './user';
export * from './molecule';
export * from './library';
export * from './experiment';
export * from './submission';
export * from './result';
export * from './notification';

/**
 * Generic interface for successful API responses.
 * Used for consistent handling of API responses across the application.
 */
export interface ApiResponse<T = any> {
  /** Whether the API request was successful */
  success: boolean;
  /** The data returned from the API */
  data: T;
  /** An optional message describing the result */
  message?: string;
}

/**
 * Interface for API error responses.
 * Used for consistent error handling across the application.
 */
export interface ApiError {
  /** Always false for error responses */
  success: false;
  /** Error message */
  error: string;
  /** HTTP status code */
  statusCode: number;
  /** Optional detailed error information */
  details?: Record<string, any>;
}

/**
 * Generic interface for paginated API responses.
 * Used when retrieving lists of items that support pagination.
 */
export interface PaginatedResponse<T = any> {
  /** Array of items for the current page */
  data: T[];
  /** Total number of items across all pages */
  total: number;
  /** Current page number (1-based) */
  page: number;
  /** Number of items per page */
  pageSize: number;
  /** Total number of pages */
  totalPages: number;
}