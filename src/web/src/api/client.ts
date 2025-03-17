/**
 * Core API client module for the Molecular Data Management and CRO Integration Platform.
 * Configures and exports an Axios instance for making HTTP requests to the backend API.
 * Handles authentication token management, request/response interceptors, and error handling.
 *
 * @module api/client
 */

import axios, { 
  AxiosInstance, 
  AxiosRequestConfig, 
  AxiosResponse, 
  AxiosError,
  InternalAxiosRequestConfig
} from 'axios'; // axios ^1.4.0
import { 
  ApiResponse, 
  ApiError 
} from '../types';
import { 
  getAuthToken, 
  isTokenExpired, 
  getRefreshToken, 
  setAuthToken, 
  setRefreshToken, 
  clearAuth 
} from '../utils/auth';

/**
 * Base URL for API requests.
 * Uses environment variable if available, otherwise defaults to '/api/v1'.
 */
const API_BASE_URL = process.env.REACT_APP_API_BASE_URL || '/api/v1';

/**
 * Endpoint for refreshing authentication tokens.
 */
const REFRESH_ENDPOINT = '/auth/refresh-token';

/**
 * Maximum number of retry attempts for failed requests.
 */
const MAX_RETRY_ATTEMPTS = 3;

/**
 * Flag indicating whether a token refresh operation is in progress.
 */
let isRefreshing = false;

/**
 * Array of callbacks to be executed when a new token is obtained.
 * Used to resolve pending requests that were waiting for token refresh.
 */
let refreshSubscribers: ((token: string) => void)[] = [];

/**
 * Creates and configures an Axios instance for API requests.
 * Sets base URL, default headers, and timeout.
 * 
 * @returns Configured Axios instance
 */
const createApiClient = (): AxiosInstance => {
  const client = axios.create({
    baseURL: API_BASE_URL,
    headers: {
      'Content-Type': 'application/json',
    },
    timeout: 30000, // 30 seconds
  });

  return client;
};

/**
 * Sets up request and response interceptors for the API client.
 * - Request interceptor adds authentication token to requests
 * - Response interceptor handles token refresh and request retry
 * 
 * @param client The Axios instance to configure
 */
const setupInterceptors = (client: AxiosInstance): void => {
  // Request interceptor - adds authentication token to requests
  client.interceptors.request.use(
    async (config: InternalAxiosRequestConfig) => {
      const token = getAuthToken();
      
      if (token && !isTokenExpired(token)) {
        config.headers = config.headers || {};
        config.headers['Authorization'] = `Bearer ${token}`;
      }
      
      return config;
    },
    (error) => Promise.reject(error)
  );
  
  // Response interceptor - handles errors and token refresh
  client.interceptors.response.use(
    // Pass through successful responses
    (response: AxiosResponse) => {
      return response;
    },
    // Handle error responses
    async (error: AxiosError) => {
      // Cannot retry if there's no config
      if (!error.config) {
        return handleApiError(error);
      }
      
      const originalRequest = { ...error.config };
      const originalRequestWithRetry = originalRequest as { _retry?: boolean };
      
      // Handle 401 Unauthorized errors (token expired)
      if (
        error.response?.status === 401 &&
        !originalRequestWithRetry._retry
      ) {
        originalRequestWithRetry._retry = true;
        
        if (isRefreshing) {
          // If a token refresh is already in progress, wait for the new token
          try {
            const newToken = await new Promise<string>((resolve, reject) => {
              subscribeToTokenRefresh((token: string) => {
                resolve(token);
              });
            });
            
            // Create a new request with the updated token
            const newRequest = { ...originalRequest };
            if (!newRequest.headers) {
              newRequest.headers = {};
            }
            newRequest.headers['Authorization'] = `Bearer ${newToken}`;
            
            return client(newRequest);
          } catch (refreshError) {
            // If the waiting for refresh fails, redirect to login
            clearAuth();
            return Promise.reject(refreshError);
          }
        }
        
        // Start the token refresh process
        isRefreshing = true;
        
        try {
          // Attempt to refresh the token
          const newToken = await refreshAuthToken();
          
          // Create a new request with the updated token
          const newRequest = { ...originalRequest };
          if (!newRequest.headers) {
            newRequest.headers = {};
          }
          newRequest.headers['Authorization'] = `Bearer ${newToken}`;
          
          // Notify all waiting requests that a new token is available
          onRefreshed(newToken);
          
          // Reset refreshing state
          isRefreshing = false;
          
          // Retry the original request with the new token
          return client(newRequest);
        } catch (refreshError) {
          // If token refresh fails, clear auth and reject
          isRefreshing = false;
          clearAuth();
          return Promise.reject(refreshError);
        }
      }
      
      // For server errors, attempt retry with exponential backoff
      if (error.response?.status && error.response.status >= 500) {
        return retryRequest(error, 0);
      }
      
      // For all other errors, standardize the error format
      return handleApiError(error);
    }
  );
};

/**
 * Refreshes the authentication token using the refresh token.
 * Makes a request to the refresh endpoint with the current refresh token.
 * 
 * @returns Promise resolving to the new authentication token
 * @throws Error if refresh token is not available or refresh fails
 */
const refreshAuthToken = async (): Promise<string> => {
  const refreshToken = getRefreshToken();
  
  if (!refreshToken) {
    return Promise.reject(new Error('No refresh token available'));
  }
  
  try {
    const response = await axios.post<ApiResponse<{ access_token: string; refresh_token: string }>>(
      `${API_BASE_URL}${REFRESH_ENDPOINT}`,
      { refresh_token: refreshToken }
    );
    
    const { access_token, refresh_token } = response.data.data;
    
    // Store the new tokens
    setAuthToken(access_token);
    setRefreshToken(refresh_token);
    
    return access_token;
  } catch (error) {
    // Clear auth on refresh failure
    clearAuth();
    return Promise.reject(new Error('Failed to refresh authentication token'));
  }
};

/**
 * Notifies subscribers when a new token is obtained.
 * Calls all waiting callbacks with the new token and clears the subscribers list.
 * 
 * @param token The new authentication token
 */
const onRefreshed = (token: string): void => {
  // Call all waiting subscribers with the new token
  refreshSubscribers.forEach(callback => callback(token));
  // Clear subscribers after notifying them
  refreshSubscribers = [];
};

/**
 * Adds a callback to be executed when a new token is obtained.
 * Used by requests that are waiting for token refresh.
 * 
 * @param callback Function to call with the new token
 */
const subscribeToTokenRefresh = (callback: (token: string) => void): void => {
  refreshSubscribers.push(callback);
};

/**
 * Transforms Axios errors into a consistent error format.
 * Extracts status code, error message, and details from the response.
 * 
 * @param error Axios error object
 * @returns Promise that always rejects with formatted error
 */
const handleApiError = (error: AxiosError): Promise<never> => {
  // Extract status code or default to 500
  const statusCode = error.response?.status || 500;
  
  // Extract error message from response or use default
  let errorMessage = 'An unexpected error occurred';
  let errorDetails = undefined;
  
  if (error.response?.data) {
    const responseData = error.response.data as any;
    
    if (responseData.error) {
      errorMessage = responseData.error;
    } else if (responseData.message) {
      errorMessage = responseData.message;
    }
    
    if (responseData.details) {
      errorDetails = responseData.details;
    }
  }
  
  // Create a standardized error response
  const apiError: ApiError = {
    success: false,
    error: errorMessage,
    statusCode,
    details: errorDetails
  };
  
  return Promise.reject(apiError);
};

/**
 * Retries a failed request with exponential backoff.
 * Waits longer between each retry attempt.
 * 
 * @param error The original error
 * @param retryCount Current retry attempt count
 * @returns Promise resolving to the response if retry succeeds
 */
const retryRequest = async (error: AxiosError, retryCount: number): Promise<AxiosResponse> => {
  // Check if we've reached the maximum retry attempts
  if (retryCount >= MAX_RETRY_ATTEMPTS) {
    return handleApiError(error);
  }
  
  // Cannot retry if there's no config
  if (!error.config) {
    return handleApiError(error);
  }
  
  // Don't retry for certain error types
  const status = error.response?.status;
  if (status && status >= 400 && status < 500 && status !== 401) {
    return handleApiError(error);
  }
  
  // Calculate delay using exponential backoff
  const delay = Math.pow(2, retryCount) * 1000; // 1s, 2s, 4s, ...
  
  // Wait for the calculated delay
  await new Promise(resolve => setTimeout(resolve, delay));
  
  try {
    // Retry the request with the original config
    const response = await axios(error.config);
    return response;
  } catch (retryError) {
    // If retry fails, recursively retry with incremented count
    return retryRequest(retryError as AxiosError, retryCount + 1);
  }
};

// Create and configure the API client
const apiClient = createApiClient();
setupInterceptors(apiClient);

// Export the configured client and error handler
export { apiClient, handleApiError };