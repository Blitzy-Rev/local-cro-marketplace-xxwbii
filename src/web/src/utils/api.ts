/**
 * Utility module that provides higher-level API request functions and helpers for the
 * Molecular Data Management and CRO Integration Platform frontend.
 * 
 * This module builds upon the core API client to offer simplified request methods,
 * response handling, error processing, and caching strategies.
 * 
 * @module utils/api
 */

import { AxiosRequestConfig, AxiosResponse } from 'axios'; // axios ^1.4.0
import { apiClient, handleApiError } from '../api/client';
import { ApiResponse, ApiError, PaginatedResponse } from '../types';
import { getLocalStorageItem, setLocalStorageItem } from './storage';

// Constants for caching
const API_CACHE_PREFIX = "api_cache_";
const DEFAULT_CACHE_DURATION = 5 * 60 * 1000; // 5 minutes in milliseconds

/**
 * Generic function to fetch data from the API with type safety
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the requested data
 */
export const fetchData = async <T = any>(
  url: string,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    const response = await apiClient.get<ApiResponse<T>>(url, config);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetch paginated data from the API with support for page, limit, and filters
 * 
 * @template T Type of the items in the paginated response
 * @param url The API endpoint URL
 * @param params Object containing pagination parameters and filters
 * @param config Optional Axios request configuration
 * @returns Promise resolving to paginated data
 */
export const fetchPaginatedData = async <T = any>(
  url: string,
  params: {
    page?: number;
    pageSize?: number;
    [key: string]: any;
  } = {},
  config?: AxiosRequestConfig
): Promise<PaginatedResponse<T>> => {
  try {
    const queryParams = { ...params };
    
    // Ensure page and pageSize are included
    if (!queryParams.page) {
      queryParams.page = 1;
    }
    if (!queryParams.pageSize) {
      queryParams.pageSize = 20;
    }
    
    const response = await apiClient.get<ApiResponse<PaginatedResponse<T>>>(
      url,
      { ...config, params: queryParams }
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Send data to the API using a POST request
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param data The data to send in the request body
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the response data
 */
export const postData = async <T = any>(
  url: string,
  data: any,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    const response = await apiClient.post<ApiResponse<T>>(url, data, config);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Update data on the API using a PUT request
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param data The data to send in the request body
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the response data
 */
export const putData = async <T = any>(
  url: string,
  data: any,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    const response = await apiClient.put<ApiResponse<T>>(url, data, config);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Partially update data on the API using a PATCH request
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param data The data to send in the request body
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the response data
 */
export const patchData = async <T = any>(
  url: string,
  data: any,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    const response = await apiClient.patch<ApiResponse<T>>(url, data, config);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Delete data on the API using a DELETE request
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the response data
 */
export const deleteData = async <T = any>(
  url: string,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    const response = await apiClient.delete<ApiResponse<T>>(url, config);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Upload a file to the API using FormData
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param file The file to upload
 * @param fieldName The name of the field to use for the file (default: 'file')
 * @param additionalData Additional data to include in the FormData
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the response data
 */
export const uploadFile = async <T = any>(
  url: string,
  file: File,
  fieldName: string = 'file',
  additionalData: Record<string, any> = {},
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    // Create FormData instance
    const formData = new FormData();
    
    // Append the file
    formData.append(fieldName, file);
    
    // Append any additional data
    Object.entries(additionalData).forEach(([key, value]) => {
      formData.append(key, typeof value === 'object' ? JSON.stringify(value) : String(value));
    });
    
    // Set content type to multipart/form-data
    const uploadConfig: AxiosRequestConfig = {
      ...config,
      headers: {
        ...config?.headers,
        'Content-Type': 'multipart/form-data',
      },
    };
    
    const response = await apiClient.post<ApiResponse<T>>(url, formData, uploadConfig);
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Download a file from the API with proper handling of binary data
 * 
 * @param url The API endpoint URL for the file
 * @param filename The name to save the file as
 * @param config Optional Axios request configuration
 * @returns Promise resolving to true if download was successful
 */
export const downloadFile = async (
  url: string,
  filename: string,
  config?: AxiosRequestConfig
): Promise<boolean> => {
  try {
    // Configure request to receive binary data
    const downloadConfig: AxiosRequestConfig = {
      ...config,
      responseType: 'blob',
    };
    
    // Make the request
    const response = await apiClient.get<Blob>(url, downloadConfig);
    
    // Create a download link and trigger download
    const blob = new Blob([response.data], { 
      type: response.headers['content-type'] || 'application/octet-stream' 
    });
    const downloadUrl = window.URL.createObjectURL(blob);
    const link = document.createElement('a');
    
    link.href = downloadUrl;
    link.download = filename;
    document.body.appendChild(link);
    link.click();
    
    // Clean up
    window.URL.revokeObjectURL(downloadUrl);
    document.body.removeChild(link);
    
    return true;
  } catch (error) {
    handleApiError(error);
    return false;
  }
};

/**
 * Interface for cached data structure
 */
interface CachedData<T> {
  data: T;
  timestamp: number;
  expiration: number;
}

/**
 * Fetch data from the API with local caching support
 * 
 * @template T Type of the expected response data
 * @param url The API endpoint URL
 * @param cacheDuration Duration in milliseconds to keep the cache valid (default: 5 minutes)
 * @param config Optional Axios request configuration
 * @returns Promise resolving to the data, either from cache or fresh from API
 */
export const getCachedData = async <T = any>(
  url: string,
  cacheDuration: number = DEFAULT_CACHE_DURATION,
  config?: AxiosRequestConfig
): Promise<T> => {
  try {
    // Generate cache key based on URL and query parameters
    const cacheKey = generateCacheKey(url, config);
    
    // Check if there's valid cached data
    const cachedItem = getLocalStorageItem<CachedData<T>>(cacheKey);
    const now = Date.now();
    
    if (cachedItem && cachedItem.expiration > now) {
      // Return cached data if still valid
      return cachedItem.data;
    }
    
    // If no valid cache, fetch fresh data
    const response = await apiClient.get<ApiResponse<T>>(url, config);
    const data = response.data.data;
    
    // Store the fresh data in cache
    const cacheData: CachedData<T> = {
      data,
      timestamp: now,
      expiration: now + cacheDuration,
    };
    
    setLocalStorageItem(cacheKey, cacheData);
    
    return data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Invalidate cached data for a specific URL or pattern
 * 
 * @param urlPattern String or RegExp to match cache keys against
 */
export const invalidateCache = (urlPattern: string | RegExp): void => {
  try {
    // Get all localStorage keys
    const keys: string[] = [];
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && key.startsWith(API_CACHE_PREFIX)) {
        keys.push(key);
      }
    }
    
    // Remove matching cache entries
    keys.forEach(key => {
      // Remove the prefix to get the original URL
      const originalKey = key.substring(API_CACHE_PREFIX.length);
      
      // Check if the key matches the pattern
      if (typeof urlPattern === 'string') {
        if (originalKey.includes(urlPattern)) {
          localStorage.removeItem(key);
        }
      } else if (urlPattern instanceof RegExp) {
        if (urlPattern.test(originalKey)) {
          localStorage.removeItem(key);
        }
      }
    });
  } catch (error) {
    console.error('Error invalidating cache:', error);
  }
};

/**
 * Build a URL query string from an object of parameters
 * 
 * @param params Object containing query parameters
 * @returns Formatted query string
 */
export const buildQueryString = (params: Record<string, any>): string => {
  if (!params || Object.keys(params).length === 0) {
    return '';
  }
  
  // Filter out null and undefined values
  const filteredParams: Record<string, any> = {};
  Object.entries(params).forEach(([key, value]) => {
    if (value !== null && value !== undefined) {
      filteredParams[key] = value;
    }
  });
  
  // Convert to URLSearchParams for proper encoding
  const searchParams = new URLSearchParams();
  
  Object.entries(filteredParams).forEach(([key, value]) => {
    // Handle arrays
    if (Array.isArray(value)) {
      value.forEach(item => searchParams.append(`${key}[]`, String(item)));
    }
    // Handle objects
    else if (typeof value === 'object') {
      searchParams.append(key, JSON.stringify(value));
    }
    // Handle primitives
    else {
      searchParams.append(key, String(value));
    }
  });
  
  return searchParams.toString();
};

/**
 * Generate a unique cache key for an API request
 * 
 * @param url The API endpoint URL
 * @param config Optional Axios request configuration containing query parameters
 * @returns Unique cache key
 */
export const generateCacheKey = (
  url: string,
  config?: AxiosRequestConfig
): string => {
  let key = url;
  
  // Include query parameters in the key if present
  if (config?.params) {
    const queryString = buildQueryString(config.params);
    key = `${url}?${queryString}`;
  }
  
  // Hash the key to create a fixed-length identifier
  // Simple string hash function
  const hash = Array.from(key)
    .reduce((hash, char) => ((hash << 5) - hash) + char.charCodeAt(0), 0)
    .toString(36);
  
  return `${API_CACHE_PREFIX}${hash}`;
};