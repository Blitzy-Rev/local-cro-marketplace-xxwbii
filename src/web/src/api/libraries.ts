/**
 * API client module for library-related operations in the Molecular Data Management and CRO Integration Platform.
 * Provides functions for fetching, creating, updating, and deleting molecule libraries, as well as
 * managing molecules within libraries through add/remove operations.
 *
 * @module api/libraries
 */

import { AxiosResponse } from 'axios'; // axios ^1.4.0
import { apiClient, handleApiError } from './client';
import { 
  ApiResponse, 
  Library, 
  LibraryDetail, 
  LibraryCreate, 
  LibraryUpdate, 
  LibraryFilter, 
  LibraryListResponse, 
  LibraryMoleculeOperation,
  LibraryOperationType,
  LibraryExportFormat
} from '../types';

/**
 * API endpoints for library operations
 */
const API_ENDPOINTS = {
  LIBRARIES: '/libraries',
  MY_LIBRARIES: '/libraries/me',
  LIBRARY_MOLECULES: '/libraries/{id}/molecules',
  LIBRARY_EXPORT: '/libraries/{id}/export'
};

/**
 * Fetches a list of libraries with optional filtering, sorting, and pagination
 * 
 * @param filter - Optional filter criteria for libraries
 * @param page - Page number for pagination (default: 1)
 * @param pageSize - Number of items per page (default: 10)
 * @returns Promise resolving to paginated list of libraries
 */
export const getLibraries = async (
  filter?: LibraryFilter,
  page: number = 1,
  pageSize: number = 10
): Promise<AxiosResponse<ApiResponse<LibraryListResponse>>> => {
  try {
    // Construct query parameters
    const params = {
      ...filter,
      page,
      page_size: pageSize
    };
    
    return await apiClient.get(API_ENDPOINTS.LIBRARIES, { params });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches libraries created by the current user
 * 
 * @param page - Page number for pagination (default: 1)
 * @param pageSize - Number of items per page (default: 10)
 * @returns Promise resolving to paginated list of user's libraries
 */
export const getMyLibraries = async (
  page: number = 1,
  pageSize: number = 10
): Promise<AxiosResponse<ApiResponse<LibraryListResponse>>> => {
  try {
    const params = {
      page,
      page_size: pageSize
    };
    
    return await apiClient.get(API_ENDPOINTS.MY_LIBRARIES, { params });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches a single library by ID
 * 
 * @param id - Library ID
 * @returns Promise resolving to library data
 */
export const getLibraryById = async (
  id: string
): Promise<AxiosResponse<ApiResponse<Library>>> => {
  try {
    return await apiClient.get(`${API_ENDPOINTS.LIBRARIES}/${id}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches detailed library information including contained molecules
 * 
 * @param id - Library ID
 * @returns Promise resolving to detailed library data
 */
export const getLibraryDetailById = async (
  id: string
): Promise<AxiosResponse<ApiResponse<LibraryDetail>>> => {
  try {
    return await apiClient.get(`${API_ENDPOINTS.LIBRARIES}/${id}`, {
      params: { detail: true }
    });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Creates a new library
 * 
 * @param libraryData - Data for creating the new library
 * @returns Promise resolving to created library data
 */
export const createLibrary = async (
  libraryData: LibraryCreate
): Promise<AxiosResponse<ApiResponse<Library>>> => {
  try {
    return await apiClient.post(API_ENDPOINTS.LIBRARIES, libraryData);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates an existing library
 * 
 * @param id - Library ID
 * @param libraryData - Updated library data
 * @returns Promise resolving to updated library data
 */
export const updateLibrary = async (
  id: string,
  libraryData: LibraryUpdate
): Promise<AxiosResponse<ApiResponse<Library>>> => {
  try {
    return await apiClient.put(`${API_ENDPOINTS.LIBRARIES}/${id}`, libraryData);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes a library
 * 
 * @param id - Library ID
 * @returns Promise resolving to success response
 */
export const deleteLibrary = async (
  id: string
): Promise<AxiosResponse<ApiResponse<{ success: boolean }>>> => {
  try {
    return await apiClient.delete(`${API_ENDPOINTS.LIBRARIES}/${id}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Adds molecules to a library
 * 
 * @param libraryId - Library ID
 * @param moleculeIds - Array of molecule IDs to add
 * @returns Promise resolving to operation result
 */
export const addMoleculesToLibrary = async (
  libraryId: string,
  moleculeIds: string[]
): Promise<AxiosResponse<ApiResponse<{ success: boolean, added_count: number }>>> => {
  try {
    const operation: LibraryMoleculeOperation = {
      operation: LibraryOperationType.ADD,
      molecule_ids: moleculeIds
    };
    
    return await updateLibraryMolecules(libraryId, operation) as Promise<AxiosResponse<ApiResponse<{ success: boolean, added_count: number }>>>;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Removes molecules from a library
 * 
 * @param libraryId - Library ID
 * @param moleculeIds - Array of molecule IDs to remove
 * @returns Promise resolving to operation result
 */
export const removeMoleculesFromLibrary = async (
  libraryId: string,
  moleculeIds: string[]
): Promise<AxiosResponse<ApiResponse<{ success: boolean, removed_count: number }>>> => {
  try {
    const operation: LibraryMoleculeOperation = {
      operation: LibraryOperationType.REMOVE,
      molecule_ids: moleculeIds
    };
    
    return await updateLibraryMolecules(libraryId, operation) as Promise<AxiosResponse<ApiResponse<{ success: boolean, removed_count: number }>>>;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Performs add or remove operations on library molecules
 * 
 * @param libraryId - Library ID
 * @param operation - Operation details (add or remove molecules)
 * @returns Promise resolving to operation result
 */
export const updateLibraryMolecules = async (
  libraryId: string,
  operation: LibraryMoleculeOperation
): Promise<AxiosResponse<ApiResponse<{ success: boolean, count: number }>>> => {
  try {
    const url = API_ENDPOINTS.LIBRARY_MOLECULES.replace('{id}', libraryId);
    return await apiClient.post(url, operation);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Exports a library in the specified format
 * 
 * @param libraryId - Library ID
 * @param format - Export format (csv, sdf)
 * @returns Promise resolving to file blob
 */
export const exportLibrary = async (
  libraryId: string,
  format: LibraryExportFormat
): Promise<Blob> => {
  try {
    const url = getLibraryExportUrl(libraryId, format);
    const response = await apiClient.get(url, { responseType: 'blob' });
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Generates a URL for library export
 * 
 * @param libraryId - Library ID
 * @param format - Export format (csv, sdf)
 * @returns URL for library export
 */
export const getLibraryExportUrl = (
  libraryId: string,
  format: LibraryExportFormat
): string => {
  const url = API_ENDPOINTS.LIBRARY_EXPORT.replace('{id}', libraryId);
  return `${url}?format=${format}`;
};