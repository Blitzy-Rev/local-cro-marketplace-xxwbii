/**
 * API client module for molecule-related operations in the Molecular Data Management and CRO Integration Platform.
 * Provides functions for fetching, creating, updating, and deleting molecules, as well as specialized
 * operations like flagging, similarity search, and substructure search.
 */

import { AxiosResponse } from 'axios'; // axios ^1.4.0
import { apiClient, handleApiError } from './client';
import {
  Molecule,
  MoleculeCreate,
  MoleculeUpdate,
  MoleculeFilter,
  MoleculeListResponse,
  PropertyFilter,
  MoleculeSimilaritySearchParams,
  MoleculeSubstructureSearchParams,
  ApiResponse
} from '../types';

/**
 * API endpoints for molecule-related operations
 */
const API_ENDPOINTS = {
  MOLECULES: '/molecules',
  MY_MOLECULES: '/molecules/me',
  PROPERTIES: '/molecules/properties',
  SIMILARITY_SEARCH: '/molecules/search/similarity',
  SUBSTRUCTURE_SEARCH: '/molecules/search/substructure',
  BULK_OPERATIONS: '/molecules/bulk'
};

/**
 * Fetches a list of molecules with optional filtering, sorting, and pagination
 * 
 * @param filter - Filter criteria for molecules
 * @param page - Page number for pagination (1-indexed)
 * @param pageSize - Number of items per page
 * @returns Promise resolving to paginated list of molecules
 */
export const getMolecules = async (
  filter: MoleculeFilter = {},
  page = 1,
  pageSize = 10
): Promise<AxiosResponse<MoleculeListResponse>> => {
  try {
    // Calculate skip value based on page and pageSize
    const skip = (page - 1) * pageSize;
    
    // Build query parameters
    const queryParams: Record<string, any> = {
      page,
      size: pageSize
    };
    
    // Add filter parameters if provided
    if (filter.smiles_pattern) {
      queryParams.smiles_pattern = filter.smiles_pattern;
    }
    
    if (filter.flag_status) {
      queryParams.flag_status = filter.flag_status;
    }
    
    if (filter.created_by) {
      queryParams.created_by = filter.created_by;
    }
    
    if (filter.created_after) {
      queryParams.created_after = filter.created_after;
    }
    
    if (filter.created_before) {
      queryParams.created_before = filter.created_before;
    }
    
    if (filter.library_ids && filter.library_ids.length > 0) {
      queryParams.library_ids = filter.library_ids.join(',');
    }
    
    if (filter.experiment_ids && filter.experiment_ids.length > 0) {
      queryParams.experiment_ids = filter.experiment_ids.join(',');
    }
    
    if (filter.property_filters && filter.property_filters.length > 0) {
      queryParams.property_filters = JSON.stringify(filter.property_filters);
    }
    
    if (filter.sort_by) {
      queryParams.sort_by = filter.sort_by;
      if (filter.sort_desc !== undefined) {
        queryParams.sort_desc = filter.sort_desc;
      }
    }
    
    return await apiClient.get(API_ENDPOINTS.MOLECULES, { params: queryParams });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches molecules created by the current user
 * 
 * @param page - Page number for pagination (1-indexed)
 * @param pageSize - Number of items per page
 * @returns Promise resolving to paginated list of user's molecules
 */
export const getMyMolecules = async (
  page = 1,
  pageSize = 10
): Promise<AxiosResponse<MoleculeListResponse>> => {
  try {
    // Calculate skip value based on page and pageSize
    const skip = (page - 1) * pageSize;
    
    return await apiClient.get(API_ENDPOINTS.MY_MOLECULES, {
      params: {
        page,
        size: pageSize
      }
    });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches a single molecule by ID
 * 
 * @param id - Molecule ID
 * @returns Promise resolving to molecule data
 */
export const getMoleculeById = async (
  id: string
): Promise<AxiosResponse<ApiResponse<Molecule>>> => {
  try {
    return await apiClient.get(`${API_ENDPOINTS.MOLECULES}/${id}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Creates a new molecule
 * 
 * @param moleculeData - Data for creating a new molecule
 * @returns Promise resolving to created molecule data
 */
export const createMolecule = async (
  moleculeData: MoleculeCreate
): Promise<AxiosResponse<ApiResponse<Molecule>>> => {
  try {
    return await apiClient.post(API_ENDPOINTS.MOLECULES, moleculeData);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates an existing molecule
 * 
 * @param id - Molecule ID
 * @param moleculeData - Data for updating the molecule
 * @returns Promise resolving to updated molecule data
 */
export const updateMolecule = async (
  id: string,
  moleculeData: MoleculeUpdate
): Promise<AxiosResponse<ApiResponse<Molecule>>> => {
  try {
    return await apiClient.put(`${API_ENDPOINTS.MOLECULES}/${id}`, moleculeData);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes a molecule
 * 
 * @param id - Molecule ID
 * @returns Promise resolving to success response
 */
export const deleteMolecule = async (
  id: string
): Promise<AxiosResponse<ApiResponse<{ success: boolean }>>> => {
  try {
    return await apiClient.delete(`${API_ENDPOINTS.MOLECULES}/${id}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates the flag status of a molecule
 * 
 * @param id - Molecule ID
 * @param flagStatus - New flag status
 * @returns Promise resolving to updated molecule data
 */
export const updateMoleculeFlag = async (
  id: string,
  flagStatus: string
): Promise<AxiosResponse<ApiResponse<Molecule>>> => {
  try {
    return await apiClient.put(`${API_ENDPOINTS.MOLECULES}/${id}/flag`, { flag_status: flagStatus });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Performs bulk operations on multiple molecules
 * 
 * @param operationData - Data describing the operation and target molecules
 * @returns Promise resolving to bulk operation results
 */
export const bulkUpdateMolecules = async (
  operationData: {
    operation_type: string;
    molecule_ids: string[];
    flag_status?: string;
    library_id?: string;
    experiment_id?: string;
  }
): Promise<AxiosResponse<ApiResponse<{ success_count: number, failures: any[] }>>> => {
  try {
    return await apiClient.post(API_ENDPOINTS.BULK_OPERATIONS, operationData);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Gets min and max values for a specific property
 * 
 * @param propertyName - Name of the property
 * @returns Promise resolving to property range data
 */
export const getPropertyRanges = async (
  propertyName: string
): Promise<AxiosResponse<ApiResponse<{ min_value: number, max_value: number }>>> => {
  try {
    return await apiClient.get(`${API_ENDPOINTS.PROPERTIES}/ranges/${propertyName}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Gets list of all available property names
 * 
 * @returns Promise resolving to list of property names
 */
export const getAvailableProperties = async (): Promise<AxiosResponse<ApiResponse<string[]>>> => {
  try {
    return await apiClient.get(`${API_ENDPOINTS.PROPERTIES}/available`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Searches for molecules similar to a reference molecule
 * 
 * @param searchParams - Parameters for similarity search
 * @returns Promise resolving to search results with similarity scores
 */
export const searchSimilarMolecules = async (
  searchParams: MoleculeSimilaritySearchParams
): Promise<AxiosResponse<ApiResponse<{ molecules: Molecule[], similarity_scores: Record<string, number> }>>> => {
  try {
    return await apiClient.post(API_ENDPOINTS.SIMILARITY_SEARCH, searchParams);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Searches for molecules containing a specific substructure
 * 
 * @param searchParams - Parameters for substructure search
 * @returns Promise resolving to search results with matching molecules
 */
export const searchSubstructureMolecules = async (
  searchParams: MoleculeSubstructureSearchParams
): Promise<AxiosResponse<ApiResponse<{ molecules: Molecule[] }>>> => {
  try {
    return await apiClient.post(API_ENDPOINTS.SUBSTRUCTURE_SEARCH, searchParams);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Generates a URL for a molecule structure image
 * 
 * @param id - Molecule ID
 * @param options - Image options (width, height, format)
 * @returns URL for the molecule image
 */
export const getMoleculeImageUrl = (
  id: string,
  options: { width?: number, height?: number, format?: string } = {}
): string => {
  let url = `${apiClient.defaults.baseURL}${API_ENDPOINTS.MOLECULES}/${id}/image`;
  
  const queryParams = [];
  if (options.width) {
    queryParams.push(`width=${options.width}`);
  }
  if (options.height) {
    queryParams.push(`height=${options.height}`);
  }
  if (options.format) {
    queryParams.push(`format=${options.format}`);
  }
  
  if (queryParams.length > 0) {
    url = `${url}?${queryParams.join('&')}`;
  }
  
  return url;
};