/**
 * API client module for result-related operations in the Molecular Data Management and CRO Integration Platform.
 * Provides functions for creating, retrieving, updating, and managing experimental results
 * uploaded by CRO users, including file uploads, structured data management, and result approval workflows.
 * @module api/results
 */

import { apiClient, handleApiError } from './client';
import { 
  Result, 
  ResultDetailed, 
  ResultCreate, 
  ResultUpdate, 
  ResultFilter, 
  ResultListResponse, 
  ResultApproval, 
  ResultDataCreate, 
  ResultFileCreate 
} from '../types/result';
import { ApiResponse } from '../types';

/** Base URL for results endpoints */
const BASE_URL = '/results';

/**
 * Fetches a paginated list of results with optional filtering
 * @param filters - Optional filter criteria for results
 * @returns Promise resolving to paginated list of results
 */
export const getResults = async (filters?: ResultFilter): Promise<ApiResponse<ResultListResponse>> => {
  try {
    // Construct query parameters from filters
    const params = filters ? { ...filters } : {};
    
    const response = await apiClient.get(`${BASE_URL}`, { params });
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches a single result by ID with detailed information
 * @param resultId - Unique identifier of the result
 * @returns Promise resolving to detailed result data
 */
export const getResult = async (resultId: string): Promise<ApiResponse<ResultDetailed>> => {
  try {
    const response = await apiClient.get(`${BASE_URL}/${resultId}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches results for a specific submission
 * @param submissionId - Unique identifier of the submission
 * @param options - Optional pagination parameters
 * @returns Promise resolving to list of results for the submission
 */
export const getResultsBySubmission = async (
  submissionId: string,
  options: { page?: number, pageSize?: number } = {}
): Promise<ApiResponse<ResultListResponse>> => {
  try {
    // Construct query parameters from options
    const params = { ...options };
    
    const response = await apiClient.get(`${BASE_URL}/submission/${submissionId}`, { params });
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Creates a new result for a submission
 * @param resultData - Data for creating the result
 * @returns Promise resolving to the created result
 */
export const createResult = async (resultData: ResultCreate): Promise<ApiResponse<ResultDetailed>> => {
  try {
    const response = await apiClient.post(`${BASE_URL}`, resultData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates an existing result
 * @param resultId - Unique identifier of the result to update
 * @param resultData - Updated result data
 * @returns Promise resolving to the updated result
 */
export const updateResult = async (
  resultId: string,
  resultData: ResultUpdate
): Promise<ApiResponse<ResultDetailed>> => {
  try {
    const response = await apiClient.put(`${BASE_URL}/${resultId}`, resultData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Approves or rejects a result
 * @param resultId - Unique identifier of the result
 * @param approvalData - Approval data including decision and optional notes
 * @returns Promise resolving to the updated result
 */
export const approveRejectResult = async (
  resultId: string,
  approvalData: ResultApproval
): Promise<ApiResponse<ResultDetailed>> => {
  try {
    const response = await apiClient.post(`${BASE_URL}/${resultId}/approve`, approvalData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Approves or rejects multiple results in a batch
 * @param batchData - Data containing result IDs and approval decision
 * @returns Promise resolving to success message with count of processed results
 */
export const batchApproveRejectResults = async (
  batchData: { result_ids: string[], approved: boolean, notes?: string }
): Promise<ApiResponse<{ message: string }>> => {
  try {
    const response = await apiClient.post(`${BASE_URL}/batch-approve`, batchData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Uploads a file for a result
 * @param resultId - Unique identifier of the result
 * @param file - File to upload
 * @param fileData - Metadata for the file
 * @returns Promise resolving to file metadata and database record
 */
export const uploadResultFile = async (
  resultId: string,
  file: File,
  fileData: ResultFileCreate
): Promise<ApiResponse<any>> => {
  try {
    // Create FormData object for file upload
    const formData = new FormData();
    formData.append('file', file);
    formData.append('file_name', fileData.file_name);
    formData.append('file_type', fileData.file_type);
    
    const response = await apiClient.post(`${BASE_URL}/${resultId}/files`, formData, {
      headers: {
        'Content-Type': 'multipart/form-data'
      }
    });
    
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Downloads a file associated with a result
 * @param fileId - Unique identifier of the file
 * @returns Promise resolving to file blob for download
 */
export const getResultFile = async (fileId: string): Promise<Blob> => {
  try {
    const response = await apiClient.get(`${BASE_URL}/files/${fileId}`, {
      responseType: 'blob'
    });
    
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes a file associated with a result
 * @param fileId - Unique identifier of the file
 * @returns Promise resolving to success message
 */
export const deleteResultFile = async (fileId: string): Promise<ApiResponse<{ message: string }>> => {
  try {
    const response = await apiClient.delete(`${BASE_URL}/files/${fileId}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Adds structured data to a result for a specific molecule
 * @param resultId - Unique identifier of the result
 * @param dataPoint - Data point to add
 * @returns Promise resolving to created data record
 */
export const addResultData = async (
  resultId: string,
  dataPoint: ResultDataCreate
): Promise<ApiResponse<any>> => {
  try {
    const response = await apiClient.post(`${BASE_URL}/${resultId}/data`, dataPoint);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Gets structured data for a result, optionally filtered by molecule
 * @param resultId - Unique identifier of the result
 * @param options - Optional filter options
 * @returns Promise resolving to list of structured data records
 */
export const getResultData = async (
  resultId: string,
  options: { molecule_id?: string } = {}
): Promise<ApiResponse<any[]>> => {
  try {
    // Construct query parameters from options
    const params = { ...options };
    
    const response = await apiClient.get(`${BASE_URL}/${resultId}/data`, { params });
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Exports result data in various formats (CSV, XLSX, JSON)
 * @param options - Export options including result ID and format
 * @returns Promise resolving to exported data blob
 */
export const exportResult = async (
  options: { result_id: string, format: 'csv' | 'xlsx' | 'json' }
): Promise<Blob> => {
  try {
    const response = await apiClient.post(`${BASE_URL}/export`, options, {
      responseType: 'blob'
    });
    
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Generates analysis of result data with summary statistics
 * @param resultId - Unique identifier of the result
 * @returns Promise resolving to analysis of result data
 */
export const analyzeResult = async (resultId: string): Promise<ApiResponse<any>> => {
  try {
    const response = await apiClient.get(`${BASE_URL}/${resultId}/analyze`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};