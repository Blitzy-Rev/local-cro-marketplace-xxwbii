/**
 * API client module for CSV file upload, mapping, and processing operations.
 * Provides functions to interact with the backend CSV endpoints for the
 * Molecular Data Management and CRO Integration Platform.
 * 
 * @module api/csv
 */

import axios from 'axios'; // axios ^1.4.0
import { apiClient, handleApiError } from './client';
import { ApiResponse } from '../types';

/**
 * Base path for CSV API endpoints
 */
const CSV_API_PATH = '/csv';

/**
 * Uploads a CSV file to the server for processing.
 * 
 * @param file - The CSV file to upload
 * @returns Promise resolving to an object containing the file ID, headers, and row count
 */
export const uploadCSV = async (file: File): Promise<{fileId: string, headers: string[], rowCount: number}> => {
  try {
    // Create a FormData object to send the file
    const formData = new FormData();
    formData.append('file', file);
    
    // Set content type to multipart/form-data (handled automatically by axios)
    const response = await apiClient.post<ApiResponse<{fileId: string, headers: string[], rowCount: number}>>(
      `${CSV_API_PATH}/upload`,
      formData,
      {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      }
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Maps CSV columns to system properties and starts processing.
 * 
 * @param params - Object containing fileId and column mapping information
 * @returns Promise resolving to an object containing the job ID and status
 */
export const mapCSVColumns = async (
  { fileId, mapping }: { fileId: string, mapping: Array<{csvColumn: string, systemProperty: string}> }
): Promise<{jobId: string, status: string}> => {
  try {
    const response = await apiClient.post<ApiResponse<{jobId: string, status: string}>>(
      `${CSV_API_PATH}/map`,
      {
        fileId,
        mapping
      }
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Gets the status of a CSV processing job.
 * 
 * @param jobId - The ID of the CSV processing job
 * @returns Promise resolving to an object containing the status, progress, and summary
 */
export const getCSVProcessingStatus = async (
  jobId: string
): Promise<{status: string, progress: number, summary: any}> => {
  try {
    const response = await apiClient.get<ApiResponse<{status: string, progress: number, summary: any}>>(
      `${CSV_API_PATH}/status/${jobId}`
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Cancels an ongoing CSV processing job.
 * 
 * @param jobId - The ID of the CSV processing job to cancel
 * @returns Promise resolving to an object containing a cancellation message
 */
export const cancelCSVProcessing = async (
  jobId: string
): Promise<{message: string}> => {
  try {
    const response = await apiClient.post<ApiResponse<{message: string}>>(
      `${CSV_API_PATH}/cancel/${jobId}`
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Gets list of available system properties for CSV mapping.
 * These properties define the possible targets for mapping CSV columns.
 * 
 * @returns Promise resolving to an object containing an array of available properties
 */
export const getAvailableProperties = async (): Promise<{
  properties: Array<{name: string, display_name: string, data_type: string, unit?: string}>
}> => {
  try {
    const response = await apiClient.get<ApiResponse<{
      properties: Array<{name: string, display_name: string, data_type: string, unit?: string}>
    }>>(
      `${CSV_API_PATH}/properties`
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Generates a detailed report for a CSV import job.
 * This is useful for obtaining comprehensive information about a completed import.
 * 
 * @param jobId - The ID of the CSV processing job
 * @returns Promise resolving to an object containing a report generation message
 */
export const generateCSVImportReport = async (
  jobId: string
): Promise<{message: string}> => {
  try {
    const response = await apiClient.post<ApiResponse<{message: string}>>(
      `${CSV_API_PATH}/report/${jobId}`
    );
    
    return response.data.data;
  } catch (error) {
    return handleApiError(error);
  }
};