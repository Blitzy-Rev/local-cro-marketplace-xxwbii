/**
 * API client module for submission-related operations in the Molecular Data Management and CRO Integration Platform.
 * This module provides functions for creating, retrieving, updating, and managing submissions to 
 * Contract Research Organizations (CROs), including quote handling and status transitions.
 * 
 * @packageDocumentation
 */

import { apiClient, handleApiError } from './client'; // axios ^1.4.0
import {
  Submission,
  SubmissionDetailed,
  SubmissionCreate,
  SubmissionUpdate,
  SubmissionFilter,
  SubmissionListResponse,
  QuoteProvide,
  QuoteResponse,
  SubmissionStatusUpdate,
  SubmissionStatus
} from '../types/submission';
import { ApiResponse } from '../types';

/**
 * Base URL for submission API endpoints
 */
const BASE_URL = '/submissions';

/**
 * Fetches a paginated list of submissions with optional filtering
 * 
 * @param filters - Optional filters to apply to the submission list
 * @returns Promise resolving to paginated list of submissions
 */
export const getSubmissions = async (filters: SubmissionFilter = {}): Promise<ApiResponse<SubmissionListResponse>> => {
  try {
    // Build query parameters from filters
    const queryParams = new URLSearchParams();
    if (filters.experiment_id) queryParams.append('experiment_id', filters.experiment_id);
    if (filters.cro_id) queryParams.append('cro_id', filters.cro_id.toString());
    if (filters.status) queryParams.append('status', filters.status);
    if (filters.submitted_after) queryParams.append('submitted_after', filters.submitted_after);
    if (filters.submitted_before) queryParams.append('submitted_before', filters.submitted_before);
    if (filters.page) queryParams.append('page', filters.page.toString());
    if (filters.page_size) queryParams.append('page_size', filters.page_size.toString());
    if (filters.sort_by) queryParams.append('sort_by', filters.sort_by);
    if (filters.sort_desc !== undefined) queryParams.append('sort_desc', filters.sort_desc.toString());

    const query = queryParams.toString() ? `?${queryParams.toString()}` : '';
    
    const response = await apiClient.get<ApiResponse<SubmissionListResponse>>(`${BASE_URL}${query}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches a single submission by ID with detailed information
 * 
 * @param submissionId - ID of the submission to fetch
 * @returns Promise resolving to detailed submission data
 */
export const getSubmission = async (submissionId: string): Promise<ApiResponse<SubmissionDetailed>> => {
  try {
    const response = await apiClient.get<ApiResponse<SubmissionDetailed>>(`${BASE_URL}/${submissionId}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Creates a new submission for an experiment to a CRO
 * 
 * @param submissionData - Data for the new submission
 * @returns Promise resolving to the created submission
 */
export const createSubmission = async (submissionData: SubmissionCreate): Promise<ApiResponse<Submission>> => {
  try {
    const response = await apiClient.post<ApiResponse<Submission>>(BASE_URL, submissionData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates an existing submission
 * 
 * @param submissionId - ID of the submission to update
 * @param submissionData - New data for the submission
 * @returns Promise resolving to the updated submission
 */
export const updateSubmission = async (
  submissionId: string,
  submissionData: SubmissionUpdate
): Promise<ApiResponse<Submission>> => {
  try {
    const response = await apiClient.put<ApiResponse<Submission>>(`${BASE_URL}/${submissionId}`, submissionData);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates the status of a submission
 * 
 * @param submissionId - ID of the submission to update
 * @param statusData - New status information
 * @returns Promise resolving to the updated submission
 */
export const updateSubmissionStatus = async (
  submissionId: string,
  statusData: SubmissionStatusUpdate
): Promise<ApiResponse<Submission>> => {
  try {
    const response = await apiClient.put<ApiResponse<Submission>>(
      `${BASE_URL}/${submissionId}/status`,
      statusData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Allows a CRO to provide a quote for a submission
 * 
 * @param submissionId - ID of the submission to quote
 * @param quoteData - Quote information including price and turnaround time
 * @returns Promise resolving to the updated submission with quote information
 */
export const provideQuote = async (
  submissionId: string,
  quoteData: QuoteProvide
): Promise<ApiResponse<Submission>> => {
  try {
    const response = await apiClient.post<ApiResponse<Submission>>(
      `${BASE_URL}/${submissionId}/quote`,
      quoteData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Allows a pharma user to approve or reject a quote
 * 
 * @param submissionId - ID of the submission with the quote
 * @param responseData - Response information including approval status
 * @returns Promise resolving to the updated submission
 */
export const respondToQuote = async (
  submissionId: string,
  responseData: QuoteResponse
): Promise<ApiResponse<Submission>> => {
  try {
    const response = await apiClient.post<ApiResponse<Submission>>(
      `${BASE_URL}/${submissionId}/quote/response`,
      responseData
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches submissions for a specific experiment
 * 
 * @param experimentId - ID of the experiment
 * @param filters - Optional filters to apply
 * @returns Promise resolving to list of submissions for the experiment
 */
export const getSubmissionsByExperiment = async (
  experimentId: string,
  filters: SubmissionFilter = {}
): Promise<ApiResponse<SubmissionListResponse>> => {
  try {
    // Build query parameters from filters
    const queryParams = new URLSearchParams();
    if (filters.status) queryParams.append('status', filters.status);
    if (filters.submitted_after) queryParams.append('submitted_after', filters.submitted_after);
    if (filters.submitted_before) queryParams.append('submitted_before', filters.submitted_before);
    if (filters.page) queryParams.append('page', filters.page.toString());
    if (filters.page_size) queryParams.append('page_size', filters.page_size.toString());
    if (filters.sort_by) queryParams.append('sort_by', filters.sort_by);
    if (filters.sort_desc !== undefined) queryParams.append('sort_desc', filters.sort_desc.toString());

    const query = queryParams.toString() ? `?${queryParams.toString()}` : '';
    
    const response = await apiClient.get<ApiResponse<SubmissionListResponse>>(
      `${BASE_URL}/experiment/${experimentId}${query}`
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches submissions assigned to the current CRO user
 * 
 * @param filters - Optional filters to apply
 * @returns Promise resolving to list of submissions for the CRO
 */
export const getSubmissionsByCRO = async (
  filters: SubmissionFilter = {}
): Promise<ApiResponse<SubmissionListResponse>> => {
  try {
    // Build query parameters from filters
    const queryParams = new URLSearchParams();
    if (filters.status) queryParams.append('status', filters.status);
    if (filters.submitted_after) queryParams.append('submitted_after', filters.submitted_after);
    if (filters.submitted_before) queryParams.append('submitted_before', filters.submitted_before);
    if (filters.page) queryParams.append('page', filters.page.toString());
    if (filters.page_size) queryParams.append('page_size', filters.page_size.toString());
    if (filters.sort_by) queryParams.append('sort_by', filters.sort_by);
    if (filters.sort_desc !== undefined) queryParams.append('sort_desc', filters.sort_desc.toString());

    const query = queryParams.toString() ? `?${queryParams.toString()}` : '';
    
    const response = await apiClient.get<ApiResponse<SubmissionListResponse>>(`${BASE_URL}/cro${query}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches submissions with a specific status
 * 
 * @param status - Status to filter by
 * @param filters - Additional optional filters
 * @returns Promise resolving to list of submissions with the specified status
 */
export const getSubmissionsByStatus = async (
  status: string,
  filters: SubmissionFilter = {}
): Promise<ApiResponse<SubmissionListResponse>> => {
  try {
    // Build query parameters from filters
    const queryParams = new URLSearchParams();
    if (filters.experiment_id) queryParams.append('experiment_id', filters.experiment_id);
    if (filters.cro_id) queryParams.append('cro_id', filters.cro_id.toString());
    if (filters.submitted_after) queryParams.append('submitted_after', filters.submitted_after);
    if (filters.submitted_before) queryParams.append('submitted_before', filters.submitted_before);
    if (filters.page) queryParams.append('page', filters.page.toString());
    if (filters.page_size) queryParams.append('page_size', filters.page_size.toString());
    if (filters.sort_by) queryParams.append('sort_by', filters.sort_by);
    if (filters.sort_desc !== undefined) queryParams.append('sort_desc', filters.sort_desc.toString());

    const query = queryParams.toString() ? `?${queryParams.toString()}` : '';
    
    const response = await apiClient.get<ApiResponse<SubmissionListResponse>>(
      `${BASE_URL}/status/${status}${query}`
    );
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches results for a specific submission
 * 
 * @param submissionId - ID of the submission
 * @param options - Pagination options
 * @returns Promise resolving to list of results for the submission
 */
export const getSubmissionResults = async (
  submissionId: string,
  options: { page?: number, pageSize?: number } = {}
): Promise<ApiResponse<any>> => {
  try {
    // Build query parameters from pagination options
    const queryParams = new URLSearchParams();
    if (options.page) queryParams.append('page', options.page.toString());
    if (options.pageSize) queryParams.append('page_size', options.pageSize.toString());

    const query = queryParams.toString() ? `?${queryParams.toString()}` : '';
    
    const response = await apiClient.get<ApiResponse<any>>(`${BASE_URL}/${submissionId}/results${query}`);
    return response.data;
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Cancels a submission
 * 
 * @param submissionId - ID of the submission to cancel
 * @param notes - Optional notes about the cancellation reason
 * @returns Promise resolving to the updated submission
 */
export const cancelSubmission = async (
  submissionId: string,
  notes: string = 'Cancelled by user'
): Promise<ApiResponse<Submission>> => {
  const statusUpdate: SubmissionStatusUpdate = {
    status: SubmissionStatus.CANCELLED,
    notes
  };
  
  return updateSubmissionStatus(submissionId, statusUpdate);
};