/**
 * API client module for experiment-related operations in the Molecular Data Management and CRO Integration Platform.
 * This module provides functions for fetching, creating, updating, and managing experiments,
 * including adding molecules to experiments and updating experiment statuses.
 * 
 * @module api/experiments
 */

import { AxiosResponse } from 'axios'; // axios ^1.4.0
import { apiClient, handleApiError } from './client';
import { 
  Experiment,
  ExperimentDetail,
  ExperimentCreate,
  ExperimentUpdate,
  ExperimentFilter,
  ExperimentListResponse,
  ExperimentTypeListResponse,
  ExperimentMoleculeOperation,
  ExperimentStatusUpdate,
  ExperimentStatus,
  ApiResponse
} from '../types';

/**
 * API endpoint definitions for experiment-related operations
 */
const API_ENDPOINTS = {
  EXPERIMENTS: '/experiments',
  MY_EXPERIMENTS: '/experiments/me',
  EXPERIMENT_TYPES: '/experiments/types',
  EXPERIMENT_MOLECULES: '/experiments/{id}/molecules',
  EXPERIMENT_STATUS: '/experiments/{id}/status'
};

/**
 * Fetches a list of experiments with optional filtering, sorting, and pagination
 * 
 * @param filter - Filter criteria for experiments
 * @param page - Page number for pagination (default: 1)
 * @param pageSize - Number of items per page (default: 10)
 * @returns Promise resolving to paginated list of experiments
 */
export const getExperiments = async (
  filter: ExperimentFilter = {},
  page: number = 1,
  pageSize: number = 10
): Promise<AxiosResponse<ExperimentListResponse>> => {
  try {
    const params = {
      ...filter,
      page,
      page_size: pageSize
    };

    return await apiClient.get<ExperimentListResponse>(API_ENDPOINTS.EXPERIMENTS, { params });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches experiments created by the current user
 * 
 * @param filter - Filter criteria for experiments
 * @param page - Page number for pagination (default: 1)
 * @param pageSize - Number of items per page (default: 10)
 * @returns Promise resolving to paginated list of user's experiments
 */
export const getMyExperiments = async (
  filter: ExperimentFilter = {},
  page: number = 1,
  pageSize: number = 10
): Promise<AxiosResponse<ExperimentListResponse>> => {
  try {
    const params = {
      ...filter,
      page,
      page_size: pageSize
    };

    return await apiClient.get<ExperimentListResponse>(API_ENDPOINTS.MY_EXPERIMENTS, { params });
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches a single experiment by ID with detailed information
 * 
 * @param id - Experiment ID
 * @returns Promise resolving to detailed experiment data
 */
export const getExperimentById = async (
  id: string
): Promise<AxiosResponse<ApiResponse<ExperimentDetail>>> => {
  try {
    return await apiClient.get<ApiResponse<ExperimentDetail>>(`${API_ENDPOINTS.EXPERIMENTS}/${id}`);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Fetches available experiment types
 * 
 * @returns Promise resolving to list of experiment types
 */
export const getExperimentTypes = async (): Promise<AxiosResponse<ExperimentTypeListResponse>> => {
  try {
    return await apiClient.get<ExperimentTypeListResponse>(API_ENDPOINTS.EXPERIMENT_TYPES);
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Creates a new experiment
 * 
 * @param experimentData - Data for creating a new experiment
 * @returns Promise resolving to created experiment data
 */
export const createExperiment = async (
  experimentData: ExperimentCreate
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  try {
    return await apiClient.post<ApiResponse<Experiment>>(
      API_ENDPOINTS.EXPERIMENTS,
      experimentData
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates an existing experiment
 * 
 * @param id - Experiment ID
 * @param experimentData - Data for updating the experiment
 * @returns Promise resolving to updated experiment data
 */
export const updateExperiment = async (
  id: string,
  experimentData: ExperimentUpdate
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  try {
    return await apiClient.put<ApiResponse<Experiment>>(
      `${API_ENDPOINTS.EXPERIMENTS}/${id}`,
      experimentData
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Deletes an experiment
 * 
 * @param id - Experiment ID
 * @returns Promise resolving to success response
 */
export const deleteExperiment = async (
  id: string
): Promise<AxiosResponse<ApiResponse<{ success: boolean }>>> => {
  try {
    return await apiClient.delete<ApiResponse<{ success: boolean }>>(
      `${API_ENDPOINTS.EXPERIMENTS}/${id}`
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Adds molecules to an existing experiment
 * 
 * @param id - Experiment ID
 * @param moleculeData - Molecule IDs to add to the experiment
 * @returns Promise resolving to operation result
 */
export const addMoleculesToExperiment = async (
  id: string,
  moleculeData: ExperimentMoleculeOperation
): Promise<AxiosResponse<ApiResponse<{ success: boolean, added_count: number }>>> => {
  try {
    const url = API_ENDPOINTS.EXPERIMENT_MOLECULES.replace('{id}', id);
    return await apiClient.post<ApiResponse<{ success: boolean, added_count: number }>>(
      url,
      moleculeData
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Removes molecules from an existing experiment
 * 
 * @param id - Experiment ID
 * @param moleculeData - Molecule IDs to remove from the experiment
 * @returns Promise resolving to operation result
 */
export const removeMoleculesFromExperiment = async (
  id: string,
  moleculeData: ExperimentMoleculeOperation
): Promise<AxiosResponse<ApiResponse<{ success: boolean, removed_count: number }>>> => {
  try {
    const url = API_ENDPOINTS.EXPERIMENT_MOLECULES.replace('{id}', id);
    return await apiClient.delete<ApiResponse<{ success: boolean, removed_count: number }>>(
      url,
      { data: moleculeData }
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates the status of an experiment
 * 
 * @param id - Experiment ID
 * @param statusData - New status and optional notes
 * @returns Promise resolving to updated experiment data
 */
export const updateExperimentStatus = async (
  id: string,
  statusData: ExperimentStatusUpdate
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  try {
    const url = API_ENDPOINTS.EXPERIMENT_STATUS.replace('{id}', id);
    return await apiClient.put<ApiResponse<Experiment>>(
      url,
      statusData
    );
  } catch (error) {
    return handleApiError(error);
  }
};

/**
 * Updates experiment status to QUEUED
 * 
 * @param id - Experiment ID
 * @param notes - Optional notes about the status change
 * @returns Promise resolving to updated experiment data
 */
export const queueExperiment = async (
  id: string,
  notes?: string
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  const statusUpdate: ExperimentStatusUpdate = {
    status: ExperimentStatus.QUEUED,
    notes
  };
  return updateExperimentStatus(id, statusUpdate);
};

/**
 * Updates experiment status to SUBMITTED
 * 
 * @param id - Experiment ID
 * @param notes - Optional notes about the status change
 * @returns Promise resolving to updated experiment data
 */
export const submitExperiment = async (
  id: string,
  notes?: string
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  const statusUpdate: ExperimentStatusUpdate = {
    status: ExperimentStatus.SUBMITTED,
    notes
  };
  return updateExperimentStatus(id, statusUpdate);
};

/**
 * Updates experiment status to CANCELLED
 * 
 * @param id - Experiment ID
 * @param notes - Optional notes about the status change
 * @returns Promise resolving to updated experiment data
 */
export const cancelExperiment = async (
  id: string,
  notes?: string
): Promise<AxiosResponse<ApiResponse<Experiment>>> => {
  const statusUpdate: ExperimentStatusUpdate = {
    status: ExperimentStatus.CANCELLED,
    notes
  };
  return updateExperimentStatus(id, statusUpdate);
};