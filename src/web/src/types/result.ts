/**
 * Result-related TypeScript definitions for the Molecular Data Management and CRO Integration Platform.
 * These types support the experimental result upload, management, and retrieval functionality.
 * @packageDocumentation
 */

import { UUID } from 'uuid';
import { Molecule } from './molecule';
import { User } from './user';

/**
 * Interface for a simplified submission reference to avoid circular dependency.
 * Contains only essential submission data needed for result relationships.
 */
export interface SubmissionReference {
  /** Unique identifier for the submission */
  id: UUID;
  /** ID of the experiment this submission is for */
  experiment_id: UUID;
  /** Current status of the submission */
  status: string;
  /** ID of the CRO assigned to this submission */
  cro_id: number;
}

/**
 * Enumeration of possible result status values.
 * Tracks the lifecycle of experimental results from upload to approval/rejection.
 */
export enum ResultStatus {
  /** Result has been created but no files uploaded yet */
  PENDING = 'pending',
  /** Result files have been uploaded but not yet reviewed */
  UPLOADED = 'uploaded',
  /** Result has been reviewed and approved by the pharma user */
  APPROVED = 'approved',
  /** Result has been reviewed and rejected by the pharma user */
  REJECTED = 'rejected'
}

/**
 * Interface representing a file associated with an experimental result.
 * Includes metadata about the file and its storage location.
 */
export interface ResultFile {
  /** Unique identifier for the file */
  id: UUID;
  /** ID of the result this file belongs to */
  result_id: UUID;
  /** Original name of the uploaded file */
  file_name: string;
  /** Server path where the file is stored */
  file_path: string;
  /** Size of the file in bytes */
  file_size: number;
  /** MIME type or format of the file */
  file_type: string;
  /** Timestamp when the file was uploaded */
  uploaded_at: string;
}

/**
 * Interface for uploading a new result file.
 * Contains the minimum required fields for file upload.
 */
export interface ResultFileCreate {
  /** The file to be uploaded */
  file: File;
  /** Name of the file */
  file_name: string;
  /** MIME type or format of the file */
  file_type: string;
}

/**
 * Interface representing a structured data point from an experimental result.
 * Links specific result values to molecules and includes measurement units.
 */
export interface ResultData {
  /** Unique identifier for the data point */
  id: UUID;
  /** ID of the result this data point belongs to */
  result_id: UUID;
  /** ID of the molecule this data point is for */
  molecule_id: UUID;
  /** Name or type of the measured property */
  data_name: string;
  /** Value of the measured property */
  data_value: number | string;
  /** Unit of measurement for the value */
  data_unit?: string;
  /** Reference to the associated molecule */
  molecule?: Molecule;
}

/**
 * Interface for creating a new result data point.
 * Contains the fields needed to add structured data to a result.
 */
export interface ResultDataCreate {
  /** ID of the molecule this data point is for */
  molecule_id: UUID;
  /** Name or type of the measured property */
  data_name: string;
  /** Value of the measured property */
  data_value: number | string;
  /** Unit of measurement for the value */
  data_unit?: string;
}

/**
 * Interface representing an experimental result with basic information.
 * Includes metadata and references but not detailed data points.
 */
export interface Result {
  /** Unique identifier for the result */
  id: UUID;
  /** ID of the submission this result is for */
  submission_id: UUID;
  /** Current status of the result */
  status: ResultStatus;
  /** Any notes or comments about the result */
  notes?: string;
  /** Timestamp when the result was uploaded */
  uploaded_at: string;
  /** Timestamp when the result was approved, if applicable */
  approved_at?: string;
  /** Reference to the associated submission */
  submission?: SubmissionReference;
  /** Associated result files */
  files?: ResultFile[];
  /** Number of files associated with this result */
  file_count?: number;
}

/**
 * Interface representing an experimental result with detailed information.
 * Includes all data points and relationships.
 */
export interface ResultDetailed extends Result {
  /** Structured data points from the experimental results */
  data_points: ResultData[];
}

/**
 * Interface for creating a new experimental result.
 * Contains all fields needed to upload a new result set.
 */
export interface ResultCreate {
  /** ID of the submission this result is for */
  submission_id: UUID;
  /** Any notes or comments about the result */
  notes?: string;
  /** Files to upload with the result */
  files?: ResultFileCreate[];
  /** Structured data points to associate with the result */
  data?: ResultDataCreate[];
}

/**
 * Interface for updating an existing result.
 * Contains fields that can be modified after upload.
 */
export interface ResultUpdate {
  /** New status to set for the result */
  status?: ResultStatus;
  /** Updated notes or comments */
  notes?: string;
}

/**
 * Interface for approving or rejecting a result.
 * Used in the result review workflow.
 */
export interface ResultApproval {
  /** Whether the result is approved */
  approved: boolean;
  /** Notes explaining the approval decision */
  notes?: string;
}

/**
 * Interface for filtering results by various criteria.
 * Used in search and list views.
 */
export interface ResultFilter {
  /** Filter by submission ID */
  submission_id?: UUID;
  /** Filter by result status */
  status?: ResultStatus;
  /** Filter for results uploaded after this date */
  uploaded_after?: string;
  /** Filter for results uploaded before this date */
  uploaded_before?: string;
  /** Page number for pagination */
  page?: number;
  /** Page size for pagination */
  page_size?: number;
  /** Field to sort by */
  sort_by?: string;
  /** Whether to sort in descending order */
  sort_desc?: boolean;
}

/**
 * Interface for paginated list of results response.
 * Used in API responses.
 */
export interface ResultListResponse {
  /** Array of result items */
  items: Result[];
  /** Total number of results matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
  /** Current filters applied */
  filters: ResultFilter;
}

/**
 * Interface for result state in Redux store.
 * Manages the state of results in the application.
 */
export interface ResultState {
  /** List of results */
  results: Result[];
  /** Currently selected or viewed result */
  currentResult?: ResultDetailed | null;
  /** Whether data is being loaded */
  loading: boolean;
  /** Error message if there was a problem */
  error: string | null;
  /** Total number of results matching current filter */
  totalResults: number;
  /** Current page number */
  currentPage: number;
  /** Number of items per page */
  pageSize: number;
  /** Current filter settings */
  filter: ResultFilter;
}