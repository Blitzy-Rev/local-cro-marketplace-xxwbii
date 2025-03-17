/**
 * Submission-related TypeScript definitions for the Molecular Data Management and CRO Integration Platform.
 * These types represent submissions to Contract Research Organizations (CROs), including submission details,
 * quotes, and status tracking for the frontend application.
 * @packageDocumentation
 */

import { Experiment } from './experiment';
import { User } from './user';
import { Result } from './result';

/**
 * Enumeration of possible submission status values.
 * Tracks the lifecycle of submissions to Contract Research Organizations.
 */
export enum SubmissionStatus {
  /** Initial state when submission is created but not yet sent to CRO */
  PENDING = 'pending',
  /** Submission was rejected by the CRO */
  REJECTED = 'rejected',
  /** CRO has provided a quote for the submission */
  QUOTE_PROVIDED = 'quote_provided',
  /** Quote was rejected by the pharma user */
  QUOTE_REJECTED = 'quote_rejected',
  /** Quote was approved by the pharma user and work can begin */
  APPROVED = 'approved',
  /** Work is currently being performed by the CRO */
  IN_PROGRESS = 'in_progress',
  /** All work has been completed by the CRO */
  COMPLETED = 'completed',
  /** Submission was cancelled before completion */
  CANCELLED = 'cancelled'
}

/**
 * Interface representing a detail for a submission.
 * Contains custom fields for specific experimental requirements.
 */
export interface SubmissionDetail {
  /** Unique identifier for the detail */
  id: string;
  /** ID of the submission this detail belongs to */
  submission_id: string;
  /** Name of the detail field */
  detail_name: string;
  /** Value of the detail field */
  detail_value: string;
}

/**
 * Interface for creating a new submission detail.
 * Used when adding custom fields to a submission.
 */
export interface SubmissionDetailCreate {
  /** Name of the detail field */
  detail_name: string;
  /** Value of the detail field */
  detail_value: string;
}

/**
 * Interface representing a submission to a CRO.
 * Contains the core data structure for experiment submissions.
 */
export interface Submission {
  /** Unique identifier for the submission */
  id: string;
  /** ID of the experiment this submission is for */
  experiment_id: string;
  /** ID of the CRO user assigned to this submission */
  cro_id: number;
  /** Current status of the submission */
  status: SubmissionStatus;
  /** Timestamp when the submission was created */
  submitted_at: string;
  /** Timestamp when the submission was last updated */
  updated_at: string;
  /** Price quoted by the CRO (if a quote has been provided) */
  price?: number;
  /** Estimated turnaround time in days */
  turnaround_days?: number;
  /** Additional notes about the submission */
  notes?: string;
  /** Related experiment data */
  experiment?: Experiment;
  /** CRO user information */
  cro?: User;
  /** Custom details for this submission */
  details?: SubmissionDetail[];
}

/**
 * Interface representing a submission with detailed information including results.
 * Extends the base Submission interface with result data.
 */
export interface SubmissionDetailed extends Submission {
  /** Results uploaded for this submission */
  results?: Result[];
}

/**
 * Interface for creating a new submission.
 * Contains the fields required to submit an experiment to a CRO.
 */
export interface SubmissionCreate {
  /** ID of the experiment to submit */
  experiment_id: string;
  /** ID of the CRO to assign the submission to */
  cro_id: number;
  /** Custom details for the submission */
  details?: SubmissionDetailCreate[];
  /** Additional notes about the submission */
  notes?: string;
}

/**
 * Interface for updating an existing submission.
 * Contains fields that can be modified after creation.
 */
export interface SubmissionUpdate {
  /** ID of the CRO to assign the submission to */
  cro_id?: number;
  /** New status for the submission */
  status?: SubmissionStatus;
  /** Updated custom details */
  details?: SubmissionDetailCreate[];
  /** Updated notes */
  notes?: string;
}

/**
 * Interface for filtering submissions by various criteria.
 * Used in search and list views.
 */
export interface SubmissionFilter {
  /** Filter by experiment ID */
  experiment_id?: string;
  /** Filter by CRO ID */
  cro_id?: number;
  /** Filter by submission status */
  status?: SubmissionStatus;
  /** Filter for submissions created after this date */
  submitted_after?: string;
  /** Filter for submissions created before this date */
  submitted_before?: string;
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
 * Interface for paginated list of submissions response.
 * Used in API responses.
 */
export interface SubmissionListResponse {
  /** Array of submission items */
  items: Submission[];
  /** Total number of submissions matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
  /** Current filters applied */
  filters: SubmissionFilter;
}

/**
 * Interface for CRO to provide a quote for a submission.
 * Used by CRO users to respond to submissions with pricing.
 */
export interface QuoteProvide {
  /** Quoted price for the requested work */
  price: number;
  /** Estimated turnaround time in business days */
  turnaround_days: number;
  /** Additional notes about the quote */
  notes?: string;
}

/**
 * Interface for pharma user to respond to a quote.
 * Used to approve or reject quotes from CROs.
 */
export interface QuoteResponse {
  /** Whether the quote is approved */
  approved: boolean;
  /** Optional notes about the approval/rejection */
  notes?: string;
}

/**
 * Interface for updating the status of a submission.
 * Used for status changes during the submission lifecycle.
 */
export interface SubmissionStatusUpdate {
  /** New status for the submission */
  status: SubmissionStatus;
  /** Notes about the status change */
  notes?: string;
}