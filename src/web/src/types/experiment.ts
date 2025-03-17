/**
 * Experiment-related TypeScript definitions for the Molecular Data Management and CRO Integration Platform.
 * These types represent experiments, their parameters, statuses, and relationships with molecules and CROs.
 * @packageDocumentation
 */

import { Molecule } from './molecule';
import { User } from './user';

/**
 * Enumeration of experiment status values.
 * Represents the possible states of an experiment throughout its lifecycle.
 */
export enum ExperimentStatus {
  /** Initial state when experiment is being defined */
  DRAFT = 'draft',
  /** Experiment is queued for submission to CRO */
  QUEUED = 'queued',
  /** Experiment has been submitted to CRO */
  SUBMITTED = 'submitted',
  /** Experiment was rejected by CRO */
  REJECTED = 'rejected',
  /** CRO has provided a quote, awaiting approval */
  QUOTE_PENDING = 'quote_pending',
  /** Quote was rejected by the pharma user */
  QUOTE_REJECTED = 'quote_rejected',
  /** Experiment is being processed by CRO */
  IN_PROGRESS = 'in_progress',
  /** Experiment is complete, results in preparation */
  RESULTS_PENDING = 'results_pending',
  /** Results have been uploaded and are available */
  RESULTS_AVAILABLE = 'results_available',
  /** Results were rejected and need revision */
  RESULTS_REJECTED = 'results_rejected',
  /** Experiment is fully completed and reviewed */
  COMPLETED = 'completed',
  /** Experiment was cancelled */
  CANCELLED = 'cancelled'
}

/**
 * Interface for experiment type definition.
 * Represents the categories and types of experiments available in the system.
 */
export interface ExperimentType {
  /** Unique identifier for the experiment type */
  id: string;
  /** Name of the experiment type */
  name: string;
  /** Detailed description of the experiment type */
  description: string;
  /** Category grouping for the experiment type */
  category: string;
}

/**
 * Interface for experiment parameters.
 * Represents configurable parameters for experiments.
 */
export interface ExperimentParameter {
  /** ID of the experiment this parameter belongs to */
  experiment_id: string;
  /** Name of the parameter */
  parameter_name: string;
  /** Value of the parameter */
  parameter_value: string | number;
}

/**
 * Interface for the relationship between experiments and molecules.
 * Represents molecules included in an experiment.
 */
export interface ExperimentMolecule {
  /** ID of the experiment */
  experiment_id: string;
  /** ID of the molecule */
  molecule_id: string;
  /** Timestamp when the molecule was added to the experiment */
  added_at: string;
  /** Optional molecule data for expanded information */
  molecule?: Molecule;
}

/**
 * Interface for experiment data with basic information.
 * Core data structure representing an experiment in the system.
 */
export interface Experiment {
  /** Unique identifier for the experiment */
  id: string;
  /** Name of the experiment */
  name: string;
  /** ID of the experiment type */
  type_id: string;
  /** Current status of the experiment */
  status: ExperimentStatus;
  /** ID of the user who created the experiment */
  created_by: number;
  /** Timestamp when the experiment was created */
  created_at: string;
  /** Timestamp when the experiment was last updated */
  updated_at?: string;
  /** Experiment type information */
  type?: ExperimentType;
  /** Creator user information */
  creator?: User;
  /** Parameters for this experiment */
  parameters?: ExperimentParameter[];
  /** Count of molecules in this experiment */
  molecule_count?: number;
}

/**
 * Interface for experiment data with detailed information.
 * Provides comprehensive information about an experiment and its relationships.
 */
export interface ExperimentDetail extends Experiment {
  /** Molecules included in this experiment */
  molecules: ExperimentMolecule[];
  /** Submissions related to this experiment */
  submissions: any[]; // This would be replaced with a Submission interface
}

/**
 * Interface for creating a new experiment.
 * Contains the required fields to create an experiment.
 */
export interface ExperimentCreate {
  /** Name of the experiment */
  name: string;
  /** ID of the experiment type */
  type_id: string;
  /** Parameters for the experiment */
  parameters: Partial<ExperimentParameter>[];
  /** IDs of molecules to include in the experiment */
  molecule_ids: string[];
}

/**
 * Interface for updating an existing experiment.
 * Contains fields that can be modified after creation.
 */
export interface ExperimentUpdate {
  /** Updated name of the experiment */
  name?: string;
  /** Updated status of the experiment */
  status?: ExperimentStatus;
  /** Updated parameters for the experiment */
  parameters?: Partial<ExperimentParameter>[];
}

/**
 * Interface for adding molecules to an experiment.
 * Used in API requests to add molecules to an experiment.
 */
export interface ExperimentMoleculeOperation {
  /** IDs of molecules to add to the experiment */
  molecule_ids: string[];
}

/**
 * Interface for filtering experiments by various criteria.
 * Used in API requests to filter experiments.
 */
export interface ExperimentFilter {
  /** Filter by experiment name (partial match) */
  name?: string;
  /** Filter by experiment type ID */
  type_id?: string;
  /** Filter by experiment status */
  status?: ExperimentStatus;
  /** Filter by creator user ID */
  created_by?: number;
  /** Filter by creation date (after this date) */
  created_after?: string;
  /** Filter by creation date (before this date) */
  created_before?: string;
  /** Filter by experiments containing a specific molecule */
  contains_molecule_id?: string;
  /** Page number for pagination */
  page?: number;
  /** Number of items per page */
  page_size?: number;
  /** Field to sort by */
  sort_by?: string;
  /** Sort in descending order if true */
  sort_desc?: boolean;
}

/**
 * Interface for paginated list of experiments in API responses.
 * Used for returning large sets of experiment data with pagination.
 */
export interface ExperimentListResponse {
  /** Array of experiment items */
  items: Experiment[];
  /** Total number of experiments matching the query */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
  /** Applied filters */
  filters: ExperimentFilter;
}

/**
 * Interface for list of experiment types in API responses.
 * Used for returning available experiment types.
 */
export interface ExperimentTypeListResponse {
  /** Array of experiment type items */
  items: ExperimentType[];
  /** Total number of experiment types */
  total: number;
}

/**
 * Interface for updating the status of an experiment.
 * Used in API requests to update experiment status.
 */
export interface ExperimentStatusUpdate {
  /** New status for the experiment */
  status: ExperimentStatus;
  /** Optional notes about the status change */
  notes?: string;
}

/**
 * Interface for experiment state in Redux store.
 * Represents the complete state of the experiment module in the application.
 */
export interface ExperimentState {
  /** Array of experiments */
  experiments: Experiment[];
  /** Array of experiment types */
  experimentTypes: ExperimentType[];
  /** Currently selected experiment */
  currentExperiment: ExperimentDetail | null;
  /** Loading state flag */
  loading: boolean;
  /** Error message if any */
  error: string | null;
  /** Total number of experiments */
  totalExperiments: number;
  /** Current page number */
  currentPage: number;
  /** Number of items per page */
  pageSize: number;
  /** Current filter criteria */
  filter: ExperimentFilter;
}