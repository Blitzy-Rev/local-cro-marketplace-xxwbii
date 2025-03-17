import { UUID } from 'uuid';
import { Molecule } from './molecule';
import { User } from './user';

/**
 * Enum for library molecule operations
 * Defines operations that can be performed on molecules within a library
 */
export enum LibraryOperationType {
  /** Add molecules to a library */
  ADD = 'add',
  /** Remove molecules from a library */
  REMOVE = 'remove'
}

/**
 * Interface for library data with basic information
 * Represents a collection of molecules created by a user
 */
export interface Library {
  /** Unique identifier for the library */
  id: UUID;
  /** Name of the library */
  name: string;
  /** Optional description of the library's purpose or contents */
  description?: string;
  /** ID of the user who created this library */
  created_by: number; // Using number to match User.id type
  /** Date and time when the library was created */
  created_at: string;
  /** Date and time when the library was last updated */
  updated_at?: string;
  /** Number of molecules in this library */
  molecule_count: number;
  /** User who created this library */
  creator?: {
    id: number;
    email: string;
  };
}

/**
 * Interface for detailed library data including contained molecules
 * Extends the basic Library interface with molecule information
 */
export interface LibraryDetail extends Library {
  /** Molecules contained in this library */
  molecules: Molecule[];
}

/**
 * Interface for creating a new library
 * Contains the minimum required fields to create a library
 */
export interface LibraryCreate {
  /** Name of the library */
  name: string;
  /** Optional description of the library's purpose or contents */
  description?: string;
  /** Optional list of molecule IDs to initially include in the library */
  molecule_ids?: UUID[];
}

/**
 * Interface for updating an existing library
 * Contains fields that can be modified after creation
 */
export interface LibraryUpdate {
  /** Updated name of the library */
  name?: string;
  /** Updated description of the library */
  description?: string;
}

/**
 * Interface for filtering libraries by various criteria
 * Used for searching and organizing libraries
 */
export interface LibraryFilter {
  /** Filter by library name (partial match) */
  name?: string;
  /** Filter by creator user ID */
  created_by?: number; // Using number to match User.id type
  /** Filter by libraries created after this date */
  created_after?: string;
  /** Filter by libraries created before this date */
  created_before?: string;
  /** Filter by libraries containing a specific molecule */
  contains_molecule_id?: UUID;
  /** Page number for pagination */
  page?: number;
  /** Number of items per page */
  page_size?: number;
  /** Field to sort results by */
  sort_by?: string;
  /** Whether to sort in descending order */
  sort_desc?: boolean;
}

/**
 * Interface for paginated list of libraries in API responses
 * Used for returning large sets of library data with pagination
 */
export interface LibraryListResponse {
  /** Array of library items */
  items: Library[];
  /** Total number of libraries matching the filter */
  total: number;
  /** Current page number */
  page: number;
  /** Number of items per page */
  size: number;
}

/**
 * Interface for adding or removing molecules from a library
 * Used for batch operations on library contents
 */
export interface LibraryMoleculeOperation {
  /** Type of operation to perform (add or remove) */
  operation: LibraryOperationType;
  /** List of molecule IDs to operate on */
  molecule_ids: UUID[];
}

/**
 * Interface for library state in Redux store
 * Represents the complete state of the library module in the application
 */
export interface LibraryState {
  /** List of libraries */
  libraries: Library[];
  /** Currently selected library detail */
  currentLibrary: LibraryDetail | null;
  /** Whether library data is being loaded */
  loading: boolean;
  /** Error message if loading failed */
  error: string | null;
  /** Total number of libraries matching the current filter */
  totalLibraries: number;
  /** Current page number */
  currentPage: number;
  /** Number of items per page */
  pageSize: number;
}

/**
 * Type for supported library export formats
 * Defines the available file formats for exporting library data
 */
export type LibraryExportFormat = 'csv' | 'sdf';