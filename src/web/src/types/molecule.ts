import { UUID } from 'uuid';

/**
 * Enum for molecule flag status values
 * Used for prioritizing and categorizing molecules
 */
export enum FlagStatus {
  IMPORTANT = 'important',
  FAVORITE = 'favorite',
  REVIEW = 'review',
  ARCHIVED = 'archived'
}

/**
 * Enum for property filter operators
 * Used in the filtering system to allow complex property-based queries
 */
export enum FilterOperator {
  EQUALS = 'equals',
  NOT_EQUALS = 'not_equals',
  GREATER_THAN = 'greater_than',
  LESS_THAN = 'less_than',
  BETWEEN = 'between',
  CONTAINS = 'contains'
}

/**
 * Interface for molecular property data
 * Represents individual properties of a molecule with their values and metadata
 */
export interface MoleculeProperty {
  id: UUID;
  molecule_id: UUID;
  property_name: string;
  property_value: number | string | boolean;
  property_unit?: string;
  is_calculated: boolean;
}

/**
 * Interface for molecule data with all properties
 * Core data structure representing a chemical molecule in the system
 */
export interface Molecule {
  id: UUID;
  smiles: string;
  flag_status?: FlagStatus;
  created_by: UUID;
  created_at: string;
  updated_at?: string;
  properties: MoleculeProperty[];
  creator?: {
    id: UUID;
    email: string;
    username?: string;
  };
  image_url?: string;
}

/**
 * Interface for creating a new molecule
 * Contains the minimum required fields to create a molecule
 */
export interface MoleculeCreate {
  smiles: string;
  flag_status?: FlagStatus;
  properties?: Partial<MoleculeProperty>[];
}

/**
 * Interface for updating an existing molecule
 * Contains fields that can be modified after creation
 */
export interface MoleculeUpdate {
  flag_status?: FlagStatus;
  properties?: Partial<MoleculeProperty>[];
}

/**
 * Interface for filtering molecules by property values
 * Enables complex filtering based on molecular properties
 */
export interface PropertyFilter {
  name: string;
  value?: string | number | boolean;
  min_value?: number;
  max_value?: number;
  operator: FilterOperator;
}

/**
 * Interface for filtering molecules by various criteria
 * Comprehensive filtering options for molecule searches
 */
export interface MoleculeFilter {
  smiles_pattern?: string;
  flag_status?: FlagStatus;
  created_by?: UUID;
  property_filters?: PropertyFilter[];
  created_after?: string;
  created_before?: string;
  library_ids?: UUID[];
  experiment_ids?: UUID[];
  sort_by?: string;
  sort_desc?: boolean;
}

/**
 * Interface for paginated list of molecules in API responses
 * Used for returning large sets of molecule data with pagination
 */
export interface MoleculeListResponse {
  items: Molecule[];
  total: number;
  page: number;
  size: number;
}

/**
 * Interface for detailed molecule data with related entities
 * Provides comprehensive information about a molecule and its relationships
 */
export interface MoleculeDetailResponse {
  molecule: Molecule;
  libraries: {
    id: UUID;
    name: string;
  }[];
  experiments: {
    id: UUID;
    name: string;
    status: string;
    date?: string;
  }[];
}

/**
 * Interface for similarity search parameters
 * Used for finding molecules with similar structures
 */
export interface MoleculeSimilaritySearchParams {
  smiles: string;
  threshold?: number;
  limit?: number;
}

/**
 * Interface for substructure search parameters
 * Used for finding molecules containing a specific substructure
 */
export interface MoleculeSubstructureSearchParams {
  smarts: string;
  limit?: number;
}

/**
 * Interface for molecule state in Redux store
 * Represents the complete state of the molecule module in the application
 */
export interface MoleculeState {
  molecules: Molecule[];
  filter: MoleculeFilter;
  selectedMolecules: UUID[];
  loading: boolean;
  error: string | null;
  totalMolecules: number;
  currentPage: number;
  pageSize: number;
}