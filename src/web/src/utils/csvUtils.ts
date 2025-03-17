import Papa from 'papaparse'; // v5.4.0
import { validateCSVHeaders, validateCSVFileSize, REQUIRED_CSV_COLUMNS, SMILES_REGEX } from '../utils/validation';
import { MoleculeProperty } from '../types/molecule';

// Maximum number of rows to show in CSV preview
const MAX_PREVIEW_ROWS = 10;

// Prefix for custom properties
const CUSTOM_PROPERTY_PREFIX = 'custom:';

/**
 * Validates a CSV file for format, size, and required columns
 * @param file - The CSV file to validate
 * @returns Object containing validation result and error message if invalid
 */
export async function validateCSVFile(file: File): Promise<{isValid: boolean, error: string | null}> {
  // Check if file is null or undefined
  if (!file) {
    return { isValid: false, error: 'No file provided' };
  }
  
  // Validate file type
  if (file.type !== 'text/csv') {
    return { isValid: false, error: 'Invalid file format. Only CSV files are accepted.' };
  }
  
  // Validate file size
  const sizeValidation = validateCSVFileSize(file.size);
  if (!sizeValidation.isValid) {
    return { isValid: false, error: sizeValidation.error };
  }
  
  // Parse CSV headers
  return new Promise<{isValid: boolean, error: string | null}>((resolve) => {
    Papa.parse(file, {
      preview: 1, // Only parse the first row to get headers
      header: false,
      skipEmptyLines: true,
      complete: (results) => {
        if (results.data.length === 0) {
          resolve({ isValid: false, error: 'CSV file is empty' });
          return;
        }
        
        const headers = results.data[0] as string[];
        const headersValidation = validateCSVHeaders(headers);
        
        if (!headersValidation.isValid) {
          resolve({ 
            isValid: false, 
            error: `CSV file is missing required columns: ${headersValidation.missingColumns.join(', ')}` 
          });
          return;
        }
        
        resolve({ isValid: true, error: null });
      },
      error: (error) => {
        resolve({ isValid: false, error: `Error parsing CSV file: ${error.message}` });
      }
    });
  });
}

/**
 * Parses a CSV file and returns a preview of the data
 * @param file - The CSV file to parse
 * @param maxRows - Maximum number of rows to include in preview
 * @returns Object containing headers and preview data rows
 */
export async function getCSVPreview(
  file: File, 
  maxRows: number = MAX_PREVIEW_ROWS
): Promise<{headers: string[], data: any[]}> {
  return new Promise<{headers: string[], data: any[]}>((resolve, reject) => {
    Papa.parse(file, {
      preview: maxRows + 1, // +1 for header row
      header: true,
      skipEmptyLines: true,
      complete: (results) => {
        resolve({
          headers: results.meta.fields || [],
          data: results.data
        });
      },
      error: (error) => {
        reject(new Error(`Error parsing CSV file: ${error.message}`));
      }
    });
  });
}

/**
 * Creates mapping options for CSV columns to system properties
 * @param headers - CSV column headers
 * @param availableProperties - Available system properties
 * @returns Array of mapping options for each CSV column
 */
export function createColumnMappingOptions(
  headers: string[], 
  availableProperties: {property_name: string, description?: string}[]
): {
  csvColumn: string;
  mappingOptions: {value: string, label: string}[];
}[] {
  return headers.map(header => {
    // Create system property options
    const systemOptions = availableProperties.map(prop => ({
      value: prop.property_name,
      label: `${prop.property_name}${prop.description ? ` (${prop.description})` : ''}`
    }));
    
    // Add custom property option
    const customOption = {
      value: `${CUSTOM_PROPERTY_PREFIX}${header}`,
      label: `Custom: ${header}`
    };
    
    return {
      csvColumn: header,
      mappingOptions: [
        ...systemOptions,
        customOption
      ]
    };
  });
}

/**
 * Automatically generates mapping configuration based on header name similarity
 * @param headers - CSV column headers
 * @param availableProperties - Available system properties
 * @returns Array of column mappings with best-guess system properties
 */
export function generateMappingFromHeaders(
  headers: string[],
  availableProperties: {property_name: string, description?: string}[]
): {csvColumn: string, systemProperty: string}[] {
  const mapping: {csvColumn: string, systemProperty: string}[] = [];
  
  // Special handling for SMILES column - it's required
  const smilesHeaderIndex = headers.findIndex(
    h => h.toUpperCase() === 'SMILES' || h.toUpperCase() === 'STRUCTURE'
  );
  
  if (smilesHeaderIndex >= 0) {
    const smilesProperty = availableProperties.find(p => 
      p.property_name.toUpperCase() === 'SMILES'
    );
    
    if (smilesProperty) {
      mapping.push({
        csvColumn: headers[smilesHeaderIndex],
        systemProperty: smilesProperty.property_name
      });
    }
  }
  
  // Map remaining headers based on similarity
  headers.forEach(header => {
    // Skip if already mapped (like SMILES)
    if (mapping.some(m => m.csvColumn === header)) {
      return;
    }
    
    const similarProperty = findSimilarPropertyName(header, availableProperties);
    
    if (similarProperty && similarProperty.score > 0.7) {
      mapping.push({
        csvColumn: header,
        systemProperty: similarProperty.property.property_name
      });
    } else {
      // If no good match found, map to custom property
      mapping.push({
        csvColumn: header,
        systemProperty: `${CUSTOM_PROPERTY_PREFIX}${header}`
      });
    }
  });
  
  return mapping;
}

/**
 * Validates a column mapping configuration for required properties
 * @param mapping - Column mapping configuration
 * @returns Validation result with errors if invalid
 */
export function validateMappingConfiguration(
  mapping: {csvColumn: string, systemProperty: string}[]
): {
  isValid: boolean;
  errors: {column: string, message: string}[];
} {
  const errors: {column: string, message: string}[] = [];
  const mappedSystemProperties = new Set<string>();
  
  // Check if SMILES column is mapped
  const smilesMapping = mapping.find(m => 
    m.systemProperty.toUpperCase() === 'SMILES'
  );
  
  if (!smilesMapping) {
    errors.push({
      column: 'mapping',
      message: 'SMILES column must be mapped'
    });
  }
  
  // Check for duplicate system property mappings
  mapping.forEach(map => {
    // Skip custom properties
    if (map.systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX)) {
      // Validate custom property name
      const customName = map.systemProperty.substring(CUSTOM_PROPERTY_PREFIX.length);
      if (!customName || customName.trim() === '') {
        errors.push({
          column: map.csvColumn,
          message: 'Custom property name cannot be empty'
        });
      }
      return;
    }
    
    // Check for duplicates of system properties
    if (mappedSystemProperties.has(map.systemProperty)) {
      errors.push({
        column: map.csvColumn,
        message: `Duplicate mapping for system property "${map.systemProperty}"`
      });
    } else {
      mappedSystemProperties.add(map.systemProperty);
    }
  });
  
  return {
    isValid: errors.length === 0,
    errors
  };
}

/**
 * Parses a CSV file and transforms data according to mapping configuration
 * @param file - The CSV file to parse
 * @param mapping - Column mapping configuration
 * @returns Parsed molecules with properties according to mapping
 */
export async function parseCSVData(
  file: File,
  mapping: {csvColumn: string, systemProperty: string}[]
): Promise<{
  molecules: Array<{
    smiles: string;
    properties: MoleculeProperty[];
  }>;
}> {
  return new Promise<{
    molecules: Array<{
      smiles: string;
      properties: MoleculeProperty[];
    }>;
  }>((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      skipEmptyLines: true,
      complete: (results) => {
        try {
          const smilesMapping = mapping.find(m => 
            m.systemProperty.toUpperCase() === 'SMILES'
          );
          
          if (!smilesMapping) {
            throw new Error('SMILES column mapping is required');
          }
          
          const molecules = results.data
            .map((row: any) => {
              const smiles = row[smilesMapping.csvColumn];
              
              // Skip rows with invalid SMILES
              if (!isValidSMILES(smiles)) {
                return null;
              }
              
              // Create molecule properties based on mapping
              const properties: MoleculeProperty[] = mapping
                .filter(m => m.systemProperty !== 'SMILES') // Skip SMILES as it's the main identifier
                .map(m => {
                  const value = row[m.csvColumn];
                  const isCustom = m.systemProperty.startsWith(CUSTOM_PROPERTY_PREFIX);
                  const propertyName = isCustom 
                    ? m.systemProperty.substring(CUSTOM_PROPERTY_PREFIX.length)
                    : m.systemProperty;
                  
                  return {
                    id: '', // Will be assigned by the backend
                    molecule_id: '', // Will be assigned by the backend
                    property_name: propertyName,
                    property_value: value,
                    is_calculated: false // CSV imported properties are not calculated
                  };
                });
              
              return {
                smiles,
                properties
              };
            })
            .filter(Boolean); // Remove null entries (invalid SMILES)
          
          resolve({ molecules });
        } catch (error) {
          reject(error);
        }
      },
      error: (error) => {
        reject(new Error(`Error parsing CSV file: ${error.message}`));
      }
    });
  });
}

/**
 * Finds the most similar property name from available properties
 * @param header - CSV column header
 * @param availableProperties - Available system properties
 * @returns Most similar property with similarity score, or null if no good match
 */
export function findSimilarPropertyName(
  header: string,
  availableProperties: {property_name: string, description?: string}[]
): {property: any, score: number} | null {
  // Normalize header for comparison
  const normalizedHeader = header.toLowerCase().replace(/[_\s-]/g, '');
  
  let bestMatch = null;
  let bestScore = 0;
  
  availableProperties.forEach(property => {
    const score = calculateStringSimilarity(
      normalizedHeader,
      property.property_name.toLowerCase().replace(/[_\s-]/g, '')
    );
    
    if (score > bestScore) {
      bestScore = score;
      bestMatch = property;
    }
  });
  
  if (bestMatch) {
    return { property: bestMatch, score: bestScore };
  }
  
  return null;
}

/**
 * Calculates similarity score between two strings
 * @param str1 - First string
 * @param str2 - Second string
 * @returns Similarity score between 0 and 1
 */
export function calculateStringSimilarity(str1: string, str2: string): number {
  // Normalize strings
  const s1 = str1.toLowerCase().replace(/[_\s-]/g, '');
  const s2 = str2.toLowerCase().replace(/[_\s-]/g, '');
  
  // Simple matching for exact or contained strings
  if (s1 === s2) return 1;
  if (s1.includes(s2) || s2.includes(s1)) return 0.9;
  
  // Calculate Levenshtein distance for string similarity
  const m = s1.length;
  const n = s2.length;
  
  // If one string is empty, distance is the length of the other
  if (m === 0) return 0;
  if (n === 0) return 0;
  
  // Create two work vectors of integer distances
  let v0 = Array(n + 1).fill(0);
  let v1 = Array(n + 1).fill(0);
  
  // Initialize v0 (edit distances for empty string)
  for (let i = 0; i <= n; i++) {
    v0[i] = i;
  }
  
  // Calculate edit distance
  for (let i = 0; i < m; i++) {
    // First element of v1 is A[i+1][0]
    v1[0] = i + 1;
    
    // Calculate row for v1
    for (let j = 0; j < n; j++) {
      const cost = s1[i] === s2[j] ? 0 : 1;
      v1[j + 1] = Math.min(
        v1[j] + 1,
        v0[j + 1] + 1,
        v0[j] + cost
      );
    }
    
    // Copy current row to previous row for next iteration
    for (let j = 0; j <= n; j++) {
      v0[j] = v1[j];
    }
  }
  
  const distance = v1[n];
  const maxLength = Math.max(m, n);
  
  // Convert distance to similarity score (0-1 range)
  return 1 - distance / maxLength;
}

/**
 * Validates a SMILES string using regex pattern
 * @param smiles - SMILES string to validate
 * @returns Whether the SMILES string is valid
 */
export function isValidSMILES(smiles: string): boolean {
  if (!smiles || typeof smiles !== 'string' || smiles.trim() === '') {
    return false;
  }
  
  return SMILES_REGEX.test(smiles);
}

export { CUSTOM_PROPERTY_PREFIX };