import { Molecule, MoleculeProperty } from '../types/molecule';
import { SMILES_REGEX, PROPERTY_RANGES } from './validation';
import axios from 'axios'; // version ^1.4.0

// Constants
export const API_BASE_URL = process.env.REACT_APP_API_URL || '/api/v1';
export const DEFAULT_MOLECULE_PROPERTIES = ["molecular_weight", "logp", "h_bond_donors", "h_bond_acceptors"];
export const PROPERTY_DISPLAY_NAMES = {
  molecular_weight: "Molecular Weight",
  logp: "LogP",
  h_bond_donors: "H-Bond Donors",
  h_bond_acceptors: "H-Bond Acceptors",
  solubility: "Solubility",
  activity: "Activity"
};
export const PROPERTY_UNITS = {
  molecular_weight: "g/mol",
  logp: "",
  h_bond_donors: "",
  h_bond_acceptors: "",
  solubility: "mg/mL",
  activity: "%"
};

/**
 * Validates a SMILES string using regex pattern and basic structure checks
 * @param smiles SMILES string to validate
 * @returns True if the SMILES string is valid, false otherwise
 */
export function isValidSMILES(smiles: string): boolean {
  // Check if smiles is null, undefined, or empty string
  if (!smiles || smiles.trim() === '') {
    return false;
  }

  // Check if smiles matches SMILES_REGEX pattern
  if (!SMILES_REGEX.test(smiles)) {
    return false;
  }

  // Check for balanced parentheses and brackets
  const openParenCount = (smiles.match(/\(/g) || []).length;
  const closeParenCount = (smiles.match(/\)/g) || []).length;
  const openBracketCount = (smiles.match(/\[/g) || []).length;
  const closeBracketCount = (smiles.match(/\]/g) || []).length;

  if (openParenCount !== closeParenCount || openBracketCount !== closeBracketCount) {
    return false;
  }

  return true;
}

/**
 * Generates a URL for rendering a molecule image from a SMILES string
 * @param smiles SMILES string of the molecule
 * @param width Width of the image in pixels
 * @param height Height of the image in pixels
 * @param format Image format (svg or png)
 * @returns URL to the molecule image
 */
export function generateMoleculeImageUrl(
  smiles: string,
  width: number = 300,
  height: number = 200,
  format: string = 'svg'
): string {
  // Check if smiles is valid using isValidSMILES
  if (!isValidSMILES(smiles)) {
    return `/assets/images/molecule-placeholder.${format}`;
  }

  // Encode SMILES string for URL
  const encodedSmiles = encodeURIComponent(smiles);
  
  // Construct URL to backend API endpoint for molecule rendering
  return `${API_BASE_URL}/molecules/render?smiles=${encodedSmiles}&width=${width}&height=${height}&format=${format}`;
}

/**
 * Extracts a specific property value from a molecule object
 * @param molecule Molecule object
 * @param propertyName Name of the property to extract
 * @returns Property value or null if not found
 */
export function getMoleculePropertyByName(
  molecule: Molecule,
  propertyName: string
): any | null {
  // Check if molecule and molecule.properties exist
  if (!molecule || !molecule.properties) {
    return null;
  }

  // Find property in molecule.properties array where property_name matches propertyName
  const property = molecule.properties.find(
    prop => prop.property_name === propertyName
  );

  // If property found, return property_value
  return property ? property.property_value : null;
}

/**
 * Formats a molecular property value with appropriate precision and units
 * @param propertyName Name of the property
 * @param value Property value
 * @param includeUnits Whether to include units in formatted value
 * @returns Formatted property value with units if applicable
 */
export function formatPropertyValue(
  propertyName: string,
  value: any,
  includeUnits: boolean = true
): string {
  // Check if value is null or undefined
  if (value === null || value === undefined) {
    return 'N/A';
  }

  let formattedValue = '';
  
  // Determine appropriate precision based on property type
  if (typeof value === 'number') {
    const precision = propertyName === 'molecular_weight' ? 2 :
                      propertyName === 'logp' ? 2 :
                      propertyName === 'solubility' ? 2 :
                      propertyName === 'activity' ? 1 :
                      0;
    
    formattedValue = value.toFixed(precision);
  } else {
    formattedValue = String(value);
  }
  
  // If includeUnits is true, append unit from PROPERTY_UNITS if available
  if (includeUnits) {
    const unit = PROPERTY_UNITS[propertyName as keyof typeof PROPERTY_UNITS];
    if (unit) {
      formattedValue += ` ${unit}`;
    }
  }
  
  return formattedValue;
}

/**
 * Gets a user-friendly display name for a molecular property
 * @param propertyName Internal property name
 * @returns Display name for the property
 */
export function getPropertyDisplayName(propertyName: string): string {
  // Look up propertyName in PROPERTY_DISPLAY_NAMES
  const displayName = PROPERTY_DISPLAY_NAMES[propertyName as keyof typeof PROPERTY_DISPLAY_NAMES];
  
  if (displayName) {
    return displayName;
  }
  
  // Convert propertyName to title case as fallback
  return propertyName
    .split('_')
    .map(word => word.charAt(0).toUpperCase() + word.slice(1))
    .join(' ');
}

/**
 * Gets the unit for a molecular property
 * @param propertyName Internal property name
 * @returns Unit for the property or empty string if no unit
 */
export function getPropertyUnit(propertyName: string): string {
  return PROPERTY_UNITS[propertyName as keyof typeof PROPERTY_UNITS] || '';
}

/**
 * Sorts an array of molecules based on a specific property
 * @param molecules Array of molecules to sort
 * @param propertyName Property to sort by
 * @param ascending Whether to sort in ascending order
 * @returns Sorted array of molecules
 */
export function sortMoleculesByProperty(
  molecules: Molecule[],
  propertyName: string,
  ascending: boolean = true
): Molecule[] {
  // Create a copy of the molecules array to avoid mutating the original
  return [...molecules].sort((a, b) => {
    const valueA = getMoleculePropertyByName(a, propertyName);
    const valueB = getMoleculePropertyByName(b, propertyName);
    
    // Handle cases where property doesn't exist in some molecules
    if (valueA === null && valueB === null) return 0;
    if (valueA === null) return ascending ? 1 : -1;
    if (valueB === null) return ascending ? -1 : 1;
    
    // Apply ascending or descending sort based on parameter
    if (typeof valueA === 'number' && typeof valueB === 'number') {
      return ascending ? valueA - valueB : valueB - valueA;
    }
    
    const strA = String(valueA);
    const strB = String(valueB);
    return ascending ? strA.localeCompare(strB) : strB.localeCompare(strA);
  });
}

/**
 * Filters an array of molecules based on property value ranges
 * @param molecules Array of molecules to filter
 * @param propertyName Property to filter by
 * @param minValue Minimum value (inclusive)
 * @param maxValue Maximum value (inclusive)
 * @returns Filtered array of molecules
 */
export function filterMoleculesByProperty(
  molecules: Molecule[],
  propertyName: string,
  minValue: number,
  maxValue: number
): Molecule[] {
  return molecules.filter(molecule => {
    const value = getMoleculePropertyByName(molecule, propertyName);
    
    // Include molecules where property value is between minValue and maxValue (inclusive)
    if (value === null || value === undefined) {
      return false;
    }
    
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      return numValue >= minValue && numValue <= maxValue;
    }
    
    return false;
  });
}

/**
 * Truncates a SMILES string for display purposes
 * @param smiles SMILES string to truncate
 * @param maxLength Maximum length before truncation
 * @returns Truncated SMILES string with ellipsis if needed
 */
export function truncateSMILES(smiles: string, maxLength: number = 20): string {
  // Check if smiles is null, undefined, or empty string
  if (!smiles || smiles.trim() === '') {
    return '';
  }
  
  // If too long, truncate and add ellipsis
  if (smiles.length > maxLength) {
    return smiles.substring(0, maxLength - 3) + '...';
  }
  
  return smiles;
}

/**
 * Attempts to get a common name for a molecule, falling back to truncated SMILES
 * @param molecule Molecule object
 * @returns Common name or truncated SMILES
 */
export function getMoleculeCommonName(molecule: Molecule): string {
  // Check if molecule has a 'name' property
  const name = getMoleculePropertyByName(molecule, 'name');
  if (name) {
    return String(name);
  }
  
  // Otherwise, return truncated SMILES using truncateSMILES function
  return truncateSMILES(molecule.smiles);
}

/**
 * Calculates statistics for a specific property across an array of molecules
 * @param molecules Array of molecules
 * @param propertyName Property to calculate statistics for
 * @returns Statistics object with min, max, avg, and count
 */
export function calculateMoleculePropertyStatistics(
  molecules: Molecule[],
  propertyName: string
): { min: number; max: number; avg: number; count: number } {
  // Initialize statistics object with null values
  const stats = {
    min: Number.MAX_VALUE,
    max: Number.MIN_VALUE,
    avg: 0,
    count: 0
  };
  
  // Filter molecules to those that have the specified property
  const values = molecules
    .map(molecule => getMoleculePropertyByName(molecule, propertyName))
    .filter(value => value !== null)
    .map(value => typeof value === 'string' ? parseFloat(value) : value)
    .filter(value => typeof value === 'number' && !isNaN(value)) as number[];
  
  // If no molecules have the property, return statistics with zeros
  if (values.length > 0) {
    stats.min = Math.min(...values);
    stats.max = Math.max(...values);
    stats.avg = values.reduce((sum, val) => sum + val, 0) / values.length;
    stats.count = values.length;
  } else {
    stats.min = 0;
    stats.max = 0;
  }
  
  return stats;
}

/**
 * Extracts properties from a molecule into a simple key-value object
 * @param molecule Molecule object
 * @param propertyNames Array of property names to extract
 * @returns Object with property names as keys and values as values
 */
export function extractPropertiesFromMolecule(
  molecule: Molecule,
  propertyNames: string[] = DEFAULT_MOLECULE_PROPERTIES
): Record<string, any> {
  // Initialize empty result object
  const result: Record<string, any> = {};
  
  // For each property name in propertyNames:
  propertyNames.forEach(propName => {
    // Get property value using getMoleculePropertyByName
    result[propName] = getMoleculePropertyByName(molecule, propName);
  });
  
  return result;
}

/**
 * Compares two SMILES strings for chemical equivalence (not just string equality)
 * @param smiles1 First SMILES string
 * @param smiles2 Second SMILES string
 * @returns True if the molecules are chemically equivalent
 */
export function compareSmiles(smiles1: string, smiles2: string): boolean {
  // Validate both SMILES strings using isValidSMILES
  if (!isValidSMILES(smiles1) || !isValidSMILES(smiles2)) {
    return false;
  }
  
  // Normalize both SMILES strings (remove whitespace, standardize format)
  const normalized1 = smiles1.trim().replace(/\s+/g, '');
  const normalized2 = smiles2.trim().replace(/\s+/g, '');
  
  // Compare normalized strings
  return normalized1 === normalized2;
}