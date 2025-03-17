/**
 * Utility module that provides a consistent interface for browser storage operations
 * in the Molecular Data Management and CRO Integration Platform frontend.
 * It includes functions for safely storing, retrieving, and removing data from
 * localStorage with proper error handling and type safety.
 */

// Prefix for all storage keys to avoid collisions with other applications
const STORAGE_PREFIX = "mdm_";

/**
 * Safely stores data in localStorage with proper error handling and prefixing.
 * This function handles stringification of the data and applies a consistent
 * prefix to the storage key.
 * 
 * @param key The key to store the data under (without prefix)
 * @param value The value to store (will be JSON stringified)
 * @returns True if storage was successful, false otherwise
 */
export const setLocalStorageItem = (key: string, value: any): boolean => {
  const prefixedKey = `${STORAGE_PREFIX}${key}`;
  try {
    const serializedValue = JSON.stringify(value);
    localStorage.setItem(prefixedKey, serializedValue);
    return true;
  } catch (error) {
    console.error(`Error storing item in localStorage: ${error}`);
    return false;
  }
};

/**
 * Safely retrieves data from localStorage with proper error handling and type casting.
 * This function handles parsing of the JSON data and returns a typed result.
 * 
 * @template T The expected type of the stored value
 * @param key The key to retrieve data from (without prefix)
 * @returns The retrieved value with the specified type or null if not found or on error
 */
export const getLocalStorageItem = <T>(key: string): T | null => {
  const prefixedKey = `${STORAGE_PREFIX}${key}`;
  try {
    const item = localStorage.getItem(prefixedKey);
    if (item === null || item === undefined) {
      return null;
    }
    
    return JSON.parse(item) as T;
  } catch (error) {
    console.error(`Error retrieving or parsing item from localStorage: ${error}`);
    return null;
  }
};

/**
 * Safely removes data from localStorage with proper error handling.
 * 
 * @param key The key to remove (without prefix)
 * @returns True if removal was successful, false otherwise
 */
export const removeLocalStorageItem = (key: string): boolean => {
  const prefixedKey = `${STORAGE_PREFIX}${key}`;
  try {
    localStorage.removeItem(prefixedKey);
    return true;
  } catch (error) {
    console.error(`Error removing item from localStorage: ${error}`);
    return false;
  }
};

/**
 * Clears all items from localStorage that match the application prefix.
 * This is useful for logout operations or clearing application state.
 * 
 * @returns True if clearing was successful, false otherwise
 */
export const clearLocalStorage = (): boolean => {
  try {
    // First collect all keys to avoid issues with removing during iteration
    const keysToRemove: string[] = [];
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && key.startsWith(STORAGE_PREFIX)) {
        keysToRemove.push(key);
      }
    }
    
    // Then remove each key
    keysToRemove.forEach(key => localStorage.removeItem(key));
    return true;
  } catch (error) {
    console.error(`Error clearing localStorage: ${error}`);
    return false;
  }
};

/**
 * Gets all keys in localStorage that match the application prefix.
 * Returns the keys without the prefix for easier use.
 * 
 * @returns Array of keys without the prefix
 */
export const getStorageKeys = (): string[] => {
  const keys: string[] = [];
  try {
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && key.startsWith(STORAGE_PREFIX)) {
        keys.push(key.replace(STORAGE_PREFIX, ""));
      }
    }
    return keys;
  } catch (error) {
    console.error(`Error getting localStorage keys: ${error}`);
    return [];
  }
};

/**
 * Checks if localStorage is available in the current browser environment.
 * This can be used to determine if the storage functionality can be used.
 * 
 * @returns True if localStorage is available, false otherwise
 */
export const isStorageAvailable = (): boolean => {
  try {
    const testKey = `${STORAGE_PREFIX}test`;
    localStorage.setItem(testKey, "test");
    localStorage.removeItem(testKey);
    return true;
  } catch (error) {
    return false;
  }
};

/**
 * Calculates the approximate storage usage in bytes for all keys
 * that match the application prefix.
 * 
 * @returns Approximate storage usage in bytes
 */
export const getStorageUsage = (): number => {
  let totalSize = 0;
  try {
    for (let i = 0; i < localStorage.length; i++) {
      const key = localStorage.key(i);
      if (key && key.startsWith(STORAGE_PREFIX)) {
        const value = localStorage.getItem(key) || "";
        // Each character in a string is 2 bytes in JavaScript
        totalSize += (key.length + value.length) * 2;
      }
    }
    return totalSize;
  } catch (error) {
    console.error(`Error calculating storage usage: ${error}`);
    return 0;
  }
};