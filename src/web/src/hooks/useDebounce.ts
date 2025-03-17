import { useState, useEffect, useRef } from 'react'; // react version 18.2+

/**
 * A hook that returns a debounced version of the provided value,
 * which only updates after the specified delay has passed without new changes.
 * 
 * This hook is useful for scenarios where you want to delay processing until
 * a user has stopped changing a value, such as:
 * - Preventing excessive API calls during typing in search fields
 * - Delaying expensive calculations in real-time filtering
 * - Improving performance for components that update frequently
 * 
 * @template T The type of the value being debounced
 * @param value The value to be debounced
 * @param delay The delay in milliseconds before the value updates
 * @returns The debounced value that only updates after the specified delay
 */
function useDebounce<T>(value: T, delay: number): T {
  // Initialize state with the initial value
  const [debouncedValue, setDebouncedValue] = useState<T>(value);
  
  // Use a ref to store the timeout ID so we can clear it if needed
  const timeoutRef = useRef<NodeJS.Timeout | null>(null);
  
  useEffect(() => {
    // Clear any existing timeout to cancel pending updates
    if (timeoutRef.current) {
      clearTimeout(timeoutRef.current);
    }
    
    // Set a new timeout that will update the debounced value after the delay
    timeoutRef.current = setTimeout(() => {
      setDebouncedValue(value);
    }, delay);
    
    // Cleanup function to clear the timeout when the component unmounts
    // or when the value or delay changes
    return () => {
      if (timeoutRef.current) {
        clearTimeout(timeoutRef.current);
      }
    };
  }, [value, delay]); // Re-run effect when value or delay changes
  
  // Return the current debounced value
  return debouncedValue;
}

export default useDebounce;