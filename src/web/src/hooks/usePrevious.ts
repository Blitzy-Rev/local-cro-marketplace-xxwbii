import { useRef, useEffect } from 'react'; // react 18.2+

/**
 * A custom React hook that stores and returns the previous value of a variable across renders.
 * This is useful for comparing current and previous values in effect hooks or for tracking
 * changes in props or state.
 * 
 * @template T - The type of the value being tracked
 * @param value - The value to track across renders
 * @returns The previous value of the provided variable (undefined on first render)
 * 
 * @example
 * ```tsx
 * function Counter() {
 *   const [count, setCount] = useState(0);
 *   const previousCount = usePrevious(count);
 *   
 *   // Now you can compare current count with previous count
 *   useEffect(() => {
 *     if (previousCount !== undefined && count > previousCount) {
 *       console.log('Count increased');
 *     }
 *   }, [count, previousCount]);
 *   
 *   return (
 *     <div>
 *       <p>Current: {count}, Previous: {previousCount ?? 'none'}</p>
 *       <button onClick={() => setCount(count + 1)}>Increment</button>
 *     </div>
 *   );
 * }
 * ```
 */
function usePrevious<T>(value: T): T | undefined {
  // Create a ref using useRef to store the previous value, initially undefined
  const ref = useRef<T | undefined>(undefined);
  
  // Set up an effect that runs after every render
  useEffect(() => {
    // In the effect, update the ref's current value to the current value passed to the hook
    ref.current = value;
  }, [value]); // Only re-run if the value changes
  
  // Return the ref's current value, which is the previous value from the last render
  return ref.current;
}

export default usePrevious;