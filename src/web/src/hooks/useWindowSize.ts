import { useState, useEffect, useCallback } from 'react'; // react version 18.2+
import { MOBILE_BREAKPOINT, TABLET_MIN_BREAKPOINT, TABLET_MAX_BREAKPOINT, DESKTOP_BREAKPOINT } from '../theme/responsive';

/**
 * Interface for window size information
 */
interface WindowSize {
  /** Current window width in pixels */
  width: number | undefined;
  /** Current window height in pixels */
  height: number | undefined;
  /** Whether the current viewport is mobile size (< 768px) */
  isMobile: boolean;
  /** Whether the current viewport is tablet size (768px-1199px) */
  isTablet: boolean;
  /** Whether the current viewport is desktop size (≥ 1200px) */
  isDesktop: boolean;
}

/**
 * A custom React hook that tracks and returns the current window dimensions and device type.
 * This hook enables responsive design by allowing components to adapt their layout and behavior
 * based on the current viewport size.
 * 
 * @returns An object containing width, height, and device type (mobile, tablet, desktop)
 */
function useWindowSize(): WindowSize {
  // Initialize with undefined to handle server-side rendering
  const [windowSize, setWindowSize] = useState<WindowSize>({
    width: undefined,
    height: undefined,
    isMobile: false,
    isTablet: false,
    isDesktop: false,
  });
  
  // Memoize the resize handler to prevent unnecessary re-creation
  const handleResize = useCallback(() => {
    const width = window.innerWidth;
    const height = window.innerHeight;
    
    // Determine device type based on breakpoint constants
    const isMobile = width <= MOBILE_BREAKPOINT;
    const isTablet = width >= TABLET_MIN_BREAKPOINT && width <= TABLET_MAX_BREAKPOINT;
    const isDesktop = width >= DESKTOP_BREAKPOINT;
    
    setWindowSize({
      width,
      height,
      isMobile,
      isTablet,
      isDesktop,
    });
  }, []);
  
  useEffect(() => {
    // Only run on client-side to prevent SSR mismatch
    if (typeof window !== 'undefined') {
      // Initialize with current dimensions
      handleResize();
      
      // Add event listener for window resize
      window.addEventListener('resize', handleResize);
      
      // Clean up event listener on component unmount to prevent memory leaks
      return () => {
        window.removeEventListener('resize', handleResize);
      };
    }
    
    return undefined;
  }, [handleResize]);
  
  return windowSize;
}

export default useWindowSize;