import { Breakpoints } from '@mui/material/styles'; // @mui/material version 5.13+

/**
 * Maximum width for mobile devices in pixels
 * Used to identify mobile viewport sizes for responsive design
 */
export const MOBILE_BREAKPOINT = 767;

/**
 * Minimum width for tablet devices in pixels
 * Used as the starting point for tablet-specific layouts
 */
export const TABLET_MIN_BREAKPOINT = 768;

/**
 * Maximum width for tablet devices in pixels
 * Used as the ending point for tablet-specific layouts
 */
export const TABLET_MAX_BREAKPOINT = 1199;

/**
 * Minimum width for desktop devices in pixels
 * Used as the starting point for full desktop layouts
 */
export const DESKTOP_BREAKPOINT = 1200;

/**
 * Creates a breakpoints configuration object for Material-UI theme
 * 
 * Configures the responsive breakpoints used throughout the application
 * to ensure consistent responsive behavior across all components
 * 
 * @returns {Breakpoints} Configured breakpoints object for responsive design
 */
function createBreakpoints(): Breakpoints {
  return {
    values: {
      xs: 0,                    // Extra small screens (mobile phones)
      sm: TABLET_MIN_BREAKPOINT, // Small screens (tablets)
      md: 1024,                 // Medium screens (small laptops)
      lg: DESKTOP_BREAKPOINT,   // Large screens (desktops)
      xl: 1536,                 // Extra large screens (large desktops)
    },
    unit: 'px',                 // Unit for breakpoint values
    step: 5,                    // Step size for breakpoint calculations
  } as Breakpoints;
}

/**
 * Breakpoints configuration for Material-UI theme
 * Contains the configuration for all responsive breakpoints used in the application theme
 */
const breakpoints = createBreakpoints();

export default { breakpoints };