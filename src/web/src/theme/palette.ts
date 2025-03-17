import { alpha, PaletteOptions } from '@mui/material/styles'; // v5.13+

// PRIMARY COLOR PALETTE
// Used for main UI elements, buttons, and primary actions
const PRIMARY = {
  lighter: '#D1E9FC',
  light: '#76B0F1',
  main: '#1976D2', // Main primary color
  dark: '#103996',
  darker: '#061B64',
  contrastText: '#FFFFFF',
};

// SECONDARY COLOR PALETTE
// Used for secondary UI elements and actions
const SECONDARY = {
  lighter: '#E6F4F1',
  light: '#8ECDC4',
  main: '#00A699', // Main secondary color
  dark: '#007A6D',
  darker: '#004D45',
  contrastText: '#FFFFFF',
};

// INFO COLOR PALETTE
// Used for informational elements and status indicators
const INFO = {
  lighter: '#D0F2FF',
  light: '#74CAFF',
  main: '#0288D1', // Main info color
  dark: '#0157A0',
  darker: '#013A6B',
  contrastText: '#FFFFFF',
};

// SUCCESS COLOR PALETTE
// Used for positive status indicators and successful actions
const SUCCESS = {
  lighter: '#E9FCD4',
  light: '#AAF27F',
  main: '#2E7D32', // Main success color
  dark: '#1B5E20',
  darker: '#114114',
  contrastText: '#FFFFFF',
};

// WARNING COLOR PALETTE
// Used for cautionary status indicators and warnings
const WARNING = {
  lighter: '#FFF7CD',
  light: '#FFE16A',
  main: '#ED6C02', // Main warning color
  dark: '#B76E00',
  darker: '#7A4100',
  contrastText: '#212B36',
};

// ERROR COLOR PALETTE
// Used for negative status indicators and errors
const ERROR = {
  lighter: '#FFE7D9',
  light: '#FFA48D',
  main: '#D32F2F', // Main error color
  dark: '#B71C1C',
  darker: '#7A0C2E',
  contrastText: '#FFFFFF',
};

// GREY COLOR PALETTE
// Used for neutral UI elements, backgrounds, and text
const GREY = {
  0: '#FFFFFF', // White
  50: '#FAFAFA', // Very light grey
  100: '#F5F5F5', // Light grey
  200: '#EEEEEE', // Light-medium grey
  300: '#E0E0E0', // Medium grey
  400: '#BDBDBD', // Medium-dark grey
  500: '#9E9E9E', // Dark grey
  600: '#757575', // Darker grey
  700: '#616161', // Very dark grey
  800: '#424242', // Almost black grey
  900: '#212121', // Nearly black
  500_8: alpha('#9E9E9E', 0.08), // Dark grey with 8% opacity
  500_12: alpha('#9E9E9E', 0.12), // Dark grey with 12% opacity
  500_16: alpha('#9E9E9E', 0.16), // Dark grey with 16% opacity
  500_24: alpha('#9E9E9E', 0.24), // Dark grey with 24% opacity
  500_32: alpha('#9E9E9E', 0.32), // Dark grey with 32% opacity
  500_48: alpha('#9E9E9E', 0.48), // Dark grey with 48% opacity
  500_56: alpha('#9E9E9E', 0.56), // Dark grey with 56% opacity
  500_80: alpha('#9E9E9E', 0.8), // Dark grey with 80% opacity
};

// COMMON COLORS
// Basic colors used throughout the application
const COMMON = {
  white: '#FFFFFF', // Pure white
  black: '#000000', // Pure black
};

/**
 * Creates a linear gradient string from two colors
 * @param {string} color1 - First color
 * @param {string} color2 - Second color
 * @returns {string} CSS linear gradient value
 */
function createGradient(color1: string, color2: string): string {
  return `linear-gradient(to bottom, ${color1}, ${color2})`;
}

// GRADIENT DEFINITIONS
// Used for special UI elements that require gradient backgrounds
const GRADIENTS = {
  primary: createGradient(PRIMARY.light, PRIMARY.main),
  secondary: createGradient(SECONDARY.light, SECONDARY.main),
  info: createGradient(INFO.light, INFO.main),
  success: createGradient(SUCCESS.light, SUCCESS.main),
  warning: createGradient(WARNING.light, WARNING.main),
  error: createGradient(ERROR.light, ERROR.main),
};

/**
 * Creates a complete palette configuration with light and dark mode support
 * @param {string} mode - Theme mode ('light' or 'dark')
 * @returns {PaletteOptions} Complete palette configuration
 */
export function createPalette(mode: 'light' | 'dark'): PaletteOptions {
  // Determine if dark mode is active
  const isDark = mode === 'dark';

  // Configure text colors based on mode
  const textColors = {
    primary: isDark ? '#FFFFFF' : '#212B36',
    secondary: isDark ? GREY[500_80] : GREY[600],
    disabled: isDark ? GREY[500_80] : GREY[500],
  };

  // Configure background colors based on mode
  const backgroundColors = {
    default: isDark ? '#121212' : '#FFFFFF',
    paper: isDark ? '#1E1E1E' : '#FFFFFF',
    neutral: isDark ? GREY[900] : GREY[200],
  };

  // Configure component colors based on mode
  const componentColors = {
    border: isDark ? GREY[500_16] : GREY[500_32],
    divider: isDark ? GREY[500_16] : GREY[500_16],
  };

  // Return the complete palette configuration with WCAG 2.1 AA compliance considerations
  return {
    common: COMMON,
    primary: PRIMARY,
    secondary: SECONDARY,
    info: INFO,
    success: SUCCESS,
    warning: WARNING,
    error: ERROR,
    grey: GREY,
    gradients: GRADIENTS,
    text: {
      primary: textColors.primary,
      secondary: textColors.secondary,
      disabled: textColors.disabled,
    },
    background: {
      paper: backgroundColors.paper,
      default: backgroundColors.default,
      neutral: backgroundColors.neutral,
    },
    action: {
      active: isDark ? GREY[500] : GREY[600],
      hover: isDark ? GREY[500_8] : GREY[500_8],
      selected: isDark ? GREY[500_16] : GREY[500_16],
      disabled: isDark ? GREY[500_80] : GREY[500_32],
      disabledBackground: isDark ? GREY[500_24] : GREY[500_12],
      focus: isDark ? GREY[500_24] : GREY[500_16],
      hoverOpacity: 0.08,
      disabledOpacity: 0.48,
    },
    divider: componentColors.divider,
    mode,
  };
}

// Default palette (light mode)
const palette = createPalette('light');

export default { palette };