import { createTheme, ThemeOptions, Theme, PaletteMode } from '@mui/material/styles'; // v5.13+

// Internal imports
import palette, { createPalette } from './palette';
import typography, { responsiveFontSizes } from './typography';
import shadows, { customShadows } from './shadows';
import components, { createComponentOverrides } from './overrides';
import breakpoints from './responsive';
import { getLocalStorageItem, setLocalStorageItem } from '../utils/storage';

// Storage key for theme mode preference
const THEME_MODE_KEY = 'theme_mode';

/**
 * Retrieves the user's preferred theme mode from localStorage
 * @returns The stored theme mode ('light' or 'dark') or null if not found
 */
export function getStoredThemeMode(): PaletteMode | null {
  const storedMode = getLocalStorageItem<string>(THEME_MODE_KEY);
  if (storedMode === 'light' || storedMode === 'dark') {
    return storedMode as PaletteMode;
  }
  return null;
}

/**
 * Stores the user's preferred theme mode in localStorage
 * @param mode The theme mode to store ('light' or 'dark')
 * @returns True if storage was successful, false otherwise
 */
export function storeThemeMode(mode: PaletteMode): boolean {
  return setLocalStorageItem(THEME_MODE_KEY, mode);
}

/**
 * Creates a custom Material-UI theme with the specified mode
 * @param mode The theme mode ('light' or 'dark')
 * @returns Complete Material-UI theme object
 */
export function createCustomTheme(mode: PaletteMode): Theme {
  // Create theme palette based on mode
  const themePalette = createPalette(mode);
  
  // Apply responsive font sizes to typography
  const themeTypography = responsiveFontSizes(typography.typography);
  
  // Create component overrides based on mode
  const themeComponents = createComponentOverrides(mode);
  
  // Combine all theme elements into theme options
  const themeOptions: ThemeOptions = {
    palette: themePalette,
    typography: themeTypography,
    shadows: shadows,
    // Configure common shape settings
    shape: { 
      borderRadius: 8 
    },
    // Apply responsive breakpoints configuration
    breakpoints: breakpoints.breakpoints,
    // Apply component style overrides
    components: themeComponents,
    // Configure common transitions
    transitions: {
      duration: {
        shortest: 150,
        shorter: 200,
        short: 250,
        standard: 300,
        complex: 375,
        enteringScreen: 225,
        leavingScreen: 195,
      },
    },
    // Apply custom shadows
    customShadows,
  };
  
  // Create and return the final theme
  return createTheme(themeOptions);
}

// Create default theme (light mode)
const defaultTheme = createCustomTheme('light');

export default defaultTheme;