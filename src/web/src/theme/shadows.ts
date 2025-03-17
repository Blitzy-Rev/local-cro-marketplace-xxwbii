import { alpha } from '@mui/material/styles'; // v5.13+
import { GREY } from './palette';

// Base shadow colors
const LIGHT_MODE_SHADOW_COLOR = GREY[500];
const DARK_MODE_SHADOW_COLOR = '#000000';

// Transparent shadow variations with different opacity levels
const transparent1 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.01);
const transparent2 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.02);
const transparent4 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.04);
const transparent8 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.08);
const transparent12 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.12);
const transparent16 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.16);
const transparent24 = alpha(LIGHT_MODE_SHADOW_COLOR, 0.24);

/**
 * Creates a shadow string with specified color and opacity values
 * @param {string} color - Base color for the shadow
 * @param {number} opacity - Opacity level for the shadow
 * @returns {string} Formatted CSS shadow value
 */
export function createShadow(color: string, opacity: number): string {
  const transparentColor = alpha(color, opacity);
  return `0px 8px 16px 0px ${transparentColor}`;
}

// Array of 25 shadow definitions for Material-UI elevation levels (0-24)
const shadows = [
  'none',
  `0px 2px 1px -1px ${transparent2}, 0px 1px 1px 0px ${transparent4}, 0px 1px 3px 0px ${transparent8}`,
  `0px 3px 1px -2px ${transparent2}, 0px 2px 2px 0px ${transparent4}, 0px 1px 5px 0px ${transparent8}`,
  `0px 3px 3px -2px ${transparent2}, 0px 3px 4px 0px ${transparent4}, 0px 1px 8px 0px ${transparent8}`,
  `0px 2px 4px -1px ${transparent2}, 0px 4px 5px 0px ${transparent4}, 0px 1px 10px 0px ${transparent8}`,
  `0px 3px 5px -1px ${transparent2}, 0px 5px 8px 0px ${transparent4}, 0px 1px 14px 0px ${transparent8}`,
  `0px 3px 5px -1px ${transparent2}, 0px 6px 10px 0px ${transparent4}, 0px 1px 18px 0px ${transparent8}`,
  `0px 4px 5px -2px ${transparent2}, 0px 7px 10px 1px ${transparent4}, 0px 2px 16px 1px ${transparent8}`,
  `0px 5px 5px -3px ${transparent2}, 0px 8px 10px 1px ${transparent4}, 0px 3px 14px 2px ${transparent8}`,
  `0px 5px 6px -3px ${transparent2}, 0px 9px 12px 1px ${transparent4}, 0px 3px 16px 2px ${transparent8}`,
  `0px 6px 6px -3px ${transparent2}, 0px 10px 14px 1px ${transparent4}, 0px 4px 18px 3px ${transparent8}`,
  `0px 6px 7px -4px ${transparent2}, 0px 11px 15px 1px ${transparent4}, 0px 4px 20px 3px ${transparent8}`,
  `0px 7px 8px -4px ${transparent2}, 0px 12px 17px 2px ${transparent4}, 0px 5px 22px 4px ${transparent8}`,
  `0px 7px 8px -4px ${transparent2}, 0px 13px 19px 2px ${transparent4}, 0px 5px 24px 4px ${transparent8}`,
  `0px 7px 9px -4px ${transparent2}, 0px 14px 21px 2px ${transparent4}, 0px 5px 26px 4px ${transparent8}`,
  `0px 8px 9px -5px ${transparent2}, 0px 15px 22px 2px ${transparent4}, 0px 6px 28px 5px ${transparent8}`,
  `0px 8px 10px -5px ${transparent2}, 0px 16px 24px 2px ${transparent4}, 0px 6px 30px 5px ${transparent8}`,
  `0px 8px 11px -5px ${transparent2}, 0px 17px 26px 2px ${transparent4}, 0px 6px 32px 5px ${transparent8}`,
  `0px 9px 11px -5px ${transparent2}, 0px 18px 28px 2px ${transparent4}, 0px 7px 34px 6px ${transparent8}`,
  `0px 9px 12px -6px ${transparent2}, 0px 19px 29px 2px ${transparent4}, 0px 7px 36px 6px ${transparent8}`,
  `0px 10px 13px -6px ${transparent2}, 0px 20px 31px 3px ${transparent4}, 0px 8px 38px 7px ${transparent8}`,
  `0px 10px 13px -6px ${transparent2}, 0px 21px 33px 3px ${transparent4}, 0px 8px 40px 7px ${transparent8}`,
  `0px 10px 14px -6px ${transparent2}, 0px 22px 35px 3px ${transparent4}, 0px 8px 42px 7px ${transparent8}`,
  `0px 11px 14px -7px ${transparent2}, 0px 23px 36px 3px ${transparent4}, 0px 9px 44px 8px ${transparent8}`,
  `0px 11px 15px -7px ${transparent2}, 0px 24px 38px 3px ${transparent4}, 0px 9px 46px 8px ${transparent8}`,
];

// Custom shadow definitions for specific components and states
export const customShadows = {
  // Z-axis shadows for different elevation levels
  z1: `0px 1px 2px 0px ${transparent8}`,
  z8: `0px 8px 16px 0px ${transparent16}`,
  z12: `0px 12px 24px -4px ${transparent12}`,
  z16: `0px 16px 32px -4px ${transparent16}`,
  z20: `0px 20px 40px -4px ${transparent24}`,
  
  // Component-specific shadows
  card: `0px 2px 8px 0px ${transparent8}`,
  dropdown: `0px 4px 16px 0px ${transparent16}`,
  dialog: `0px 24px 48px -12px ${transparent24}`,
  
  // Colored shadows for different state indicators
  primary: `0px 8px 16px 0px ${alpha('#1976D2', 0.24)}`,
  secondary: `0px 8px 16px 0px ${alpha('#00A699', 0.24)}`,
  info: `0px 8px 16px 0px ${alpha('#0288D1', 0.24)}`,
  success: `0px 8px 16px 0px ${alpha('#2E7D32', 0.24)}`,
  warning: `0px 8px 16px 0px ${alpha('#ED6C02', 0.24)}`,
  error: `0px 8px 16px 0px ${alpha('#D32F2F', 0.24)}`,
};

export default shadows;