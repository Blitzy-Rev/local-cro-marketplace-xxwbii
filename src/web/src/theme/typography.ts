import { TypographyOptions } from '@mui/material/styles'; // @mui/material/styles v5.13+

// Font families
const FONT_PRIMARY = "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif";
const FONT_SECONDARY = "'Roboto Slab', serif";
const FONT_MONO = "'Roboto Mono', monospace";

// Base font size in pixels for rem calculations
const BASE_FONT_SIZE = 16;

/**
 * Converts pixel values to rem units for responsive typography
 * @param size Font size in pixels
 * @returns Font size in rem units
 */
function pxToRem(size: number): string {
  return `${size / BASE_FONT_SIZE}rem`;
}

/**
 * Adjusts font sizes based on screen size breakpoints
 * @param typography Typography configuration
 * @returns Updated typography with responsive font sizes
 */
function responsiveFontSizes(typography: TypographyOptions): TypographyOptions {
  // Create a deep copy to avoid modifying the original
  const typographyCopy = JSON.parse(JSON.stringify(typography)) as TypographyOptions;
  
  // Define sizes for different breakpoints
  const fontSizes = {
    h1: { xs: 32, sm: 36, md: 38, lg: 40 },
    h2: { xs: 24, sm: 28, md: 30, lg: 32 },
    h3: { xs: 20, sm: 22, md: 23, lg: 24 },
    h4: { xs: 18, sm: 19, md: 19, lg: 20 },
    h5: { xs: 16, sm: 17, md: 17, lg: 18 },
    h6: { xs: 14, sm: 15, md: 15, lg: 16 },
  };
  
  // Adjust heading sizes for smaller screens
  // Note: This simplifies the responsive approach for the typography system
  // In practice, these would be applied using Material-UI's breakpoint system
  if (typographyCopy.h1) {
    typographyCopy.h1.fontSize = pxToRem(fontSizes.h1.xs);
  }
  
  if (typographyCopy.h2) {
    typographyCopy.h2.fontSize = pxToRem(fontSizes.h2.xs);
  }
  
  if (typographyCopy.h3) {
    typographyCopy.h3.fontSize = pxToRem(fontSizes.h3.xs);
  }
  
  if (typographyCopy.h4) {
    typographyCopy.h4.fontSize = pxToRem(fontSizes.h4.xs);
  }
  
  if (typographyCopy.h5) {
    typographyCopy.h5.fontSize = pxToRem(fontSizes.h5.xs);
  }
  
  if (typographyCopy.h6) {
    typographyCopy.h6.fontSize = pxToRem(fontSizes.h6.xs);
  }
  
  return typographyCopy;
}

// Typography configuration
const typography: TypographyOptions = {
  fontFamily: FONT_PRIMARY,
  fontWeightRegular: 400,
  fontWeightMedium: 500,
  fontWeightBold: 700,
  h1: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 700,
    lineHeight: 1.2,
    fontSize: pxToRem(40),
  },
  h2: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 700,
    lineHeight: 1.3,
    fontSize: pxToRem(32),
  },
  h3: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 600,
    lineHeight: 1.4,
    fontSize: pxToRem(24),
  },
  h4: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 600,
    lineHeight: 1.4,
    fontSize: pxToRem(20),
  },
  h5: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 600,
    lineHeight: 1.5,
    fontSize: pxToRem(18),
  },
  h6: {
    fontFamily: FONT_SECONDARY,
    fontWeight: 600,
    lineHeight: 1.6,
    fontSize: pxToRem(16),
  },
  subtitle1: {
    fontWeight: 500,
    lineHeight: 1.5,
    fontSize: pxToRem(16),
  },
  subtitle2: {
    fontWeight: 500,
    lineHeight: 1.6,
    fontSize: pxToRem(14),
  },
  body1: {
    fontWeight: 400,
    lineHeight: 1.5,
    fontSize: pxToRem(16),
  },
  body2: {
    fontWeight: 400,
    lineHeight: 1.6,
    fontSize: pxToRem(14),
  },
  caption: {
    fontWeight: 400,
    lineHeight: 1.5,
    fontSize: pxToRem(12),
  },
  overline: {
    fontWeight: 600,
    lineHeight: 1.5,
    fontSize: pxToRem(12),
    textTransform: 'uppercase',
    letterSpacing: '1px',
  },
  button: {
    fontWeight: 500,
    lineHeight: 1.5,
    fontSize: pxToRem(14),
    textTransform: 'none',
  },
  code: {
    fontFamily: FONT_MONO,
    fontSize: pxToRem(14),
  },
};

export { responsiveFontSizes };
export default { typography };