import { alpha, Components } from '@mui/material/styles'; // v5.13+
import { GREY, PRIMARY, SECONDARY } from './palette';
import { customShadows } from './shadows';

/**
 * Creates component style overrides based on the current theme mode
 * @param {string} mode - Theme mode ('light' or 'dark')
 * @returns {Components} Component style overrides configuration
 */
export function createComponentOverrides(mode: 'light' | 'dark'): Components {
  const isDark = mode === 'dark';
  
  // Common colors based on theme mode
  const backgroundDefault = isDark ? '#121212' : GREY[0];
  const backgroundPaper = isDark ? '#1E1E1E' : GREY[0];
  const backgroundNeutral = isDark ? GREY[900] : GREY[200];
  
  // Text colors based on theme mode
  const textPrimary = isDark ? GREY[0] : GREY[900];
  const textSecondary = isDark ? GREY[500_80] : GREY[600];
  const textDisabled = isDark ? GREY[600] : GREY[500];
  
  // Border and divider colors
  const borderColor = isDark ? GREY[500_16] : GREY[500_32];
  const dividerColor = isDark ? GREY[500_16] : GREY[500_16];
  
  // Action colors
  const actionActive = isDark ? GREY[500] : GREY[600];
  const actionHover = isDark ? GREY[500_8] : GREY[500_8];
  const actionSelected = isDark ? GREY[500_16] : GREY[500_16];
  const actionDisabled = isDark ? GREY[500_80] : GREY[500_32];
  const actionDisabledBackground = isDark ? GREY[500_24] : GREY[500_12];
  const actionFocus = isDark ? GREY[500_24] : GREY[500_16];
  
  return {
    MuiCssBaseline: {
      styleOverrides: {
        '*': {
          boxSizing: 'border-box',
        },
        html: {
          margin: 0,
          padding: 0,
          width: '100%',
          height: '100%',
          WebkitOverflowScrolling: 'touch',
        },
        body: {
          margin: 0,
          padding: 0,
          width: '100%',
          height: '100%',
        },
        '#root': {
          width: '100%',
          height: '100%',
        },
        input: {
          '&[type=number]': {
            MozAppearance: 'textfield',
            '&::-webkit-outer-spin-button': {
              margin: 0,
              WebkitAppearance: 'none',
            },
            '&::-webkit-inner-spin-button': {
              margin: 0,
              WebkitAppearance: 'none',
            },
          },
        },
        img: {
          display: 'block',
          maxWidth: '100%',
        },
        a: {
          textDecoration: 'none',
          color: PRIMARY.main,
        },
      },
    },
    MuiBackdrop: {
      styleOverrides: {
        root: {
          backgroundColor: alpha(GREY[900], 0.6),
        },
        invisible: {
          background: 'transparent',
        },
      },
    },
    MuiButton: {
      styleOverrides: {
        root: {
          textTransform: 'none',
          borderRadius: 8,
          fontWeight: 600,
          boxShadow: 'none',
          '&:hover': {
            boxShadow: 'none',
          },
          '&:active': {
            boxShadow: 'none',
          },
          '&.Mui-disabled': {
            boxShadow: 'none',
          },
        },
        sizeLarge: {
          height: 48,
          fontSize: 15,
        },
        sizeMedium: {
          height: 40,
          fontSize: 14,
        },
        sizeSmall: {
          height: 32,
          fontSize: 13,
          padding: '0 8px',
        },
        containedInherit: {
          color: textPrimary,
          backgroundColor: backgroundNeutral,
          '&:hover': {
            backgroundColor: isDark ? GREY[800] : GREY[300],
          },
        },
        containedPrimary: {
          boxShadow: customShadows.primary,
          '&:hover': {
            backgroundColor: PRIMARY.dark,
          },
          '&:active': {
            backgroundColor: PRIMARY.darker,
          },
        },
        containedSecondary: {
          boxShadow: customShadows.secondary,
          '&:hover': {
            backgroundColor: SECONDARY.dark,
          },
          '&:active': {
            backgroundColor: SECONDARY.darker,
          },
        },
        outlined: {
          borderColor: borderColor,
          '&:hover': {
            backgroundColor: actionHover,
            borderColor: isDark ? GREY[500] : GREY[700],
          },
        },
        outlinedInherit: {
          color: textPrimary,
          '&:hover': {
            backgroundColor: actionHover,
          },
        },
        outlinedPrimary: {
          borderColor: PRIMARY.main,
          '&:hover': {
            backgroundColor: alpha(PRIMARY.main, 0.08),
            borderColor: PRIMARY.dark,
          },
        },
        outlinedSecondary: {
          borderColor: SECONDARY.main,
          '&:hover': {
            backgroundColor: alpha(SECONDARY.main, 0.08),
            borderColor: SECONDARY.dark,
          },
        },
        text: {
          '&:hover': {
            backgroundColor: actionHover,
          },
        },
        textInherit: {
          color: textPrimary,
        },
        textPrimary: {
          '&:hover': {
            backgroundColor: alpha(PRIMARY.main, 0.08),
          },
        },
        textSecondary: {
          '&:hover': {
            backgroundColor: alpha(SECONDARY.main, 0.08),
          },
        },
      },
    },
    MuiCard: {
      styleOverrides: {
        root: {
          boxShadow: customShadows.card,
          borderRadius: 12,
          position: 'relative',
          zIndex: 0,
          backgroundColor: backgroundPaper,
          transition: 'all .2s ease-in-out',
          '&:hover': {
            boxShadow: customShadows.z8,
            transform: 'translateY(-4px)',
          },
        },
      },
    },
    MuiCardHeader: {
      defaultProps: {
        titleTypographyProps: { variant: 'h6' },
        subheaderTypographyProps: { variant: 'body2', marginTop: 0.5 },
      },
      styleOverrides: {
        root: {
          padding: 24,
        },
      },
    },
    MuiCardContent: {
      styleOverrides: {
        root: {
          padding: 24,
          '&:last-child': {
            paddingBottom: 24,
          },
        },
      },
    },
    MuiCardActions: {
      styleOverrides: {
        root: {
          padding: '16px 24px',
        },
      },
    },
    MuiCheckbox: {
      styleOverrides: {
        root: {
          padding: 8,
          '&:hover': {
            backgroundColor: actionHover,
          },
          '&.Mui-checked': {
            '&:hover': {
              backgroundColor: alpha(PRIMARY.main, 0.08),
            },
          },
        },
      },
    },
    MuiChip: {
      styleOverrides: {
        root: {
          borderRadius: 8,
          fontWeight: 500,
          fontSize: 13,
        },
        filled: {
          '&.MuiChip-colorDefault': {
            backgroundColor: isDark ? GREY[800] : GREY[300],
            color: isDark ? GREY[0] : GREY[800],
          },
        },
        outlined: {
          borderColor: borderColor,
          '&.MuiChip-colorDefault': {
            color: textPrimary,
          },
        },
      },
    },
    MuiDialog: {
      styleOverrides: {
        paper: {
          boxShadow: customShadows.dialog,
          borderRadius: 12,
          backgroundColor: backgroundPaper,
        },
      },
    },
    MuiDialogTitle: {
      styleOverrides: {
        root: {
          padding: 24,
          fontSize: 18,
          fontWeight: 600,
        },
      },
    },
    MuiDialogContent: {
      styleOverrides: {
        root: {
          padding: 24,
          borderTop: `1px solid ${dividerColor}`,
          borderBottom: `1px solid ${dividerColor}`,
        },
      },
    },
    MuiDialogActions: {
      styleOverrides: {
        root: {
          padding: 16,
          '& > :not(:first-of-type)': {
            marginLeft: 8,
          },
        },
      },
    },
    MuiDivider: {
      styleOverrides: {
        root: {
          borderColor: dividerColor,
        },
      },
    },
    MuiDrawer: {
      styleOverrides: {
        paper: {
          backgroundColor: backgroundPaper,
          boxShadow: customShadows.z16,
        },
      },
    },
    MuiIconButton: {
      styleOverrides: {
        root: {
          borderRadius: 8,
          padding: 8,
          '&:hover': {
            backgroundColor: actionHover,
          },
          '& svg': {
            fontSize: 24,
          },
        },
        sizeSmall: {
          padding: 4,
          '& svg': {
            fontSize: 20,
          },
        },
      },
    },
    MuiInputBase: {
      styleOverrides: {
        root: {
          '&.Mui-disabled': {
            '& svg': { color: textDisabled },
          },
        },
        input: {
          '&::placeholder': {
            opacity: 1,
            color: textDisabled,
          },
        },
      },
    },
    MuiFilledInput: {
      styleOverrides: {
        root: {
          backgroundColor: isDark ? GREY[900] : GREY[200],
          borderRadius: 8,
          '&:hover': {
            backgroundColor: isDark ? GREY[800] : GREY[300],
          },
          '&.Mui-focused': {
            backgroundColor: isDark ? GREY[900] : GREY[200],
          },
          '&.Mui-disabled': {
            backgroundColor: actionDisabledBackground,
          },
        },
        underline: {
          '&:before, &:after': {
            display: 'none',
          },
        },
      },
    },
    MuiOutlinedInput: {
      styleOverrides: {
        root: {
          borderRadius: 8,
          '& .MuiOutlinedInput-notchedOutline': {
            borderColor: borderColor,
          },
          '&:hover .MuiOutlinedInput-notchedOutline': {
            borderColor: PRIMARY.light,
          },
          '&.Mui-focused .MuiOutlinedInput-notchedOutline': {
            borderWidth: 1,
            borderColor: PRIMARY.main,
          },
          '&.Mui-disabled .MuiOutlinedInput-notchedOutline': {
            borderColor: actionDisabledBackground,
          },
        },
      },
    },
    MuiLink: {
      defaultProps: {
        underline: 'hover',
      },
      styleOverrides: {
        root: {
          color: PRIMARY.main,
          '&:hover': {
            color: PRIMARY.dark,
          },
        },
      },
    },
    MuiLinearProgress: {
      styleOverrides: {
        root: {
          borderRadius: 4,
          overflow: 'hidden',
        },
        bar: {
          borderRadius: 4,
        },
        colorPrimary: {
          backgroundColor: isDark ? PRIMARY.darker : PRIMARY.lighter,
        },
      },
    },
    MuiList: {
      styleOverrides: {
        root: {
          padding: 0,
        },
      },
    },
    MuiListItemButton: {
      styleOverrides: {
        root: {
          padding: '10px 12px',
          borderRadius: 8,
          '&.Mui-selected': {
            backgroundColor: isDark ? PRIMARY.darker : PRIMARY.lighter,
            color: isDark ? PRIMARY.lighter : PRIMARY.darker,
            '&:hover': {
              backgroundColor: isDark ? PRIMARY.darker : PRIMARY.lighter,
            },
          },
          '&:hover': {
            backgroundColor: actionHover,
          },
        },
      },
    },
    MuiListItemIcon: {
      styleOverrides: {
        root: {
          minWidth: 'auto',
          marginRight: 16,
        },
      },
    },
    MuiListItemText: {
      styleOverrides: {
        root: {
          margin: 0,
        },
      },
    },
    MuiMenu: {
      styleOverrides: {
        paper: {
          boxShadow: customShadows.dropdown,
          borderRadius: 8,
        },
      },
    },
    MuiMenuItem: {
      styleOverrides: {
        root: {
          padding: '10px 12px',
          '&.Mui-selected': {
            backgroundColor: isDark ? PRIMARY.darker : PRIMARY.lighter,
            color: isDark ? PRIMARY.lighter : PRIMARY.darker,
            '&:hover': {
              backgroundColor: isDark ? alpha(PRIMARY.darker, 0.8) : alpha(PRIMARY.lighter, 0.8),
            },
          },
          '&:hover': {
            backgroundColor: actionHover,
          },
        },
      },
    },
    MuiPaper: {
      defaultProps: {
        elevation: 0,
      },
      styleOverrides: {
        root: {
          backgroundImage: 'none',
          backgroundColor: backgroundPaper,
        },
        outlined: {
          borderColor: borderColor,
        },
      },
    },
    MuiPopover: {
      styleOverrides: {
        paper: {
          boxShadow: customShadows.dropdown,
          borderRadius: 8,
        },
      },
    },
    MuiRadio: {
      styleOverrides: {
        root: {
          padding: 8,
          '&:hover': {
            backgroundColor: actionHover,
          },
          '&.Mui-checked': {
            '&:hover': {
              backgroundColor: alpha(PRIMARY.main, 0.08),
            },
          },
        },
      },
    },
    MuiSelect: {
      styleOverrides: {
        icon: {
          right: 10,
          width: 18,
          height: 18,
          top: 'calc(50% - 9px)',
        },
      },
    },
    MuiSlider: {
      styleOverrides: {
        root: {
          height: 8,
          '& .MuiSlider-thumb': {
            width: 20,
            height: 20,
            '&::before': {
              boxShadow: '0 2px 12px 0 rgba(0,0,0,0.1)',
            },
            '&:hover, &.Mui-focusVisible': {
              boxShadow: '0px 0px 0px 8px rgba(25, 118, 210, 0.16)',
            },
          },
          '& .MuiSlider-valueLabel': {
            borderRadius: 8,
            backgroundColor: isDark ? GREY[900] : GREY[50],
            color: isDark ? GREY[0] : GREY[900],
            fontWeight: 500,
            fontSize: 13,
            lineHeight: 1.4,
            padding: '4px 8px',
            boxShadow: customShadows.z8,
          },
          '& .MuiSlider-track': {
            border: 'none',
            height: 8,
            borderRadius: 4,
          },
          '& .MuiSlider-rail': {
            opacity: 0.5,
            height: 8,
            borderRadius: 4,
          },
          '& .MuiSlider-mark': {
            backgroundColor: GREY[500],
            width: 8,
            height: 8,
            borderRadius: '50%',
          },
        },
      },
    },
    MuiSnackbar: {
      styleOverrides: {
        root: {
          '& .MuiPaper-root': {
            boxShadow: customShadows.z8,
            borderRadius: 8,
          },
        },
        anchorOriginTopRight: {
          top: 16,
          right: 16,
        },
        anchorOriginBottomRight: {
          bottom: 16,
          right: 16,
        },
        anchorOriginTopLeft: {
          top: 16,
          left: 16,
        },
        anchorOriginBottomLeft: {
          bottom: 16,
          left: 16,
        },
      },
    },
    MuiSwitch: {
      styleOverrides: {
        root: {
          width: 40,
          height: 24,
          padding: 0,
          '& .MuiSwitch-switchBase': {
            padding: 1,
            '&.Mui-checked': {
              transform: 'translateX(16px)',
              '& + .MuiSwitch-track': {
                opacity: 1,
              },
            },
          },
          '& .MuiSwitch-thumb': {
            width: 22,
            height: 22,
            borderRadius: 11,
          },
          '& .MuiSwitch-track': {
            opacity: 1,
            borderRadius: 12,
          },
        },
      },
    },
    MuiTableCell: {
      styleOverrides: {
        root: {
          padding: 12,
          borderBottom: `1px solid ${borderColor}`,
        },
        head: {
          backgroundColor: isDark ? GREY[900] : GREY[100],
          color: textPrimary,
          fontWeight: 600,
        },
      },
    },
    MuiTableRow: {
      styleOverrides: {
        root: {
          '&.Mui-selected': {
            backgroundColor: isDark ? PRIMARY.darker : PRIMARY.lighter,
            '&:hover': {
              backgroundColor: isDark ? alpha(PRIMARY.darker, 0.8) : alpha(PRIMARY.lighter, 0.8),
            },
          },
          '&:hover': {
            backgroundColor: actionHover,
          },
        },
      },
    },
    MuiTab: {
      styleOverrides: {
        root: {
          textTransform: 'none',
          fontWeight: 600,
          fontSize: 14,
          padding: '12px 16px',
          minWidth: 100,
          '&.Mui-selected': {
            color: PRIMARY.main,
          },
        },
      },
    },
    MuiTabs: {
      styleOverrides: {
        root: {
          minHeight: 44,
        },
        indicator: {
          backgroundColor: PRIMARY.main,
          height: 3,
        },
      },
    },
    MuiTextField: {
      defaultProps: {
        variant: 'outlined',
      },
    },
    MuiTooltip: {
      styleOverrides: {
        tooltip: {
          backgroundColor: isDark ? GREY[800] : GREY[900],
          color: isDark ? GREY[0] : GREY[0],
          fontSize: 13,
          fontWeight: 500,
          padding: '8px 12px',
          borderRadius: 6,
          boxShadow: customShadows.z8,
          maxWidth: 300,
          wordWrap: 'break-word',
        },
        arrow: {
          color: isDark ? GREY[800] : GREY[900],
        },
      },
    },
    MuiTypography: {
      styleOverrides: {
        root: {
          color: textPrimary,
        },
        paragraph: {
          marginBottom: 16,
        },
        gutterBottom: {
          marginBottom: 8,
        },
        body1: {
          fontSize: 16,
          lineHeight: 1.5,
        },
        body2: {
          fontSize: 14,
          lineHeight: 1.5,
        },
      },
    },
  };
}

// Default light mode component overrides
const components = createComponentOverrides('light');

export default { components };