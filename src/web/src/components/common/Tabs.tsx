import React from 'react'; // React 18.2+
import { Tabs as MuiTabs, Tab as MuiTab, Box } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+

// Internal imports
import theme from '../../theme'; // Access theme configuration for consistent styling

// Interface for tab panel props
interface TabPanelProps extends React.HTMLAttributes<HTMLDivElement> {
  children: React.ReactNode;
  value: number;
  index: number;
}

// Interface for tab item definition
interface TabItem {
  label: string;
  icon?: React.ReactNode;
  disabled?: boolean;
}

// Interface for Tabs component props
interface TabsProps {
  tabs: TabItem[];
  value: number;
  onChange: (event: React.SyntheticEvent, newValue: number) => void;
  orientation?: 'horizontal' | 'vertical';
  variant?: 'standard' | 'scrollable' | 'fullWidth';
  children?: React.ReactNode;
}

// Custom styled Tabs component with theme styling
const StyledTabs = styled(MuiTabs)(({ theme }) => ({
  borderBottom: '1px solid',
  borderColor: theme.palette.divider,
  '& .MuiTabs-indicator': {
    backgroundColor: theme.palette.primary.main,
    height: '3px',
  },
}));

// Custom styled Tab component with theme styling
const StyledTab = styled(MuiTab)(({ theme }) => ({
  textTransform: 'none',
  fontWeight: theme.typography.fontWeightMedium,
  fontSize: theme.typography.pxToRem(15),
  marginRight: theme.spacing(1),
  color: theme.palette.text.secondary,
  '&.Mui-selected': {
    color: theme.palette.primary.main,
    fontWeight: theme.typography.fontWeightBold,
  },
  '&:hover': {
    color: theme.palette.primary.main,
    opacity: '0.8',
  },
}));

/**
 * Helper function to generate accessibility props for tabs
 * @param index - Tab index
 * @returns Object with id and aria-controls attributes
 */
function a11yProps(index: number) {
  return {
    id: `tab-${index}`,
    'aria-controls': `tabpanel-${index}`,
  };
}

/**
 * Component that renders the content for a specific tab
 */
export const TabPanel: React.FC<TabPanelProps> = (props) => {
  const { children, value, index, ...other } = props;

  if (value !== index) {
    return null;
  }

  return (
    <Box
      role="tabpanel"
      id={`tabpanel-${index}`}
      aria-labelledby={`tab-${index}`}
      p={3}
      {...other}
    >
      {children}
    </Box>
  );
};

/**
 * Main tabs component that provides a tabbed interface
 */
const Tabs: React.FC<TabsProps> = ({
  tabs,
  value,
  onChange,
  orientation = 'horizontal',
  variant = 'standard',
  children,
}) => {
  return (
    <Box sx={{ width: '100%' }}>
      <StyledTabs 
        value={value} 
        onChange={onChange} 
        orientation={orientation}
        variant={variant}
        aria-label="tabs"
      >
        {tabs.map((tab, index) => (
          <StyledTab
            key={index}
            label={tab.label}
            icon={tab.icon}
            disabled={tab.disabled}
            {...a11yProps(index)}
          />
        ))}
      </StyledTabs>
      {children}
    </Box>
  );
};

export default Tabs;