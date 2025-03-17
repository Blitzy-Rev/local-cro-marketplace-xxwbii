import React, { useState, useEffect, useCallback } from 'react'; // react 18.2+
import { useLocation, useNavigate, NavLink } from 'react-router-dom'; // react-router-dom 6.4+
import {
  Drawer,
  List,
  ListItem,
  ListItemIcon,
  ListItemText,
  Divider,
  Box,
  Tooltip,
  IconButton,
  Collapse,
} from '@mui/material'; // @mui/material 5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles 5.13+
import {
  DashboardIcon,
  MoleculeIcon,
  LibraryBooksIcon,
  ScienceIcon,
  SendIcon,
  ReceiptIcon,
  MessageIcon,
  ChevronLeftIcon,
  ChevronRightIcon,
  ExpandLess,
  ExpandMore,
  PeopleIcon,
  SettingsIcon,
} from '@mui/icons-material'; // @mui/icons-material 5.13+
import useAuth from '../../features/auth/hooks/useAuth';
import { UserRole } from '../../types/user';
import theme from '../../theme';
import useWindowSize from '../../hooks/useWindowSize';

/**
 * Props interface for the Sidebar component
 */
interface SidebarProps {
  /**
   * Whether the sidebar is open or closed
   */
  open: boolean;
  /**
   * Function to toggle the sidebar open/closed state
   */
  onToggle: () => void;
  /**
   * Function to close the sidebar (used for mobile view)
   */
  onClose: () => void;
}

/**
 * Interface for navigation menu items
 */
interface NavItem {
  /**
   * Unique identifier for the menu item
   */
  id: string;
  /**
   * Display text for the menu item
   */
  label: string;
  /**
   * Icon component to display with the menu item
   */
  icon: React.ReactNode;
  /**
   * Navigation path for the menu item
   */
  path: string;
  /**
   * User roles that can access this menu item
   */
  roles: UserRole[];
  /**
   * Submenu items (optional)
   */
  children?: NavItem[];
}

/**
 * Styled Drawer component for the sidebar
 */
const SidebarDrawer = styled(Drawer, { shouldForwardProp: (prop) => prop !== 'open' })(
  ({ theme, open }) => ({
    width: 240,
    flexShrink: 0,
    whiteSpace: 'nowrap',
    boxSizing: 'border-box',
    '& .MuiDrawer-paper': {
      width: 240,
      transition: theme.transitions.create('width', {
        easing: theme.transitions.easing.sharp,
        duration: theme.transitions.duration.enteringScreen,
      }),
      overflowX: 'hidden',
      backgroundColor: theme.palette.background.paper,
      borderRight: '1px solid',
    },
    '&.closed': {
      width: theme.spacing(7),
      '& .MuiDrawer-paper': {
        width: theme.spacing(7),
        transition: theme.transitions.create('width', {
          easing: theme.transitions.easing.sharp,
          duration: theme.transitions.duration.leavingScreen,
        }),
      },
    },
  })
);

/**
 * Styled Box component for the drawer header
 */
const DrawerHeader = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'flex-end',
  padding: theme.spacing(0, 1),
  minHeight: '64px',
}));

/**
 * Styled ListItem component for navigation items
 */
const NavItemStyled = styled(ListItem)(({ theme }) => ({
  padding: theme.spacing(1, 2),
  cursor: 'pointer',
  '&.active': {
    backgroundColor: theme.palette.action.selected,
  },
  '&:hover': {
    backgroundColor: theme.palette.action.hover,
  },
}));

/**
 * Styled List component for submenus
 */
const SubMenuList = styled(List)(({ theme }) => ({
  paddingLeft: theme.spacing(4),
}));

/**
 * Navigation sidebar component with role-based menu items
 */
const Sidebar: React.FC<SidebarProps> = ({ open, onToggle, onClose }) => {
  // Get user data and role checking function from useAuth hook
  const { user, checkRole } = useAuth();
  // Get current route using useLocation hook
  const location = useLocation();
  // Get navigation function using useNavigate hook
  const navigate = useNavigate();
  // Get screen size using useWindowSize hook
  const windowSize = useWindowSize();

  // State to track expanded submenus
  const [expandedMenus, setExpandedMenus] = useState<Set<string>>(new Set());

  /**
   * Toggles the sidebar open/closed state
   */
  const handleDrawerToggle = useCallback(() => {
    onToggle();
  }, [onToggle]);

  /**
   * Toggles the expansion state of a submenu
   * @param menuId - Unique identifier for the menu
   */
  const handleSubMenuToggle = useCallback((menuId: string) => {
    setExpandedMenus((prevExpandedMenus) => {
      const newExpandedMenus = new Set(prevExpandedMenus);
      if (newExpandedMenus.has(menuId)) {
        newExpandedMenus.delete(menuId);
      } else {
        newExpandedMenus.add(menuId);
      }
      return newExpandedMenus;
    });
  }, []);

  /**
   * Checks if a submenu is currently expanded
   * @param menuId - Unique identifier for the menu
   * @returns True if the submenu is expanded, false otherwise
   */
  const isSubmenuExpanded = useCallback((menuId: string) => {
    return expandedMenus.has(menuId);
  }, [expandedMenus]);

  /**
   * Checks if a route is currently active
   * @param path - Path to check
   * @returns True if the route is active, false otherwise
   */
  const isActiveRoute = useCallback((path: string) => {
    return location.pathname.startsWith(path);
  }, [location.pathname]);

  /**
   * Navigates to a specified route
   * @param path - Path to navigate to
   */
  const navigateTo = useCallback((path: string) => {
    navigate(path);
    if (windowSize.isMobile) {
      onClose();
    }
  }, [navigate, windowSize, onClose]);

  // Close submenus when sidebar collapses
  useEffect(() => {
    if (!open) {
      setExpandedMenus(new Set());
    }
  }, [open]);

  // Define navigation items with role-based access control
  const navItems: NavItem[] = [
    {
      id: 'dashboard',
      label: 'Dashboard',
      icon: <DashboardIcon />,
      path: '/app/dashboard',
      roles: [UserRole.PHARMA, UserRole.CRO, UserRole.ADMIN],
    },
    {
      id: 'molecules',
      label: 'Molecules',
      icon: <MoleculeIcon />,
      path: '/app/molecules',
      roles: [UserRole.PHARMA],
      children: [
        {
          id: 'molecules-list',
          label: 'All Molecules',
          path: '/app/molecules',
          roles: [UserRole.PHARMA],
        },
        {
          id: 'molecules-upload',
          label: 'Upload CSV',
          path: '/app/molecules/upload',
          roles: [UserRole.PHARMA],
        },
      ],
    },
    {
      id: 'libraries',
      label: 'Libraries',
      icon: <LibraryBooksIcon />,
      path: '/app/libraries',
      roles: [UserRole.PHARMA],
      children: [
        {
          id: 'libraries-list',
          label: 'All Libraries',
          path: '/app/libraries',
          roles: [UserRole.PHARMA],
        },
        {
          id: 'libraries-create',
          label: 'Create Library',
          path: '/app/libraries/create',
          roles: [UserRole.PHARMA],
        },
      ],
    },
    {
      id: 'experiments',
      label: 'Experiments',
      icon: <ScienceIcon />,
      path: '/app/experiments',
      roles: [UserRole.PHARMA],
      children: [
        {
          id: 'experiments-list',
          label: 'All Experiments',
          path: '/app/experiments',
          roles: [UserRole.PHARMA],
        },
        {
          id: 'experiments-create',
          label: 'Create Experiment',
          path: '/app/experiments/create',
          roles: [UserRole.PHARMA],
        },
      ],
    },
    {
      id: 'submissions',
      label: 'Submissions',
      icon: <SendIcon />,
      path: '/app/submissions',
      roles: [UserRole.PHARMA],
    },
    {
      id: 'cro-submissions',
      label: 'Submissions',
      icon: <ReceiptIcon />,
      path: '/cro/submissions',
      roles: [UserRole.CRO],
    },
    {
      id: 'results',
      label: 'Results',
      icon: <ReceiptIcon />,
      path: '/app/results',
      roles: [UserRole.PHARMA],
    },
    {
      id: 'communications',
      label: 'Communications',
      icon: <MessageIcon />,
      path: '/app/communications',
      roles: [UserRole.PHARMA, UserRole.CRO],
    },
    {
      id: 'admin-users',
      label: 'User Management',
      icon: <PeopleIcon />,
      path: '/admin/users',
      roles: [UserRole.ADMIN],
    },
    {
      id: 'admin-system',
      label: 'System Monitoring',
      icon: <SettingsIcon />,
      path: '/admin/system',
      roles: [UserRole.ADMIN],
    },
  ];

  return (
    <SidebarDrawer
      variant="permanent"
      open={open}
      className={open ? '' : 'closed'}
    >
      <DrawerHeader>
        <IconButton onClick={handleDrawerToggle}>
          {theme.direction === 'rtl' ? <ChevronRightIcon /> : <ChevronLeftIcon />}
        </IconButton>
      </DrawerHeader>
      <Divider />
      <List>
        {navItems.map((navItem) => {
          if (user && navItem.roles.some(role => checkRole(role))) {
            if (navItem.children) {
              return (
                <React.Fragment key={navItem.id}>
                  <NavItemStyled button onClick={() => handleSubMenuToggle(navItem.id)}>
                    <ListItemIcon>
                      {React.cloneElement(navItem.icon as React.ReactElement, {})}
                    </ListItemIcon>
                    <ListItemText primary={navItem.label} />
                    {isSubmenuExpanded(navItem.id) ? <ExpandLess /> : <ExpandMore />}
                  </NavItemStyled>
                  <Collapse in={isSubmenuExpanded(navItem.id)} timeout="auto" unmountOnExit>
                    <SubMenuList component="div" disablePadding>
                      {navItem.children.map((child) => {
                        if (user && child.roles.some(role => checkRole(role))) {
                          return (
                            <Tooltip title={!open ? child.label : ''} placement="right" key={child.id}>
                              <NavItemStyled
                                button
                                key={child.id}
                                component={NavLink}
                                to={child.path}
                                className={isActiveRoute(child.path) ? 'active' : ''}
                                onClick={() => navigateTo(child.path)}
                              >
                                <ListItemIcon>
                                  {React.cloneElement(child.icon as React.ReactElement, {})}
                                </ListItemIcon>
                                <ListItemText primary={child.label} />
                              </NavItemStyled>
                            </Tooltip>
                          );
                        }
                        return null;
                      })}
                    </SubMenuList>
                  </Collapse>
                </React.Fragment>
              );
            } else {
              return (
                <Tooltip title={!open ? navItem.label : ''} placement="right" key={navItem.id}>
                  <NavItemStyled
                    button
                    key={navItem.id}
                    component={NavLink}
                    to={navItem.path}
                    className={isActiveRoute(navItem.path) ? 'active' : ''}
                    onClick={() => navigateTo(navItem.path)}
                  >
                    <ListItemIcon>
                      {React.cloneElement(navItem.icon as React.ReactElement, {})}
                    </ListItemIcon>
                    <ListItemText primary={navItem.label} />
                  </NavItemStyled>
                </Tooltip>
              );
            }
          }
          return null;
        })}
      </List>
    </SidebarDrawer>
  );
};

export default Sidebar;