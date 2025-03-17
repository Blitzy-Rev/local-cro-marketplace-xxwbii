import React, { useState, useEffect, useCallback } from 'react'; // React v18.2+
import { useDispatch, useSelector } from 'react-redux'; // react-redux v8.0+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.4+
import {
  AppBar,
  Toolbar,
  IconButton,
  Typography,
  Avatar,
  Menu,
  MenuItem,
  Badge,
  Box,
  Tooltip,
  Divider,
  ListItemIcon,
  ListItemText,
  Switch,
} from '@mui/material'; // @mui/material v5.13+
import { styled, PaletteMode } from '@mui/material/styles'; // @mui/material/styles v5.13+
import {
  MenuIcon,
  NotificationsIcon,
  AccountCircle,
  Logout,
  Settings,
  DarkMode,
  LightMode,
} from '@mui/icons-material'; // @mui/icons-material v5.13+

import useAuth from '../../features/auth/hooks/useAuth';
import { UserRole } from '../../types/user';
import theme, { createCustomTheme, getStoredThemeMode, storeThemeMode } from '../../theme';
import Notifications from './Notifications';
import useWindowSize from '../../hooks/useWindowSize';
import { selectSidebarOpen, setSidebarOpen, setThemeMode, selectThemeMode } from '../../store/ui/uiSlice';
import { selectUnreadCount, fetchNotifications } from '../../store/notifications/notificationsSlice';

/**
 * Props interface for the Header component
 */
interface HeaderProps {}

/**
 * Styled AppBar component for consistent header styling
 */
const StyledAppBar = styled(AppBar)(({ theme }) => ({
  position: 'fixed',
  zIndex: theme.zIndex.drawer + 1,
  boxShadow: theme.shadows[2],
  backgroundColor: theme.palette.background.paper,
  color: theme.palette.text.primary,
}));

/**
 * Styled Box component for the logo container
 */
const LogoContainer = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  marginRight: theme.spacing(2),
}));

/**
 * Styled Typography component for the logo text
 */
const LogoText = styled(Typography)(({ theme }) => ({
  fontWeight: '600',
  fontSize: '1.2rem',
  display: { xs: 'none', sm: 'block' },
}));

/**
 * Styled Box component for header actions
 */
const HeaderActions = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  marginLeft: 'auto',
}));

/**
 * Main header component with navigation controls, user profile, and notifications
 */
const Header: React.FC<HeaderProps> = () => {
  // Component state for menu anchors
  const [anchorEl, setAnchorEl] = useState<HTMLElement | null>(null);
  const [notificationsAnchorEl, setNotificationsAnchorEl] = useState<HTMLElement | null>(null);

  // Authentication state and functions
  const { user, logout } = useAuth();

  // Redux state and dispatch
  const sidebarOpen = useSelector(selectSidebarOpen);
  const themeMode = useSelector(selectThemeMode);
  const unreadCount = useSelector(selectUnreadCount);
  const dispatch = useDispatch();

  // Navigation hook
  const navigate = useNavigate();

  // Responsive design hook
  const { isMobile } = useWindowSize();

  /**
   * Toggles the sidebar open/closed state
   */
  const handleToggleSidebar = useCallback(() => {
    dispatch(setSidebarOpen(!sidebarOpen));
  }, [dispatch, sidebarOpen]);

  /**
   * Opens the user profile menu
   */
  const handleProfileMenuOpen = useCallback((event: React.MouseEvent<HTMLElement>) => {
    setAnchorEl(event.currentTarget);
  }, []);

  /**
   * Closes the user profile menu
   */
  const handleProfileMenuClose = useCallback(() => {
    setAnchorEl(null);
  }, []);

  /**
   * Opens the notifications panel
   */
  const handleNotificationsOpen = useCallback((event: React.MouseEvent<HTMLElement>) => {
    setNotificationsAnchorEl(event.currentTarget);
    dispatch(fetchNotifications({ page: 1, pageSize: 5 }));
  }, [dispatch]);

  /**
   * Closes the notifications panel
   */
  const handleNotificationsClose = useCallback(() => {
    setNotificationsAnchorEl(null);
  }, []);

  /**
   * Handles user logout
   */
  const handleLogout = useCallback(() => {
    handleProfileMenuClose();
    logout();
    navigate('/login');
  }, [logout, navigate, handleProfileMenuClose]);

  /**
   * Navigates to user settings page
   */
  const handleSettings = useCallback(() => {
    handleProfileMenuClose();
    navigate('/settings');
  }, [navigate, handleProfileMenuClose]);

  /**
   * Toggles between light and dark theme modes
   */
  const handleThemeToggle = useCallback(() => {
    const newThemeMode: PaletteMode = themeMode === 'light' ? 'dark' : 'light';
    dispatch(setThemeMode(newThemeMode));
    storeThemeMode(newThemeMode);
  }, [dispatch, themeMode]);

  // Initialize theme mode from storage
  useEffect(() => {
    const storedTheme = getStoredThemeMode();
    if (storedTheme) {
      dispatch(setThemeMode(storedTheme));
    }
  }, [dispatch]);

  return (
    <StyledAppBar>
      <Toolbar>
        <IconButton
          edge="start"
          color="inherit"
          aria-label="menu"
          onClick={handleToggleSidebar}
          sx={{ mr: 2 }}
        >
          <MenuIcon />
        </IconButton>
        <LogoContainer>
          <LogoText variant="h6">
            Molecular Data Platform
          </LogoText>
        </LogoContainer>
        <HeaderActions>
          <Tooltip title="Notifications">
            <IconButton color="inherit" onClick={handleNotificationsOpen}>
              <Badge badgeContent={unreadCount} color="error">
                <NotificationsIcon />
              </Badge>
            </IconButton>
          </Tooltip>
          <IconButton
            edge="end"
            aria-label="account of current user"
            aria-controls="profile-menu"
            aria-haspopup="true"
            onClick={handleProfileMenuOpen}
            color="inherit"
          >
            <Avatar>{user?.email.charAt(0).toUpperCase()}</Avatar>
          </IconButton>
          <Menu
            id="profile-menu"
            anchorEl={anchorEl}
            anchorOrigin={{
              vertical: 'bottom',
              horizontal: 'right',
            }}
            keepMounted
            transformOrigin={{
              vertical: 'top',
              horizontal: 'right',
            }}
            open={Boolean(anchorEl)}
            onClose={handleProfileMenuClose}
          >
            <MenuItem disabled>
              <Typography variant="body2">
                {user?.email}
              </Typography>
            </MenuItem>
            <Divider />
            <MenuItem onClick={handleThemeToggle}>
              <ListItemIcon>
                {themeMode === 'light' ? <DarkMode fontSize="small" /> : <LightMode fontSize="small" />}
              </ListItemIcon>
              <ListItemText>
                {themeMode === 'light' ? 'Dark Mode' : 'Light Mode'}
              </ListItemText>
              <Switch checked={themeMode === 'dark'} onChange={handleThemeToggle} />
            </MenuItem>
            <MenuItem onClick={handleSettings}>
              <ListItemIcon>
                <Settings fontSize="small" />
              </ListItemIcon>
              <ListItemText>Settings</ListItemText>
            </MenuItem>
            <MenuItem onClick={handleLogout}>
              <ListItemIcon>
                <Logout fontSize="small" />
              </ListItemIcon>
              <ListItemText>Logout</ListItemText>
            </MenuItem>
          </Menu>
        </HeaderActions>
      </Toolbar>
      <Notifications anchorEl={notificationsAnchorEl} onClose={handleNotificationsClose} />
    </StyledAppBar>
  );
};

export default Header;