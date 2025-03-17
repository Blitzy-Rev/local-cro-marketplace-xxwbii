import React, { useEffect } from 'react'; // React v18.2+
import { useDispatch, useSelector } from 'react-redux'; // react-redux v8.0+
import { Navigate, Outlet, useLocation } from 'react-router-dom'; // react-router-dom v6.4+
import { Box, CssBaseline, styled } from '@mui/material'; // @mui/material v5.13+
import { ThemeProvider } from '@mui/material/styles'; // @mui/material/styles v5.13+

import Header from './components/Header';
import Sidebar from './components/Sidebar';
import Footer from './components/Footer';
import useAuth from '../features/auth/hooks/useAuth';
import theme from '../theme';
import useWindowSize from '../hooks/useWindowSize';
import { selectSidebarOpen, setSidebarOpen } from '../store/ui/uiSlice';

/**
 * @file Main layout component for the Molecular Data Management and CRO Integration Platform that provides the common structure for authenticated pages. It includes a responsive header, sidebar navigation, and footer, adapting to different screen sizes and user roles.
 */

/**
 * @interface MainLayoutProps
 * @description Props interface for the MainLayout component
 */
interface MainLayoutProps {}

/**
 * @component MainContainer
 * @description Styled Box component for the main container
 */
const MainContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  minHeight: '100vh',
  backgroundColor: theme.palette.background.default,
});

/**
 * @component ContentWrapper
 * @description Styled Box component for the content wrapper
 */
const ContentWrapper = styled(Box)({
  display: 'flex',
  flex: '1',
});

/**
 * @component MainContent
 * @description Styled Box component for the main content area
 */
const MainContent = styled(Box)<{ sidebarOpen: boolean }>(({ theme, sidebarOpen }) => ({
  flexGrow: '1',
  padding: theme.spacing(3),
  marginTop: '64px',
  transition: theme.transitions.create('margin', {
    easing: theme.transitions.easing.sharp,
    duration: theme.transitions.duration.leavingScreen,
  }),
  marginLeft: sidebarOpen ? '240px' : theme.spacing(7),
  width: sidebarOpen ? `calc(100% - 240px)` : `calc(100% - ${theme.spacing(7)})`,
}));

/**
 * @component MainLayout
 * @description Main layout component with header, sidebar, and footer for authenticated pages
 */
const MainLayout: React.FC<MainLayoutProps> = () => {
  // Get authentication state and user data using useAuth hook
  const { isAuthenticated } = useAuth();

  // Get sidebar open state from Redux store using useSelector and selectSidebarOpen
  const sidebarOpen = useSelector(selectSidebarOpen);

  // Get dispatch function from Redux store using useDispatch
  const dispatch = useDispatch();

  // Get current location using useLocation hook
  const location = useLocation();

  // Get window size using useWindowSize hook
  const { isMobile } = useWindowSize();

  /**
   * @function handleToggleSidebar
   * @description Toggles the sidebar open/closed state
   */
  const handleToggleSidebar = () => {
    // Dispatch setSidebarOpen action with the opposite of current sidebar state
    dispatch(setSidebarOpen(!sidebarOpen));
  };

  /**
   * @function handleCloseSidebar
   * @description Closes the sidebar (used for mobile view)
   */
  const handleCloseSidebar = () => {
    // Dispatch setSidebarOpen action with false to close the sidebar
    dispatch(setSidebarOpen(false));
  };

  // Use useEffect to close sidebar automatically on mobile when route changes
  useEffect(() => {
    if (isMobile) {
      handleCloseSidebar();
    }
  }, [location, isMobile, handleCloseSidebar]);

  // Use useEffect to close sidebar automatically when screen size changes to mobile
  useEffect(() => {
    if (isMobile) {
      handleCloseSidebar();
    }
  }, [isMobile, handleCloseSidebar]);

  // If user is not authenticated, redirect to login page using Navigate component
  if (!isAuthenticated) {
    return <Navigate to="/login" />;
  }

  return (
    <ThemeProvider theme={theme}>
      <CssBaseline />
      <MainContainer>
        <Header handleToggleSidebar={handleToggleSidebar} />
        <ContentWrapper>
          <Sidebar open={sidebarOpen} onToggle={handleToggleSidebar} onClose={handleCloseSidebar} />
          <MainContent sidebarOpen={sidebarOpen}>
            <Outlet />
          </MainContent>
        </ContentWrapper>
        <Footer />
      </MainContainer>
    </ThemeProvider>
  );
};

export default MainLayout;