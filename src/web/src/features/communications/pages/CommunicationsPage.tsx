import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import {
  Box,
  Paper,
  Typography,
  Divider,
  Button,
  IconButton,
  useTheme,
  useMediaQuery,
} from '@mui/material'; // @mui/material 5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles 5.13+
import { Add, ArrowBack } from '@mui/icons-material'; // @mui/icons-material 5.13+

// Internal imports
import useCommunications from '../hooks/useCommunications';
import MessageList from '../components/MessageList';
import MessageDetail from '../components/MessageDetail';
import useWindowSize from '../../../hooks/useWindowSize';
import useAuth from '../../auth/hooks/useAuth';

/**
 * Interface defining the props for the CommunicationsPage component
 */
interface CommunicationsPageProps { }

/**
 * Styled component for the main page container
 */
const PageContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: 'calc(100vh - 64px - 48px)', // Adjust height based on header and footer
  padding: '24px',
  backgroundColor: 'background.default',
});

/**
 * Styled component for the content container
 */
const ContentContainer = styled(Box)({
  display: 'flex',
  flex: '1', // Take up remaining vertical space
  overflow: 'hidden', // Prevent content from overflowing
  borderRadius: '8px',
  boxShadow: '1', // Use a predefined shadow from the theme
  backgroundColor: 'background.paper',
});

/**
 * Styled component for the list container
 */
const ListContainer = styled(Box)({
  width: '320px', // Fixed width for the list
  borderRight: '1px solid',
  borderColor: 'divider', // Use the theme's divider color
  display: 'flex',
  flexDirection: 'column',
});

/**
 * Styled component for the detail container
 */
const DetailContainer = styled(Box)({
  flex: '1', // Take up remaining horizontal space
  display: 'flex',
  flexDirection: 'column',
  overflow: 'hidden', // Prevent content from overflowing
});

/**
 * Styled component for the list header
 */
const ListHeader = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  padding: '16px',
  borderBottom: '1px solid',
  borderColor: 'divider', // Use the theme's divider color
});

/**
 * Styled component for the empty state
 */
const EmptyState = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  height: '100%',
  padding: '24px',
  textAlign: 'center',
});

/**
 * Main component for the communications page
 */
const CommunicationsPage: React.FC<CommunicationsPageProps> = () => {
  // Get communications data and functions from useCommunications hook
  const { currentConversation, setCurrentConversation } = useCommunications();

  // Get current user data from useAuth hook
  const { user } = useAuth();

  // Get theme and media query for responsive design
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('md'));

  // Set up state for mobile view and drawer open status
  const [showDetail, setShowDetail] = useState<boolean>(false);

  /**
   * Handles creating a new conversation
   */
  const handleCreateConversation = useCallback(() => {
    // Open a dialog or navigate to create conversation page
    // This is a placeholder for future implementation
    console.log('Create new conversation');
  }, []);

  /**
   * Handles selecting a conversation from the list
   * @param conversation
   */
  const handleSelectConversation = useCallback((conversation: any) => {
    // Set the selected conversation using setCurrentConversation
    setCurrentConversation(conversation);

    // If in mobile view, set showDetail to true to display the message detail view
    if (isMobile) {
      setShowDetail(true);
    }
  }, [setCurrentConversation, isMobile]);

  /**
   * Handles back button click in mobile view
   */
  const handleBackClick = useCallback(() => {
    // Set showDetail to false to return to the conversation list view
    setShowDetail(false);

    // Clear the current conversation selection
    setCurrentConversation(null);
  }, [setCurrentConversation]);

  // Use useEffect to handle responsive layout changes based on screen size
  useEffect(() => {
    // If not mobile and detail is showing, hide detail
    if (!isMobile && showDetail) {
      setShowDetail(false);
    }

    // If not mobile and no conversation is selected, clear selection
    if (!isMobile && !currentConversation) {
      setCurrentConversation(null);
    }
  }, [isMobile, currentConversation, setCurrentConversation, showDetail]);

  // Render the main container with appropriate layout based on screen size
  return (
    <PageContainer>
      {/* Page title */}
      <Typography variant="h4" gutterBottom>
        Communications
      </Typography>

      {/* Main content container */}
      <ContentContainer>
        {/* For desktop view, render a two-column layout with MessageList and MessageDetail */}
        {!isMobile ? (
          <>
            {/* List container with MessageList */}
            <ListContainer>
              <ListHeader>
                <Typography variant="h6">Conversations</Typography>
                <Button variant="contained" size="small" onClick={handleCreateConversation}>
                  New
                </Button>
              </ListHeader>
              <MessageList onSelectMessage={handleSelectConversation} />
            </ListContainer>

            {/* Detail container with MessageDetail or EmptyState */}
            <DetailContainer>
              {currentConversation ? (
                <MessageDetail />
              ) : (
                <EmptyState>
                  <Typography variant="subtitle1">Select a conversation to view details.</Typography>
                </EmptyState>
              )}
            </DetailContainer>
          </>
        ) : (
          <>
            {/* For mobile view, conditionally render either MessageList or MessageDetail based on selected conversation */}
            {!showDetail ? (
              <MessageList onSelectMessage={handleSelectConversation} />
            ) : (
              <MessageDetail showBackButton onBackClick={handleBackClick} />
            )}
          </>
        )}
      </ContentContainer>
    </PageContainer>
  );
};

export default CommunicationsPage;