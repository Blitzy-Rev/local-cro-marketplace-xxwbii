import React, { useState, useEffect } from 'react'; // React v18.2+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.10+
import { useSelector } from 'react-redux'; // react-redux v8.0+
import { 
  Box, 
  Typography, 
  Skeleton, 
  Divider, 
  List, 
  ListItem, 
  ListItemText, 
  ListItemIcon,
  useTheme,
  useMediaQuery
} from '@mui/material'; // @mui/material v5.13+
import { 
  Message as MessageIcon, 
  Notifications as NotificationsIcon 
} from '@mui/icons-material'; // @mui/icons-material v5.13+

import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Badge from '../../../components/common/Badge';
import useCommunications from '../../communications/hooks/useCommunications';
import { selectUnreadCount } from '../../../store/notifications/notificationsSlice';
import { formatDate } from '../../../utils/formatters';

// Define the props for the CROCommunications component
interface CROCommunicationsProps {
  onOpenMessages?: () => void;
  maxItems?: number;
  className?: string;
  style?: React.CSSProperties;
}

/**
 * Dashboard widget component that displays a summary of CRO communications
 * @param {CROCommunicationsProps} props - The props for the component
 * @returns {JSX.Element} The rendered component
 */
const CROCommunications: React.FC<CROCommunicationsProps> = ({
  onOpenMessages,
  maxItems = 3,
  className = '',
  style = {},
}) => {
  // Initialize navigation hook for routing
  const navigate = useNavigate();

  // Initialize theme and media query for responsive design
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // Initialize communications hook to fetch recent conversations
  const { conversations, loading } = useCommunications();

  // Get unread notification count from Redux store
  const unreadCount = useSelector(selectUnreadCount);

  // Define a function to handle navigation to the messages page
  const handleNavigateToMessages = () => {
    if (onOpenMessages) {
      onOpenMessages();
    } else {
      navigate('/communications');
    }
  };

  // Define a function to render a list item for a single message
  const renderMessageItem = (conversation: any) => (
    <ListItem key={conversation.id} alignItems="flex-start">
      <ListItemIcon>
        <MessageIcon color="primary" />
      </ListItemIcon>
      <ListItemText
        primary={
          <Typography variant="subtitle2" component="span" style={{ fontWeight: 500 }}>
            {conversation.senderName}
          </Typography>
        }
        secondary={
          <React.Fragment>
            <Typography
              sx={{ display: 'inline' }}
              component="span"
              variant="body2"
              color="text.primary"
            >
              {conversation.subject}
            </Typography>
            {` — ${formatDate(conversation.lastMessageDate)}`}
          </React.Fragment>
        }
      />
      {conversation.unread && <Badge label="New" color="primary" />}
    </ListItem>
  );

  // Define a function to render a skeleton placeholder list item for loading state
  const renderSkeletonItem = () => (
    <ListItem key="skeleton" alignItems="flex-start">
      <ListItemIcon>
        <Skeleton variant="circular" width={24} height={24} />
      </ListItemIcon>
      <ListItemText
        primary={<Skeleton variant="text" width={120} />}
        secondary={<Skeleton variant="text" width={80} />}
      />
    </ListItem>
  );

  // Apply responsive layout based on screen size
  const MessageListContainer = styled(Box)({
    width: '100%',
    height: '100%',
    display: 'flex',
    flexDirection: 'column',
    overflow: 'hidden',
  });

  const EmptyStateContainer = styled(Box)({
    display: 'flex',
    flexDirection: 'column',
    alignItems: 'center',
    justifyContent: 'center',
    padding: '32px',
    color: theme.palette.text.secondary,
  });

  // Fetch conversations when maxItems changes
  useEffect(() => {
    if (!loading) {
    }
  }, [maxItems, loading]);

  return (
    <Card className={className} style={style}>
      <Box display="flex" justifyContent="space-between" alignItems="center" mb={2} p={2}>
        <Typography variant="h6" component="div">
          CRO Communications
          {unreadCount > 0 && (
            <Badge label={unreadCount.toString()} color="primary" style={{ marginLeft: 8 }} />
          )}
        </Typography>
      </Box>
      <Divider />
      <MessageListContainer>
        <List>
          {loading ? (
            Array.from(new Array(maxItems)).map((_, index) => renderSkeletonItem())
          ) : conversations && conversations.length > 0 ? (
            conversations.slice(0, maxItems).map(conversation => renderMessageItem(conversation))
          ) : (
            <EmptyStateContainer>
              <NotificationsIcon sx={{ fontSize: 40, mb: 1 }} />
              <Typography variant="body2">No recent communications</Typography>
            </EmptyStateContainer>
          )}
        </List>
      </MessageListContainer>
      <Divider />
      <Box p={2} display="flex" justifyContent="flex-end">
        <Button 
          variant="contained" 
          color="primary" 
          onClick={handleNavigateToMessages}
          size={isMobile ? 'small' : 'medium'}
        >
          Open Messages
        </Button>
      </Box>
    </Card>
  );
};

export default CROCommunications;