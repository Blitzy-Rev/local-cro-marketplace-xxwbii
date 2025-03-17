import React from 'react'; // react 18.2+
import { useState, useEffect, useCallback } from 'react'; // react 18.2+
import {
  Box,
  Typography,
  Divider,
  Avatar,
  Badge,
  CircularProgress,
  List,
  ListItem,
  ListItemAvatar,
  ListItemText,
  ListItemButton,
} from '@mui/material'; // @mui/material 5.13+
import { styled } from '@mui/material/styles'; // @mui/material/styles 5.13+
import { Search } from '@mui/icons-material'; // @mui/icons-material 5.13+

// Internal imports
import useCommunications from '../hooks/useCommunications';
import Card from '../../../components/common/Card';
import { formatRelativeTime, formatDate } from '../../../utils/formatters';
import { MessageListProps } from '../../../types';

// Styled components for consistent styling
const ListContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: '100%',
  overflow: 'hidden',
});

const SearchContainer = styled(Box)(({ theme }) => ({
  padding: theme.spacing(2),
  borderBottom: `1px solid ${theme.palette.divider}`,
}));

const ConversationList = styled(List)({
  flex: '1',
  overflow: 'auto',
  padding: 0,
});

const ConversationItem = styled(ListItemButton)(({ theme }) => ({
  borderBottom: `1px solid ${theme.palette.divider}`,
  padding: theme.spacing(1.5, 2),
  '&.selected': {
    backgroundColor: theme.palette.action.selected,
  },
}));

const EmptyState = styled(Box)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  height: '100%',
  padding: theme.spacing(4),
  color: theme.palette.text.secondary,
}));

const MessagePreview = styled(Typography)({
  color: 'text.secondary',
  whiteSpace: 'nowrap',
  overflow: 'hidden',
  textOverflow: 'ellipsis',
  maxWidth: '100%',
  fontSize: '0.875rem',
});

const TimeStamp = styled(Typography)({
  color: 'text.secondary',
  fontSize: '0.75rem',
  marginLeft: 'auto',
  whiteSpace: 'nowrap',
});

/**
 * Component that displays a list of conversations or messages
 */
const MessageList: React.FC<MessageListProps> = (props) => {
  // Extract props including onSelectMessage, className, and mode
  const { onSelectMessage, className = '', mode = 'conversations' } = props;

  // Get communications data and functions from useCommunications hook
  const { conversations, currentConversation, setCurrentConversation, markAsRead, loading } = useCommunications();

  // Set up state for search query and loading state
  const [searchQuery, setSearchQuery] = useState('');

  /**
   * Handles selecting a conversation from the list
   * @param conversation
   */
  const handleConversationSelect = useCallback((conversation: any) => {
    // Set the selected conversation using setCurrentConversation
    setCurrentConversation(conversation);

    // If onSelectMessage callback is provided, call it with the conversation
    if (onSelectMessage) {
      onSelectMessage(conversation);
    }

    // If conversation has unread messages, mark them as read
    // This is now handled by useEffect
  }, [setCurrentConversation, onSelectMessage]);

  /**
   * Handles changes to the search input
   * @param e
   */
  const handleSearchChange = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    // Update searchQuery state with the input value
    setSearchQuery(e.target.value);
  }, []);

  /**
   * Filters conversations based on search query
   */
  const filterConversations = useCallback(() => {
    // If search query is empty, return all conversations
    if (!searchQuery) {
      return conversations || [];
    }

    // Filter conversations by matching subject or participant names with search query
    const filtered = (conversations || []).filter((conversation: any) => {
      const subject = conversation.subject || '';
      const participants = conversation.participants || [];
      const participantNames = participants.map((p: any) => p.name).join(' ');
      const searchString = `${subject} ${participantNames}`.toLowerCase();
      return searchString.includes(searchQuery.toLowerCase());
    });

    // Return filtered conversations array
    return filtered;
  }, [conversations, searchQuery]);

  // Use useEffect to mark messages as read when a conversation is selected
  useEffect(() => {
    if (currentConversation && currentConversation.unread_count > 0) {
      const messageIds = currentConversation.messages.filter((message: any) => !message.read).map((message: any) => message.id);
      if (messageIds.length > 0) {
        markAsRead({ conversationId: currentConversation.id, messageIds });
      }
    }
  }, [currentConversation, markAsRead]);

  // Render the component
  return (
    <ListContainer className={className}>
      {/* Search input */}
      <SearchContainer>
        <Box sx={{ display: 'flex', alignItems: 'center', border: '1px solid #ccc', borderRadius: '4px', padding: '4px' }}>
          <Search sx={{ mr: 1, color: 'text.secondary' }} />
          <input
            type="text"
            placeholder="Search conversations..."
            value={searchQuery}
            onChange={handleSearchChange}
            style={{ border: 'none', outline: 'none', flex: 1, padding: '4px' }}
          />
        </Box>
      </SearchContainer>

      {/* Conversation list */}
      <ConversationList>
        {loading ? (
          // Display loading indicator when loading conversations
          <Box sx={{ display: 'flex', justifyContent: 'center', alignItems: 'center', height: '100%' }}>
            <CircularProgress />
          </Box>
        ) : filterConversations().length === 0 ? (
          // Display empty state when no conversations exist
          <EmptyState>
            <Typography variant="subtitle1">No conversations found.</Typography>
          </EmptyState>
        ) : (
          // Render each conversation as a ListItem with appropriate styling
          filterConversations().map((conversation: any) => (
            <ConversationItem
              key={conversation.id}
              selected={currentConversation?.id === conversation.id}
              onClick={() => handleConversationSelect(conversation)}
            >
              <ListItemAvatar>
                <Badge
                  color="primary"
                  badgeContent={conversation.unread_count}
                  max={99}
                  invisible={conversation.unread_count === 0}
                  overlap="circular"
                >
                  <Avatar>{conversation.participants[0].name.charAt(0)}</Avatar>
                </Badge>
              </ListItemAvatar>
              <ListItemText
                primary={conversation.subject}
                secondary={
                  <MessagePreview>
                    {conversation.messages[0]?.text}
                  </MessagePreview>
                }
              />
              <TimeStamp>
                {formatRelativeTime(conversation.messages[0]?.created_at)}
              </TimeStamp>
            </ConversationItem>
          ))
        )}
      </ConversationList>
    </ListContainer>
  );
};

export default MessageList;