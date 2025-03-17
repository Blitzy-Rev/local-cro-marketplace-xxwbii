import React, { useState, useEffect, useRef, useCallback } from 'react'; // React 18.2+
import { Box, Typography, Divider, Avatar, CircularProgress, IconButton, Tooltip, Paper } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { Download, ArrowBack } from '@mui/icons-material'; // v5.13+
import useCommunications from '../hooks/useCommunications';
import MessageComposer from './MessageComposer';
import Card from '../../../components/common/Card';
import { formatDate, formatRelativeTime } from '../../../utils/formatters';
import { MessageDetailProps, MessageGroup } from '../../../types/notification';

/**
 * Styled component for the main detail container
 */
const DetailContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: '100%',
  overflow: 'hidden',
});

/**
 * Styled component for the header section
 */
const Header = styled(Box)(({ theme }) => ({
  padding: theme.spacing(2),
  borderBottom: `1px solid ${theme.palette.divider}`,
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'space-between',
}));

/**
 * Styled component for the messages container
 */
const MessagesContainer = styled(Box)({
  flex: 1,
  overflow: 'auto',
  padding: '16px',
  display: 'flex',
  flexDirection: 'column',
});

/**
 * Styled component for individual message bubbles
 */
const MessageBubble = styled(Paper, {
  shouldForwardProp: (prop) => prop !== 'isCurrentUser',
})<{ isCurrentUser: boolean }>(({ theme, isCurrentUser }) => ({
  padding: theme.spacing(1.5),
  borderRadius: theme.spacing(1.5),
  maxWidth: '70%',
  wordBreak: 'break-word',
  marginBottom: theme.spacing(1),
  backgroundColor: isCurrentUser ? theme.palette.primary.light : theme.palette.background.paper,
  color: isCurrentUser ? theme.palette.primary.contrastText : theme.palette.text.primary,
  alignSelf: isCurrentUser ? 'flex-end' : 'flex-start',
}));

/**
 * Styled component for grouping messages by date
 */
const MessageGroupStyle = styled(Box)({
  marginBottom: '16px',
  width: '100%',
});

/**
 * Styled component for date separators
 */
const DateSeparator = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  margin: '16px 0',
  width: '100%',
});

/**
 * Styled component for the date line
 */
const DateLine = styled(Divider)({
  flex: 1,
});

/**
 * Styled component for the date text
 */
const DateText = styled(Typography)(({ theme }) => ({
  margin: `0 ${theme.spacing(2)}`,
  color: theme.palette.text.secondary,
  fontSize: '0.75rem',
}));

/**
 * Styled component for the message header
 */
const MessageHeader = styled(Box)({
  display: 'flex',
  alignItems: 'center',
  marginBottom: '8px',
});

/**
 * Styled component for the sender's name
 */
const SenderName = styled(Typography)({
  fontWeight: 500,
  marginLeft: '8px',
});

/**
 * Styled component for the timestamp
 */
const TimeStamp = styled(Typography)(({ theme }) => ({
  color: theme.palette.text.secondary,
  fontSize: '0.75rem',
  marginLeft: '8px',
}));

/**
 * Styled component for the attachments container
 */
const AttachmentsContainer = styled(Box)({
  marginTop: '8px',
  display: 'flex',
  flexWrap: 'wrap',
  gap: '8px',
});

/**
 * Styled component for individual attachment items
 */
const AttachmentItem = styled(Box)(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  backgroundColor: theme.palette.action.hover,
  borderRadius: '4px',
  padding: '4px 8px',
  fontSize: '0.875rem',
}));

/**
 * Styled component for the empty state
 */
const EmptyState = styled(Box)(({ theme }) => ({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  height: '100%',
  padding: theme.spacing(4),
  color: theme.palette.text.secondary,
}));

/**
 * Styled component for the loading container
 */
const LoadingContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'center',
  padding: '16px',
});

/**
 * Interface defining the props for the MessageDetail component
 */
interface MessageDetailProps {
  /** Callback function when the back button is clicked */
  onBackClick?: () => void;
  /** Additional CSS class for styling */
  className?: string;
  /** Whether to show the back button (for mobile view) */
  showBackButton?: boolean;
}

/**
 * Component that displays the details of a selected conversation and its messages
 * @param props - MessageDetailProps
 * @returns JSX.Element
 */
const MessageDetail: React.FC<MessageDetailProps> = (props) => {
  // Extract props including onBackClick and className
  const { onBackClick, className = '', showBackButton = false } = props;

  // Get communications data and functions from useCommunications hook
  const { currentConversation, messages, loading, fetchMessages, markAsRead } = useCommunications();

  // Set up state for pagination and scroll position
  const [page, setPage] = useState<number>(1);
  const [hasMore, setHasMore] = useState<boolean>(true);
  const [loadingMore, setLoadingMore] = useState<boolean>(false);
  const [scrollPosition, setScrollPosition] = useState<number>(0);

  // Create a ref for the message container element
  const messageContainerRef = useRef<HTMLDivElement>(null);

  /**
   * Handles loading more messages when scrolling up
   */
  const handleLoadMore = useCallback(async () => {
    // Check if there are more messages to load
    if (!hasMore || loadingMore || !currentConversation) {
      return;
    }

    // Set loading state to true
    setLoadingMore(true);

    // Increment page number
    const newPage = page + 1;

    // Call fetchMessages with conversation ID and new page number
    try {
      const result = await fetchMessages(currentConversation.id, newPage);
      if (result && result.messages.length > 0) {
        setPage(newPage);
        setHasMore(result.page < result.pages);
        // Save current scroll position before new messages are added
        if (messageContainerRef.current) {
          setScrollPosition(messageContainerRef.current.scrollHeight);
        }
      } else {
        setHasMore(false);
      }
    } finally {
      // Set loading state to false after messages are loaded
      setLoadingMore(false);
    }
  }, [hasMore, loadingMore, page, currentConversation, fetchMessages]);

  /**
   * Handles downloading file attachments
   * @param fileUrl - string
   * @param fileName - string
   */
  const handleDownload = useCallback((fileUrl: string, fileName: string) => {
    // Create a temporary anchor element
    const link = document.createElement('a');

    // Set the href attribute to the file URL
    link.href = fileUrl;

    // Set the download attribute to the file name
    link.download = fileName;

    // Set the target attribute to '_blank'
    link.target = '_blank';

    // Programmatically click the anchor element
    link.click();

    // Remove the anchor element from the DOM
    link.remove();
  }, []);

  /**
   * Scrolls the message container to the bottom
   * @param smooth - boolean
   */
  const scrollToBottom = useCallback((smooth: boolean = true) => {
    // Check if messageContainerRef.current exists
    if (messageContainerRef.current) {
      // Get the scrollHeight of the container
      const { scrollHeight } = messageContainerRef.current;

      // Scroll to the bottom with smooth scrolling if specified
      messageContainerRef.current.scrollTo({
        top: scrollHeight,
        behavior: smooth ? 'smooth' : 'instant',
      });
    }
  }, []);

  /**
   * Handles events after a message is successfully sent
   */
  const handleMessageSent = useCallback(() => {
    // Scroll the message container to the bottom
    scrollToBottom();
  }, [scrollToBottom]);

  /**
   * Handles scroll events to detect when to load more messages
   * @param e - React.UIEvent<HTMLDivElement>
   */
  const handleScroll = useCallback((e: React.UIEvent<HTMLDivElement>) => {
    // Get the scroll position from the event target
    const { scrollTop } = e.target as HTMLDivElement;

    // Check if the user has scrolled near the top of the container
    if (scrollTop <= 50 && !loadingMore) {
      // If near the top and not currently loading, call handleLoadMore
      handleLoadMore();
    }
  }, [loadingMore, handleLoadMore]);

  useEffect(() => {
    if (messages && messageContainerRef.current) {
      scrollToBottom(false);
    }
  }, [messages, scrollToBottom]);

  useEffect(() => {
    if (currentConversation && messages && messages.length > 0) {
      const unreadMessageIds = messages.filter(message => !message.read).map(message => message.id);
      if (unreadMessageIds.length > 0) {
        markAsRead({ conversationId: currentConversation.id, messageIds: unreadMessageIds });
      }
    }
  }, [currentConversation, messages, markAsRead]);

  useEffect(() => {
    if (currentConversation) {
      setPage(1);
      setHasMore(true);
      fetchMessages(currentConversation.id, 1);
    }
  }, [currentConversation, fetchMessages]);

  // Render a container with conversation details and messages
  return (
    <DetailContainer className={className}>
      {/* If showBackButton is true, render back button in Header */}
      {showBackButton && (
        <Header>
          <IconButton onClick={onBackClick} aria-label="back">
            <ArrowBack />
          </IconButton>
        </Header>
      )}

      {/* If no currentConversation, render EmptyState with message */}
      {!currentConversation ? (
        <EmptyState>
          <Typography variant="subtitle1">No conversation selected.</Typography>
        </EmptyState>
      ) : (
        <>
          {/* Render Header with conversation subject and participants */}
          <Header>
            <Typography variant="h6">{currentConversation.subject}</Typography>
            <Typography variant="body2">Participants: {currentConversation.participants?.length}</Typography>
          </Header>

          {/* Render MessagesContainer with ref and onScroll handler */}
          <MessagesContainer ref={messageContainerRef} onScroll={handleScroll}>
            {/* If loadingMore, render LoadingContainer with CircularProgress at top */}
            {loadingMore && (
              <LoadingContainer>
                <CircularProgress size={20} />
              </LoadingContainer>
            )}

            {/* Group messages by date using formatDate */}
            {messages.reduce<MessageGroup[]>((groups, message) => {
              const date = formatDate(message.created_at);
              const existingGroup = groups.find((group) => group.date === date);

              if (existingGroup) {
                existingGroup.messages.push(message);
              } else {
                groups.push({ date, messages: [message] });
              }

              return groups;
            }, []).map((group) => (
              <MessageGroupStyle key={group.date}>
                {/* Render DateSeparator with formatted date */}
                <DateSeparator>
                  <DateLine />
                  <DateText>{formatDate(group.date)}</DateText>
                  <DateLine />
                </DateSeparator>

                {/* For each message in group, render MessageBubble with appropriate styling */}
                {group.messages.map((message) => (
                  <MessageBubble key={message.id} isCurrentUser={message.sender?.id === currentConversation.current_user_id}>
                    {/* Determine if message is from current user for alignment */}
                    <MessageHeader>
                      {/* Render message header with sender Avatar and name */}
                      <Avatar alt={message.sender?.email} src={message.sender?.image_url} />
                      <SenderName>{message.sender?.email}</SenderName>
                      <TimeStamp>{formatRelativeTime(message.created_at)}</TimeStamp>
                    </MessageHeader>

                    {/* Render message text content */}
                    <Typography variant="body2">{message.text}</Typography>

                    {/* If message has attachments, render AttachmentsContainer */}
                    {message.attachments && message.attachments.length > 0 && (
                      <AttachmentsContainer>
                        {/* For each attachment, render AttachmentItem with download button */}
                        {message.attachments.map((attachment) => (
                          <AttachmentItem key={attachment.id}>
                            <Typography variant="body2">{attachment.file_name}</Typography>
                            <Tooltip title="Download">
                              <IconButton
                                color="primary"
                                onClick={() => handleDownload(attachment.file_path, attachment.file_name)}
                                aria-label="download attachment"
                              >
                                <Download />
                              </IconButton>
                            </Tooltip>
                          </AttachmentItem>
                        ))}
                      </AttachmentsContainer>
                    )}
                  </MessageBubble>
                ))}
              </MessageGroupStyle>
            ))}
          </MessagesContainer>

          {/* Render MessageComposer at bottom for replying */}
          <MessageComposer onMessageSent={handleMessageSent} placeholder="Reply to message..." disabled={loading} />
        </>
      )}
    </DetailContainer>
  );
};

export default MessageDetail;