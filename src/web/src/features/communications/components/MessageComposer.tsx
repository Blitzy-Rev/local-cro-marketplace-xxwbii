import React, { useState, useRef, useCallback, useEffect } from 'react'; // React 18.2+
import { Box, TextField, IconButton, Tooltip, Paper, Collapse } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { Send, AttachFile, Close } from '@mui/icons-material'; // v5.13+
import useCommunications from '../hooks/useCommunications';
import Button from '../../../components/common/Button';
import FileUpload from '../../../components/common/FileUpload';
import { useToast } from '../../../hooks/useToast';

/**
 * Interface defining the props for the MessageComposer component
 */
interface MessageComposerProps {
  /** Callback function when a message is successfully sent */
  onMessageSent?: () => void;
  /** Additional CSS class for styling */
  className?: string;
    /** Placeholder text for the message input */
  placeholder?: string;
  /** Whether the composer is disabled */
  disabled?: boolean;
}

/**
 * Component for composing and sending messages
 * @param props - MessageComposerProps
 * @returns JSX.Element
 */
const MessageComposer: React.FC<MessageComposerProps> = (props) => {
  // Extract props including onMessageSent and className
  const { onMessageSent, className, placeholder, disabled } = props;

  // Get communications data and functions from useCommunications hook
  const { currentConversation, sendMessage, loading: communicationsLoading } = useCommunications();

  // Set up state for message text, attachments, and loading state
  const [messageText, setMessageText] = useState<string>('');
  const [attachments, setAttachments] = useState<File[]>([]);
  const [loading, setLoading] = useState<boolean>(false);

  // Create a ref for the text input element
  const textInputRef = useRef<HTMLInputElement>(null);

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  /**
   * Handles changes to the message text input
   * @param e - React.ChangeEvent<HTMLInputElement>
   */
  const handleTextChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    // Update messageText state with the input value
    setMessageText(e.target.value);
  };

  /**
   * Handles file selection for attachments
   * @param files - File[]
   */
  const handleFilesSelected = (files: File[]) => {
    // Update attachments state with the selected files
    setAttachments(files);
  };

  /**
   * Handles removing a file from attachments
   * @param index - number
   */
  const handleFileRemove = (index: number) => {
    // Create a copy of the current attachments array
    const newAttachments = [...attachments];
    // Remove the file at the specified index
    newAttachments.splice(index, 1);
    // Update attachments state with the modified array
    setAttachments(newAttachments);
  };

  /**
   * Handles sending the composed message
   */
  const handleSendMessage = useCallback(async () => {
    // Check if there is a current conversation
    if (!currentConversation) {
      showToast({ type: 'warning', message: 'No conversation selected.' });
      return;
    }

    // Check if message text is not empty or there are attachments
    if (messageText.trim() === '' && attachments.length === 0) {
      showToast({ type: 'warning', message: 'Message text cannot be empty.' });
      return;
    }

    try {
      // Set loading state to true
      setLoading(true);

      // Call sendMessage with conversation ID, message text, and attachments
      await sendMessage({
        conversationId: currentConversation.id,
        text: messageText,
        attachments: attachments,
      });

      // Clear message text and attachments after successful send
      setMessageText('');
      setAttachments([]);

      // Call onMessageSent callback if provided
      if (onMessageSent) {
        onMessageSent();
      }

      // Focus the text input for the next message
      if (textInputRef.current) {
        textInputRef.current.focus();
      }
    } catch (error: any) {
      // Handle errors with appropriate error messages
      showToast({ type: 'error', message: `Failed to send message: ${error.message}` });
    } finally {
      // Set loading state to false
      setLoading(false);
    }
  }, [currentConversation, messageText, attachments, sendMessage, onMessageSent, showToast]);

  /**
   * Handles keyboard shortcuts for sending messages
   * @param e - React.KeyboardEvent<HTMLDivElement>
   */
  const handleKeyDown = (e: React.KeyboardEvent<HTMLDivElement>) => {
    // Check if Ctrl+Enter or Cmd+Enter was pressed
    if ((e.ctrlKey || e.metaKey) && e.key === 'Enter') {
      // If the shortcut was pressed, prevent default behavior
      e.preventDefault();
      // Call handleSendMessage to send the message
      handleSendMessage();
    }
  };

  // Use useEffect to focus the text input when the component mounts
  useEffect(() => {
    if (textInputRef.current) {
      textInputRef.current.focus();
    }
  }, []);

  // Render a container with text input and action buttons
  return (
    <Box className={className}>
      {/* Render TextField for message input with appropriate event handlers */}
      <TextField
        fullWidth
        multiline
        rows={3}
        placeholder={placeholder || "Type your message here..."}
        variant="outlined"
        value={messageText}
        onChange={handleTextChange}
        onKeyDown={handleKeyDown}
        inputRef={textInputRef}
        disabled={disabled || communicationsLoading}
      />

      <Box mt={1} display="flex" alignItems="center">
        {/* Render attachment button with FileUpload component */}
        <FileUpload
          multiple
          disabled={disabled || communicationsLoading}
          showPreview={false}
          onFilesSelected={handleFilesSelected}
          onFileValidationFail={(error) => showToast({ type: 'error', message: error })}
          label=""
          buttonText={<Tooltip title="Attach File"><IconButton color="primary" component="span" disabled={disabled || communicationsLoading}><AttachFile /></IconButton></Tooltip>}
        />

        <Box flexGrow={1} />

        {/* Render send button with loading state */}
        <Button
          variant="contained"
          color="primary"
          onClick={handleSendMessage}
          disabled={disabled || loading || communicationsLoading}
          startIcon={<Send />}
        >
          Send
        </Button>
      </Box>

      {/* If attachments exist, render attachment preview with remove option */}
      <Collapse in={attachments.length > 0} timeout="auto" unmountOnExit>
        <Paper elevation={1} sx={{ mt: 2, p: 1 }}>
          {attachments.map((file, index) => (
            <Box key={index} display="flex" alignItems="center" justifyContent="space-between" sx={{ mb: 1 }}>
              <Box display="flex" alignItems="center">
                <AttachFile sx={{ mr: 1 }} />
                <Typography variant="body2">{file.name}</Typography>
              </Box>
              <IconButton onClick={() => handleFileRemove(index)} aria-label="remove attachment">
                <Close />
              </IconButton>
            </Box>
          ))}
        </Paper>
      </Collapse>
    </Box>
  );
};

export default MessageComposer;