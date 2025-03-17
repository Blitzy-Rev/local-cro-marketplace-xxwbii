import { useState, useEffect, useCallback, useRef } from 'react'; // react 18.2+
import { useQuery, useMutation, useQueryClient } from 'react-query'; // react-query 4.28+
import {
  apiClient,
  handleApiError,
} from '../../../api/client';
import { useToast } from '../../../hooks/useToast';
import {
  Conversation,
  Message,
  User,
  Attachment,
  CreateConversationData,
} from '../../../types';

const BASE_URL = '/communications';

/**
 * Custom hook for managing communications between users
 */
const useCommunications = () => {
  const queryClient = useQueryClient();
  const { showToast } = useToast();
  const [currentConversation, setCurrentConversation] = useState<Conversation | null>(null);

  const {
    data: conversations,
    isLoading: conversationsLoading,
    error: conversationsError,
    refetch: fetchConversations,
  } = useQuery<Conversation[]>('conversations', () => fetchConversations());

  const {
    data: messagesData,
    isLoading: messagesLoading,
    error: messagesError,
    refetch: fetchMessages,
  } = useQuery<{ messages: Message[]; total: number; page: number; pages: number }>(
    ['messages', currentConversation?.id],
    () => fetchMessages(currentConversation?.id),
    {
      enabled: !!currentConversation?.id,
    }
  );

  const messages = messagesData?.messages || [];

  const { mutate: sendMessage, isLoading: sendingMessage } = useMutation<Message, any, { conversationId: number; text: string; attachments?: File[] }>(
    ({ conversationId, text, attachments = [] }) => sendMessage({ conversationId, text, attachments }),
    {
      onSuccess: () => {
        queryClient.invalidateQueries(['messages', currentConversation?.id]);
      },
      onError: (error) => {
        showToast({ type: 'error', message: `Failed to send message: ${error.message}` });
      },
    }
  );

  const { mutate: createConversation, isLoading: creatingConversation } = useMutation<Conversation, any, CreateConversationData>(
    (data) => createConversation(data),
    {
      onSuccess: () => {
        queryClient.invalidateQueries('conversations');
      },
      onError: (error) => {
        showToast({ type: 'error', message: `Failed to create conversation: ${error.message}` });
      },
    }
  );

  const { mutate: markAsRead, isLoading: markingAsRead } = useMutation<
    { success: boolean; count: number },
    any,
    { conversationId: number; messageIds: number[] }
  >(
    ({ conversationId, messageIds }) => markAsRead({ conversationId, messageIds }),
    {
      onSuccess: () => {
        queryClient.invalidateQueries(['messages', currentConversation?.id]);
        queryClient.invalidateQueries('conversations');
      },
      onError: (error) => {
        showToast({ type: 'error', message: `Failed to mark as read: ${error.message}` });
      },
    }
  );

  /**
   * Fetches all conversations for the current user
   */
  const fetchConversations = useCallback(async (): Promise<Conversation[]> => {
    try {
      const response = await apiClient.get<Conversation[]>(`${BASE_URL}`);
      return response.data;
    } catch (error: any) {
      handleApiError(error);
      showToast({ type: 'error', message: `Failed to fetch conversations: ${error.message}` });
      return [];
    }
  }, [showToast]);

  /**
   * Fetches messages for a specific conversation
   * @param conversationId
   * @param page
   * @param pageSize
   */
  const fetchMessages = useCallback(
    async (conversationId: number | undefined, page: number = 1, pageSize: number = 20): Promise<{ messages: Message[]; total: number; page: number; pages: number }> => {
      if (!conversationId) {
        return { messages: [], total: 0, page: 1, pages: 1 };
      }

      try {
        const params = new URLSearchParams();
        params.append('page', page.toString());
        params.append('pageSize', pageSize.toString());

        const response = await apiClient.get<{ messages: Message[]; total: number; page: number; pages: number }>(
          `${BASE_URL}/${conversationId}/messages?${params.toString()}`
        );
        return response.data;
      } catch (error: any) {
        handleApiError(error);
        showToast({ type: 'error', message: `Failed to fetch messages: ${error.message}` });
        return { messages: [], total: 0, page: 1, pages: 1 };
      }
    },
    [showToast]
  );

  /**
   * Sends a new message in a conversation with optional file attachments
   * @param conversationId
   * @param text
   * @param attachments
   */
  const sendMessage = useCallback(
    async ({ conversationId, text, attachments = [] }: { conversationId: number; text: string; attachments?: File[] }): Promise<Message> => {
      if (!conversationId || !text) {
        throw new Error('Conversation ID and text are required');
      }

      try {
        const formData = new FormData();
        formData.append('text', text);

        attachments.forEach((file) => {
          formData.append('attachments', file);
        });

        const response = await apiClient.post<Message>(`${BASE_URL}/${conversationId}/messages`, formData, {
          headers: {
            'Content-Type': 'multipart/form-data',
          },
        });

        queryClient.invalidateQueries(['messages', conversationId]);
        return response.data;
      } catch (error: any) {
        handleApiError(error);
        showToast({ type: 'error', message: `Failed to send message: ${error.message}` });
        throw error;
      }
    },
    [showToast, queryClient]
  );

  /**
   * Creates a new conversation with another user
   * @param data
   */
  const createConversation = useCallback(
    async (data: CreateConversationData): Promise<Conversation> => {
      if (!data.recipientId && !data.experimentId && !data.submissionId) {
        throw new Error('Recipient ID, experiment ID, or submission ID is required');
      }

      if (!data.subject) {
        throw new Error('Subject is required');
      }

      try {
        const response = await apiClient.post<Conversation>(`${BASE_URL}`, data);
        queryClient.invalidateQueries('conversations');
        return response.data;
      } catch (error: any) {
        handleApiError(error);
        showToast({ type: 'error', message: `Failed to create conversation: ${error.message}` });
        throw error;
      }
    },
    [showToast, queryClient]
  );

  /**
   * Marks messages in a conversation as read
   * @param conversationId
   * @param messageIds
   */
  const markAsRead = useCallback(
    async ({ conversationId, messageIds }: { conversationId: number; messageIds: number[] }): Promise<{ success: boolean; count: number }> => {
      if (!conversationId || !messageIds || messageIds.length === 0) {
        throw new Error('Conversation ID and message IDs are required');
      }

      try {
        const response = await apiClient.put<{ success: boolean; count: number }>(`${BASE_URL}/${conversationId}/messages/read`, {
          messageIds,
        });
        queryClient.invalidateQueries(['messages', conversationId]);
        queryClient.invalidateQueries('conversations');
        return response.data;
      } catch (error: any) {
        handleApiError(error);
        showToast({ type: 'error', message: `Failed to mark as read: ${error.message}` });
        throw error;
      }
    },
    [showToast, queryClient]
  );

  return {
    conversations,
    currentConversation,
    messages,
    loading: conversationsLoading || messagesLoading || sendingMessage || creatingConversation || markingAsRead,
    error: conversationsError || messagesError,
    fetchConversations,
    fetchMessages,
    sendMessage,
    createConversation,
    markAsRead,
    setCurrentConversation,
  };
};

export default useCommunications;