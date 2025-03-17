import React from 'react'; // react ^18.2.0
import { MemoryRouter, Routes, Route } from 'react-router-dom'; // react-router-dom ^6.11.2
import { act } from '@testing-library/react'; // @testing-library/react ^14.0.0
import { vi } from 'vitest'; // vitest ^0.31.0

import ResultDetail from '../components/ResultDetail';
import { customRender, screen, waitFor, fireEvent, userEvent, setupApiMock, rest } from '../../../__tests__/utils/test-utils';
import { server } from '../../../__tests__/mocks/server';
import { ResultStatus } from '../../../types/result';
import useResults from '../hooks/useResults';

/**
 * Creates a mock result detail object for testing
 * @param overrides - Optional overrides for the mock result
 * @returns Mock result detail object
 */
const mockResultDetail = (overrides = {}) => ({
  id: '1',
  submission_id: '1',
  status: ResultStatus.UPLOADED,
  uploaded_at: '2023-07-20T10:00:00Z',
  data_points: [],
  files: [],
  ...overrides,
});

/**
 * Sets up mocks for the useResults hook and API endpoints
 * @param options - Optional overrides for the mock functions
 * @returns Mock functions and handlers
 */
const setupMocks = (options = {}) => {
  const fetchResultDetailMock = vi.fn().mockResolvedValue(mockResultDetail());
  const approveResultMock = vi.fn().mockResolvedValue({});
  const rejectResultMock = vi.fn().mockResolvedValue({});

  vi.spyOn(useResults, 'default').mockReturnValue({
    fetchResultDetail: fetchResultDetailMock,
    approveResult: approveResultMock,
    rejectResult: rejectResultMock,
    exportResult: vi.fn(),
    analyzeResult: vi.fn(),
    requestAdditionalData: vi.fn(),
    currentResult: mockResultDetail(),
    loading: false,
    error: null,
    clearError: vi.fn(),
  } as any);

  const handlers = [
    rest.get('/api/v1/results/1', (req, res, ctx) => {
      return res(ctx.status(200), ctx.json(mockResultDetail()));
    }),
    rest.post('/api/v1/results/1/approve', (req, res, ctx) => {
      return res(ctx.status(200), ctx.json({}));
    }),
    rest.post('/api/v1/results/1/reject', (req, res, ctx) => {
      return res(ctx.status(200), ctx.json({}));
    }),
  ];

  setupApiMock(handlers);

  return {
    fetchResultDetailMock,
    approveResultMock,
    rejectResultMock,
    handlers
  };
};

/**
 * Helper function to render the ResultDetail component with proper routing
 * @param props - Optional props for the ResultDetail component
 * @returns Rendered component and utilities
 */
const renderResultDetail = (props = {}) => {
  return customRender(
    <MemoryRouter initialEntries={['/results/1']}>
      <Routes>
        <Route path="/results/:id" element={<ResultDetail {...props} />} />
      </Routes>
    </MemoryRouter>
  );
};

describe('ResultDetail Component', () => {
  it('renders loading state initially', () => {
    vi.spyOn(useResults, 'default').mockReturnValue({
      fetchResultDetail: vi.fn(),
      approveResult: vi.fn(),
      rejectResult: vi.fn(),
      exportResult: vi.fn(),
      analyzeResult: vi.fn(),
      requestAdditionalData: vi.fn(),
      currentResult: null,
      loading: true,
      error: null,
      clearError: vi.fn(),
    } as any);
    renderResultDetail();
    expect(screen.getByRole('progressbar')).toBeInTheDocument();
  });

  it('renders error state when fetch fails', async () => {
    vi.spyOn(useResults, 'default').mockReturnValue({
      fetchResultDetail: vi.fn().mockRejectedValue(new Error('Fetch failed')),
      approveResult: vi.fn(),
      rejectResult: vi.fn(),
      exportResult: vi.fn(),
      analyzeResult: vi.fn(),
      requestAdditionalData: vi.fn(),
      currentResult: null,
      loading: false,
      error: 'Fetch failed',
      clearError: vi.fn(),
    } as any);
    renderResultDetail();
    await waitFor(() => {
      expect(screen.getByText('Fetch failed')).toBeInTheDocument();
    });
  });

  it('renders result details when fetch succeeds', async () => {
    const { fetchResultDetailMock } = setupMocks();
    renderResultDetail();
    await waitFor(() => {
      expect(fetchResultDetailMock).toHaveBeenCalled();
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
      expect(screen.getByText('Files')).toBeInTheDocument();
      expect(screen.getByText('Molecule Results')).toBeInTheDocument();
    });
  });

  it('allows switching between tabs', async () => {
    setupMocks();
    renderResultDetail();
    await waitFor(() => {
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
    });

    await act(async () => {
      fireEvent.click(screen.getByText('Molecule Results'));
    });
    expect(screen.getByText('Molecule Results')).toBeInTheDocument();

    await act(async () => {
      fireEvent.click(screen.getByText('Notes'));
    });
    expect(screen.getByText('Notes')).toBeInTheDocument();
  });

  it('handles approving a result', async () => {
    const { approveResultMock, fetchResultDetailMock } = setupMocks();
    renderResultDetail();
    await waitFor(() => {
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
    });

    await act(async () => {
      fireEvent.click(screen.getByText('Approve'));
    });

    const dialog = screen.getByRole('dialog');
    const confirmButton = within(dialog).getByRole('button', { name: 'Approve' });

    await act(async () => {
      fireEvent.click(confirmButton);
    });

    expect(approveResultMock).toHaveBeenCalledWith('1', '');
    expect(fetchResultDetailMock).toHaveBeenCalledTimes(2);
  });

  it('handles rejecting a result', async () => {
    const { rejectResultMock, fetchResultDetailMock } = setupMocks();
    renderResultDetail();
    await waitFor(() => {
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
    });

    await act(async () => {
      fireEvent.click(screen.getByText('Reject'));
    });

    const dialog = screen.getByRole('dialog');
    const confirmButton = within(dialog).getByRole('button', { name: 'Reject' });

    await act(async () => {
      fireEvent.click(confirmButton);
    });

    expect(rejectResultMock).toHaveBeenCalledWith('1', '');
    expect(fetchResultDetailMock).toHaveBeenCalledTimes(2);
  });

  it('disables approve/reject buttons for already processed results', async () => {
    setupMocks({ status: ResultStatus.APPROVED });
    renderResultDetail();
    await waitFor(() => {
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
    });

    expect(screen.getByText('Approve')).toBeDisabled();
    expect(screen.getByText('Reject')).toBeDisabled();
  });

  it('navigates back when back button is clicked', async () => {
    const onBackMock = vi.fn();
    setupMocks();
    renderResultDetail({ onBack: onBackMock });
    await waitFor(() => {
      expect(screen.getByText('Result Summary')).toBeInTheDocument();
    });

    await act(async () => {
      fireEvent.click(screen.getByText('Back to Results'));
    });

    expect(onBackMock).toHaveBeenCalled();
  });
});