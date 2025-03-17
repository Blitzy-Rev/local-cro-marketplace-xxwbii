# src/web/src/features/submissions/__tests__/SubmissionDetail.test.tsx
```typescript
import React from 'react'; // React core library for component testing // IE2
import { render, screen, fireEvent } from '@testing-library/react'; // React Testing Library for component testing // IE2
import { vi } from 'vitest'; // Vitest for mocking and assertions // IE2
import SubmissionDetail from '../components/SubmissionDetail';
import { customRender, screen as customScreen, waitFor, fireEvent as customFireEvent, userEvent, within } from '../../../__tests__/utils/test-utils';
import { server, rest } from '../../../__tests__/mocks/server';
import useSubmissions from '../hooks/useSubmissions';
import useAuth from '../../auth/hooks/useAuth';
import { SubmissionStatus } from '../../../types/submission';

// Mock the useAuth hook
vi.mock('../../auth/hooks/useAuth');

// Mock the useSubmissions hook
vi.mock('../hooks/useSubmissions');

// Helper function to create mock submission data for tests
const mockSubmissionData = (overrides: { status?: SubmissionStatus; price?: number; turnaround_days?: number } = {}) => {
  const baseSubmission = {
    id: '1',
    experiment_id: '1',
    cro_id: '2',
    status: SubmissionStatus.PENDING,
    submitted_at: '2023-06-01T10:00:00Z',
    updated_at: '2023-06-01T10:00:00Z',
    price: null,
    turnaround_days: null,
    notes: 'Test submission',
    experiment: {
      id: '1',
      name: 'Test Experiment',
      type: 'Binding Assay',
      status: 'Submitted',
      created_at: '2023-05-30T10:00:00Z',
    },
    cro: {
      id: '2',
      email: 'cro@example.com',
      role: 'cro',
    },
    details: [],
    results: [],
  };

  return { ...baseSubmission, ...overrides };
};

// Setup and teardown for the test suite
beforeAll(() => server.listen());
afterEach(() => server.resetHandlers());
afterAll(() => server.close());

describe('SubmissionDetail Component', () => {
  // Helper function to set up the component for testing with mocks
  const setup = (options: { submissionId: string; submission?: any; userRole?: string; isAuthenticated?: boolean } = { submissionId: '1' }) => {
    const { submissionId, submission, userRole = 'pharma', isAuthenticated = true } = options;

    // Mock useAuth hook to return specified user role and authentication status
    (useAuth as vi.Mock).mockReturnValue({
      user: { id: '1', email: 'test@example.com', role: userRole },
      isAuthenticated: isAuthenticated,
    });

    // Mock useSubmissions hook to return specified submission data and mock functions
    const mockCancelSubmission = vi.fn();
    const mockRespondToQuote = vi.fn();
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: false,
      error: null,
      currentSubmission: submission,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: mockCancelSubmission,
      respondToQuote: mockRespondToQuote,
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component with customRender
    customRender(<SubmissionDetail submissionId={submissionId} />);

    // Return rendered component and mock functions for assertions
    return { mockCancelSubmission, mockRespondToQuote };
  };

  it('renders loading state when fetching submission', () => {
    // Mock useSubmissions to return loading: true
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: true,
      error: null,
      currentSubmission: null,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: vi.fn(),
      respondToQuote: vi.fn(),
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component
    customRender(<SubmissionDetail submissionId="1" />);

    // Assert that loading indicator is displayed
    expect(screen.getByText('Loading submission details...')).toBeInTheDocument();
  });

  it('renders error state when fetch fails', () => {
    // Mock useSubmissions to return error: 'Failed to fetch submission'
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: false,
      error: 'Failed to fetch submission',
      currentSubmission: null,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: vi.fn(),
      respondToQuote: vi.fn(),
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component
    customRender(<SubmissionDetail submissionId="1" />);

    // Assert that error message is displayed
    expect(screen.getByText('Failed to fetch submission')).toBeInTheDocument();
  });

  it('renders submission details correctly', () => {
    // Create mock submission data
    const mockSubmission = mockSubmissionData();

    // Mock useSubmissions to return the mock submission
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: false,
      error: null,
      currentSubmission: mockSubmission,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: vi.fn(),
      respondToQuote: vi.fn(),
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component
    customRender(<SubmissionDetail submissionId="1" />);

    // Assert that submission details are displayed correctly
    expect(screen.getByText('Test Experiment - Binding Assay')).toBeInTheDocument();
    expect(screen.getByText('Experiment Details')).toBeInTheDocument();
    expect(screen.getByText('CRO Details')).toBeInTheDocument();
    expect(screen.getByText('Submission Details')).toBeInTheDocument();
    expect(screen.getByText('Timeline')).toBeInTheDocument();

    // Assert that experiment details are displayed
    expect(screen.getByText('Experiment Type:')).toBeInTheDocument();
    expect(screen.getByText('Created By:')).toBeInTheDocument();

    // Assert that CRO details are displayed
    expect(screen.getByText('CRO:')).toBeInTheDocument();

    // Assert that status badge is displayed with correct status
    expect(screen.getByText('Submitted')).toBeInTheDocument();
  });

  it('shows quote approval UI for pharma user when status is QUOTE_PROVIDED', () => {
    // Create mock submission with status QUOTE_PROVIDED
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.QUOTE_PROVIDED });

    // Mock useAuth to return role: 'pharma'
    (useAuth as vi.Mock).mockReturnValue({
      user: { id: '1', email: 'test@example.com', role: 'pharma' },
      isAuthenticated: true,
    });

    // Render SubmissionDetail component
    setup({ submissionId: '1', submission: mockSubmission, userRole: 'pharma' });

    // Assert that quote approval UI is displayed
    expect(screen.getByText('Quote Details')).toBeInTheDocument();

    // Assert that approve and reject buttons are displayed
    expect(screen.getByText('Approve Quote')).toBeInTheDocument();
    expect(screen.getByText('Reject Quote')).toBeInTheDocument();
  });

  it('does not show quote approval UI for CRO user', () => {
    // Create mock submission with status QUOTE_PROVIDED
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.QUOTE_PROVIDED });

    // Mock useAuth to return role: 'cro'
    (useAuth as vi.Mock).mockReturnValue({
      user: { id: '2', email: 'cro@example.com', role: 'cro' },
      isAuthenticated: true,
    });

    // Render SubmissionDetail component
    setup({ submissionId: '1', submission: mockSubmission, userRole: 'cro' });

    // Assert that quote approval UI is not displayed
    expect(screen.queryByText('Quote Details')).toBeNull();
  });

  it('shows cancel button for pharma user when status is not COMPLETED or CANCELLED', () => {
    // Create mock submission with status IN_PROGRESS
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.IN_PROGRESS });

    // Mock useAuth to return role: 'pharma'
    (useAuth as vi.Mock).mockReturnValue({
      user: { id: '1', email: 'test@example.com', role: 'pharma' },
      isAuthenticated: true,
    });

    // Render SubmissionDetail component
    setup({ submissionId: '1', submission: mockSubmission, userRole: 'pharma' });

    // Assert that cancel button is displayed
    expect(screen.getByText('Cancel Submission')).toBeInTheDocument();
  });

  it('does not show cancel button for completed or cancelled submissions', () => {
    // Create mock submission with status COMPLETED
    const mockSubmissionCompleted = mockSubmissionData({ status: SubmissionStatus.COMPLETED });

    // Mock useAuth to return role: 'pharma'
    (useAuth as vi.Mock).mockReturnValue({
      user: { id: '1', email: 'test@example.com', role: 'pharma' },
      isAuthenticated: true,
    });

    // Render SubmissionDetail component
    setup({ submissionId: '1', submission: mockSubmissionCompleted, userRole: 'pharma' });

    // Assert that cancel button is not displayed
    expect(screen.queryByText('Cancel Submission')).toBeNull();

    // Update mock submission with status CANCELLED
    const mockSubmissionCancelled = mockSubmissionData({ status: SubmissionStatus.CANCELLED });

    // Re-render component
    setup({ submissionId: '1', submission: mockSubmissionCancelled, userRole: 'pharma' });

    // Assert that cancel button is not displayed
    expect(screen.queryByText('Cancel Submission')).toBeNull();
  });

  it('shows view results button when results are available', () => {
    // Create mock submission with results array
    const mockSubmission = mockSubmissionData({ results: [{ id: '1' }] });

    // Mock useSubmissions to return the mock submission
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: false,
      error: null,
      currentSubmission: mockSubmission,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: vi.fn(),
      respondToQuote: vi.fn(),
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component
    customRender(<SubmissionDetail submissionId="1" />);

    // Assert that view results button is displayed
    expect(screen.getByText('View Results')).toBeInTheDocument();
  });

  it('calls cancelSubmission when cancel button is clicked', async () => {
    // Create mock submission with status PENDING
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.PENDING });

    // Mock useSubmissions with mock cancelSubmission function
    const { mockCancelSubmission } = setup({ submissionId: '1', submission: mockSubmission });

    // Click cancel button
    fireEvent.click(screen.getByText('Cancel Submission'));

    // Enter cancellation reason in dialog
    fireEvent.change(screen.getByLabelText('Cancellation Reason'), { target: { value: 'Test reason' } });

    // Click confirm button
    fireEvent.click(within(screen.getByRole('dialog')).getByText('Cancel Submission'));

    // Assert that cancelSubmission was called with correct parameters
    await waitFor(() => {
      expect(mockCancelSubmission).toHaveBeenCalledWith('1', 'Test reason');
    });
  });

  it('calls respondToQuote when quote is approved', async () => {
    // Create mock submission with status QUOTE_PROVIDED
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.QUOTE_PROVIDED });

    // Mock useSubmissions with mock respondToQuote function
    const { mockRespondToQuote } = setup({ submissionId: '1', submission: mockSubmission, userRole: 'pharma' });

    // Click approve quote button
    fireEvent.click(screen.getByText('Approve Quote'));

    // Assert that respondToQuote was called with correct parameters (approved: true)
    await waitFor(() => {
      expect(mockRespondToQuote).toHaveBeenCalledWith('1', { approved: true });
    });
  });

  it('calls respondToQuote when quote is rejected', async () => {
    // Create mock submission with status QUOTE_PROVIDED
    const mockSubmission = mockSubmissionData({ status: SubmissionStatus.QUOTE_PROVIDED });

    // Mock useSubmissions with mock respondToQuote function
    const { mockRespondToQuote } = setup({ submissionId: '1', submission: mockSubmission, userRole: 'pharma' });

    // Click reject quote button
    fireEvent.click(screen.getByText('Reject Quote'));

    // Enter rejection reason
    fireEvent.change(screen.getByLabelText('Rejection Notes'), { target: { value: 'Test rejection reason' } });

    // Click confirm button
    fireEvent.click(within(screen.getByRole('dialog')).getByText('Reject'));

    // Assert that respondToQuote was called with correct parameters (approved: false)
    await waitFor(() => {
      expect(mockRespondToQuote).toHaveBeenCalledWith('1', { approved: false, notes: 'Test rejection reason' });
    });
  });

  it('calls onBack when back button is clicked', async () => {
    // Create mock submission
    const mockSubmission = mockSubmissionData();

    // Create mock onBack function
    const mockOnBack = vi.fn();

    // Mock useSubmissions to return the mock submission
    (useSubmissions as vi.Mock).mockReturnValue({
      loading: false,
      error: null,
      currentSubmission: mockSubmission,
      fetchSubmissionDetails: vi.fn(),
      cancelSubmission: vi.fn(),
      respondToQuote: vi.fn(),
      updateSubmissionStatus: vi.fn(),
    });

    // Render SubmissionDetail component with onBack prop
    customRender(<SubmissionDetail submissionId="1" onBack={mockOnBack} />);

    // Click back button
    fireEvent.click(screen.getByText('Back to Submissions'));

    // Assert that onBack was called
    await waitFor(() => {
      expect(mockOnBack).toHaveBeenCalled();
    });
  });
});