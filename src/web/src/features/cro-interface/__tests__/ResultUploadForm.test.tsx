import React from 'react'; // React v18.2+
import { render, screen, waitFor, fireEvent } from '@testing-library/react'; // @testing-library/react ^14.0.0
import { vi, describe, it, expect, beforeEach, afterEach } from 'vitest'; // vitest ^0.30.1
import ResultUploadForm from '../components/ResultUploadForm';
import { customRender, setupApiMock, rest } from '../../../__tests__/utils/test-utils';
import { SubmissionStatus } from '../../../types/submission';

// Mock submission data for testing
const mockSubmission = {
  id: 1,
  experiment_id: 101,
  cro_id: 201,
  status: SubmissionStatus.COMPLETED,
  submitted_at: '2023-08-01T10:00:00Z',
  updated_at: '2023-08-05T14:00:00Z',
  experiment: {
    id: 101,
    name: 'Binding Assay',
    type: 'Binding',
    molecules: [
      { id: 'mol1', smiles: 'CCO' },
      { id: 'mol2', smiles: 'c1ccccc1' },
    ],
  },
};

// Mock callback function for successful form submission
const mockOnSuccess = vi.fn();

// Mock callback function for form cancellation
const mockOnCancel = vi.fn();

// Mock file object for testing file upload
const mockFile = new File(['test content'], 'test-result.xlsx', { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });

// Mock API handlers for testing API interactions
const apiMocks = [
  rest.post('/api/results', (req, res, ctx) => {
    return res(ctx.status(201), ctx.json({ id: 123, message: 'Result created successfully' }));
  }),
  rest.post('/api/results/:id/data', (req, res, ctx) => {
    return res(ctx.status(201), ctx.json({ message: 'Result data added successfully' }));
  }),
];

describe('ResultUploadForm', () => {
  beforeEach(() => {
    // Reset all mocks before each test
    vi.clearAllMocks();

    // Set up API mocks for result upload endpoints
    setupApiMock(apiMocks);
  });

  afterEach(() => {
    // Clean up any resources or mocks after each test
  });

  it('renders the form correctly', () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Check that form title is displayed
    expect(screen.getByText('Upload Results')).toBeInTheDocument();

    // Check that file upload section is displayed
    expect(screen.getByText('Result Files')).toBeInTheDocument();

    // Check that structured data entry section is displayed
    expect(screen.getByText('Structured Data Entry')).toBeInTheDocument();

    // Check that notes field is displayed
    expect(screen.getByLabelText('Notes')).toBeInTheDocument();

    // Check that submit and cancel buttons are displayed
    expect(screen.getByRole('button', { name: 'Submit Results' })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Cancel' })).toBeInTheDocument();
  });

  it('handles file upload correctly', async () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Upload a mock file using the FileUpload component
    const fileInput = screen.getByLabelText('Result Files');
    fireEvent.change(fileInput, { target: { files: [mockFile] } });

    // Check that the file is displayed in the file list
    await waitFor(() => {
      expect(screen.getByText('test-result.xlsx')).toBeInTheDocument();
    });
  });

  it('handles structured data entry correctly', async () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Enter data values for molecule properties
    const solubilityInput = screen.getByPlaceholderText('Enter Solubility');
    fireEvent.change(solubilityInput, { target: { value: '1.23' } });

    // Check that the entered values are reflected in the form
    await waitFor(() => {
      expect(solubilityInput).toHaveValue('1.23');
    });
  });

  it('validates required fields', async () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Submit the form without entering required data
    const submitButton = screen.getByRole('button', { name: 'Submit Results' });
    fireEvent.click(submitButton);

    // Check that validation error messages are displayed
    await waitFor(() => {
      expect(screen.getByText('At least one result file is required')).toBeInTheDocument();
      expect(screen.getByText('Notes are required')).toBeInTheDocument();
    });
  });

  it('submits the form with valid data', async () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Upload a mock file
    const fileInput = screen.getByLabelText('Result Files');
    fireEvent.change(fileInput, { target: { files: [mockFile] } });

    // Enter structured data for molecules
    const solubilityInput = screen.getByPlaceholderText('Enter Solubility');
    fireEvent.change(solubilityInput, { target: { value: '1.23' } });

    // Enter notes
    const notesInput = screen.getByLabelText('Notes');
    fireEvent.change(notesInput, { target: { value: 'Test notes' } });

    // Submit the form
    const submitButton = screen.getByRole('button', { name: 'Submit Results' });
    fireEvent.click(submitButton);

    // Check that API endpoints are called with correct data
    await waitFor(() => {
      expect(mockOnSuccess).toHaveBeenCalled();
    });
  });

  it('handles API errors correctly', async () => {
    // Set up API mocks to return error responses
    setupApiMock([
      rest.post('/api/results', (req, res, ctx) => {
        return res(ctx.status(500), ctx.json({ message: 'Internal Server Error' }));
      }),
    ]);

    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Fill in form data and submit
    const fileInput = screen.getByLabelText('Result Files');
    fireEvent.change(fileInput, { target: { files: [mockFile] } });
    const notesInput = screen.getByLabelText('Notes');
    fireEvent.change(notesInput, { target: { value: 'Test notes' } });
    const submitButton = screen.getByRole('button', { name: 'Submit Results' });
    fireEvent.click(submitButton);

    // Check that error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Internal Server Error')).toBeInTheDocument();
    });

    // Check that onSuccess callback is not called
    expect(mockOnSuccess).not.toHaveBeenCalled();
  });

  it('calls onCancel when cancel button is clicked', () => {
    // Render ResultUploadForm with mock submission
    customRender(<ResultUploadForm submission={mockSubmission} onSuccess={mockOnSuccess} onCancel={mockOnCancel} />);

    // Click the cancel button
    const cancelButton = screen.getByRole('button', { name: 'Cancel' });
    fireEvent.click(cancelButton);

    // Check that onCancel callback is called
    expect(mockOnCancel).toHaveBeenCalled();
  });
});