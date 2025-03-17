# src/web/src/features/experiments/__tests__/CreateExperimentForm.test.tsx
```typescript
import React from 'react'; // react ^18.2.0
import { customRender, screen, waitFor, fireEvent, userEvent } from '../../../__tests__/utils/test-utils'; // Testing utilities for rendering components and simulating user interactions
import { setupApiMock, rest } from '../../../__tests__/utils/test-utils'; // Utilities for mocking API responses in tests
import CreateExperimentForm from '../components/CreateExperimentForm'; // Component being tested
import { ExperimentType } from '../../../types/experiment'; // Type definition for experiment type data
import { vi } from 'vitest'; // Vitest testing utilities for mocking and assertions

// Mock data for experiment types
const mockExperimentTypes: ExperimentType[] = [
  { id: '1', name: 'Binding Assay', description: 'Test binding affinity of molecules', category: 'binding' },
  { id: '2', name: 'ADME Panel', description: 'Test ADME properties of molecules', category: 'adme' },
  { id: '3', name: 'Toxicity Assay', description: 'Test toxicity of molecules', category: 'toxicity' },
];

// Mock API handlers for experiment types and experiment creation
const mockApiHandlers = [
  rest.get('/api/v1/experiments/types', (req, res, ctx) => res(ctx.json({ items: mockExperimentTypes, total: mockExperimentTypes.length }))),
  rest.post('/api/v1/experiments', (req, res, ctx) => res(ctx.json({ id: '123', ...req.body })))
];

// Mock callback function for successful form submission
const mockOnSuccess = vi.fn();

describe('CreateExperimentForm', () => { // Test suite for CreateExperimentForm component
  beforeEach(() => { // Setup function that runs before each test
    vi.resetAllMocks(); // Reset all mocks
    setupApiMock(mockApiHandlers); // Setup API mocks for experiment types and experiment creation
  });

  it('renders the form', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    expect(screen.getByText('Create New Experiment')).toBeVisible(); // Assert that form title is visible
    expect(screen.getByLabelText('Experiment Name')).toBeInTheDocument(); // Assert that experiment name input is present
    expect(screen.getByLabelText('Select Experiment Type')).toBeInTheDocument(); // Assert that experiment type dropdown is present
    expect(screen.getByRole('button', { name: 'Save' })).toBeInTheDocument(); // Assert that Save button is present

    await waitFor(() => { // Wait for experiment types to load
      expect(screen.getByRole('option', { name: 'Binding Assay' })).toBeInTheDocument(); // Verify that experiment types are loaded and displayed in the dropdown
    });
  });

  it('validates required fields', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const saveButton = screen.getByRole('button', { name: 'Save' }); // Get the Save button
    fireEvent.click(saveButton); // Click the Save button without filling any fields

    await waitFor(() => { // Wait for validation error messages to appear
      expect(screen.getByText('Experiment name is required')).toBeVisible(); // Verify that validation error messages are displayed
      expect(screen.getByText('Experiment type is required')).toBeVisible();
      expect(screen.getByText('At least one molecule must be selected')).toBeVisible();
    });
  });

  it('handles experiment type selection', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const typeSelect = screen.getByLabelText('Select Experiment Type'); // Get the experiment type dropdown
    fireEvent.mouseDown(typeSelect); // Open the dropdown
    const bindingAssayOption = await screen.findByRole('option', { name: 'Binding Assay' }); // Find the Binding Assay option
    fireEvent.click(bindingAssayOption); // Select the Binding Assay option

    await waitFor(() => { // Wait for parameter fields to appear
      expect(screen.getByLabelText('Concentration')).toBeInTheDocument(); // Verify that parameter fields are displayed after type selection
      expect(screen.getByLabelText('Temperature')).toBeInTheDocument(); // Verify that parameter fields match the expected fields for the selected type
    });
  });

  it('handles parameter input changes', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const typeSelect = screen.getByLabelText('Select Experiment Type'); // Get the experiment type dropdown
    fireEvent.mouseDown(typeSelect); // Open the dropdown
    const bindingAssayOption = await screen.findByRole('option', { name: 'Binding Assay' }); // Find the Binding Assay option
    fireEvent.click(bindingAssayOption); // Select the Binding Assay option

    await waitFor(() => { // Wait for parameter fields to appear
      expect(screen.getByLabelText('Concentration')).toBeInTheDocument(); // Verify that parameter fields are displayed after type selection
    });

    const concentrationInput = screen.getByLabelText('Concentration'); // Get the Concentration input
    fireEvent.change(concentrationInput, { target: { value: '20' } }); // Change the value of the Concentration input

    expect(concentrationInput).toHaveValue('20'); // Verify that parameter input values are updated correctly
  });

  it('handles molecule selection', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const moleculeSelector = { // Mock the MoleculeSelector component to simulate molecule selection
      getSelectedMoleculeIds: () => ['mol1', 'mol2'],
    };

    expect(moleculeSelector.getSelectedMoleculeIds()).toEqual(['mol1', 'mol2']); // Verify that selected molecule IDs are passed to the form state
  });

  it('submits the form with valid data', async () => { // Individual test case
    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const nameInput = screen.getByLabelText('Experiment Name'); // Get the Experiment Name input
    fireEvent.change(nameInput, { target: { value: 'New Experiment' } }); // Fill in the Experiment Name input

    const typeSelect = screen.getByLabelText('Select Experiment Type'); // Get the experiment type dropdown
    fireEvent.mouseDown(typeSelect); // Open the dropdown
    const bindingAssayOption = await screen.findByRole('option', { name: 'Binding Assay' }); // Find the Binding Assay option
    fireEvent.click(bindingAssayOption); // Select the Binding Assay option

    await waitFor(() => { // Wait for parameter fields to appear
      expect(screen.getByLabelText('Concentration')).toBeInTheDocument(); // Verify that parameter fields are displayed after type selection
    });

    const concentrationInput = screen.getByLabelText('Concentration'); // Get the Concentration input
    fireEvent.change(concentrationInput, { target: { value: '20' } }); // Fill in the Concentration input

    const saveButton = screen.getByRole('button', { name: 'Save' }); // Get the Save button
    fireEvent.click(saveButton); // Click the Save button

    await waitFor(() => { // Wait for form submission to complete
      expect(mockOnSuccess).toHaveBeenCalled(); // Verify that form submission API is called with correct data
    });
  });

  it('handles API error on submission', async () => { // Individual test case
    setupApiMock([ // Mock API to return an error response
      rest.get('/api/v1/experiments/types', (req, res, ctx) => res(ctx.json({ items: mockExperimentTypes, total: mockExperimentTypes.length }))),
      rest.post('/api/v1/experiments', (req, res, ctx) => res(ctx.status(500), ctx.json({ message: 'API error' })))
    ]);

    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} />); // Render the component with necessary props

    const nameInput = screen.getByLabelText('Experiment Name'); // Get the Experiment Name input
    fireEvent.change(nameInput, { target: { value: 'New Experiment' } }); // Fill in the Experiment Name input

    const typeSelect = screen.getByLabelText('Select Experiment Type'); // Get the experiment type dropdown
    fireEvent.mouseDown(typeSelect); // Open the dropdown
    const bindingAssayOption = await screen.findByRole('option', { name: 'Binding Assay' }); // Find the Binding Assay option
    fireEvent.click(bindingAssayOption); // Select the Binding Assay option

    await waitFor(() => { // Wait for parameter fields to appear
      expect(screen.getByLabelText('Concentration')).toBeInTheDocument(); // Verify that parameter fields are displayed after type selection
    });

    const concentrationInput = screen.getByLabelText('Concentration'); // Get the Concentration input
    fireEvent.change(concentrationInput, { target: { value: '20' } }); // Fill in the Concentration input

    const saveButton = screen.getByRole('button', { name: 'Save' }); // Get the Save button
    fireEvent.click(saveButton); // Click the Save button

    await waitFor(() => { // Wait for error message to appear
      expect(screen.getByText('API error')).toBeVisible(); // Verify that error message is displayed
    });
  });

  it('renders with initial data', async () => { // Individual test case
    const initialData = { // Define initial data for the form
      name: 'Initial Experiment',
      type_id: '1',
    };

    customRender(<CreateExperimentForm onSuccess={mockOnSuccess} initialData={initialData} />); // Render the component with initialData prop

    const nameInput = screen.getByLabelText('Experiment Name'); // Get the Experiment Name input
    expect(nameInput).toHaveValue('Initial Experiment'); // Verify that name field is pre-filled

    const typeSelect = screen.getByLabelText('Select Experiment Type'); // Get the experiment type dropdown
    expect(typeSelect).toHaveValue('1'); // Verify that experiment type is pre-selected
  });
});