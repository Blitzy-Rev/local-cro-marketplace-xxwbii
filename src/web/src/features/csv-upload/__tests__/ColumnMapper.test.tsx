# src/web/src/features/csv-upload/__tests__/ColumnMapper.test.tsx
```tsx
import React from 'react'; // react ^18.2.0
import { render, screen, fireEvent, waitFor, within } from '../../../__tests__/utils/test-utils'; // Testing utilities for rendering components and interacting with them
import userEvent from '../../../__tests__/utils/test-utils'; // User event simulation for testing user interactions
import ColumnMapper from '../components/ColumnMapper'; // Component under test
import * as useCSVImportHook from '../hooks/useCSVImport'; // Hook for mocking in tests
import * as csvUtils from '../../../utils/csvUtils'; // Utilities for mocking in tests
import { vi } from 'vitest'; // Mocking and test utilities

// Mock data for useCSVImport hook
const mockUseCSVImportData = {
  headers: ['SMILES', 'MW', 'LogP', 'Activity'],
  preview: [
    { SMILES: 'CCO', MW: '46.07', LogP: '-0.14', Activity: '78.5' },
    { SMILES: 'CCCCO', MW: '74.12', LogP: '0.88', Activity: '65.2' },
  ],
  rowCount: 2,
  currentMapping: [],
  mappingOptions: [
    { csvColumn: 'SMILES', mappingOptions: [{ value: 'smiles', label: 'SMILES', isCustom: false }, { value: 'custom:', label: 'Custom Property', isCustom: true }] },
    { csvColumn: 'MW', mappingOptions: [{ value: 'molecular_weight', label: 'Molecular Weight', isCustom: false }, { value: 'custom:', label: 'Custom Property', isCustom: true }] },
    { csvColumn: 'LogP', mappingOptions: [{ value: 'logp', label: 'LogP', isCustom: false }, { value: 'custom:', label: 'Custom Property', isCustom: true }] },
    { csvColumn: 'Activity', mappingOptions: [{ value: 'activity', label: 'Activity', isCustom: false }, { value: 'custom:', label: 'Custom Property', isCustom: true }] },
  ],
  availableProperties: [
    { name: 'smiles', display_name: 'SMILES', data_type: 'string' },
    { name: 'molecular_weight', display_name: 'Molecular Weight', data_type: 'number', unit: 'g/mol' },
    { name: 'logp', display_name: 'LogP', data_type: 'number' },
    { name: 'activity', display_name: 'Activity', data_type: 'number', unit: '%' },
  ],
  handleUpdateMapping: vi.fn(),
  handleGenerateMapping: vi.fn(),
  handleSubmitMapping: vi.fn().mockResolvedValue({}),
  isLoading: false,
  errors: []
};

// Mock implementation of validateMappingConfiguration function
const mockValidateMappingConfiguration = { isValid: true, errors: [] };

// Mock callback function for when mapping is complete
const mockOnMappingComplete = vi.fn();

// Mock callback function for when mapping validity changes
const mockOnMappingChange = vi.fn();

/**
 * Setup function to create mocks and render the component
 * @param {object} props
 * @returns {object} Rendered component and mocks
 */
const setup = (props = {}) => {
  // Create mock functions for onMappingComplete and onMappingChange callbacks
  const onMappingComplete = props.onMappingComplete || mockOnMappingComplete;
  const onMappingChange = props.onMappingChange || mockOnMappingChange;

  // Create mock data for useCSVImport hook
  const useCSVImportMock = props.useCSVImport || mockUseCSVImportData;

  // Mock useCSVImport hook to return mock data
  vi.spyOn(useCSVImportHook, 'useCSVImport').mockReturnValue(useCSVImportMock);

  // Mock validateMappingConfiguration function
  vi.spyOn(csvUtils, 'validateMappingConfiguration').mockReturnValue(mockValidateMappingConfiguration);

  // Render ColumnMapper component with mock props
  render(<ColumnMapper onMappingComplete={onMappingComplete} onMappingChange={onMappingChange} />);

  // Return rendered component and mock functions
  return {
    onMappingComplete,
    onMappingChange,
  };
};

describe('ColumnMapper', () => {
  it('renders correctly with headers and preview data', () => {
    // Setup component with mock data
    setup();

    // Check that component title is rendered
    expect(screen.getByText('Map CSV Columns')).toBeInTheDocument();

    // Check that all CSV headers are displayed
    mockUseCSVImportData.headers.forEach(header => {
      expect(screen.getByText(header)).toBeInTheDocument();
    });

    // Check that dropdown selectors are rendered for each header
    mockUseCSVImportData.headers.forEach(header => {
      expect(screen.getByLabelText(`Select Property for ${header}`)).toBeInTheDocument();
    });
  });

  it('allows selecting properties for columns', async () => {
    // Setup component with mock data
    const { onMappingChange } = setup();

    // Find dropdown for a specific header
    const smilesDropdown = screen.getByLabelText('Select Property for SMILES');

    // Open dropdown and select a property
    fireEvent.mouseDown(smilesDropdown);
    const listbox = await waitFor(() => screen.getByRole('listbox'));
    fireEvent.click(within(listbox).getByText('SMILES'));

    // Verify that handleUpdateMapping was called with correct parameters
    expect(useCSVImportHook.useCSVImport().handleUpdateMapping).toHaveBeenCalledWith('SMILES', 'smiles', false);
    expect(onMappingChange).toHaveBeenCalledWith(true);
  });

  it('allows entering custom property names', async () => {
    // Setup component with mock data
    setup();

    // Find dropdown for a specific header
    const smilesDropdown = screen.getByLabelText('Select Property for SMILES');

    // Open dropdown and select custom property option
    fireEvent.mouseDown(smilesDropdown);
    const listbox = await waitFor(() => screen.getByRole('listbox'));
    fireEvent.click(within(listbox).getByText('Custom Property'));

    // Verify that custom property input field appears
    const customNameInput = screen.getByLabelText('Custom Property Name');
    expect(customNameInput).toBeInTheDocument();

    // Enter custom property name
    fireEvent.change(customNameInput, { target: { value: 'My Custom Property' } });

    // Verify that handleUpdateMapping was called with correct parameters
    expect(useCSVImportHook.useCSVImport().handleUpdateMapping).toHaveBeenCalledWith('SMILES', 'My Custom Property', true);
  });

  it('validates mapping configuration', () => {
    // Setup component with mock data
    setup();

    // Mock validateMappingConfiguration to return validation errors
    const mockErrors = [{ column: 'SMILES', message: 'SMILES column must be mapped' }];
    vi.spyOn(csvUtils, 'validateMappingConfiguration').mockReturnValue({ isValid: false, errors: mockErrors });

    // Verify that error messages are displayed
    mockErrors.forEach(error => {
      expect(screen.getByText(error.message)).toBeInTheDocument();
    });

    // Verify that onMappingChange was called with false
    expect(mockOnMappingChange).toHaveBeenCalledWith(false);
  });

  it('allows auto-mapping columns', () => {
    // Setup component with mock data
    setup();

    // Find auto-map button
    const autoMapButton = screen.getByText('Auto Map');

    // Click auto-map button
    fireEvent.click(autoMapButton);

    // Verify that handleGenerateMapping was called
    expect(useCSVImportHook.useCSVImport().handleGenerateMapping).toHaveBeenCalled();
  });

  it('submits mapping when valid', async () => {
    // Setup component with mock data
    const { onMappingComplete } = setup();

    // Mock validateMappingConfiguration to return no errors
    vi.spyOn(csvUtils, 'validateMappingConfiguration').mockReturnValue({ isValid: true, errors: [] });

    // Find submit button
    const importMoleculesButton = screen.getByText('Import Molecules');

    // Click submit button
    fireEvent.click(importMoleculesButton);

    // Verify that handleSubmitMapping was called
    expect(useCSVImportHook.useCSVImport().handleSubmitMapping).toHaveBeenCalled();

    // Verify that onMappingComplete was called
    await waitFor(() => expect(onMappingComplete).toHaveBeenCalled());
  });

  it('shows error when submitting invalid mapping', () => {
    // Setup component with mock data
    setup();

    // Mock validateMappingConfiguration to return validation errors
    const mockErrors = [{ column: 'SMILES', message: 'SMILES column must be mapped' }];
    vi.spyOn(csvUtils, 'validateMappingConfiguration').mockReturnValue({ isValid: false, errors: mockErrors });

    // Find submit button
    const importMoleculesButton = screen.getByText('Import Molecules');

    // Click submit button
    fireEvent.click(importMoleculesButton);

    // Verify that error toast is shown
    expect(useCSVImportHook.useCSVImport().handleSubmitMapping).not.toHaveBeenCalled();
  });

  it('displays preview table with current mapping', () => {
    // Setup component with mock data
    setup();

    // Check that preview table is rendered
    expect(screen.getByRole('table')).toBeInTheDocument();

    // Verify that preview table has correct headers based on mapping
    mockUseCSVImportData.headers.forEach(header => {
      expect(screen.getByText(header)).toBeInTheDocument();
    });
  });

  it('handles loading state', () => {
    // Setup component with isLoading set to true
    setup({ useCSVImport: { ...mockUseCSVImportData, isLoading: true } });

    // Verify that loading indicator is displayed
    expect(screen.getByRole('progressbar')).toBeInTheDocument();

    // Verify that submit button is disabled during loading
    const importMoleculesButton = screen.getByText('Import Molecules');
    expect(importMoleculesButton).toBeDisabled();
  });

  it('handles errors from useCSVImport', () => {
    // Setup component with errors in useCSVImport
    setup({ useCSVImport: { ...mockUseCSVImportData, errors: ['Test error'] } });

    // Verify that error alert is displayed with correct message
    expect(screen.getByText('Test error')).toBeInTheDocument();
  });
});