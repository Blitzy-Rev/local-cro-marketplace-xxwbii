# src/web/src/features/molecules/__tests__/MoleculeTable.test.tsx
```typescript
import React from 'react'; // React v18.2+
import { render, screen, waitFor, fireEvent, userEvent, within } from '../../../__tests__/utils/test-utils'; // Testing utilities
import MoleculeTable from '../components/MoleculeTable'; // Component under test
import { Molecule, FlagStatus } from '../../../types/molecule'; // Type definitions
import { getMoleculePropertyByName } from '../../../utils/molecularUtils'; // Utility function

/**
 * @function generateMockMolecules
 * @description Generates an array of mock molecule objects for testing
 * @param {number} count - The number of mock molecules to generate
 * @returns {Molecule[]} - An array of mock molecule objects
 */
const generateMockMolecules = (count: number): Molecule[] => {
  const molecules: Molecule[] = []; // Create an empty array to hold mock molecules
  for (let i = 0; i < count; i++) { // Loop count times to create the specified number of molecules
    const id = `mol-${i + 1}`; // For each molecule, generate a unique ID
    const smiles = `CCO${i}`; // Assign a valid SMILES string
    const properties = [ // Create properties array with common molecular properties
      { id: `prop-${i + 1}-1`, molecule_id: id, property_name: 'molecular_weight', property_value: 46.07 + i, is_calculated: false },
      { id: `prop-${i + 1}-2`, molecule_id: id, property_name: 'logp', property_value: -0.14 + i, is_calculated: true },
    ];
    const flag_status = i % 3 === 0 ? FlagStatus.IMPORTANT : i % 3 === 1 ? FlagStatus.FAVORITE : null; // Assign random flag_status or null
    molecules.push({ id, smiles, properties, flag_status, created_by: 'user-1', created_at: new Date().toISOString() }); // Add the molecule to the array
  }
  return molecules; // Return the array of mock molecules
};

/**
 * @testsuite MoleculeTable
 * @description Test suite for the MoleculeTable component
 */
describe('MoleculeTable', () => {
  /**
   * @test renders without crashing
   * @description Verifies that the component renders without errors
   */
  it('renders without crashing', () => {
    // Render MoleculeTable with minimal props
    render(<MoleculeTable molecules={[]} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={false} />);
    // Verify that the component is in the document
    expect(screen.getByText('No data available')).toBeInTheDocument();
  });

  /**
   * @test displays loading state correctly
   * @description Verifies that loading indicator is shown when loading prop is true
   */
  it('displays loading state correctly', () => {
    // Render MoleculeTable with loading=true
    render(<MoleculeTable molecules={[]} loading={true} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={false} />);
    // Verify that loading indicator is displayed
    expect(screen.getByRole('progressbar')).toBeInTheDocument();

    // Re-render with loading=false
    render(<MoleculeTable molecules={[]} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={false} />);
    // Verify that loading indicator is no longer displayed
    expect(screen.queryByRole('progressbar')).not.toBeInTheDocument();
  });

  /**
   * @test renders molecule data correctly
   * @description Verifies that molecule data is displayed correctly in the table
   */
  it('renders molecule data correctly', () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(2);

    // Render MoleculeTable with mock data
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={['molecular_weight', 'logp']} selectable={false} showActions={false} />);

    // Verify that molecule SMILES are displayed
    expect(screen.getByText('CCO0')).toBeInTheDocument();
    expect(screen.getByText('CCO1')).toBeInTheDocument();

    // Verify that molecule structures are rendered
    expect(screen.getAllByRole('img', { name: /molecular structure/i }).length).toBe(2);

    // Verify that molecule properties are displayed correctly
    expect(screen.getByText('46.07')).toBeInTheDocument();
    expect(screen.getByText('-0.14')).toBeInTheDocument();
    expect(screen.getByText('47.07')).toBeInTheDocument();
    expect(screen.getByText('-0.13')).toBeInTheDocument();
  });

  /**
   * @test handles selection correctly
   * @description Verifies that molecule selection works correctly
   */
  it('handles selection correctly', async () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(2);

    // Create mock selection change handler
    const onSelectionChange = jest.fn();

    // Render MoleculeTable with selectable=true
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={onSelectionChange} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={true} showActions={false} />);

    // Click on selection checkbox for a molecule
    const checkbox = screen.getAllByRole('checkbox')[1];
    await userEvent.click(checkbox);

    // Verify that selection change handler was called with correct molecule ID
    expect(onSelectionChange).toHaveBeenCalledWith([mockMolecules[0].id]);

    // Verify that selected row is visually indicated
    expect(checkbox).toBeChecked();
  });

  /**
   * @test handles view details action correctly
   * @description Verifies that view details button works correctly
   */
  it('handles view details action correctly', async () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(1);

    // Create mock view details handler
    const onViewDetails = jest.fn();

    // Render MoleculeTable with showActions=true
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={onViewDetails} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={true} />);

    // Click on view details button for a molecule
    const viewDetailsButton = screen.getByRole('button', { name: /view details/i });
    await userEvent.click(viewDetailsButton);

    // Verify that view details handler was called with correct molecule
    expect(onViewDetails).toHaveBeenCalledWith(mockMolecules[0]);
  });

  /**
   * @test handles add to queue action correctly
   * @description Verifies that add to queue button works correctly
   */
  it('handles add to queue action correctly', async () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(1);

    // Create mock add to queue handler
    const onAddToQueue = jest.fn();

    // Render MoleculeTable with showActions=true
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={onAddToQueue} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={true} />);

    // Click on add to queue button for a molecule
    const addToQueueButton = screen.getByRole('button', { name: /add to experiment queue/i });
    await userEvent.click(addToQueueButton);

    // Verify that add to queue handler was called with correct molecule
    expect(onAddToQueue).toHaveBeenCalledWith(mockMolecules[0]);
  });

  /**
   * @test handles flag change action correctly
   * @description Verifies that flag status change works correctly
   */
  it('handles flag change action correctly', async () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(1);

    // Create mock flag change handler
    const onFlagChange = jest.fn();

    // Render MoleculeTable with showActions=true
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={onFlagChange} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={true} />);

    // Click on flag button for a molecule
    const flagButton = screen.getByRole('button', { name: /flag as important/i });
    await userEvent.click(flagButton);

    // Verify that flag change handler was called with correct molecule and flag status
    expect(onFlagChange).toHaveBeenCalledWith(mockMolecules[0], FlagStatus.IMPORTANT);
  });

  /**
   * @test displays pagination correctly
   * @description Verifies that pagination controls work correctly
   */
  it('displays pagination correctly', async () => {
    // Create mock molecule data
    const mockMolecules = generateMockMolecules(25);

    // Create mock pagination handlers
    const onPageChange = jest.fn();
    const onPageSizeChange = jest.fn();

    // Render MoleculeTable with pagination props
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 25, onPageChange: onPageChange, onPageSizeChange: onPageSizeChange }} displayProperties={[]} selectable={false} showActions={false} />);

    // Verify that pagination controls are displayed
    expect(screen.getByRole('button', { name: /next page/i })).toBeInTheDocument();

    // Click on next page button
    const nextPageButton = screen.getByRole('button', { name: /next page/i });
    await userEvent.click(nextPageButton);

    // Verify that page change handler was called with correct page number
    expect(onPageChange).toHaveBeenCalledWith(2);
  });

  /**
   * @test displays custom properties correctly
   * @description Verifies that custom property columns are displayed correctly
   */
  it('displays custom properties correctly', () => {
    // Create mock molecule data with custom properties
    const mockMolecules = generateMockMolecules(1);

    // Render MoleculeTable with displayProperties prop
    render(<MoleculeTable molecules={mockMolecules} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={['molecular_weight']} selectable={false} showActions={false} />);

    // Verify that specified property columns are displayed
    expect(screen.getByText(/molecular weight/i)).toBeInTheDocument();

    // Verify that non-specified properties are not displayed
    expect(screen.queryByText(/logp/i)).not.toBeInTheDocument();
  });

  /**
   * @test handles empty data correctly
   * @description Verifies that empty state is handled correctly
   */
  it('handles empty data correctly', () => {
    // Render MoleculeTable with empty molecules array
    render(<MoleculeTable molecules={[]} loading={false} error={null} selectedMolecules={[]} onSelectionChange={() => {}} onViewDetails={() => {}} onAddToQueue={() => {}} onFlagChange={() => {}} pagination={{ page: 1, pageSize: 10, totalItems: 0, onPageChange: () => {}, onPageSizeChange: () => {} }} displayProperties={[]} selectable={false} showActions={false} />);

    // Verify that empty state message is displayed
    expect(screen.getByText(/no data available/i)).toBeInTheDocument();

    // Verify that table headers are still displayed
    expect(screen.getByText(/structure/i)).toBeInTheDocument();
    expect(screen.getByText(/smiles/i)).toBeInTheDocument();
  });
});