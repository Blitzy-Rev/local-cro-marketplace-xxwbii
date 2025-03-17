import React from 'react'; // react ^18.2.0
import { render, screen, waitFor, fireEvent, userEvent } from '../../__tests__/utils/test-utils'; // @testing-library/react ^14.0.0
import MoleculeViewer from '../MoleculeViewer';
import { mockMolecule, validSmiles, invalidSmiles } from './mockData';

/**
 * Helper function to render the MoleculeViewer component with default or custom props
 * @param {object} props = {}
 * @returns {object} Rendered component with testing utilities
 */
const renderMoleculeViewer = (props = {}) => {
  // LD1: Merge default props with provided props
  const defaultProps = {
    smiles: validSmiles,
    width: 300,
    height: 200,
    interactive: true,
    showControls: true,
    showProperties: true,
  };
  const mergedProps = { ...defaultProps, ...props };

  // LD1: Render the MoleculeViewer component with customRender
  const rendered = render(<MoleculeViewer {...mergedProps} />);

  // LD1: Return the rendered component
  return rendered;
};

describe('MoleculeViewer', () => {
  it('renders without crashing', () => {
    // LD1: Render MoleculeViewer with valid SMILES
    renderMoleculeViewer();

    // LD1: Verify that the component is in the document
    expect(screen.getByRole('region')).toBeInTheDocument();
  });

  it('renders molecule structure correctly', async () => {
    // LD1: Render MoleculeViewer with valid SMILES
    renderMoleculeViewer();

    // LD1: Wait for the molecule image to load
    const moleculeImage = await screen.findByAltText(`Molecular structure of ${validSmiles}`);

    // LD1: Verify that the molecule image is displayed
    expect(moleculeImage).toBeInTheDocument();

    // LD1: Image should have correct alt text
    expect(moleculeImage).toHaveAttribute('alt', `Molecular structure of ${validSmiles}`);
  });

  it('displays error message for invalid SMILES', () => {
    // LD1: Render MoleculeViewer with invalid SMILES
    renderMoleculeViewer({ smiles: invalidSmiles });

    // LD1: Verify that an error message is displayed
    expect(screen.getByText('Invalid SMILES structure')).toBeInTheDocument();
  });

  it('shows controls when showControls is true', () => {
    // LD1: Render MoleculeViewer with showControls set to true
    renderMoleculeViewer({ showControls: true });

    // LD1: Verify that control buttons are displayed
    expect(screen.getByLabelText('Zoom in')).toBeInTheDocument();
    expect(screen.getByLabelText('Zoom out')).toBeInTheDocument();
    expect(screen.getByLabelText('View molecule details')).toBeInTheDocument();
  });

  it('hides controls when showControls is false', () => {
    // LD1: Render MoleculeViewer with showControls set to false
    renderMoleculeViewer({ showControls: false });

    // LD1: Verify that control buttons are not in the document
    expect(screen.queryByLabelText('Zoom in')).not.toBeInTheDocument();
    expect(screen.queryByLabelText('Zoom out')).not.toBeInTheDocument();
    expect(screen.queryByLabelText('View molecule details')).not.toBeInTheDocument();
  });

  it('shows properties when showProperties is true', () => {
    // LD1: Create a mock molecule object with properties
    const molecule = mockMolecule;

    // LD1: Render MoleculeViewer with the molecule and showProperties set to true
    renderMoleculeViewer({ molecule, showProperties: true });

    // LD1: Verify that property chips are displayed
    expect(screen.getByText('Molecular Weight: 46.07 g/mol')).toBeInTheDocument();
    expect(screen.getByText('LogP: -0.14')).toBeInTheDocument();

    // LD1: Property values should be correctly formatted
    expect(screen.getByText('H-Bond Donors: 1')).toBeInTheDocument();
  });

  it('hides properties when showProperties is false', () => {
    // LD1: Create a mock molecule object with properties
    const molecule = mockMolecule;

    // LD1: Render MoleculeViewer with the molecule and showProperties set to false
    renderMoleculeViewer({ molecule, showProperties: false });

    // LD1: Verify that property chips are not displayed
    expect(screen.queryByText('Molecular Weight: 46.07 g/mol')).not.toBeInTheDocument();
    expect(screen.queryByText('LogP: -0.14')).not.toBeInTheDocument();
  });

  it('filters properties based on propertiesToShow', () => {
    // LD1: Create a mock molecule object with multiple properties
    const molecule = mockMolecule;

    // LD1: Render MoleculeViewer with the molecule and specific propertiesToShow
    renderMoleculeViewer({ molecule, showProperties: true, propertiesToShow: ['molecular_weight'] });

    // LD1: Verify that only specified properties are displayed
    expect(screen.getByText('Molecular Weight: 46.07 g/mol')).toBeInTheDocument();

    // LD1: Non-specified properties should not be in the document
    expect(screen.queryByText('LogP: -0.14')).not.toBeInTheDocument();
  });

  it('calls onClick when clicked and interactive is true', () => {
    // LD1: Create a mock onClick function
    const onClick = jest.fn();

    // LD1: Render MoleculeViewer with the mock function and interactive set to true
    renderMoleculeViewer({ interactive: true, onClick });

    // LD1: Simulate a click on the component
    fireEvent.click(screen.getByRole('region'));

    // LD1: Verify that the mock function was called
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it('does not call onClick when clicked and interactive is false', () => {
    // LD1: Create a mock onClick function
    const onClick = jest.fn();

    // LD1: Render MoleculeViewer with the mock function and interactive set to false
    renderMoleculeViewer({ interactive: false, onClick });

    // LD1: Simulate a click on the component
    fireEvent.click(screen.getByRole('region'));

    // LD1: Verify that the mock function was not called
    expect(onClick).not.toHaveBeenCalled();
  });

  it('zooms in when zoom in button is clicked', async () => {
    // LD1: Render MoleculeViewer with showControls set to true
    renderMoleculeViewer({ showControls: true });

    // LD1: Find and click the zoom in button
    const zoomInButton = screen.getByLabelText('Zoom in');
    fireEvent.click(zoomInButton);

    // LD1: Verify that the zoom level increases
    const moleculeContainer = screen.getByRole('region');
    await waitFor(() => {
      expect(moleculeContainer).toHaveStyle('width: 330px');
    });
  });

  it('zooms out when zoom out button is clicked', async () => {
    // LD1: Render MoleculeViewer with showControls set to true
    renderMoleculeViewer({ showControls: true });

    // LD1: Find and click the zoom out button
    const zoomOutButton = screen.getByLabelText('Zoom out');
    fireEvent.click(zoomOutButton);

    // LD1: Verify that the zoom level decreases
    const moleculeContainer = screen.getByRole('region');
    await waitFor(() => {
      expect(moleculeContainer).toHaveStyle('width: 270px');
    });
  });

  it('applies custom width and height', () => {
    // LD1: Render MoleculeViewer with custom width and height
    renderMoleculeViewer({ width: 400, height: 300 });

    // LD1: Verify that the container has the correct dimensions
    const moleculeContainer = screen.getByRole('region');
    expect(moleculeContainer).toHaveStyle('width: 400px');
    expect(moleculeContainer).toHaveStyle('height: 300px');
  });

  it('applies custom className and style', () => {
    // LD1: Render MoleculeViewer with custom className and style
    const customClassName = 'custom-class';
    const customStyle = { border: '1px solid red' };
    renderMoleculeViewer({ className: customClassName, style: customStyle });

    // LD1: Verify that the custom styling is applied
    const moleculeContainer = screen.getByRole('region');
    expect(moleculeContainer).toHaveClass(customClassName);
    expect(moleculeContainer).toHaveStyle(customStyle);
  });
});