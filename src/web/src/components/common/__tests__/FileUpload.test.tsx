# src/web/src/components/common/__tests__/FileUpload.test.tsx
```typescript
import React from 'react'; // React v18.2+
import { FileUpload } from '../FileUpload';
import { customRender, screen, waitFor, fireEvent, userEvent } from '../../../__tests__/utils/test-utils';
import '@testing-library/jest-dom/extend-expect'; // v5.16.5

// Mock functions for testing callbacks
const onFilesSelected = jest.fn();
const onFileValidationFail = jest.fn();

// Helper function to create a mock File object
const createMockFile = (name: string, type: string, size: number): File => {
  // LD1: Create a new File object with empty array buffer
  const file = new File([''], name, { type });

  // LD1: Set the name, type, and size properties
  Object.defineProperty(file, 'size', { value: size });
  Object.defineProperty(file, 'type', { value: type });

  // LD1: Return the mock File object
  return file;
};

// Helper function to create a mock DataTransfer object for drag and drop testing
const createDataTransfer = (files: File[]): DataTransfer => {
  // LD1: Create a new object with files property
  const dataTransfer: any = {};

  // LD1: Add methods like getData, setData, and getFiles
  dataTransfer.files = files;
  dataTransfer.items = files.map(file => ({
    kind: 'file',
    type: file.type,
    getAsFile: () => file
  }));
  dataTransfer.getData = jest.fn();
  dataTransfer.setData = jest.fn();
  dataTransfer.clearData = jest.fn();
  dataTransfer.setDragImage = jest.fn();
  dataTransfer.addElement = jest.fn();
  dataTransfer.removeElement = jest.fn();

  // LD1: Return the mock DataTransfer object
  return dataTransfer as DataTransfer;
};

describe('FileUpload Component', () => {
  // Reset mocks before each test
  beforeEach(() => {
    onFilesSelected.mockClear();
    onFileValidationFail.mockClear();
  });

  // Define test file data
  const validFile = createMockFile('test.csv', 'text/csv', 1024);
  const invalidFile = createMockFile('test.txt', 'text/plain', 1024);
  const largeFile = createMockFile('test.csv', 'text/csv', 20485760); // 20MB

  it('renders with default props', () => {
    // LD1: Render the FileUpload component with default props
    customRender(<FileUpload onFilesSelected={onFilesSelected} />);

    // LD1: Assert that the component renders without errors
    expect(screen.getByText('Drag and drop files here, or click to select files')).toBeInTheDocument();

    // LD1: Assert that the upload area is visible
    expect(screen.getByRole('button', { name: 'Select Files' })).toBeVisible();
  });

  it('allows file selection via button click', async () => {
    // LD1: Render the FileUpload component with a mock onFilesSelected callback
    customRender(<FileUpload onFilesSelected={onFilesSelected} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the selected files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile]);
    });
  });

  it('handles drag and drop operations', async () => {
    // LD1: Render the FileUpload component with a mock onFilesSelected callback
    customRender(<FileUpload onFilesSelected={onFilesSelected} />);

    // LD1: Get the upload area element
    const uploadArea = screen.getByText('Drag and drop files here, or click to select files');

    // LD1: Create a mock DataTransfer object with the test file
    const dataTransfer = createDataTransfer([validFile]);

    // LD1: Simulate a drag enter event
    fireEvent.dragEnter(uploadArea, { dataTransfer });

    // LD1: Assert that the component shows active state during drag over
    expect(uploadArea).toHaveClass('MuiPaper-root');

    // LD1: Simulate a drop event
    fireEvent.drop(uploadArea, { dataTransfer });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the dropped files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile]);
    });
  });

  it('validates file types', async () => {
    // LD1: Render the FileUpload component with a mock onFilesSelected callback and accepted file types
    customRender(<FileUpload accept={['text/csv']} onFilesSelected={onFilesSelected} onFileValidationFail={onFileValidationFail} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event with an invalid file type
    fireEvent.change(fileInput, { target: { files: [invalidFile] } });

    // LD1: Wait for the onFileValidationFail callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFileValidationFail callback is called with the correct error message
      expect(onFileValidationFail).toHaveBeenCalledWith('File type not accepted: test.txt');
    });

    // LD1: Simulate a file selection event with a valid file type
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the selected files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile]);
    });
  });

  it('validates file size', async () => {
    // LD1: Render the FileUpload component with a mock onFilesSelected callback and a maximum file size
    customRender(<FileUpload maxSize={10485760} onFilesSelected={onFilesSelected} onFileValidationFail={onFileValidationFail} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event with a file that exceeds the maximum size
    fireEvent.change(fileInput, { target: { files: [largeFile] } });

    // LD1: Wait for the onFileValidationFail callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFileValidationFail callback is called with the correct error message
      expect(onFileValidationFail).toHaveBeenCalledWith('File is too large: test.csv');
    });

    // LD1: Simulate a file selection event with a file that is within the maximum size
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the selected files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile]);
    });
  });

  it('supports custom validation', async () => {
    // LD1: Define a mock custom validation function
    const customValidation = jest.fn().mockResolvedValue({ isValid: false, error: 'Custom validation failed' });

    // LD1: Render the FileUpload component with a mock onFilesSelected callback and a custom validation function
    customRender(<FileUpload customValidation={customValidation} onFilesSelected={onFilesSelected} onFileValidationFail={onFileValidationFail} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the custom validation function to be called
    await waitFor(() => {
      // LD1: Assert that the custom validation function is called with the correct file
      expect(customValidation).toHaveBeenCalledWith(validFile);
    });

    // LD1: Wait for the onFileValidationFail callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFileValidationFail callback is called with the correct error message
      expect(onFileValidationFail).toHaveBeenCalledWith('Custom validation failed');
    });

    // LD1: Update the mock custom validation function to return a valid result
    customValidation.mockResolvedValue({ isValid: true });

    // LD1: Simulate a file selection event again
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the selected files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile]);
    });
  });

  it('displays file previews when enabled', async () => {
    // LD1: Render the FileUpload component with showPreview enabled
    customRender(<FileUpload showPreview onFilesSelected={onFilesSelected} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event
    fireEvent.change(fileInput, { target: { files: [validFile] } });

    // LD1: Wait for the file preview to be displayed
    await waitFor(() => {
      // LD1: Assert that the file preview is displayed with the correct file name
      expect(screen.getByText('test.csv')).toBeInTheDocument();
    });
  });

  it('handles disabled state correctly', () => {
    // LD1: Render the FileUpload component with the disabled prop set to true
    customRender(<FileUpload disabled onFilesSelected={onFilesSelected} />);

    // LD1: Get the upload area element
    const uploadArea = screen.getByText('Drag and drop files here, or click to select files');

    // LD1: Assert that the upload area has the disabled styling
    expect(uploadArea).toHaveStyle('opacity: 0.6');

    // LD1: Assert that the upload area has the disabled pointer events
    expect(uploadArea).toHaveStyle('pointer-events: none');
  });

  it('displays loading state', () => {
    // LD1: Render the FileUpload component with the isLoading prop set to true
    customRender(<FileUpload isLoading onFilesSelected={onFilesSelected} />);

    // LD1: Assert that the loading spinner is displayed
    expect(screen.getByRole('progressbar')).toBeInTheDocument();
  });

  it('displays error messages', async () => {
    // LD1: Render the FileUpload component with a mock onFileValidationFail callback
    customRender(<FileUpload onFilesSelected={onFilesSelected} onFileValidationFail={onFileValidationFail} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event with an invalid file type
    fireEvent.change(fileInput, { target: { files: [invalidFile] } });

    // LD1: Wait for the error message to be displayed
    await waitFor(() => {
      // LD1: Assert that the error message is displayed with the correct text
      expect(screen.getByText('File type not accepted: test.txt')).toBeInTheDocument();
    });
  });

  it('supports multiple file selection', async () => {
    // LD1: Render the FileUpload component with the multiple prop set to true
    customRender(<FileUpload multiple onFilesSelected={onFilesSelected} />);

    // LD1: Get the file input element
    const fileInput = screen.getByLabelText('Select CSV File') as HTMLInputElement;

    // LD1: Simulate a file selection event with multiple files
    fireEvent.change(fileInput, { target: { files: [validFile, invalidFile] } });

    // LD1: Wait for the onFilesSelected callback to be called
    await waitFor(() => {
      // LD1: Assert that the onFilesSelected callback is called with the selected files
      expect(onFilesSelected).toHaveBeenCalledWith([validFile, invalidFile]);
    });
  });
});