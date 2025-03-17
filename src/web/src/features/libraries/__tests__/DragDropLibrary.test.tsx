# src/web/src/features/libraries/__tests__/DragDropLibrary.test.tsx
```typescript
import React from 'react'; // react ^18.2.0
import { render, screen, fireEvent, waitFor } from '../../../__tests__/utils/test-utils'; // Testing utilities for rendering and interacting with components
import DragDropLibrary from '../components/DragDropLibrary'; // Component under test
import { useToast } from '../../../hooks/useToast'; // Hook for displaying toast notifications
import useDragDrop from '../hooks/useDragDrop'; // Hook for drag and drop functionality
import { DndProvider } from 'react-dnd'; // react-dnd ^16.0.1
import { HTML5Backend } from 'react-dnd-html5-backend'; // react-dnd-html5-backend ^16.0.1
import jest from 'jest'; // Testing framework

// Mock the useDragDrop hook
jest.mock('../hooks/useDragDrop', () => ({
  __esModule: true,
  default: jest.fn(() => ({
    libraryDropRef: jest.fn(),
    isOver: false,
    handleDrop: jest.fn()
  }))
}));

// Mock the useToast hook
jest.mock('../../../hooks/useToast', () => ({
  useToast: jest.fn(() => ({ showToast: jest.fn() }))
}));

describe('DragDropLibrary', () => {
  const mockLibrary = { id: 'lib-1', name: 'Test Library' };
  const onDrop = jest.fn();

  const setup = (props: any = {}) => {
    (useDragDrop as jest.Mock).mockImplementation(() => ({
      libraryDropRef: { current: {} },
      isOver: false,
      handleDrop: jest.fn()
    }));
    (useToast as jest.Mock).mockImplementation(() => ({ showToast: jest.fn() }));

    const utils = render(
      <DndProvider backend={HTML5Backend}>
        <DragDropLibrary library={mockLibrary} onDrop={onDrop} {...props} />
      </DndProvider>
    );
    return {
      ...utils,
      mockUseDragDrop: useDragDrop as jest.Mock,
      mockUseToast: useToast as jest.Mock,
    };
  };

  it('renders empty state message when no children are provided', () => {
    setup({ children: null });
    expect(screen.getByText('Drag molecules here to add to library')).toBeInTheDocument();
  });

  it('renders children when provided', () => {
    setup({ children: <div data-testid="test-child">Test Child</div> });
    expect(screen.getByTestId('test-child')).toBeInTheDocument();
  });

  it('applies drop indicator styling when drag is over', () => {
    (useDragDrop as jest.Mock).mockImplementation(() => ({
      libraryDropRef: { current: {} },
      isOver: true,
      handleDrop: jest.fn()
    }));
    setup();
    expect(screen.getByText('Drop here to add')).toBeInTheDocument();
  });

  it('calls onDrop when a molecule is dropped', async () => {
    const handleDropMock = jest.fn();
    (useDragDrop as jest.Mock).mockImplementation(() => ({
      libraryDropRef: { current: {} },
      isOver: true,
      handleDrop: handleDropMock
    }));
    setup({ children: <div data-testid="test-child">Test Child</div> });

    const dropEvent = new Event('drop', { bubbles: true, cancelable: true });
    Object.defineProperty(dropEvent, 'dataTransfer', {
      value: {
        getData: jest.fn().mockReturnValue(JSON.stringify({ moleculeId: 'mol-1' }))
      }
    });

    fireEvent(screen.getByTestId('test-child'), dropEvent);

    await waitFor(() => {
      expect(handleDropMock).toHaveBeenCalledWith(mockLibrary.id, mockLibrary.name);
    });
  });

  it('shows toast notification on successful drop', async () => {
    const showToastMock = jest.fn();
    (useToast as jest.Mock).mockImplementation(() => ({ showToast: showToastMock }));
    (useDragDrop as jest.Mock).mockImplementation(() => ({
      libraryDropRef: { current: {} },
      isOver: true,
      handleDrop: jest.fn()
    }));
    setup({ children: <div data-testid="test-child">Test Child</div> });

    const dropEvent = new Event('drop', { bubbles: true, cancelable: true });
    Object.defineProperty(dropEvent, 'dataTransfer', {
      value: {
        getData: jest.fn().mockReturnValue(JSON.stringify({ moleculeId: 'mol-1' }))
      }
    });

    fireEvent(screen.getByTestId('test-child'), dropEvent);

    await waitFor(() => {
      expect(showToastMock).toHaveBeenCalledWith({
        type: 'success',
        message: `Added molecule to ${mockLibrary.name}`
      });
    });
  });

  it('applies custom className and style when provided', () => {
    const { container } = setup({ className: 'custom-class', style: { backgroundColor: 'red' } });
    expect(container.firstChild).toHaveClass('custom-class');
    expect(container.firstChild).toHaveStyle('backgroundColor: red');
  });

  it('does not show drop indicator when showDropIndicator is false', () => {
    (useDragDrop as jest.Mock).mockImplementation(() => ({
      libraryDropRef: { current: {} },
      isOver: true,
      handleDrop: jest.fn()
    }));
    setup({ showDropIndicator: false });
    expect(screen.queryByText('Drop here to add')).not.toBeInTheDocument();
  });
});