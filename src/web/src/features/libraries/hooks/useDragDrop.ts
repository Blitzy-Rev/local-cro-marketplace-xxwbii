import { useState, useCallback, useRef, useEffect } from 'react'; // react 18.2+
import { useDrag, useDrop } from 'react-dnd'; // react-dnd 16.0.1
import { getEmptyImage } from 'react-dnd-html5-backend'; // react-dnd-html5-backend 16.0.1
import { Molecule } from '../../../types/molecule';
import { Library } from '../../../types/library';
import { useToast } from '../../../hooks/useToast';

// Constants for drag item types
const ItemTypes = {
  MOLECULE: 'MOLECULE'
};

// Interface for the data transferred during drag operations
interface DragItem {
  type: string;
  molecule: Molecule;
  sourceLibraryId: string | null;
}

// Interface for the options passed to the useDragDrop hook
interface UseDragDropOptions {
  // Callback function triggered when a molecule is dropped on a library
  onDrop: (molecule: Molecule, libraryId: string, sourceLibraryId: string | null) => void;
  // Optional callback function triggered when drag starts
  onDragStart?: (molecule: Molecule, sourceLibraryId: string | null) => void;
  // Optional callback function triggered when drag ends
  onDragEnd?: () => void;
}

// Interface for the return value of the useDragDrop hook
interface UseDragDropResult {
  // Ref to attach to draggable molecule elements
  moleculeDragRef: React.RefObject<HTMLElement>;
  // Ref to attach to library drop target elements
  libraryDropRef: React.RefObject<HTMLElement>;
  // Whether a drag operation is in progress
  isDragging: boolean;
  // Whether a draggable item is currently over the drop target
  isOver: boolean;
  // Function to handle drag start events
  handleDragStart: (molecule: Molecule, sourceLibraryId: string | null) => void;
  // Function to handle drag end events
  handleDragEnd: () => void;
  // Function to handle drop events
  handleDrop: (libraryId: string, libraryName: string) => void;
}

/**
 * A custom React hook that provides drag and drop functionality for molecules and libraries.
 * This hook enables intuitive organization of molecules into libraries through a drag-and-drop interface.
 *
 * @param options Configuration options including callbacks for drag and drop events
 * @returns Object containing refs, state, and handler functions for drag and drop functionality
 */
const useDragDrop = (options: UseDragDropOptions): UseDragDropResult => {
  // Destructure options
  const { onDrop, onDragStart, onDragEnd } = options;
  
  // State for the current drag operation
  const [currentMolecule, setCurrentMolecule] = useState<Molecule | null>(null);
  const [sourceLibraryId, setSourceLibraryId] = useState<string | null>(null);
  
  // State for the drop target
  const [targetLibraryId, setTargetLibraryId] = useState<string>('');
  const [targetLibraryName, setTargetLibraryName] = useState<string>('');
  
  // Access toast notification functionality
  const { showToast } = useToast();
  
  // Refs for drag source and drop target DOM elements
  const moleculeRef = useRef<HTMLElement>(null);
  const libraryRef = useRef<HTMLElement>(null);
  
  // Handler for starting a drag operation
  const handleDragStart = useCallback((molecule: Molecule, sourceLibId: string | null) => {
    setCurrentMolecule(molecule);
    setSourceLibraryId(sourceLibId);
    
    if (onDragStart) {
      onDragStart(molecule, sourceLibId);
    }
  }, [onDragStart]);
  
  // Handler for ending a drag operation
  const handleDragEnd = useCallback(() => {
    setCurrentMolecule(null);
    setSourceLibraryId(null);
    
    if (onDragEnd) {
      onDragEnd();
    }
  }, [onDragEnd]);
  
  // Handler for configuring a drop target
  const handleDrop = useCallback((libraryId: string, libraryName: string) => {
    setTargetLibraryId(libraryId);
    setTargetLibraryName(libraryName);
  }, []);
  
  // Configure molecule drag source
  const [{ isDragging }, drag, preview] = useDrag(() => ({
    type: ItemTypes.MOLECULE,
    item: (): DragItem | null => {
      if (!currentMolecule) return null;
      
      return {
        type: ItemTypes.MOLECULE,
        molecule: currentMolecule,
        sourceLibraryId
      };
    },
    canDrag: !!currentMolecule,
    collect: (monitor) => ({
      isDragging: monitor.isDragging()
    }),
    end: (item, monitor) => {
      const didDrop = monitor.didDrop();
      
      handleDragEnd();
      
      if (!didDrop && item) {
        showToast({
          type: 'info',
          message: 'Molecule not added to any library'
        });
      }
    }
  }), [currentMolecule, sourceLibraryId, handleDragEnd, showToast]);
  
  // Configure library drop target
  const [{ isOver }, drop] = useDrop(() => ({
    accept: ItemTypes.MOLECULE,
    drop: (item: DragItem) => {
      if (!targetLibraryId) {
        console.warn('No library ID set for drop target');
        return;
      }
      
      onDrop(item.molecule, targetLibraryId, item.sourceLibraryId);
      
      // Only show success toast if different from source library
      if (item.sourceLibraryId !== targetLibraryId) {
        showToast({
          type: 'success',
          message: `Added molecule to ${targetLibraryName || 'library'}`
        });
      }
      
      return { droppedInLibrary: true };
    },
    collect: (monitor) => ({
      isOver: monitor.isOver()
    })
  }), [targetLibraryId, targetLibraryName, onDrop, showToast]);
  
  // Connect drag ref
  useEffect(() => {
    drag(moleculeRef);
  }, [drag]);
  
  // Connect drop ref
  useEffect(() => {
    drop(libraryRef);
  }, [drop]);
  
  // Use empty image as drag preview
  useEffect(() => {
    preview(getEmptyImage(), { captureDraggingState: true });
  }, [preview]);
  
  return {
    moleculeDragRef: moleculeRef,
    libraryDropRef: libraryRef,
    isDragging,
    isOver,
    handleDragStart,
    handleDragEnd,
    handleDrop
  };
};

export default useDragDrop;