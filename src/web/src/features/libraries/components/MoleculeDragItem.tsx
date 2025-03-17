import React, { useRef, useEffect, useMemo } from 'react'; // react 18.2+
import { styled } from '@mui/material/styles'; // 5.13+
import Box from '@mui/material/Box'; // 5.13+
import Paper from '@mui/material/Paper'; // 5.13+
import Typography from '@mui/material/Typography'; // 5.13+
import Chip from '@mui/material/Chip'; // 5.13+

import { useDragDrop } from '../hooks/useDragDrop';
import { Molecule } from '../../../types/molecule';
import MoleculeViewer from '../../../components/molecular/MoleculeViewer';
import { getMoleculePropertyByName, formatPropertyValue } from '../../../utils/molecularUtils';

/**
 * Props interface for the MoleculeDragItem component
 */
interface MoleculeDragItemProps {
  /** The molecule object to be rendered as a draggable item */
  molecule: Molecule;
  /** ID of the library the molecule is currently in, or null if not in a library */
  sourceLibraryId?: string | null;
  /** Callback function triggered when drag starts */
  onDragStart?: (molecule: Molecule, sourceLibraryId: string | null) => void;
  /** Callback function triggered when drag ends */
  onDragEnd?: () => void;
  /** Optional callback function triggered when the molecule is clicked */
  onClick?: (molecule: Molecule) => void;
  /** Whether to display molecule properties */
  showProperties?: boolean;
  /** Array of property names to display */
  propertiesToShow?: string[];
  /** Optional CSS class name for styling */
  className?: string;
  /** Optional inline styles */
  style?: React.CSSProperties;
}

// Styled components
const DraggableItem = styled(Paper)(({ theme, isDragging }: { theme: any; isDragging: boolean }) => ({
  padding: '8px',
  margin: '8px',
  cursor: 'grab',
  transition: 'all 0.2s ease-in-out',
  opacity: isDragging ? 0.5 : 1,
  transform: isDragging ? 'scale(0.95)' : 'scale(1)',
  boxShadow: isDragging ? 'none' : theme.shadows[1],
  '&:hover': {
    boxShadow: theme.shadows[3],
    transform: isDragging ? 'scale(0.95)' : 'scale(1.02)',
  },
}));

const MoleculeContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
});

const MoleculeInfo = styled(Box)({
  width: '100%',
  marginTop: '8px',
});

const MoleculeName = styled(Typography)({
  fontSize: '0.875rem',
  fontWeight: '500',
  textAlign: 'center',
  overflow: 'hidden',
  textOverflow: 'ellipsis',
  whiteSpace: 'nowrap',
});

const PropertiesContainer = styled(Box)({
  display: 'flex',
  flexWrap: 'wrap',
  justifyContent: 'center',
  gap: '4px',
  marginTop: '4px',
});

const PropertyChip = styled(Chip)({
  height: '20px',
  fontSize: '0.7rem',
});

/**
 * A component that renders a molecule as a draggable item for library organization
 */
const MoleculeDragItem: React.FC<MoleculeDragItemProps> = ({
  molecule,
  sourceLibraryId = null,
  onDragStart,
  onDragEnd,
  onClick,
  showProperties = true,
  propertiesToShow = ['molecular_weight', 'logp'],
  className = '',
  style = {},
}) => {
  // Define a callback function to handle drag start events
  const handleDragStartCallback = (molecule: Molecule, sourceLibraryId: string | null) => {
    if (onDragStart) {
      onDragStart(molecule, sourceLibraryId);
    }
  };

  // Use the useDragDrop hook to enable drag and drop functionality
  const { moleculeDragRef, isDragging } = useDragDrop({
    onDrop: undefined, // This component is only a drag source, so no onDrop is needed
    onDragStart: handleDragStartCallback,
    onDragEnd: onDragEnd,
  });

  // Create a ref for the drag start callback to prevent unnecessary hook re-renders
  const dragStartRef = useRef(onDragStart);

  // Update the drag start callback ref when the prop changes
  useEffect(() => {
    dragStartRef.current = onDragStart;
  }, [onDragStart, molecule, sourceLibraryId]);

  // Handle click events
  const handleClick = () => {
    if (onClick) {
      onClick(molecule);
    }
  };

  // Memoize the molecule properties to display
  const moleculeProperties = useMemo(() => {
    if (!molecule || !molecule.properties) return [];

    return propertiesToShow.map(propName => {
      const value = getMoleculePropertyByName(molecule, propName);
      return value !== null ? {
        name: propName,
        value: formatPropertyValue(propName, value)
      } : null;
    }).filter(Boolean);
  }, [molecule, propertiesToShow]);

  return (
    <DraggableItem
      ref={moleculeDragRef}
      className={className}
      style={style}
      isDragging={isDragging}
      onClick={handleClick}
    >
      <MoleculeContainer>
        <MoleculeViewer
          smiles={molecule.smiles}
          width={150}
          height={100}
          interactive={false}
        />
        <MoleculeInfo>
          <MoleculeName variant="subtitle2">
            {molecule.smiles}
          </MoleculeName>
        </MoleculeInfo>
        {showProperties && moleculeProperties.length > 0 && (
          <PropertiesContainer>
            {moleculeProperties.map((prop) => (
              <PropertyChip
                key={prop.name}
                label={`${prop.name}: ${prop.value}`}
                size="small"
              />
            ))}
          </PropertiesContainer>
        )}
      </MoleculeContainer>
    </DraggableItem>
  );
};

export default MoleculeDragItem;