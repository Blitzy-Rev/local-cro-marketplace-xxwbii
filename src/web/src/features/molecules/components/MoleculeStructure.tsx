import React from 'react';
import { styled } from '@mui/material/styles'; // v5.13+
import Box from '@mui/material/Box'; // v5.13+

import SmileRenderer from '../../../components/molecular/SmileRenderer';
import { isValidSMILES, truncateSMILES } from '../../../utils/molecularUtils';
import { Molecule } from '../../../types/molecule';
import theme from '../../../theme';

/**
 * Props interface for the MoleculeStructure component
 */
export interface MoleculeStructureProps {
  /** SMILES string representing the molecular structure */
  smiles?: string;
  /** Molecule object containing structure and property data */
  molecule?: Molecule;
  /** Width of the rendering in pixels */
  width?: number;
  /** Height of the rendering in pixels */
  height?: number;
  /** Whether the component should have interactive behavior */
  interactive?: boolean;
  /** Whether to show a border around the structure */
  showBorder?: boolean;
  /** Additional CSS class */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
  /** Click handler function */
  onClick?: () => void;
  /** Alternative text for accessibility */
  alt?: string;
}

/**
 * Styled container for the molecular structure
 */
const StructureContainer = styled(Box, {
  shouldForwardProp: (prop) => 
    prop !== 'width' && 
    prop !== 'height' && 
    prop !== 'interactive' && 
    prop !== 'showBorder',
})<{ 
  width?: number; 
  height?: number; 
  interactive?: boolean; 
  showBorder?: boolean;
}>(({ width, height, interactive, showBorder }) => ({
  width: width ? `${width}px` : '100%',
  height: height ? `${height}px` : '200px',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  overflow: 'hidden',
  position: 'relative',
  backgroundColor: theme.palette.background.paper,
  borderRadius: '4px',
  border: showBorder ? `1px solid ${theme.palette.divider}` : 'none',
  cursor: interactive ? 'pointer' : 'default',
  transition: 'all 0.2s ease-in-out',
  '&:hover': {
    transform: interactive ? 'scale(1.05)' : 'none',
    boxShadow: interactive ? theme.shadows[2] : 'none',
  },
}));

/**
 * A component that renders a molecular structure visualization with customizable options.
 * This component serves as a specialized wrapper around the SmileRenderer component,
 * providing additional features specific to the molecules feature area of the application.
 */
const MoleculeStructure: React.FC<MoleculeStructureProps> = ({
  smiles,
  molecule,
  width = 200,
  height = 200,
  interactive = false,
  showBorder = true,
  className,
  style,
  onClick,
  alt,
}) => {
  // Determine effective SMILES string from props or molecule object
  const effectiveSmiles = smiles || (molecule?.smiles || '');
  
  // Generate appropriate alt text for accessibility
  const defaultAlt = `Molecular structure of ${truncateSMILES(effectiveSmiles, 30)}`;
  
  // Handle click events if interactive and onClick provided
  const handleClick = () => {
    if (interactive && onClick) {
      onClick();
    }
  };

  return (
    <StructureContainer
      width={width}
      height={height}
      interactive={interactive}
      showBorder={showBorder}
      className={className}
      style={style}
      onClick={handleClick}
      aria-label={alt || defaultAlt}
      role={interactive ? 'button' : 'img'}
      tabIndex={interactive ? 0 : undefined}
    >
      <SmileRenderer
        smiles={effectiveSmiles}
        width={width}
        height={height}
        interactive={interactive}
        alt={alt || defaultAlt}
      />
    </StructureContainer>
  );
};

export default MoleculeStructure;