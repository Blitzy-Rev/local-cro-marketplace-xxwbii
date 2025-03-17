import React, { useState, useMemo } from 'react';
import { styled } from '@mui/material/styles'; // v5.13+
import Box from '@mui/material/Box'; // v5.13+
import CircularProgress from '@mui/material/CircularProgress'; // v5.13+
import Typography from '@mui/material/Typography'; // v5.13+

import { isValidSMILES, generateMoleculeImageUrl, truncateSMILES } from '../../utils/molecularUtils';
import theme from '../../theme';

// Interface for the component props
export interface SmileRendererProps {
  /** SMILES string representing the molecular structure */
  smiles: string;
  /** Width of the rendering in pixels */
  width?: number;
  /** Height of the rendering in pixels */
  height?: number;
  /** Image format ('svg' or 'png') */
  format?: string;
  /** Whether the component should have interactive behavior */
  interactive?: boolean;
  /** Additional CSS class */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
  /** Click handler function */
  onClick?: () => void;
  /** Alternative text for accessibility */
  alt?: string;
  /** Whether to show a placeholder or just text for invalid SMILES */
  showPlaceholderOnError?: boolean;
}

// Styled components
const MoleculeContainer = styled(Box, {
  shouldForwardProp: (prop) => prop !== 'interactive',
})<{ interactive?: boolean }>(({ theme, interactive }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  overflow: 'hidden',
  position: 'relative',
  backgroundColor: theme.palette.background.paper,
  borderRadius: '4px',
  cursor: interactive ? 'pointer' : 'default',
  transition: 'all 0.2s ease-in-out',
  '&:hover': {
    transform: interactive ? 'scale(1.05)' : 'none',
    boxShadow: interactive ? theme.shadows[2] : 'none',
  },
}));

const MoleculeImage = styled('img')({
  maxWidth: '100%',
  maxHeight: '100%',
  objectFit: 'contain',
});

const ErrorContainer = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '8px',
  textAlign: 'center',
});

const ErrorText = styled(Typography)(({ theme }) => ({
  fontSize: '0.75rem',
  color: theme.palette.error.main,
}));

/**
 * Component for rendering molecular structures from SMILES strings.
 * 
 * This component serves as the foundation for molecular visualization throughout the
 * application, providing a visual representation of chemical structures.
 */
const SmileRenderer: React.FC<SmileRendererProps> = ({
  smiles,
  width = 200,
  height = 200,
  format = 'svg',
  interactive = false,
  className,
  style,
  onClick,
  alt,
  showPlaceholderOnError = true,
}) => {
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  
  const isSmileValid = isValidSMILES(smiles);
  
  const imageUrl = useMemo(() => {
    if (!isSmileValid) return '';
    return generateMoleculeImageUrl(smiles, width, height, format);
  }, [smiles, width, height, format, isSmileValid]);
  
  const handleImageLoad = () => {
    setLoading(false);
    setError(null);
  };
  
  const handleImageError = () => {
    setLoading(false);
    setError('Failed to load molecular structure image');
  };

  const handleContainerClick = () => {
    if (interactive && onClick) {
      onClick();
    }
  };
  
  // Generate accessible alt text
  const imageAlt = alt || `Molecular structure of ${truncateSMILES(smiles, 30)}`;
  
  if (!isSmileValid) {
    if (showPlaceholderOnError) {
      return (
        <MoleculeContainer
          width={`${width}px`}
          height={`${height}px`}
          className={className}
          style={style}
          interactive={false}
        >
          <ErrorContainer>
            <ErrorText variant="caption">Invalid SMILES: {truncateSMILES(smiles, 20)}</ErrorText>
          </ErrorContainer>
        </MoleculeContainer>
      );
    }
    return (
      <ErrorText variant="caption">Invalid SMILES: {truncateSMILES(smiles, 20)}</ErrorText>
    );
  }
  
  return (
    <MoleculeContainer
      width={`${width}px`}
      height={`${height}px`}
      className={className}
      style={style}
      interactive={interactive}
      onClick={handleContainerClick}
    >
      {loading && (
        <CircularProgress size={24} thickness={4} />
      )}
      
      {error && (
        <ErrorContainer>
          <ErrorText variant="caption">{error}</ErrorText>
        </ErrorContainer>
      )}
      
      <MoleculeImage
        src={imageUrl}
        alt={imageAlt}
        onLoad={handleImageLoad}
        onError={handleImageError}
        style={{ display: loading || error ? 'none' : 'block' }}
      />
    </MoleculeContainer>
  );
};

export default SmileRenderer;