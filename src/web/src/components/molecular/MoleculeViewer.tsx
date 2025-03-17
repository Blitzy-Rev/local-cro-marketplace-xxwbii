import React, { useState, useEffect, useMemo } from 'react';
import { styled } from '@mui/material/styles'; // v5.13+
import Box from '@mui/material/Box'; // v5.13+
import Card from '@mui/material/Card'; // v5.13+
import CardContent from '@mui/material/CardContent'; // v5.13+
import CardActions from '@mui/material/CardActions'; // v5.13+
import Typography from '@mui/material/Typography'; // v5.13+
import Tooltip from '@mui/material/Tooltip'; // v5.13+
import IconButton from '@mui/material/IconButton'; // v5.13+
import Chip from '@mui/material/Chip'; // v5.13+
import ZoomIn from '@mui/icons-material/ZoomIn'; // v5.13+
import ZoomOut from '@mui/icons-material/ZoomOut'; // v5.13+
import Fullscreen from '@mui/icons-material/Fullscreen'; // v5.13+
import FullscreenExit from '@mui/icons-material/FullscreenExit'; // v5.13+
import Info from '@mui/icons-material/Info'; // v5.13+

import SmileRenderer from './SmileRenderer';
import { isValidSMILES, getMoleculePropertyByName, formatPropertyValue, getPropertyDisplayName } from '../../utils/molecularUtils';
import { Molecule } from '../../types/molecule';
import theme from '../../theme';

/**
 * Props interface for the MoleculeViewer component
 */
export interface MoleculeViewerProps {
  /** SMILES string representing the molecular structure */
  smiles?: string;
  /** Molecule object containing SMILES and property data */
  molecule?: Molecule;
  /** Width of the viewer in pixels */
  width?: number;
  /** Height of the viewer in pixels */
  height?: number;
  /** Whether the component should be interactive */
  interactive?: boolean;
  /** Whether to show control buttons */
  showControls?: boolean;
  /** Whether to show molecule properties */
  showProperties?: boolean;
  /** Array of property names to display */
  propertiesToShow?: string[];
  /** Additional CSS class */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
  /** Click handler function */
  onClick?: () => void;
  /** Alternative text for accessibility */
  alt?: string;
}

// Styled components
const ViewerContainer = styled(Card, {
  shouldForwardProp: (prop) => !['width', 'height', 'interactive'].includes(prop as string),
})<{ width?: number; height?: number; interactive?: boolean }>(({ theme, width, height, interactive }) => ({
  width: width ? `${width}px` : '100%',
  height: height ? `${height}px` : 'auto',
  display: 'flex',
  flexDirection: 'column',
  overflow: 'hidden',
  borderRadius: '8px',
  boxShadow: theme.shadows[2],
  transition: 'all 0.2s ease-in-out',
  '&:hover': {
    transform: interactive ? 'scale(1.02)' : 'none',
    boxShadow: interactive ? theme.shadows[4] : theme.shadows[2],
  },
}));

const MoleculeContainer = styled(Box)({
  flex: '1',
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'center',
  padding: '16px',
  backgroundColor: theme.palette.background.paper,
  position: 'relative',
});

const ControlsContainer = styled(CardActions)({
  display: 'flex',
  justifyContent: 'space-between',
  padding: '8px 16px',
  backgroundColor: theme.palette.background.default,
});

const PropertiesContainer = styled(Box)({
  display: 'flex',
  flexWrap: 'wrap',
  gap: '8px',
  padding: '8px 16px',
  backgroundColor: theme.palette.background.default,
  borderTop: `1px solid ${theme.palette.divider}`,
});

const PropertyChip = styled(Chip)({
  fontSize: '0.75rem',
  height: '24px',
});

const InfoPanel = styled(Box)(({ theme }) => ({
  position: 'absolute',
  top: 0,
  left: 0,
  right: 0,
  bottom: 0,
  backgroundColor: 'rgba(0,0,0,0.8)',
  color: 'white',
  padding: theme.spacing(2),
  overflow: 'auto',
  zIndex: 10,
  display: 'flex',
  flexDirection: 'column',
}));

/**
 * A component that provides an enhanced visualization of molecular structures with additional features
 * like interactive controls, property display, and customization options.
 */
const MoleculeViewer: React.FC<MoleculeViewerProps> = ({
  smiles,
  molecule,
  width = 300,
  height = 250,
  interactive = true,
  showControls = true,
  showProperties = true,
  propertiesToShow = ['molecular_weight', 'logp', 'h_bond_donors', 'h_bond_acceptors'],
  className,
  style,
  onClick,
  alt,
}) => {
  // State for zoom, fullscreen and info display
  const [zoom, setZoom] = useState(1);
  const [isFullscreen, setIsFullscreen] = useState(false);
  const [showInfo, setShowInfo] = useState(false);
  
  // Determine which SMILES string to use (from props.smiles or props.molecule.smiles)
  const effectiveSmiles = useMemo(() => {
    if (smiles) return smiles;
    if (molecule?.smiles) return molecule.smiles;
    return '';
  }, [smiles, molecule]);
  
  // Calculate effective dimensions based on zoom and fullscreen state
  const effectiveWidth = useMemo(() => {
    const baseWidth = isFullscreen ? window.innerWidth * 0.9 : width;
    return Math.floor(baseWidth * zoom);
  }, [width, zoom, isFullscreen]);
  
  const effectiveHeight = useMemo(() => {
    const baseHeight = isFullscreen ? window.innerHeight * 0.8 : height;
    return Math.floor(baseHeight * zoom);
  }, [height, zoom, isFullscreen]);
  
  // Extract and format molecule properties
  const moleculeProperties = useMemo(() => {
    if (!molecule || !molecule.properties) return [];
    
    return propertiesToShow
      .map(propName => {
        const value = getMoleculePropertyByName(molecule, propName);
        if (value === null) return null;
        
        return {
          name: propName,
          displayName: getPropertyDisplayName(propName),
          value: formatPropertyValue(propName, value)
        };
      })
      .filter(Boolean) as Array<{name: string, displayName: string, value: string}>;
  }, [molecule, propertiesToShow]);
  
  // Handle fullscreen mode and window resize
  useEffect(() => {
    if (isFullscreen) {
      document.body.style.overflow = 'hidden';
      
      // Add event listener for ESC key to exit fullscreen
      const handleEscKey = (e: KeyboardEvent) => {
        if (e.key === 'Escape') {
          setIsFullscreen(false);
        }
      };
      
      // Add event listener for window resize to adjust dimensions
      const handleResize = () => {
        // Force rerender by updating state
        setZoom(prevZoom => {
          // Ensure staying within min/max zoom bounds
          return Math.max(0.5, Math.min(2.0, prevZoom));
        });
      };
      
      window.addEventListener('keydown', handleEscKey);
      window.addEventListener('resize', handleResize);
      
      return () => {
        document.body.style.overflow = '';
        window.removeEventListener('keydown', handleEscKey);
        window.removeEventListener('resize', handleResize);
      };
    } else {
      document.body.style.overflow = '';
    }
  }, [isFullscreen]);
  
  // Event handlers
  const handleZoomIn = (e: React.MouseEvent) => {
    e.stopPropagation();
    setZoom(prev => Math.min(prev + 0.1, 2.0));
  };
  
  const handleZoomOut = (e: React.MouseEvent) => {
    e.stopPropagation();
    setZoom(prev => Math.max(prev - 0.1, 0.5));
  };
  
  const handleToggleFullscreen = (e: React.MouseEvent) => {
    e.stopPropagation();
    setIsFullscreen(prev => !prev);
  };
  
  const handleToggleInfo = (e: React.MouseEvent) => {
    e.stopPropagation();
    setShowInfo(prev => !prev);
  };
  
  const handleClick = () => {
    if (interactive && onClick && !showInfo) {
      onClick();
    }
  };
  
  // If SMILES is invalid, return error message
  if (!isValidSMILES(effectiveSmiles)) {
    return (
      <ViewerContainer 
        width={width} 
        height={height} 
        className={className} 
        style={style}
        interactive={false}
        aria-label="Invalid molecular structure"
      >
        <MoleculeContainer>
          <Typography variant="body2" color="error">
            Invalid SMILES structure
          </Typography>
        </MoleculeContainer>
      </ViewerContainer>
    );
  }
  
  // Set appropriate ARIA attributes for accessibility
  const ariaLabel = alt || `Molecular structure of ${
    molecule?.id || effectiveSmiles
  }${moleculeProperties.length > 0 ? ` with properties: ${
    moleculeProperties.map(p => `${p.displayName} ${p.value}`).join(', ')
  }` : ''}`;
  
  return (
    <ViewerContainer 
      width={isFullscreen ? undefined : width}
      height={isFullscreen ? undefined : height}
      className={className}
      style={{
        ...style,
        position: isFullscreen ? 'fixed' : 'relative',
        top: isFullscreen ? '50%' : undefined,
        left: isFullscreen ? '50%' : undefined,
        transform: isFullscreen ? 'translate(-50%, -50%)' : undefined,
        zIndex: isFullscreen ? 1300 : undefined,
        width: isFullscreen ? '90vw' : undefined,
        height: isFullscreen ? '90vh' : undefined,
      }}
      interactive={interactive && !isFullscreen}
      role="region"
      aria-label={ariaLabel}
    >
      <MoleculeContainer 
        onClick={handleClick}
        role="img"
        aria-label={ariaLabel}
        sx={{ cursor: interactive && onClick ? 'pointer' : 'default' }}
      >
        <SmileRenderer 
          smiles={effectiveSmiles}
          width={effectiveWidth} 
          height={effectiveHeight}
          interactive={false}
          alt={alt || `Molecular structure of ${molecule?.id || effectiveSmiles}`}
        />
        
        {showInfo && molecule && (
          <InfoPanel
            onClick={(e) => e.stopPropagation()}
            role="dialog"
            aria-modal="true"
            aria-labelledby="molecule-info-title"
          >
            <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 2 }}>
              <Typography id="molecule-info-title" variant="h6">Molecule Details</Typography>
              <IconButton 
                size="small"
                onClick={(e) => {
                  e.stopPropagation();
                  setShowInfo(false);
                }}
                sx={{ color: 'white' }}
                aria-label="Close details"
              >
                <FullscreenExit fontSize="small" />
              </IconButton>
            </Box>
            
            <Typography variant="body2" sx={{ mb: 1 }}>
              ID: {molecule.id}
            </Typography>
            <Typography variant="body2" sx={{ mb: 2, wordBreak: 'break-all' }}>
              SMILES: {molecule.smiles}
            </Typography>
            
            {molecule.properties && molecule.properties.length > 0 && (
              <>
                <Typography variant="subtitle2" sx={{ mt: 2, mb: 1, fontWeight: 'bold' }}>Properties:</Typography>
                <Box sx={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: 1 }}>
                  {molecule.properties.map((prop) => (
                    <Typography key={prop.property_name} variant="body2" sx={{ display: 'flex', justifyContent: 'space-between' }}>
                      <span>{getPropertyDisplayName(prop.property_name)}:</span>
                      <span>{formatPropertyValue(prop.property_name, prop.property_value)}</span>
                    </Typography>
                  ))}
                </Box>
              </>
            )}
            
            {molecule.created_at && (
              <Typography variant="body2" sx={{ mt: 2 }}>
                Created: {new Date(molecule.created_at).toLocaleString()}
              </Typography>
            )}
            
            {molecule.creator && (
              <Typography variant="body2">
                Created by: {molecule.creator.email}
              </Typography>
            )}
          </InfoPanel>
        )}
      </MoleculeContainer>
      
      {showControls && (
        <ControlsContainer>
          <Box>
            <Tooltip title="Zoom in">
              <IconButton 
                onClick={handleZoomIn} 
                size="small" 
                aria-label="Zoom in"
                disabled={zoom >= 2.0}
              >
                <ZoomIn fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title="Zoom out">
              <IconButton 
                onClick={handleZoomOut} 
                size="small" 
                aria-label="Zoom out"
                disabled={zoom <= 0.5}
              >
                <ZoomOut fontSize="small" />
              </IconButton>
            </Tooltip>
          </Box>
          
          <Box>
            <Tooltip title="View details">
              <IconButton 
                onClick={handleToggleInfo} 
                size="small" 
                aria-label="View molecule details"
                aria-pressed={showInfo}
              >
                <Info fontSize="small" />
              </IconButton>
            </Tooltip>
            <Tooltip title={isFullscreen ? "Exit fullscreen" : "Fullscreen"}>
              <IconButton 
                onClick={handleToggleFullscreen} 
                size="small" 
                aria-label={isFullscreen ? "Exit fullscreen" : "View in fullscreen"}
              >
                {isFullscreen ? (
                  <FullscreenExit fontSize="small" />
                ) : (
                  <Fullscreen fontSize="small" />
                )}
              </IconButton>
            </Tooltip>
          </Box>
        </ControlsContainer>
      )}
      
      {showProperties && molecule && moleculeProperties.length > 0 && (
        <PropertiesContainer>
          {moleculeProperties.map((prop) => (
            <PropertyChip
              key={prop.name}
              label={`${prop.displayName}: ${prop.value}`}
              variant="outlined"
              size="small"
              title={`${prop.displayName}: ${prop.value}`}
            />
          ))}
        </PropertiesContainer>
      )}
    </ViewerContainer>
  );
};

export default MoleculeViewer;