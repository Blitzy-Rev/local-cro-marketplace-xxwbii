import React from 'react'; // react v18.2+
import { styled } from '@mui/material/styles'; // v5.13+
import { Box, Typography, Tooltip, IconButton } from '@mui/material'; // v5.13+
import { Visibility, Add, Flag, Star, StarBorder } from '@mui/icons-material'; // v5.13+

import Card from '../../../components/common/Card';
import MoleculeStructure from './MoleculeStructure';
import PropertyBadge from '../../../components/molecular/PropertyBadge';
import Button from '../../../components/common/Button';
import { Molecule, FlagStatus } from '../../../types/molecule';
import { getMoleculePropertyByName, getMoleculeCommonName, DEFAULT_MOLECULE_PROPERTIES } from '../../../utils/molecularUtils';
import theme from '../../../theme';

/**
 * Props interface for the MoleculeCard component
 */
interface MoleculeCardProps {
  /** Molecule data object */
  molecule: Molecule;
  /** Whether the card is in a selected state */
  selected?: boolean;
  /** Whether the card has interactive behavior */
  interactive?: boolean;
  /** Array of property names to display */
  displayProperties?: string[];
  /** Whether to show action buttons */
  showActions?: boolean;
  /** Click handler for the card */
  onClick?: () => void;
  /** Click handler specifically for molecule interactions */
  onMoleculeClick?: (molecule: Molecule) => void;
  /** Handler for viewing molecule details */
  onViewDetails?: (molecule: Molecule) => void;
  /** Handler for adding molecule to queue */
  onAddToQueue?: (molecule: Molecule) => void;
  /** Handler for changing flag status */
  onFlagChange?: (molecule: Molecule, flagStatus: FlagStatus | null) => void;
  /** Additional CSS class */
  className?: string;
  /** Additional inline styles */
  style?: React.CSSProperties;
  /** Width for structure visualization */
  structureWidth?: number;
  /** Height for structure visualization */
  structureHeight?: number;
}

// Styled components for the card layout
const CardContent = styled(Box)({
  display: 'flex',
  flexDirection: 'column',
  height: '100%',
});

const StructureContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'center',
  alignItems: 'center',
  marginBottom: '8px',
});

const MoleculeInfo = styled(Box)({
  marginBottom: '8px',
});

const MoleculeName = styled(Typography)({
  fontWeight: '500',
  fontSize: '1rem',
  marginBottom: '4px',
  overflow: 'hidden',
  textOverflow: 'ellipsis',
  whiteSpace: 'nowrap',
});

const MoleculeSmiles = styled(Typography)({
  fontSize: '0.75rem',
  color: theme.palette.text.secondary,
  overflow: 'hidden',
  textOverflow: 'ellipsis',
  whiteSpace: 'nowrap',
});

const PropertiesContainer = styled(Box)({
  display: 'flex',
  flexWrap: 'wrap',
  gap: '4px',
  marginBottom: '8px',
  flex: '1',
});

const ActionsContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  marginTop: 'auto',
});

/**
 * A component that displays a molecule in a card format with structure visualization,
 * properties, and actions. This component serves as a reusable building block for
 * displaying molecular information throughout the application.
 */
const MoleculeCard: React.FC<MoleculeCardProps> = ({
  molecule,
  selected = false,
  interactive = true,
  displayProperties = DEFAULT_MOLECULE_PROPERTIES,
  showActions = true,
  onClick,
  onMoleculeClick,
  onViewDetails,
  onAddToQueue,
  onFlagChange,
  className,
  style,
  structureWidth = 150,
  structureHeight = 150,
}) => {
  // Handle card click based on provided handlers
  const handleClick = () => {
    if (onClick) {
      onClick();
    } else if (onMoleculeClick) {
      onMoleculeClick(molecule);
    }
  };

  // Get molecule name from properties or fallback to SMILES
  const moleculeName = getMoleculeCommonName(molecule);
  
  // Determine if the molecule is flagged
  const isFlagged = molecule.flag_status === FlagStatus.IMPORTANT || 
                   molecule.flag_status === FlagStatus.FAVORITE;

  return (
    <Card
      selected={selected}
      interactive={interactive}
      onClick={interactive ? handleClick : undefined}
      elevation={2}
      className={className}
      style={style}
      fullHeight
    >
      <CardContent>
        <StructureContainer>
          <MoleculeStructure
            molecule={molecule}
            width={structureWidth}
            height={structureHeight}
          />
        </StructureContainer>
        
        <MoleculeInfo>
          <MoleculeName variant="subtitle1" title={moleculeName}>
            {moleculeName}
          </MoleculeName>
          <MoleculeSmiles variant="caption" title={molecule.smiles}>
            {molecule.smiles}
          </MoleculeSmiles>
        </MoleculeInfo>
        
        <PropertiesContainer>
          {displayProperties.map((propertyName) => {
            const value = getMoleculePropertyByName(molecule, propertyName);
            if (value !== null) {
              return (
                <PropertyBadge
                  key={propertyName}
                  propertyName={propertyName}
                  value={value}
                  showTooltip
                />
              );
            }
            return null;
          })}
        </PropertiesContainer>
        
        {showActions && (
          <ActionsContainer>
            <Box>
              {onViewDetails && (
                <Tooltip title="View details">
                  <IconButton 
                    size="small" 
                    onClick={(e) => {
                      e.stopPropagation();
                      onViewDetails(molecule);
                    }}
                  >
                    <Visibility fontSize="small" />
                  </IconButton>
                </Tooltip>
              )}
              
              {onAddToQueue && (
                <Tooltip title="Add to queue">
                  <IconButton 
                    size="small" 
                    onClick={(e) => {
                      e.stopPropagation();
                      onAddToQueue(molecule);
                    }}
                  >
                    <Add fontSize="small" />
                  </IconButton>
                </Tooltip>
              )}
            </Box>
            
            {onFlagChange && (
              <Tooltip title={isFlagged ? "Remove flag" : "Flag as important"}>
                <IconButton 
                  size="small" 
                  color={isFlagged ? "primary" : "default"}
                  onClick={(e) => {
                    e.stopPropagation();
                    onFlagChange(
                      molecule, 
                      isFlagged ? null : FlagStatus.IMPORTANT
                    );
                  }}
                >
                  {isFlagged ? (
                    <Star fontSize="small" />
                  ) : (
                    <StarBorder fontSize="small" />
                  )}
                </IconButton>
              </Tooltip>
            )}
          </ActionsContainer>
        )}
      </CardContent>
    </Card>
  );
};

export default MoleculeCard;