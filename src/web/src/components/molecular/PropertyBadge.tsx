import React from 'react';
import { Tooltip } from '@mui/material'; // v5.13+
import Badge from '../common/Badge';
import theme from '../../theme';
import { formatPropertyValue, getPropertyDisplayName, getPropertyUnit } from '../../utils/molecularUtils';

/**
 * Props interface for the PropertyBadge component
 */
interface PropertyBadgeProps {
  /** Name of the molecular property */
  propertyName: string;
  /** Value of the property */
  value: any;
  /** Optional unit override (if not using default unit for the property) */
  unit?: string;
  /** Whether to show tooltip with full property details */
  showTooltip?: boolean;
  /** Whether to show property name in the badge */
  showName?: boolean;
  /** Size of the badge */
  size?: 'small' | 'medium' | 'large';
  /** Custom inline styles */
  style?: React.CSSProperties;
  /** CSS class name */
  className?: string;
  /** Click handler */
  onClick?: () => void;
}

/**
 * Determines the appropriate color for a property badge based on property type and value
 * 
 * @param propertyName Name of the property
 * @param value Value of the property
 * @returns Color name from theme palette
 */
function getPropertyColor(propertyName: string, value: any): string {
  const propertyNameLower = propertyName.toLowerCase();
  
  // Handle LogP property - common medicinal chemistry metric
  if (propertyNameLower === 'logp') {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      if (numValue < 5) return 'success'; // Good drug-likeness (Lipinski's rule of 5)
      if (numValue <= 7) return 'warning'; // Borderline acceptable
      return 'error'; // Too lipophilic, likely poor drug-likeness
    }
  }
  
  // Handle Molecular Weight property - another important drug-likeness metric
  if (propertyNameLower === 'molecular_weight' || propertyNameLower === 'mw') {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      if (numValue < 500) return 'success'; // Good drug-likeness (Lipinski's rule of 5)
      if (numValue <= 800) return 'warning'; // Getting high but potentially acceptable
      return 'error'; // Exceeds typical drug-like range
    }
  }
  
  // Handle H-Bond donors and acceptors - important for drug-likeness
  if (propertyNameLower === 'h_bond_donors' || propertyNameLower === 'hbd') {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      return numValue <= 5 ? 'success' : 'warning'; // Lipinski's rule of 5
    }
  }
  
  if (propertyNameLower === 'h_bond_acceptors' || propertyNameLower === 'hba') {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      return numValue <= 10 ? 'success' : 'warning'; // Lipinski's rule of 5
    }
  }
  
  // Handle solubility - important for drug bioavailability
  if (propertyNameLower === 'solubility') {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      if (numValue > 100) return 'success'; // Highly soluble
      if (numValue > 10) return 'info'; // Moderately soluble
      if (numValue > 1) return 'warning'; // Poorly soluble
      return 'error'; // Very poorly soluble
    }
  }
  
  // Handle activity data - typically percentage values for things like inhibition
  if (propertyNameLower === 'activity' || propertyNameLower.includes('inhibition')) {
    const numValue = typeof value === 'string' ? parseFloat(value) : value;
    if (typeof numValue === 'number' && !isNaN(numValue)) {
      if (numValue > 80) return 'success'; // High activity
      if (numValue > 50) return 'info'; // Moderate activity
      if (numValue > 20) return 'warning'; // Low activity
      return 'error'; // Very low activity
    }
  }
  
  // Default color for other properties
  return 'primary';
}

/**
 * A specialized badge component for displaying molecular property values with appropriate formatting and styling
 * 
 * This component extends the base Badge component with molecular property-specific features like
 * automatic color coding based on property values, unit display, and tooltips.
 */
const PropertyBadge: React.FC<PropertyBadgeProps> = ({
  propertyName,
  value,
  unit = '',
  showTooltip = true,
  showName = false,
  size = 'small',
  style,
  className,
  onClick,
}) => {
  // Get display name for property using getPropertyDisplayName
  const displayName = getPropertyDisplayName(propertyName);
  
  // Get unit for property using getPropertyUnit if not provided in props
  const propertyUnit = unit || getPropertyUnit(propertyName);
  
  // Format property value using formatPropertyValue
  const formattedValue = formatPropertyValue(propertyName, value, Boolean(propertyUnit));
  
  // Determine badge color using getPropertyColor function
  const color = getPropertyColor(propertyName, value);
  
  // Construct label with property name and value based on showName prop
  const label = showName ? `${displayName}: ${formattedValue}` : formattedValue;
  
  // Create the base badge component
  const badge = (
    <Badge
      label={label}
      color={color}
      size={size}
      style={style}
      className={className}
      onClick={onClick}
    />
  );
  
  // If showTooltip is true, wrap Badge in Tooltip with full property details
  if (showTooltip) {
    return (
      <Tooltip title={`${displayName}: ${formattedValue}`}>
        {badge}
      </Tooltip>
    );
  }
  
  return badge;
};

export default PropertyBadge;