import React from 'react'; // React v18.2+
import Badge from '../../../components/common/Badge'; // Reusable badge component for displaying status indicators
import { ExperimentStatus as ExperimentStatusEnum } from '../../../types/experiment'; // Import experiment status enum for type checking and status values

/**
 * Props for the ExperimentStatus component
 */
interface ExperimentStatusProps {
  /** The status of the experiment */
  status: ExperimentStatusEnum;
  /** Size of the status badge */
  size?: 'small' | 'medium' | 'large';
  /** Visual style of the badge */
  variant?: 'filled' | 'outlined';
  /** Additional CSS class */
  className?: string;
  /** Inline styles */
  style?: React.CSSProperties;
}

/**
 * Determines the appropriate color and label for a given experiment status
 */
const getStatusConfig = (status: ExperimentStatusEnum): { color: string, label: string } => {
  switch (status) {
    case ExperimentStatusEnum.DRAFT:
      return { color: 'default', label: 'Draft' };
    case ExperimentStatusEnum.QUEUED:
      return { color: 'info', label: 'Queued' };
    case ExperimentStatusEnum.SUBMITTED:
      return { color: 'primary', label: 'Submitted' };
    case ExperimentStatusEnum.REJECTED:
      return { color: 'error', label: 'Rejected' };
    case ExperimentStatusEnum.QUOTE_PENDING:
      return { color: 'warning', label: 'Quote Pending' };
    case ExperimentStatusEnum.QUOTE_REJECTED:
      return { color: 'error', label: 'Quote Rejected' };
    case ExperimentStatusEnum.IN_PROGRESS:
      return { color: 'primary', label: 'In Progress' };
    case ExperimentStatusEnum.RESULTS_PENDING:
      return { color: 'warning', label: 'Results Pending' };
    case ExperimentStatusEnum.RESULTS_AVAILABLE:
      return { color: 'success', label: 'Results Available' };
    case ExperimentStatusEnum.RESULTS_REJECTED:
      return { color: 'error', label: 'Results Rejected' };
    case ExperimentStatusEnum.COMPLETED:
      return { color: 'success', label: 'Completed' };
    case ExperimentStatusEnum.CANCELLED:
      return { color: 'error', label: 'Cancelled' };
    default:
      return { color: 'default', label: 'Unknown' };
  }
};

/**
 * Component that displays an experiment status as a colored badge
 * 
 * Provides a visual indicator of where an experiment is in its workflow
 * with appropriate color coding based on status
 */
const ExperimentStatus: React.FC<ExperimentStatusProps> = ({
  status,
  size = 'small',
  variant = 'filled',
  className,
  style,
}) => {
  const { color, label } = getStatusConfig(status);
  
  return (
    <Badge
      label={label}
      color={color as 'default' | 'primary' | 'secondary' | 'success' | 'error' | 'warning' | 'info'}
      size={size}
      variant={variant}
      className={className}
      style={style}
    />
  );
};

export default ExperimentStatus;