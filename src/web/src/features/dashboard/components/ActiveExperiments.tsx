import React, { useState, useEffect } from 'react'; // React v18.2+
import { Box, Typography, Button, Skeleton, List, ListItem, ListItemText, ListItemSecondaryAction, useTheme } from '@mui/material'; // @mui/material 5.13+
import { useNavigate } from 'react-router-dom'; // react-router-dom v6.10+
import Card from '../../../components/common/Card'; // Container component for the active experiments widget
import ProgressBar from '../../../components/common/ProgressBar'; // Visual indicator for experiment progress
import ExperimentStatus from '../../experiments/components/ExperimentStatus'; // Component to display experiment status with appropriate color coding
import { useExperiments } from '../../experiments/hooks/useExperiments'; // Hook for fetching and managing experiment data
import { ExperimentStatus as ExperimentStatusEnum } from '../../../types/experiment'; // Enum representing the possible status values for experiments

/**
 * Interface defining the props for the ActiveExperiments component
 */
interface ActiveExperimentsProps {
  /** Callback function to execute when the "View All" button is clicked */
  onViewAll?: () => void;
  /** Maximum number of experiments to display (default: 5) */
  maxItems?: number;
  /** Optional CSS class name for styling */
  className?: string;
}

/**
 * Helper function to calculate the progress percentage based on experiment status
 * @param status ExperimentStatus
 * @returns Progress percentage (0-100)
 */
const getExperimentProgress = (status: ExperimentStatusEnum): number => {
  switch (status) {
    case ExperimentStatusEnum.DRAFT:
      return 5;
    case ExperimentStatusEnum.QUEUED:
      return 10;
    case ExperimentStatusEnum.SUBMITTED:
      return 25;
    case ExperimentStatusEnum.QUOTE_PENDING:
      return 30;
    case ExperimentStatusEnum.IN_PROGRESS:
      return 50;
    case ExperimentStatusEnum.RESULTS_PENDING:
      return 75;
    case ExperimentStatusEnum.RESULTS_AVAILABLE:
      return 90;
    case ExperimentStatusEnum.COMPLETED:
      return 100;
    default:
      return 0;
  }
};

/**
 * Helper function to determine if an experiment is considered active
 * @param status ExperimentStatus
 * @returns Whether the experiment is active
 */
const isActiveExperiment = (status: ExperimentStatusEnum): boolean => {
  return status !== ExperimentStatusEnum.COMPLETED && status !== ExperimentStatusEnum.CANCELLED;
};

/**
 * Component that displays a list of active experiments with their status and progress
 * @param props 
 * @returns Rendered component
 */
const ActiveExperiments: React.FC<ActiveExperimentsProps> = (props) => {
  // Initialize navigate function using useNavigate hook
  const navigate = useNavigate();
  // Initialize theme using useTheme hook
  const theme = useTheme();

  // Destructure props with default values
  const { onViewAll, maxItems = 5, className } = props;

  // Get experiments data, loading state, and fetchExperiments function from useExperiments hook
  const { experiments, loading, error, fetchExperiments } = useExperiments();

  // Set up filter to only show active experiments (not completed or cancelled)
  const [activeExperiments, setActiveExperiments] = useState<Experiment[]>([]);

  // Use useEffect to fetch experiments on component mount
  useEffect(() => {
    fetchExperiments();
  }, [fetchExperiments]);

  // Use useEffect to filter experiments and limit the number of items shown
  useEffect(() => {
    // Filter experiments to only include active ones
    const active = experiments.filter(experiment => isActiveExperiment(experiment.status));
    // Limit the number of experiments shown based on maxItems prop
    setActiveExperiments(active.slice(0, maxItems));
  }, [experiments, maxItems]);

  /**
   * Handler for when an experiment is clicked
   * @param experimentId 
   */
  const handleExperimentClick = (experimentId: string) => {
    // Navigate to the experiment detail page using the experiment ID
    navigate(`/experiments/${experimentId}`);
  };

  /**
   * Handler for when the View All button is clicked
   */
  const handleViewAll = () => {
    // If onViewAll prop is provided, call it
    if (onViewAll) {
      onViewAll();
    } else {
      // Otherwise, navigate to the experiments page
      navigate('/experiments');
    }
  };

  return (
    // Render a Card component as the container
    <Card className={className} fullHeight>
      {/* Render a header with title 'Active Experiments' and experiment count */}
      <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', padding: theme.spacing(2) }}>
        <Typography variant="h6">
          Active Experiments ({activeExperiments.length})
        </Typography>
      </Box>
      {/* If loading, render Skeleton components as placeholders */}
      {loading ? (
        <List>
          {Array.from(new Array(maxItems)).map((_, index) => (
            <ListItem key={index} alignItems="flex-start">
              <Skeleton variant="text" width={120} />
              <Skeleton variant="rectangular" width={80} height={10} sx={{ ml: 2 }} />
            </ListItem>
          ))}
        </List>
      ) : error ? (
        // If error, render an error message
        <Typography color="error" sx={{ p: 2 }}>
          Error: {error}
        </Typography>
      ) : activeExperiments.length === 0 ? (
        // If no active experiments, render a message indicating no active experiments
        <Typography sx={{ p: 2 }}>No active experiments.</Typography>
      ) : (
        // Otherwise, render a List of experiments
        <List>
          {/* For each experiment, render a ListItem with experiment name, type, status badge, and progress bar */}
          {activeExperiments.map((experiment) => (
            <ListItem
              key={experiment.id}
              alignItems="flex-start"
              button
              onClick={() => handleExperimentClick(experiment.id)}
            >
              <ListItemText
                primary={experiment.name}
                secondary={
                  <React.Fragment>
                    <Typography
                      sx={{ display: 'inline' }}
                      component="span"
                      variant="body2"
                      color="text.primary"
                    >
                      {experiment.type?.name}
                    </Typography>
                    {` — Progress: `}
                    <ProgressBar
                      value={getExperimentProgress(experiment.status)}
                      showPercentage
                      size="small"
                    />
                  </React.Fragment>
                }
              />
              <ListItemSecondaryAction>
                <ExperimentStatus status={experiment.status} />
              </ListItemSecondaryAction>
            </ListItem>
          ))}
        </List>
      )}
      {/* Add a View All button in the card footer */}
      <Box sx={{ p: 2, display: 'flex', justifyContent: 'flex-end' }}>
        <Button onClick={handleViewAll}>View All</Button>
      </Box>
    </Card>
  );
};

export default ActiveExperiments;