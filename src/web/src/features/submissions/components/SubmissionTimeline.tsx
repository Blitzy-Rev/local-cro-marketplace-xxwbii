import React, { useMemo } from 'react';
import { 
  Timeline, 
  TimelineItem, 
  TimelineSeparator, 
  TimelineConnector, 
  TimelineContent, 
  TimelineDot, 
  TimelineOppositeContent 
} from '@mui/lab';
import { Typography, Box, Paper, styled } from '@mui/material';
import { 
  CheckCircle, 
  Cancel, 
  AttachMoney, 
  Pending, 
  HourglassEmpty, 
  Assignment, 
  Science 
} from '@mui/icons-material';

import Card from '../../../components/common/Card';
import Badge from '../../../components/common/Badge';
import { Submission, SubmissionStatus } from '../../../types/submission';
import { formatDate } from '../../../utils/formatters';

/**
 * Props for the SubmissionTimeline component
 */
interface SubmissionTimelineProps {
  /** The submission object containing status history */
  submission: Submission;
  /** Additional CSS class name for styling */
  className?: string;
}

/**
 * Internal interface for timeline event data
 */
interface TimelineEvent {
  /** Status of the submission at this point in time */
  status: SubmissionStatus;
  /** Date when the status changed */
  date: Date;
  /** Description of the status change event */
  description: string;
}

/**
 * Styled container for the timeline
 */
const TimelineContainer = styled(Box)(({ theme }) => ({
  width: '100%',
  marginTop: theme.spacing(2),
  marginBottom: theme.spacing(2)
}));

/**
 * Styled content container for timeline items
 */
const TimelineItemContent = styled(Paper)(({ theme }) => ({
  padding: theme.spacing(2),
  backgroundColor: theme.palette.background.default
}));

/**
 * Returns the appropriate icon component for a given submission status
 * @param status Submission status
 * @returns Icon component representing the status
 */
const getStatusIcon = (status: SubmissionStatus): React.ReactElement => {
  switch (status) {
    case SubmissionStatus.APPROVED:
    case SubmissionStatus.COMPLETED:
      return <CheckCircle />;
    case SubmissionStatus.REJECTED:
    case SubmissionStatus.CANCELLED:
      return <Cancel />;
    case SubmissionStatus.QUOTE_PROVIDED:
      return <AttachMoney />;
    case SubmissionStatus.QUOTE_REJECTED:
      return <Cancel />;
    case SubmissionStatus.IN_PROGRESS:
      return <HourglassEmpty />;
    case SubmissionStatus.PENDING:
      return <Pending />;
    default:
      return <Assignment />;
  }
};

/**
 * Determines the appropriate color for a status dot in the timeline
 * @param status Submission status
 * @returns Color name (success, error, warning, info, etc.)
 */
const getStatusColor = (status: SubmissionStatus): string => {
  switch (status) {
    case SubmissionStatus.APPROVED:
    case SubmissionStatus.COMPLETED:
      return 'success';
    case SubmissionStatus.REJECTED:
    case SubmissionStatus.CANCELLED:
      return 'error';
    case SubmissionStatus.QUOTE_PROVIDED:
    case SubmissionStatus.QUOTE_REJECTED:
      return 'warning';
    case SubmissionStatus.IN_PROGRESS:
      return 'info';
    case SubmissionStatus.PENDING:
    default:
      return 'grey';
  }
};

/**
 * Converts submission status enum value to a user-friendly display label
 * @param status Submission status enum value
 * @returns Human-readable status label
 */
const getStatusLabel = (status: SubmissionStatus): string => {
  switch (status) {
    case SubmissionStatus.PENDING:
      return 'Submitted';
    case SubmissionStatus.REJECTED:
      return 'Rejected';
    case SubmissionStatus.QUOTE_PROVIDED:
      return 'Quote Received';
    case SubmissionStatus.QUOTE_REJECTED:
      return 'Quote Rejected';
    case SubmissionStatus.APPROVED:
      return 'Approved';
    case SubmissionStatus.IN_PROGRESS:
      return 'In Progress';
    case SubmissionStatus.COMPLETED:
      return 'Completed';
    case SubmissionStatus.CANCELLED:
      return 'Cancelled';
    default:
      return 'Unknown';
  }
};

/**
 * Generates timeline events from submission history data
 * @param submission The submission object with status history
 * @returns Array of timeline events with status, date, and description
 */
const generateTimelineEvents = (submission: Submission): TimelineEvent[] => {
  if (!submission) return [];
  
  // Initialize events array with submission creation event
  const events: TimelineEvent[] = [
    {
      status: SubmissionStatus.PENDING,
      date: new Date(submission.submitted_at),
      description: 'Submission created'
    }
  ];
  
  // Extract status history if available
  const statusHistory = (submission as any).statusHistory;
  if (statusHistory && Array.isArray(statusHistory)) {
    // Sort history by date
    const historyEvents = statusHistory
      .sort((a: any, b: any) => new Date(a.date).getTime() - new Date(b.date).getTime())
      .map((history: any) => ({
        status: history.status,
        date: new Date(history.date),
        description: history.description || `Status changed to ${getStatusLabel(history.status)}`
      }));
    
    events.push(...historyEvents);
  }
  // If no history is available but status has changed from initial
  else if (submission.status !== SubmissionStatus.PENDING) {
    events.push({
      status: submission.status,
      date: new Date(submission.updated_at || submission.submitted_at),
      description: `Status updated to ${getStatusLabel(submission.status)}`
    });
  }
  
  // Add current status as the latest event if not already included
  const lastEvent = events[events.length - 1];
  if (lastEvent && lastEvent.status !== submission.status) {
    events.push({
      status: submission.status,
      date: new Date(submission.updated_at || submission.submitted_at),
      description: `Current status: ${getStatusLabel(submission.status)}`
    });
  }
  
  // Ensure events are sorted by date
  return events.sort((a, b) => a.date.getTime() - b.date.getTime());
};

/**
 * Component that displays a timeline visualization of submission status history
 * It shows the progression of a submission through various status stages from creation to completion,
 * with timestamps and status indicators.
 */
const SubmissionTimeline: React.FC<SubmissionTimelineProps> = ({ 
  submission, 
  className = '' 
}) => {
  // If submission is null or undefined, render nothing
  if (!submission) {
    return null;
  }

  // Generate timeline events from submission data
  const timelineEvents = useMemo(() => generateTimelineEvents(submission), [submission]);

  return (
    <Card className={className}>
      <TimelineContainer>
        <Timeline position="alternate">
          {timelineEvents.map((event, index) => (
            <TimelineItem key={`${event.status}-${index}`}>
              <TimelineOppositeContent color="text.secondary">
                {formatDate(event.date)}
              </TimelineOppositeContent>
              <TimelineSeparator>
                <TimelineDot color={getStatusColor(event.status)}>
                  {getStatusIcon(event.status)}
                </TimelineDot>
                {index < timelineEvents.length - 1 && <TimelineConnector />}
              </TimelineSeparator>
              <TimelineContent>
                <TimelineItemContent elevation={1}>
                  <Typography variant="h6" component="h3">
                    {getStatusLabel(event.status)}
                  </Typography>
                  <Typography>{event.description}</Typography>
                </TimelineItemContent>
              </TimelineContent>
            </TimelineItem>
          ))}
        </Timeline>
      </TimelineContainer>
    </Card>
  );
};

export default SubmissionTimeline;