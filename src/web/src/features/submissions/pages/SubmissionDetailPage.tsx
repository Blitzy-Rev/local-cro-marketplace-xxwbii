import React, { useState, useEffect, useCallback } from 'react'; // React 18.2+
import { useParams, useNavigate, Link } from 'react-router-dom'; // react-router-dom v6.11.0
import { Box, Typography, Divider, Grid, Paper, Chip, Stack, Alert, Container } from '@mui/material'; // @mui/material v5.13.0
import { styled } from '@mui/material/styles'; // @mui/material/styles v5.13.0
import { Science, Assignment, Cancel, ArrowBack, Visibility } from '@mui/icons-material'; // @mui/icons-material v5.13.0
import Card from '../../../components/common/Card';
import Button from '../../../components/common/Button';
import Dialog from '../../../components/common/Dialog';
import Badge from '../../../components/common/Badge';
import SubmissionTimeline from '../components/SubmissionTimeline';
import QuoteApproval from '../components/QuoteApproval';
import useSubmissions from '../hooks/useSubmissions';
import useToast from '../../../hooks/useToast';
import { Submission, SubmissionDetailed, SubmissionStatus } from '../../../types/submission';
import useAuth from '../../auth/hooks/useAuth';
import { formatDate, formatCurrency } from '../../../utils/formatters';
import MainLayout from '../../../layouts/MainLayout';
import Loading from '../../../components/common/Loading';

/**
 * Interface defining the props for the SubmissionDetailPage component.
 */
interface SubmissionDetailPageProps {
  // No props are explicitly passed to this component
}

/**
 * Styled component for the page container.
 */
const PageContainer = styled(Container)(({ theme }) => ({
  maxWidth: 'lg',
  py: 3,
}));

/**
 * Styled component for the page header.
 */
const PageHeader = styled(Box)(({ theme }) => ({
  display: 'flex',
  justifyContent: 'space-between',
  alignItems: 'center',
  mb: 3,
}));

/**
 * Styled component for the content paper.
 */
const ContentPaper = styled(Paper)(({ theme }) => ({
  p: 3,
  mb: 3,
}));

/**
 * Page component that displays detailed information about a submission to a CRO.
 */
const SubmissionDetailPage: React.FC<SubmissionDetailPageProps> = () => {
  // LD1: Extract submission ID from URL parameters using useParams
  const { id } = useParams<{ id: string }>();

  // LD1: Initialize navigate function using useNavigate
  const navigate = useNavigate();

  // LD1: Set up state for loading status and error handling
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // LD1: Get submission operations from useSubmissions hook
  const { fetchSubmissionDetails } = useSubmissions();

  /**
   * LD1: Handles navigation back to submissions list.
   */
  const handleBack = useCallback(() => {
    navigate('/app/submissions');
  }, [navigate]);

  /**
   * LD1: Handles submission status changes by refreshing data.
   */
  const handleStatusChange = useCallback(() => {
    if (id) {
      fetchSubmissionDetails(id);
    }
  }, [fetchSubmissionDetails, id]);

  // LD1: Fetch submission details when component mounts or ID changes
  useEffect(() => {
    if (id) {
      setLoading(true);
      setError(null);
      fetchSubmissionDetails(id)
        .catch((err) => {
          setError(err?.message || 'Failed to fetch submission details.');
        })
        .finally(() => {
          setLoading(false);
        });
    }
  }, [fetchSubmissionDetails, id]);

  // LD1: Render MainLayout as container
  return (
    <MainLayout>
      {/* LD1: Render PageContainer with appropriate width and padding */}
      <PageContainer maxWidth="lg" sx={{ py: 3 }}>
        {/* LD1: Render PageHeader with title 'Submission Details' */}
        <PageHeader sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', mb: 3 }}>
          <Typography variant="h4" component="h1">
            Submission Details
          </Typography>
        </PageHeader>

        {/* LD1: Render loading indicator when data is being fetched */}
        {loading && <Loading message="Loading submission details..." />}

        {/* LD1: Render error message if fetch failed */}
        {error && <Alert severity="error">{error}</Alert>}

        {/* LD1: Render SubmissionDetail component with:
            - submissionId prop set to the ID from URL params
            - onBack prop set to handleBack function
            - onStatusChange prop set to handleStatusChange function */}
        {!loading && !error && (
          <SubmissionDetail
            submissionId={id}
            onBack={handleBack}
            onStatusChange={handleStatusChange}
          />
        )}
      </PageContainer>
    </MainLayout>
  );
};

export default SubmissionDetailPage;