import React, { useState, useEffect } from 'react'; // react v18.2+
import { Grid, Box, Typography, Skeleton, useTheme, useMediaQuery } from '@mui/material'; // @mui/material v5.13+
import { ScienceOutlined, LibraryBooksOutlined, AssignmentOutlined, PeopleOutlined } from '@mui/icons-material'; // @mui/icons-material v5.13+
import { useQuery } from 'react-query'; // react-query v4.28+
import Card from '../../../components/common/Card';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Component that displays key statistics based on user role
 * @returns Rendered statistics component
 */
const Stats: React.FC = () => {
  const theme = useTheme();
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));
  const { isPharma, isCRO, isAdmin } = usePermissions();

  // Define types for the statistics data
  type PharmaStats = {
    moleculeCount: number;
    libraryCount: number;
    experimentCount: number;
    resultCount: number;
  };

  type CROStats = {
    submissionCount: number;
    inProgressCount: number;
    completedCount: number;
    pendingReviewCount: number;
  };

  type AdminStats = {
    userCount: number;
    moleculeCount: number;
    experimentCount: number;
    systemMetrics: any; // Replace 'any' with a more specific type if possible
  };

  // Define the API URLs
  const API_URLS = {
    pharmaStats: '/stats/pharma',
    croStats: '/stats/cro',
    adminStats: '/stats/admin'
  };

  /**
   * Fetches statistics for Pharma users
   * @returns Promise resolving to pharma statistics
   */
  const fetchPharmaStats = async (): Promise<PharmaStats> => {
    try {
      // Simulate API calls to get molecule count, library count, active experiment count, and result count
      const moleculeCount = Math.floor(Math.random() * 1000);
      const libraryCount = Math.floor(Math.random() * 100);
      const experimentCount = Math.floor(Math.random() * 50);
      const resultCount = Math.floor(Math.random() * 200);

      return { moleculeCount, libraryCount, experimentCount, resultCount };
    } catch (error) {
      console.error('Error fetching Pharma stats:', error);
      return { moleculeCount: 0, libraryCount: 0, experimentCount: 0, resultCount: 0 };
    }
  };

  /**
   * Fetches statistics for CRO users
   * @returns Promise resolving to CRO statistics
   */
  const fetchCROStats = async (): Promise<CROStats> => {
    try {
      // Simulate API calls to get submission count, in-progress count, completed count, and pending review count
      const submissionCount = Math.floor(Math.random() * 500);
      const inProgressCount = Math.floor(Math.random() * 100);
      const completedCount = Math.floor(Math.random() * 300);
      const pendingReviewCount = Math.floor(Math.random() * 50);

      return { submissionCount, inProgressCount, completedCount, pendingReviewCount };
    } catch (error) {
      console.error('Error fetching CRO stats:', error);
      return { submissionCount: 0, inProgressCount: 0, completedCount: 0, pendingReviewCount: 0 };
    }
  };

  /**
   * Fetches statistics for Admin users
   * @returns Promise resolving to admin statistics
   */
  const fetchAdminStats = async (): Promise<AdminStats> => {
    try {
      // Simulate API call to admin statistics endpoint
      const userCount = Math.floor(Math.random() * 20);
      const moleculeCount = Math.floor(Math.random() * 1500);
      const experimentCount = Math.floor(Math.random() * 75);
      const systemMetrics = {
        cpuUsage: Math.random() * 100,
        memoryUsage: Math.random() * 100,
        diskUsage: Math.random() * 100
      };

      return { userCount, moleculeCount, experimentCount, systemMetrics };
    } catch (error) {
      console.error('Error fetching Admin stats:', error);
      return { userCount: 0, moleculeCount: 0, experimentCount: 0, systemMetrics: {} };
    }
  };

  // Fetch data based on user role
  const { data: pharmaStats, isLoading: isPharmaStatsLoading } = useQuery('pharmaStats', fetchPharmaStats, { enabled: isPharma() });
  const { data: croStats, isLoading: isCROStatsLoading } = useQuery('croStats', fetchCROStats, { enabled: isCRO() });
  const { data: adminStats, isLoading: isAdminStatsLoading } = useQuery('adminStats', fetchAdminStats, { enabled: isAdmin() });

  /**
   * Renders a single statistic card
   * @param title 
   * @param icon 
   * @param value 
   * @param isLoading 
   * @returns Rendered statistic card
   */
  const renderStatCard = (title: string, icon: React.ReactNode, value: number | string, isLoading: boolean) => (
    <Card interactive={false}>
      <Box display="flex" alignItems="center" mb={2}>
        <Box mr={1}>{icon}</Box>
        <Typography variant="subtitle1" fontWeight="bold">{title}</Typography>
      </Box>
      {isLoading ? (
        <Skeleton variant="text" width={100} height={30} />
      ) : (
        <Typography variant="h4" fontWeight="bold">{value}</Typography>
      )}
    </Card>
  );

  return (
    <Grid container spacing={isMobile ? 2 : 3}>
      {isPharma() && (
        <>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Molecules', <ScienceOutlined />, pharmaStats?.moleculeCount ?? 0, isPharmaStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Libraries', <LibraryBooksOutlined />, pharmaStats?.libraryCount ?? 0, isPharmaStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Active Experiments', <ScienceOutlined />, pharmaStats?.experimentCount ?? 0, isPharmaStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Results', <AssignmentOutlined />, pharmaStats?.resultCount ?? 0, isPharmaStatsLoading)}
          </Grid>
        </>
      )}

      {isCRO() && (
        <>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Submissions', <AssignmentOutlined />, croStats?.submissionCount ?? 0, isCROStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('In Progress', <ScienceOutlined />, croStats?.inProgressCount ?? 0, isCROStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Completed', <AssignmentOutlined />, croStats?.completedCount ?? 0, isCROStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Pending Review', <AssignmentOutlined />, croStats?.pendingReviewCount ?? 0, isCROStatsLoading)}
          </Grid>
        </>
      )}

      {isAdmin() && (
        <>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Users', <PeopleOutlined />, adminStats?.userCount ?? 0, isAdminStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Molecules', <ScienceOutlined />, adminStats?.moleculeCount ?? 0, isAdminStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('Experiments', <ScienceOutlined />, adminStats?.experimentCount ?? 0, isAdminStatsLoading)}
          </Grid>
          <Grid item xs={12} sm={6} md={3}>
            {renderStatCard('System Metrics', <AssignmentOutlined />, 'N/A', isAdminStatsLoading)}
          </Grid>
        </>
      )}
    </Grid>
  );
};

export default Stats;