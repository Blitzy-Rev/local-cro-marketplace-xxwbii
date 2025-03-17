# src/web/src/features/dashboard/components/QuickActions.tsx
```typescript
import React from 'react'; // React 18.2+
import { useNavigate } from 'react-router-dom'; // v6.10+
import { Grid, Box, Typography, useTheme, useMediaQuery } from '@mui/material'; // v5.13+
import { CloudUploadOutlined, LibraryAddOutlined, ScienceOutlined, AssignmentOutlined, UploadFileOutlined, UpdateOutlined, PeopleOutlined, SettingsOutlined, BackupOutlined } from '@mui/icons-material'; // v5.13+
import Button from '../../../components/common/Button';
import Card from '../../../components/common/Card';
import { usePermissions } from '../../../hooks/usePermissions';

/**
 * Component that displays quick action buttons based on user role
 * @returns Rendered quick actions component
 */
const QuickActions: React.FC = () => {
  // Initialize navigate function using useNavigate hook
  const navigate = useNavigate();

  // Initialize theme using useTheme hook
  const theme = useTheme();

  // Check if screen is mobile using useMediaQuery hook
  const isMobile = useMediaQuery(theme.breakpoints.down('sm'));

  // Get user permissions using usePermissions hook
  const { isPharma, isCRO, isAdmin, canAccessMolecules, canManageLibraries, canSubmitExperiments, canReviewSubmissions, canUploadResults, canManageUsers } = usePermissions();

  /**
   * Handles navigation to the CSV upload page
   */
  const handleNavigateToCSVUpload = () => {
    navigate('/app/molecules/upload');
  };

  /**
   * Handles navigation to the create library page
   */
  const handleNavigateToCreateLibrary = () => {
    navigate('/app/libraries/create');
  };

  /**
   * Handles navigation to the create experiment page
   */
  const handleNavigateToCreateExperiment = () => {
    navigate('/app/experiments/create');
  };

  /**
   * Handles navigation to the submissions page for CRO users
   */
  const handleNavigateToSubmissions = () => {
    navigate('/cro/submissions');
  };

  /**
   * Handles navigation to the result upload page for CRO users
   */
  const handleNavigateToResultUpload = () => {
    navigate('/cro/results/upload');
  };

  /**
   * Handles navigation to the user management page for admin users
   */
  const handleNavigateToUserManagement = () => {
    navigate('/admin/users');
  };

  /**
   * Handles navigation to the system configuration page for admin users
   */
  const handleNavigateToSystemConfig = () => {
    navigate('/admin/system');
  };

  /**
   * Renders a single action button with icon and label
   * @param label 
   * @param icon 
   * @param onClick 
   * @param disabled 
   * @returns Rendered action button
   */
  const renderActionButton = (
    label: string,
    icon: React.ReactNode,
    onClick: () => void,
    disabled: boolean = false
  ) => (
    <Button
      variant="outlined"
      color="primary"
      startIcon={icon}
      onClick={onClick}
      disabled={disabled}
      fullWidth={isMobile}
      style={{ marginBottom: theme.spacing(2) }}
    >
      {label}
    </Button>
  );

  // Render card with title 'QUICK ACTIONS'
  return (
    <Card title="QUICK ACTIONS">
      <Box sx={{ padding: 2 }}>
        <Grid container spacing={2} direction="column">
          {/* Render grid of action buttons based on user role and permissions */}
          {/* For Pharma users: show CSV Upload, Create Library, and Create Experiment buttons */}
          {isPharma() && canAccessMolecules() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Upload CSV',
                <CloudUploadOutlined />,
                handleNavigateToCSVUpload
              )}
            </Grid>
          )}
          {isPharma() && canManageLibraries() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Create Library',
                <LibraryAddOutlined />,
                handleNavigateToCreateLibrary
              )}
            </Grid>
          )}
          {isPharma() && canSubmitExperiments() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Create Experiment',
                <ScienceOutlined />,
                handleNavigateToCreateExperiment
              )}
            </Grid>
          )}

          {/* For CRO users: show Review Submissions, Upload Results, and Update Status buttons */}
          {isCRO() && canReviewSubmissions() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Review Submissions',
                <AssignmentOutlined />,
                handleNavigateToSubmissions
              )}
            </Grid>
          )}
          {isCRO() && canUploadResults() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Upload Results',
                <UploadFileOutlined />,
                handleNavigateToResultUpload
              )}
            </Grid>
          )}
          {isCRO() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Update Status',
                <UpdateOutlined />,
                () => console.log('Update Status clicked')
              )}
            </Grid>
          )}

          {/* For Admin users: show Add User, System Config, and Backup System buttons */}
          {isAdmin() && canManageUsers() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Manage Users',
                <PeopleOutlined />,
                handleNavigateToUserManagement
              )}
            </Grid>
          )}
          {isAdmin() && (
            <Grid item xs={12}>
              {renderActionButton(
                'System Config',
                <SettingsOutlined />,
                handleNavigateToSystemConfig
              )}
            </Grid>
          )}
          {isAdmin() && (
            <Grid item xs={12}>
              {renderActionButton(
                'Backup System',
                <BackupOutlined />,
                () => console.log('Backup System clicked')
              )}
            </Grid>
          )}
        </Grid>
      </Box>
    </Card>
  );
};

export default QuickActions;