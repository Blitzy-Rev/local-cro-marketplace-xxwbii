import React, { useState, useEffect } from 'react'; // React 18.2+
import { Box, Typography, Paper, Stepper, Step, StepLabel, StepContent, Button } from '@mui/material'; // v5.13+
import { styled } from '@mui/material/styles'; // v5.13+
import { CloudUpload, Map, CheckCircle } from '@mui/icons-material'; // v5.13+
import FileUploader from '../components/FileUploader';
import ColumnMapper from '../components/ColumnMapper';
import ImportSummary from '../components/ImportSummary';
import { useCSVImport } from '../hooks/useCSVImport';
import { useToast } from '../../../hooks/useToast';

/**
 * Styled component for the page container
 */
const PageContainer = styled(Box)({
  maxWidth: '1200px',
  margin: '0 auto',
  padding: '24px',
});

/**
 * Styled component for the page title
 */
const PageTitle = styled(Typography)({
  marginBottom: '16px',
  fontWeight: '500',
});

/**
 * Styled component for the page description
 */
const PageDescription = styled(Typography)(({ theme }) => ({
  marginBottom: '32px',
  color: theme.palette.text.secondary,
}));

/**
 * Styled component for the stepper container
 */
const StepperContainer = styled(Paper)({
  padding: '24px',
  marginBottom: '24px',
  borderRadius: '8px',
});

/**
 * Styled component for the button container
 */
const ButtonContainer = styled(Box)({
  display: 'flex',
  justifyContent: 'space-between',
  marginTop: '24px',
});

/**
 * Page component for CSV upload workflow
 * @returns JSX.Element
 */
const CSVUploadPage: React.FC = () => {
  // Extract necessary state and functions from useCSVImport hook
  const { status, summary, handleReset } = useCSVImport();

  // Initialize state for current step in the workflow
  const [activeStep, setActiveStep] = useState<number>(0);

  // Initialize state for mapping validity
  const [isMappingValid, setIsMappingValid] = useState<boolean>(false);

  // Get the showToast function from the useToast hook
  const { showToast } = useToast();

  /**
   * Advances to the mapping step after successful file upload
   * @param headers 
   * @param rowCount 
   */
  const handleUploadComplete = (headers: string[], rowCount: number) => {
    setActiveStep(1); // Set activeStep to 1 (mapping step)
    showToast({
      type: 'success',
      message: `File uploaded successfully with ${rowCount} rows and ${headers.length} columns`,
    });
  };

  /**
   * Advances to the results step after successful mapping submission
   */
  const handleMappingComplete = () => {
    setActiveStep(2); // Set activeStep to 2 (results step)
    showToast({
      type: 'success',
      message: 'Processing started',
    });
  };

  /**
   * Updates mapping validity state
   * @param isValid 
   */
  const handleMappingChange = (isValid: boolean) => {
    setIsMappingValid(isValid); // Set isMappingValid to the provided isValid value
  };

  /**
   * Resets the workflow and returns to the first step
   */
  const handleStartNew = () => {
    handleReset(); // Call handleReset from useCSVImport hook
    setActiveStep(0); // Set activeStep to 0 (upload step)
    setIsMappingValid(false);
  };

  /**
   * Goes back to the previous step
   */
  const handleBack = () => {
    setActiveStep((prevActiveStep) => prevActiveStep - 1); // Decrement activeStep by 1
  };

  // Automatically advance to results step when processing completes
  useEffect(() => {
    if (status === 'completed' && summary) {
      setActiveStep(2); // Set activeStep to 2 (results step)
      showToast({
        type: 'success',
        message: 'Processing completed',
      });
    }
  }, [status, summary, showToast]);

  // Define steps array with labels and icons for the stepper
  const steps = [
    {
      label: 'Upload',
      icon: <CloudUpload />,
    },
    {
      label: 'Map',
      icon: <Map />,
    },
    {
      label: 'Results',
      icon: <CheckCircle />,
    },
  ];

  return (
    <PageContainer>
      {/* Render PageContainer with title 'Import Molecules from CSV' */}
      <PageTitle variant="h4">Import Molecules from CSV</PageTitle>

      {/* Render PageDescription with instructions for the CSV import process */}
      <PageDescription variant="body1">
        Upload a CSV file containing molecular data. Map the columns to the
        appropriate system properties and review the import summary.
      </PageDescription>

      {/* Render StepperContainer with Material-UI Stepper component */}
      <StepperContainer>
        <Stepper activeStep={activeStep} orientation="vertical">
          {steps.map((step, index) => (
            <Step key={step.label}>
              {/* Render Step and StepLabel with appropriate icon */}
              <StepLabel icon={step.icon}>{step.label}</StepLabel>
              <StepContent>
                {/* Render appropriate component based on current step */}
                {index === 0 && (
                  <FileUploader onUploadComplete={handleUploadComplete} />
                )}
                {index === 1 && (
                  <ColumnMapper onMappingComplete={handleMappingComplete} onMappingChange={handleMappingChange} />
                )}
                {index === 2 && (
                  <ImportSummary onStartNew={handleStartNew} />
                )}
              </StepContent>
            </Step>
          ))}
        </Stepper>
      </StepperContainer>

      {/* Render ButtonContainer with navigation buttons */}
      <ButtonContainer>
        <Button disabled={activeStep === 0} onClick={handleBack}>
          Back
        </Button>
        <div>
          {activeStep === steps.length - 1 ? (
            <Button variant="contained" color="primary" onClick={handleStartNew}>
              Start New
            </Button>
          ) : (
            <Button variant="contained" color="primary" disabled={activeStep === 1 && !isMappingValid} onClick={handleSubmitMapping}>
              Next
            </Button>
          )}
        </div>
      </ButtonContainer>
    </PageContainer>
  );
};

export default CSVUploadPage;