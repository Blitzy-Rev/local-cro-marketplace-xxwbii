import React, { useState, useEffect } from 'react'; // React 18.2+
import { Box, Typography, Paper, Grid, Divider } from '@mui/material'; // @mui/material v5.13.0
import AttachMoney from '@mui/icons-material/AttachMoney'; // @mui/icons-material v5.13.0
import CalendarToday from '@mui/icons-material/CalendarToday'; // @mui/icons-material v5.13.0
import Notes from '@mui/icons-material/Notes'; // @mui/icons-material v5.13.0
import Input from '../../../components/common/Input';
import Button from '../../../components/common/Button';
import useCROInterface from '../hooks/useCROInterface';
import useToast from '../../../hooks/useToast';
import { Submission, QuoteProvide } from '../../../types/submission';

/**
 * Interface defining the props for the QuoteForm component
 */
interface QuoteFormProps {
  /** The submission data for which to provide a quote */
  submission: Submission;
  /** Callback function called after successful quote submission */
  onSubmitted: () => void;
  /** Callback function called when quote form is cancelled */
  onCancel: () => void;
}

/**
 * A form component for CRO users to provide pricing quotes for experiment submissions
 * @param {QuoteFormProps} props - The props for the component
 * @returns {JSX.Element} Rendered component
 */
const QuoteForm: React.FC<QuoteFormProps> = ({ submission, onSubmitted, onCancel }) => {
  // 1. Extract submission, onSubmitted, and onCancel from props
  
  // 2. Initialize form state with price, turnaroundDays, and notes
  const [formData, setFormData] = useState<QuoteProvide>({
    price: '',
    turnaround_days: '',
    notes: '',
  });

  // 3. Initialize validation state for form fields
  const [errors, setErrors] = useState<{ price: string | null; turnaround_days: string | null }>({
    price: null,
    turnaround_days: null,
  });

  // 4. Get handleProvideQuote function from useCROInterface hook
  const { handleProvideQuote } = useCROInterface();

  // 5. Get showToast function from useToast hook
  const { showToast } = useToast();

  /**
   * 6. Define handleChange function to update form state
   * @param {React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>} e - The event object
   */
  const handleChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setFormData({
      ...formData,
      [name]: value,
    });
    setErrors({
      ...errors,
      [name]: null, // Clear error when the user types
    });
  };

  /**
   * 7. Define validateForm function to check form validity
   * @returns {boolean} True if the form is valid, false otherwise
   */
  const validateForm = (): boolean => {
    let newErrors = { price: null, turnaround_days: null };
    let isValid = true;

    if (!formData.price || isNaN(Number(formData.price))) {
      newErrors.price = 'Price is required and must be a number';
      isValid = false;
    }

    if (!formData.turnaround_days || isNaN(Number(formData.turnaround_days))) {
      newErrors.turnaround_days = 'Turnaround days is required and must be a number';
      isValid = false;
    }

    setErrors(newErrors);
    return isValid;
  };

  /**
   * 8. Define handleSubmit function to submit the quote
   * @param {React.FormEvent} e - The form event
   */
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();

    // 9. Validate form data
    if (!validateForm()) {
      return;
    }

    // 10. Set isSubmitting to true
    setIsSubmitting(true);

    try {
      // 11. Call handleProvideQuote with submission ID and form data
      if (submission?.id) {
        await handleProvideQuote(submission.id, formData);

        // 12. On success, show success toast and call onSubmitted callback
        showToast({
          type: 'success',
          message: 'Quote submitted successfully!',
        });
        onSubmitted();
      } else {
        showToast({
          type: 'error',
          message: 'Submission ID is missing.',
        });
      }
    } catch (error: any) {
      // 13. On error, show error toast
      showToast({
        type: 'error',
        message: error?.message || 'Failed to submit quote.',
      });
    } finally {
      // 14. Set isSubmitting to false regardless of outcome
      setIsSubmitting(false);
    }
  };

  const [isSubmitting, setIsSubmitting] = useState(false);

  // 15. Render form with Material-UI components
  return (
    <Paper elevation={3} style={{ padding: '20px' }}>
      <Typography variant="h6" gutterBottom>
        Provide Quote
      </Typography>
      <form onSubmit={handleSubmit}>
        <Grid container spacing={2}>
          {/* 16. Include price input with currency symbol */}
          <Grid item xs={12}>
            <Input
              label="Price"
              name="price"
              value={formData.price}
              onChange={handleChange}
              error={errors.price}
              fullWidth
              InputProps={{
                startAdornment: <AttachMoney />,
              }}
            />
          </Grid>

          {/* 17. Include turnaround days input with calendar icon */}
          <Grid item xs={12}>
            <Input
              label="Estimated Turnaround (Days)"
              name="turnaround_days"
              value={formData.turnaround_days}
              onChange={handleChange}
              error={errors.turnaround_days}
              fullWidth
              InputProps={{
                startAdornment: <CalendarToday />,
              }}
            />
          </Grid>

          {/* 18. Include notes textarea for additional information */}
          <Grid item xs={12}>
            <Input
              label="Additional Notes"
              name="notes"
              value={formData.notes}
              onChange={handleChange}
              multiline
              rows={4}
              fullWidth
              InputProps={{
                startAdornment: <Notes />,
              }}
            />
          </Grid>
        </Grid>

        {/* 19. Divider for visual separation */}
        <Divider style={{ margin: '20px 0' }} />

        {/* 20. Box component for button container */}
        <Box display="flex" justifyContent="flex-end">
          {/* 21. Cancel button with onCancel handler */}
          <Button onClick={onCancel} style={{ marginRight: '10px' }}>
            Cancel
          </Button>

          {/* 22. Submit button with loading state and type='submit' */}
          <Button type="submit" loading={isSubmitting}>
            Submit Quote
          </Button>
        </Box>
      </form>
    </Paper>
  );
};

export default QuoteForm;