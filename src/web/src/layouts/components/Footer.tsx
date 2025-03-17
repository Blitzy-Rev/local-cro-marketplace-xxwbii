import React from 'react'; // React version 18.2+
import { Box, Typography, Link, Container } from '@mui/material'; // @mui/material version 5.13+
import { styled } from '@mui/material/styles'; // @mui/material version 5.13+
import theme from '../../theme'; // Theme configuration for consistent styling
import useWindowSize from '../../hooks/useWindowSize'; // Hook for responsive design adaptations

/**
 * Props interface for the Footer component
 */
interface FooterProps {}

// Styled components for the footer
const FooterContainer = styled(Box)(({ theme }) => ({
  backgroundColor: theme.palette.background.paper,
  color: theme.palette.text.secondary,
  padding: theme.spacing(2, 0),
  borderTop: '1px solid',
  borderColor: theme.palette.divider,
  marginTop: 'auto',
}));

const FooterContent = styled(Container)(({ theme }) => ({
  display: 'flex',
  flexDirection: {
    xs: 'column',
    sm: 'row',
  },
  justifyContent: 'space-between',
  alignItems: {
    xs: 'center',
    sm: 'flex-start',
  },
  gap: theme.spacing(2),
}));

const FooterLinks = styled(Box)(({ theme }) => ({
  display: 'flex',
  gap: theme.spacing(2),
  flexWrap: 'wrap',
  justifyContent: {
    xs: 'center',
    sm: 'flex-end',
  },
}));

const FooterLink = styled(Link)(({ theme }) => ({
  color: theme.palette.text.secondary,
  textDecoration: 'none',
  fontSize: '0.875rem',
  '&:hover': {
    color: theme.palette.primary.main,
    textDecoration: 'underline',
  },
  '&:focus': {
    outline: '2px solid',
    outlineColor: theme.palette.primary.main,
    outlineOffset: '2px',
  },
}));

const VersionText = styled(Typography)(({ theme }) => ({
  fontSize: '0.75rem',
  marginTop: theme.spacing(1),
  textAlign: {
    xs: 'center',
    sm: 'left',
  },
}));

/**
 * Footer component with copyright information, version details, and essential links
 * Adapts to different screen sizes with responsive layout changes
 * Ensures accessibility with proper contrast and focus states
 */
const Footer: React.FC<FooterProps> = () => {
  // Get window size information for responsive design
  const { isMobile } = useWindowSize();
  
  // Get current year for copyright notice
  const currentYear = new Date().getFullYear();
  
  return (
    <FooterContainer>
      <FooterContent maxWidth="lg">
        <Box>
          <Typography variant="body2">
            © {currentYear} Molecular Data Management Platform
          </Typography>
          <VersionText variant="caption">
            Version 1.0.0
          </VersionText>
        </Box>
        
        <FooterLinks>
          <FooterLink href="/help">Help</FooterLink>
          <FooterLink href="/privacy">Privacy Policy</FooterLink>
          <FooterLink href="/terms">Terms of Service</FooterLink>
        </FooterLinks>
      </FooterContent>
    </FooterContainer>
  );
};

export default Footer;