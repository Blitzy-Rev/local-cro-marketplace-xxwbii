import React from 'react'; // React v18.2+
import { Navigate, Outlet } from 'react-router-dom'; // react-router-dom v6.4+

import MainLayout from './MainLayout';
import useAuth from '../features/auth/hooks/useAuth';
import { UserRole } from '../types/user';

/**
 * @interface CROLayoutProps
 * @description Props interface for the CROLayout component
 */
interface CROLayoutProps {}

/**
 * @component CROLayout
 * @description Layout component for CRO-specific pages with role-based access control
 */
const CROLayout: React.FC<CROLayoutProps> = () => {
  // Retrieve authentication state and user role using useAuth hook
  const { user, checkRole } = useAuth();

  // Check if the user has the CRO role using checkRole function
  const isCROUser = checkRole(UserRole.CRO);

  // If user is not a CRO, redirect to dashboard using Navigate component
  if (!isCROUser) {
    return <Navigate to="/app/dashboard" />;
  }

  // If user is a CRO, render the MainLayout component
  // Render Outlet component inside MainLayout to display nested route content
  return (
    <MainLayout>
      <Outlet />
    </MainLayout>
  );
};

export default CROLayout;