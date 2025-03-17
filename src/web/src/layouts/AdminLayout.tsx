import React from 'react'; // React v18.2+
import { Navigate, Outlet } from 'react-router-dom'; // react-router-dom v6.4+

import MainLayout from './MainLayout';
import useAuth from '../features/auth/hooks/useAuth';
import { UserRole } from '../types/user';

/**
 * @interface AdminLayoutProps
 * @description Props interface for the AdminLayout component
 */
interface AdminLayoutProps {}

/**
 * @component AdminLayout
 * @description Layout component for Admin-specific pages with role-based access control
 */
const AdminLayout: React.FC<AdminLayoutProps> = () => {
  // LD1: Retrieve authentication state and user role using useAuth hook
  const { user, checkRole } = useAuth();

  // LD1: Check if the user has the Admin role using checkRole function
  const isAdmin = checkRole(UserRole.ADMIN);

  // LD1: If user is not an Admin, redirect to dashboard using Navigate component
  if (!isAdmin) {
    return <Navigate to="/app/dashboard" />;
  }

  // LD1: If user is an Admin, render the MainLayout component
  // LD1: Render Outlet component inside MainLayout to display nested route content
  return (
    <MainLayout>
      <Outlet />
    </MainLayout>
  );
};

// IE3: Export the AdminLayout component for use in routing
export default AdminLayout;