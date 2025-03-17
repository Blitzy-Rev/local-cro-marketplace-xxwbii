import React from 'react'; // React v18.2+
import { createBrowserRouter, Navigate, RouteObject } from 'react-router-dom'; // react-router-dom v6.11.1

import MainLayout from './layouts/MainLayout';
import AuthLayout from './layouts/AuthLayout';
import CROLayout from './layouts/CROLayout';
import AdminLayout from './layouts/AdminLayout';
import LoginPage from './features/auth/pages/LoginPage';
import RegisterPage from './features/auth/pages/RegisterPage';
import ForgotPasswordPage from './features/auth/pages/ForgotPasswordPage';
import ResetPasswordPage from './features/auth/pages/ResetPasswordPage';
import DashboardPage from './features/dashboard/pages/DashboardPage';
import CRODashboardPage from './features/dashboard/pages/CRODashboardPage';
import AdminDashboardPage from './features/dashboard/pages/AdminDashboardPage';
import MoleculesListPage from './features/molecules/pages/MoleculesListPage';
import MoleculeDetailPage from './features/molecules/pages/MoleculeDetailPage';
import CSVUploadPage from './features/csv-upload/pages/CSVUploadPage';
import LibrariesPage from './features/libraries/pages/LibrariesPage';
import LibraryDetailPage from './features/libraries/pages/LibraryDetailPage';
import CreateLibraryPage from './features/libraries/pages/CreateLibraryPage';
import ExperimentsPage from './features/experiments/pages/ExperimentsPage';
import ExperimentDetailPage from './features/experiments/pages/ExperimentDetailPage';
import CreateExperimentPage from './features/experiments/pages/CreateExperimentPage';
import SubmissionsPage from './features/submissions/pages/SubmissionsPage';
import SubmissionDetailPage from './features/submissions/pages/SubmissionDetailPage';
import CreateSubmissionPage from './features/submissions/pages/CreateSubmissionPage';
import ResultsPage from './features/results/pages/ResultsPage';
import ResultDetailPage from './features/results/pages/ResultDetailPage';
import CROSubmissionsPage from './features/cro-interface/pages/CROSubmissionsPage';
import SubmissionReviewPage from './features/cro-interface/pages/SubmissionReviewPage';
import ResultUploadPage from './features/cro-interface/pages/ResultUploadPage';
import UserManagementPage from './features/admin/pages/UserManagementPage';
import SystemMonitoringPage from './features/admin/pages/SystemMonitoringPage';
import CommunicationsPage from './features/communications/pages/CommunicationsPage';

// LD1: Define the route configuration object with nested routes
const routes: RouteObject[] = [
  // LD1: Root path redirects to login page
  {
    path: '/',
    element: <Navigate to="/auth/login" />,
  },
  // LD1: Authentication layout for unauthenticated pages
  {
    path: '/auth',
    element: <AuthLayout />,
    children: [
      // LD1: Login page
      {
        path: 'login',
        element: <LoginPage />,
        // LD2: Add a description for the login route
        meta: {
          description: 'Login page for user authentication',
        },
      },
      // LD1: Registration page
      {
        path: 'register',
        element: <RegisterPage />,
        // LD2: Add a description for the register route
        meta: {
          description: 'Registration page for new user signup',
        },
      },
      // LD1: Forgot password page
      {
        path: 'forgot-password',
        element: <ForgotPasswordPage />,
        // LD2: Add a description for the forgot password route
        meta: {
          description: 'Page for initiating password recovery',
        },
      },
      // LD1: Reset password page
      {
        path: 'reset-password',
        element: <ResetPasswordPage />,
        // LD2: Add a description for the reset password route
        meta: {
          description: 'Page for resetting password with token',
        },
      },
    ],
  },
  // LD1: Main layout for authenticated Pharma users
  {
    path: '/app',
    element: <MainLayout />,
    children: [
      // LD1: Default redirect to dashboard
      {
        path: '',
        element: <Navigate to="/app/dashboard" />,
        // LD2: Add a description for the default redirect route
        meta: {
          description: 'Default redirect to dashboard',
        },
      },
      // LD1: Pharma user dashboard
      {
        path: 'dashboard',
        element: <DashboardPage />,
        // LD2: Add a description for the dashboard route
        meta: {
          description: 'Pharma user dashboard',
        },
      },
      // LD1: Molecules routes
      {
        path: 'molecules',
        children: [
          // LD1: Molecule list page
          {
            path: '',
            element: <MoleculesListPage />,
            // LD2: Add a description for the molecule list route
            meta: {
              description: 'Molecule list page',
            },
          },
          // LD1: Molecule detail page with dynamic ID parameter
          {
            path: ':id',
            element: <MoleculeDetailPage />,
            // LD2: Add a description for the molecule detail route
            meta: {
              description: 'Molecule detail page with dynamic ID parameter',
            },
          },
          // LD1: CSV upload page
          {
            path: 'upload',
            element: <CSVUploadPage />,
            // LD2: Add a description for the CSV upload route
            meta: {
              description: 'CSV upload page',
            },
          },
        ],
      },
      // LD1: Libraries routes
      {
        path: 'libraries',
        children: [
          // LD1: Libraries list page
          {
            path: '',
            element: <LibrariesPage />,
            // LD2: Add a description for the libraries list route
            meta: {
              description: 'Libraries list page',
            },
          },
          // LD1: Library detail page with dynamic ID parameter
          {
            path: ':id',
            element: <LibraryDetailPage />,
            // LD2: Add a description for the library detail route
            meta: {
              description: 'Library detail page with dynamic ID parameter',
            },
          },
          // LD1: Create library page
          {
            path: 'create',
            element: <CreateLibraryPage />,
            // LD2: Add a description for the create library route
            meta: {
              description: 'Create library page',
            },
          },
        ],
      },
      // LD1: Experiments routes
      {
        path: 'experiments',
        children: [
          // LD1: Experiments list page
          {
            path: '',
            element: <ExperimentsPage />,
            // LD2: Add a description for the experiments list route
            meta: {
              description: 'Experiments list page',
            },
          },
          // LD1: Experiment detail page with dynamic ID parameter
          {
            path: ':id',
            element: <ExperimentDetailPage />,
            // LD2: Add a description for the experiment detail route
            meta: {
              description: 'Experiment detail page with dynamic ID parameter',
            },
          },
          // LD1: Create experiment page
          {
            path: 'create',
            element: <CreateExperimentPage />,
            // LD2: Add a description for the create experiment route
            meta: {
              description: 'Create experiment page',
            },
          },
        ],
      },
      // LD1: Submissions routes
      {
        path: 'submissions',
        children: [
          // LD1: Submissions list page
          {
            path: '',
            element: <SubmissionsPage />,
            // LD2: Add a description for the submissions list route
            meta: {
              description: 'Submissions list page',
            },
          },
          // LD1: Submission detail page with dynamic ID parameter
          {
            path: ':id',
            element: <SubmissionDetailPage />,
            // LD2: Add a description for the submission detail route
            meta: {
              description: 'Submission detail page with dynamic ID parameter',
            },
          },
          // LD1: Create submission page
          {
            path: 'create',
            element: <CreateSubmissionPage />,
            // LD2: Add a description for the create submission route
            meta: {
              description: 'Create submission page',
            },
          },
        ],
      },
      // LD1: Results routes
      {
        path: 'results',
        children: [
          // LD1: Results list page
          {
            path: '',
            element: <ResultsPage />,
            // LD2: Add a description for the results list route
            meta: {
              description: 'Results list page',
            },
          },
          // LD1: Result detail page with dynamic ID parameter
          {
            path: ':id',
            element: <ResultDetailPage />,
            // LD2: Add a description for the result detail route
            meta: {
              description: 'Result detail page with dynamic ID parameter',
            },
          },
        ],
      },
      // LD1: Communications page
      {
        path: 'communications',
        element: <CommunicationsPage />,
        // LD2: Add a description for the communications route
        meta: {
          description: 'Communications page',
        },
      },
    ],
  },
  // LD1: CRO layout for authenticated CRO users
  {
    path: '/cro',
    element: <CROLayout />,
    children: [
      // LD1: Default redirect to CRO dashboard
      {
        path: '',
        element: <Navigate to="/cro/dashboard" />,
        // LD2: Add a description for the default redirect route
        meta: {
          description: 'Default redirect to CRO dashboard',
        },
      },
      // LD1: CRO user dashboard
      {
        path: 'dashboard',
        element: <CRODashboardPage />,
        // LD2: Add a description for the CRO dashboard route
        meta: {
          description: 'CRO user dashboard',
        },
      },
      // LD1: CRO submissions routes
      {
        path: 'submissions',
        children: [
          // LD1: CRO submissions list page
          {
            path: '',
            element: <CROSubmissionsPage />,
            // LD2: Add a description for the CRO submissions list route
            meta: {
              description: 'CRO submissions list page',
            },
          },
          // LD1: Submission review page with dynamic ID parameter
          {
            path: ':id',
            element: <SubmissionReviewPage />,
            // LD2: Add a description for the submission review route
            meta: {
              description: 'Submission review page with dynamic ID parameter',
            },
          },
        ],
      },
      // LD1: Result upload page with dynamic submission ID parameter
      {
        path: 'results/upload/:submissionId',
        element: <ResultUploadPage />,
        // LD2: Add a description for the result upload route
        meta: {
          description: 'Result upload page with dynamic submission ID parameter',
        },
      },
      // LD1: Communications page for CRO users
      {
        path: 'communications',
        element: <CommunicationsPage />,
        // LD2: Add a description for the communications route
        meta: {
          description: 'Communications page for CRO users',
        },
      },
    ],
  },
  // LD1: Admin layout for authenticated Admin users
  {
    path: '/admin',
    element: <AdminLayout />,
    children: [
      // LD1: Default redirect to Admin dashboard
      {
        path: '',
        element: <Navigate to="/admin/dashboard" />,
        // LD2: Add a description for the default redirect route
        meta: {
          description: 'Default redirect to Admin dashboard',
        },
      },
      // LD1: Admin dashboard
      {
        path: 'dashboard',
        element: <AdminDashboardPage />,
        // LD2: Add a description for the admin dashboard route
        meta: {
          description: 'Admin dashboard',
        },
      },
      // LD1: User management page
      {
        path: 'users',
        element: <UserManagementPage />,
        // LD2: Add a description for the user management route
        meta: {
          description: 'User management page',
        },
      },
      // LD1: System monitoring page
      {
        path: 'system',
        element: <SystemMonitoringPage />,
        // LD2: Add a description for the system monitoring route
        meta: {
          description: 'System monitoring page',
        },
      },
    ],
  },
  // LD1: Catch-all route for 404 handling, redirects to root
  {
    path: '*',
    element: <Navigate to="/" />,
    // LD2: Add a description for the catch-all route
    meta: {
      description: 'Catch-all route for 404 handling, redirects to root',
    },
  },
];

// LD1: Create and configure the application router with all routes and layouts
const createRouter = () => {
  // LD1: Create a browser router with the route configuration
  const router = createBrowserRouter(routes);
  return router;
};

// IE3: Export the configured router for use in the App component
const router = createRouter();
export default router;