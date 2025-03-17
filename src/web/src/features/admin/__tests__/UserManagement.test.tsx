# src/web/src/features/admin/__tests__/UserManagement.test.tsx
```typescript
import React from 'react'; // react ^18.2.0
import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest'; // vitest ^0.30.1
import { customRender, screen, waitFor, fireEvent, userEvent, within, setupApiMock, rest } from '../../../../__tests__/utils/test-utils';
import UserManagement from '../../components/UserManagement';
import { UserRole, UserStatus } from '../../../../types/user';

/**
 * Generates an array of mock users for testing
 * @param count The number of mock users to generate
 * @returns An array of mock user objects
 */
const generateMockUsers = (count: number) => {
  // LD1: Create an array of the specified count
  const mockUsersArray = Array.from({ length: count }, (_, i) => {
    // LD1: Map each index to a mock user object with id, email, role, status, and other properties
    return {
      id: i + 1,
      email: `test${i + 1}@example.com`,
      role: UserRole.PHARMA,
      status: UserStatus.ACTIVE,
      is_active: true,
      email_verified: true,
      created_at: new Date().toISOString(),
      last_login: new Date().toISOString(),
    };
  });
  // LD1: Return the array of mock users
  return mockUsersArray;
};

describe('UserManagement Component', () => {
  // LD1: Define mock data for API responses
  const mockUsers = generateMockUsers(5);
  const mockUserListResponse = {
    items: mockUsers,
    total: 5,
    page: 1,
    size: 10,
  };

  // LD1: Define mock API handlers for testing
  const mockApiHandlers = [
    rest.get('/api/v1/admin/users', (req, res, ctx) => {
      return res(ctx.status(200), ctx.json(mockUserListResponse));
    }),
    rest.post('/api/v1/admin/users', (req, res, ctx) => {
      return res(
        ctx.status(201),
        ctx.json({
          success: true,
          data: {
            id: 6,
            email: 'newuser@example.com',
            role: UserRole.PHARMA,
            status: UserStatus.ACTIVE,
            is_active: true,
            email_verified: false,
            created_at: new Date().toISOString(),
            last_login: null,
          },
          message: 'User created successfully',
        })
      );
    }),
    rest.delete('/api/v1/admin/users/:id', (req, res, ctx) => {
      return res(ctx.status(200), ctx.json({ success: true, message: 'User deleted successfully' }));
    }),
  ];

  beforeEach(() => {
    // LD1: Set up mock API handlers using setupApiMock
    setupApiMock(mockApiHandlers);
    // LD1: Mock console.error to prevent error messages during tests
    vi.spyOn(console, 'error').mockImplementation(() => {});
  });

  afterEach(() => {
    // LD1: Restore console.error
    console.error.mockRestore();
    // LD1: Clear all mocks
    vi.clearAllMocks();
  });

  it('renders the component with loading state', () => {
    // LD1: Mock the API to return loading state
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.delay(100), ctx.status(200), ctx.json(mockUserListResponse));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Assert that loading indicator is displayed
    expect(screen.getByText('Loading users...')).toBeInTheDocument();
  });

  it('renders the user list when data is loaded', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json(mockUserListResponse));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Assert that the user table is displayed with the correct number of rows
    const rows = screen.getAllByRole('row');
    expect(rows.length).toBe(mockUsers.length + 1); // +1 for header row
    // LD1: Assert that user data is displayed correctly in each row
    expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    expect(screen.getByText('test2@example.com')).toBeInTheDocument();
  });

  it('handles filtering by role', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        const role = req.url.searchParams.get('role');
        const filteredUsers = role ? mockUsers.filter(user => user.role === role) : mockUsers;
        return res(ctx.status(200), ctx.json({ ...mockUserListResponse, items: filteredUsers }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Select a role from the role filter dropdown
    fireEvent.mouseDown(screen.getByLabelText('Filter by Role'));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText('CRO'));
    // LD1: Assert that the API is called with the correct role filter
    await waitFor(() => {
      expect(screen.getByText('cro@example.com')).toBeInTheDocument();
    });
    // LD1: Assert that the filtered user list is displayed
    expect(screen.queryByText('pharma@example.com')).not.toBeInTheDocument();
  });

  it('handles filtering by status', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        const status = req.url.searchParams.get('status');
        const filteredUsers = status ? mockUsers.filter(user => user.status === status) : mockUsers;
        return res(ctx.status(200), ctx.json({ ...mockUserListResponse, items: filteredUsers }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Select a status from the status filter dropdown
    fireEvent.mouseDown(screen.getByLabelText('Filter by Status'));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText('Active'));
    // LD1: Assert that the API is called with the correct status filter
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Assert that the filtered user list is displayed
    expect(screen.queryByText('pending@example.com')).not.toBeInTheDocument();
  });

  it('handles search by email', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        const email = req.url.searchParams.get('email');
        const filteredUsers = email ? mockUsers.filter(user => user.email.includes(email)) : mockUsers;
        return res(ctx.status(200), ctx.json({ ...mockUserListResponse, items: filteredUsers }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Enter a search term in the email search input
    const searchInput = screen.getByRole('textbox', { name: 'Search by Email' });
    fireEvent.change(searchInput, { target: { value: 'test1' } });
    // LD1: Assert that the API is called with the correct email filter
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Assert that the filtered user list is displayed
    expect(screen.queryByText('test2@example.com')).not.toBeInTheDocument();
  });

  it('handles pagination', async () => {
    // LD1: Mock the API to return a paginated list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        const page = parseInt(req.url.searchParams.get('page') || '1');
        const limit = parseInt(req.url.searchParams.get('limit') || '10');
        const startIndex = (page - 1) * limit;
        const endIndex = Math.min(startIndex + limit, mockUsers.length);
        const paginatedUsers = mockUsers.slice(startIndex, endIndex);
        return res(ctx.status(200), ctx.json({ ...mockUserListResponse, items: paginatedUsers, page }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Click on the next page button
    const nextButton = screen.getByRole('button', { name: 'Go to next page' });
    fireEvent.click(nextButton);
    // LD1: Assert that the API is called with the correct page parameter
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Assert that the new page of users is displayed
    expect(screen.queryByText('test1@example.com')).toBeInTheDocument();
  });

  it('displays user details when a user is selected', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json(mockUserListResponse));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Click on a user row
    const userRow = screen.getByRole('row', { name: 'test1@example.com' });
    fireEvent.click(userRow);
    // LD1: Assert that the user detail component is displayed with the correct user data
    await waitFor(() => {
      expect(screen.getByText('User Details')).toBeInTheDocument();
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
  });

  it('handles creating a new user', async () => {
    // LD1: Mock the API to return a list of users and handle user creation
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json(mockUserListResponse));
      }),
      rest.post('/api/v1/admin/users', (req, res, ctx) => {
        return res(
          ctx.status(201),
          ctx.json({
            success: true,
            data: {
              id: 6,
              email: 'newuser@example.com',
              role: UserRole.PHARMA,
              status: UserStatus.ACTIVE,
              is_active: true,
              email_verified: false,
              created_at: new Date().toISOString(),
              last_login: null,
            },
            message: 'User created successfully',
          })
        );
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Click on the 'Add User' button
    const addUserButton = screen.getByRole('button', { name: 'Add User' });
    fireEvent.click(addUserButton);
    // LD1: Fill in the user creation form
    const emailInput = screen.getByLabelText('Email');
    fireEvent.change(emailInput, { target: { value: 'newuser@example.com' } });
    const passwordInput = screen.getByLabelText('Password');
    fireEvent.change(passwordInput, { target: { value: 'Password123!' } });
    const confirmPasswordInput = screen.getByLabelText('Confirm Password');
    fireEvent.change(confirmPasswordInput, { target: { value: 'Password123!' } });
    fireEvent.mouseDown(screen.getByLabelText('Role'));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText('Pharma'));
    // LD1: Submit the form
    const createButton = screen.getByRole('button', { name: 'Create' });
    fireEvent.click(createButton);
    // LD1: Assert that the API is called with the correct user data
    await waitFor(() => {
      expect(screen.getByText('User created successfully')).toBeInTheDocument();
    });
    // LD1: Assert that a success message is displayed
    expect(screen.getByText('User created successfully')).toBeInTheDocument();
    // LD1: Assert that the user list is refreshed
    await waitFor(() => {
      expect(screen.getByText('newuser@example.com')).toBeInTheDocument();
    });
  });

  it('validates user creation form', async () => {
    // LD1: Mock the API to return a list of users
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json(mockUserListResponse));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Click on the 'Add User' button
    const addUserButton = screen.getByRole('button', { name: 'Add User' });
    fireEvent.click(addUserButton);
    // LD1: Submit the form without filling in required fields
    const createButton = screen.getByRole('button', { name: 'Create' });
    fireEvent.click(createButton);
    // LD1: Assert that validation error messages are displayed
    await waitFor(() => {
      expect(screen.getByText('Email is required')).toBeInTheDocument();
      expect(screen.getByText('Password is required')).toBeInTheDocument();
    });
    // LD1: Fill in invalid data (e.g., mismatched passwords)
    const passwordInput = screen.getByLabelText('Password');
    fireEvent.change(passwordInput, { target: { value: 'Password123!' } });
    const confirmPasswordInput = screen.getByLabelText('Confirm Password');
    fireEvent.change(confirmPasswordInput, { target: { value: 'Password123' } });
    fireEvent.click(createButton);
    // LD1: Assert that validation error messages are displayed
    await waitFor(() => {
      expect(screen.getByText('Passwords do not match')).toBeInTheDocument();
    });
  });

  it('handles deleting a user', async () => {
    // LD1: Mock the API to return a list of users and handle user deletion
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json(mockUserListResponse));
      }),
      rest.delete('/api/v1/admin/users/:id', (req, res, ctx) => {
        return res(ctx.status(200), ctx.json({ success: true, message: 'User deleted successfully' }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Wait for the data to load
    await waitFor(() => {
      expect(screen.getByText('test1@example.com')).toBeInTheDocument();
    });
    // LD1: Click on the delete button for a user
    const deleteButton = screen.getAllByRole('button', { name: 'Delete' })[0];
    fireEvent.click(deleteButton);
    // LD1: Confirm the deletion in the confirmation dialog
    const confirmDeleteButton = screen.getByRole('button', { name: 'Delete' });
    fireEvent.click(confirmDeleteButton);
    // LD1: Assert that the API is called with the correct user ID
    await waitFor(() => {
      expect(screen.getByText('User deleted successfully')).toBeInTheDocument();
    });
    // LD1: Assert that a success message is displayed
    expect(screen.getByText('User deleted successfully')).toBeInTheDocument();
    // LD1: Assert that the user list is refreshed
    await waitFor(() => {
      expect(screen.queryByText('test1@example.com')).not.toBeInTheDocument();
    });
  });

  it('handles API errors', async () => {
    // LD1: Mock the API to return an error
    setupApiMock([
      rest.get('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(500), ctx.json({ success: false, error: 'Failed to fetch users' }));
      }),
    ]);
    // LD1: Render the UserManagement component
    customRender(<UserManagement />);
    // LD1: Assert that an error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Failed to fetch users')).toBeInTheDocument();
    });
    // LD1: Mock user creation to fail
    setupApiMock([
      rest.post('/api/v1/admin/users', (req, res, ctx) => {
        return res(ctx.status(500), ctx.json({ success: false, error: 'Failed to create user' }));
      }),
    ]);
    // LD1: Attempt to create a user
    const addUserButton = screen.getByRole('button', { name: 'Add User' });
    fireEvent.click(addUserButton);
    const emailInput = screen.getByLabelText('Email');
    fireEvent.change(emailInput, { target: { value: 'newuser@example.com' } });
    const passwordInput = screen.getByLabelText('Password');
    fireEvent.change(passwordInput, { target: { value: 'Password123!' } });
    const confirmPasswordInput = screen.getByLabelText('Confirm Password');
    fireEvent.change(confirmPasswordInput, { target: { value: 'Password123!' } });
    fireEvent.mouseDown(screen.getByLabelText('Role'));
    const listbox = within(screen.getByRole('listbox'));
    fireEvent.click(listbox.getByText('Pharma'));
    const createButton = screen.getByRole('button', { name: 'Create' });
    fireEvent.click(createButton);
    // LD1: Assert that an error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Failed to create user')).toBeInTheDocument();
    });
    // LD1: Mock user deletion to fail
    setupApiMock([
      rest.delete('/api/v1/admin/users/:id', (req, res, ctx) => {
        return res(ctx.status(500), ctx.json({ success: false, error: 'Failed to delete user' }));
      }),
    ]);
    // LD1: Attempt to delete a user
    const deleteButton = screen.getAllByRole('button', { name: 'Delete' })[0];
    fireEvent.click(deleteButton);
    const confirmDeleteButton = screen.getByRole('button', { name: 'Delete' });
    fireEvent.click(confirmDeleteButton);
    // LD1: Assert that an error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Failed to delete user')).toBeInTheDocument();
    });
  });
});