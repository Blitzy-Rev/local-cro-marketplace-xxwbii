import { renderHook, act } from '@testing-library/react-hooks'; // @testing-library/react-hooks ^8.0.1
import { describe, it, expect, beforeEach, afterEach, jest } from 'jest'; // jest ^29.5.0
import { useAuth } from '../../features/auth/hooks/useAuth';
import { customRender, screen, waitFor, setupApiMock, rest } from '../../__tests__/utils/test-utils';
import { LoginCredentials, RegisterData } from '../../types/auth';
import { UserRole } from '../../types/user';

/**
 * Helper function to render the useAuth hook with all necessary providers
 */
const renderAuthHook = () => {
  // LD1: Use renderHook to render the useAuth hook
  const { result } = renderHook(() => useAuth(), {
    // LD2: Wrap the hook with the necessary providers using customRender's wrapper
    wrapper: customRender,
  });
  // LD3: Return the result of renderHook
  return result;
};

describe('useAuth hook', () => {
  // Before each test, reset the MSW handlers
  beforeEach(() => {
    setupApiMock([]);
  });

  // After each test, clear all mocks
  afterEach(() => {
    jest.clearAllMocks();
  });

  it('should return the initial authentication state', () => {
    // LD1: Render the useAuth hook
    const { result } = renderAuthHook();

    // LD2: Verify that isAuthenticated is false
    expect(result.current.isAuthenticated).toBe(false);

    // LD3: Verify that user is null
    expect(result.current.user).toBeNull();

    // LD4: Verify that loading is false
    expect(result.current.loading).toBe(false);

    // LD5: Verify that error is null
    expect(result.current.error).toBeNull();
  });

  it('should handle login success', async () => {
    // LD1: Set up API mock for successful login response
    setupApiMock([
      rest.post('/api/v1/auth/login', (req, res, ctx) => {
        return res(
          ctx.status(200),
          ctx.json({
            success: true,
            data: {
              user: { id: '1', email: 'test@example.com', role: 'pharma' },
              tokens: { access_token: 'mock-access-token', refresh_token: 'mock-refresh-token' },
            },
          })
        );
      }),
    ]);

    // LD2: Render the useAuth hook
    const { result } = renderAuthHook();

    // LD3: Call the login function with valid credentials
    act(() => {
      result.current.login({ email: 'test@example.com', password: 'password' } as LoginCredentials);
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.isAuthenticated).toBe(true));

    // LD6: Verify that isAuthenticated becomes true
    expect(result.current.isAuthenticated).toBe(true);

    // LD7: Verify that user contains the expected data
    expect(result.current.user).toEqual({ id: '1', email: 'test@example.com', role: 'pharma' });

    // LD8: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD9: Verify that error remains null
    expect(result.current.error).toBeNull();
  });

  it('should handle login failure', async () => {
    // LD1: Set up API mock for login failure response
    setupApiMock([
      rest.post('/api/v1/auth/login', (req, res, ctx) => {
        return res(
          ctx.status(401),
          ctx.json({ success: false, error: 'Invalid credentials' })
        );
      }),
    ]);

    // LD2: Render the useAuth hook
    const { result } = renderAuthHook();

    // LD3: Call the login function with invalid credentials
    act(() => {
      result.current.login({ email: 'test@example.com', password: 'wrongpassword' } as LoginCredentials);
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.loading).toBe(false));

    // LD6: Verify that isAuthenticated remains false
    expect(result.current.isAuthenticated).toBe(false);

    // LD7: Verify that user remains null
    expect(result.current.user).toBeNull();

    // LD8: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD9: Verify that error contains the expected error message
    expect(result.current.error).toBe('Invalid credentials');
  });

  it('should handle registration success', async () => {
    // LD1: Set up API mock for successful registration response
    setupApiMock([
      rest.post('/api/v1/auth/register', (req, res, ctx) => {
        return res(
          ctx.status(201),
          ctx.json({ success: true, message: 'Registration successful' })
        );
      }),
    ]);

    // LD2: Render the useAuth hook
    const { result } = renderAuthHook();

    // LD3: Call the register function with valid registration data
    act(() => {
      result.current.register({ email: 'test@example.com', password: 'password', confirm_password: 'password', role: UserRole.PHARMA } as RegisterData);
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.loading).toBe(false));

    // LD6: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD7: Verify that error remains null
    expect(result.current.error).toBeNull();
  });

  it('should handle registration failure', async () => {
    // LD1: Set up API mock for registration failure response
    setupApiMock([
      rest.post('/api/v1/auth/register', (req, res, ctx) => {
        return res(
          ctx.status(400),
          ctx.json({ success: false, error: 'Registration failed' })
        );
      }),
    ]);

    // LD2: Render the useAuth hook
    const { result } = renderAuthHook();

    // LD3: Call the register function with invalid registration data
    act(() => {
      result.current.register({ email: 'test@example.com', password: 'password', confirm_password: 'password', role: UserRole.PHARMA } as RegisterData);
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.loading).toBe(false));

    // LD6: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD7: Verify that error contains the expected error message
    expect(result.current.error).toBe('Registration failed');
  });

  it('should handle logout', async () => {
    // LD1: Set up API mock for successful logout response
    setupApiMock([
      rest.post('/api/v1/auth/logout', (req, res, ctx) => {
        return res(
          ctx.status(200),
          ctx.json({ success: true, message: 'Logout successful' })
        );
      }),
    ]);

    // LD2: Render the useAuth hook with authenticated state
    const { result } = renderAuthHook();
    act(() => {
      result.current.login({ email: 'test@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(result.current.isAuthenticated).toBe(true));

    // LD3: Call the logout function
    act(() => {
      result.current.logout();
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.loading).toBe(false));

    // LD6: Verify that isAuthenticated becomes false
    expect(result.current.isAuthenticated).toBe(false);

    // LD7: Verify that user becomes null
    expect(result.current.user).toBeNull();

    // LD8: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD9: Verify that error remains null
    expect(result.current.error).toBeNull();
  });

  it('should handle token refresh', async () => {
    // LD1: Set up API mock for successful token refresh response
    setupApiMock([
      rest.post('/api/v1/auth/refresh-token', (req, res, ctx) => {
        return res(
          ctx.status(200),
          ctx.json({
            success: true,
            data: {
              access_token: 'new-mock-access-token',
              refresh_token: 'new-mock-refresh-token',
              user: { id: '1', email: 'test@example.com', role: 'pharma' },
            },
          })
        );
      }),
    ]);

    // LD2: Render the useAuth hook with authenticated state
    const { result } = renderAuthHook();
    act(() => {
      result.current.login({ email: 'test@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(result.current.isAuthenticated).toBe(true));

    // LD3: Call the refreshToken function
    act(() => {
      result.current.refreshToken();
    });

    // LD4: Verify that loading becomes true during the request
    expect(result.current.loading).toBe(true);

    // LD5: Wait for the request to complete
    await waitFor(() => expect(result.current.loading).toBe(false));

    // LD6: Verify that isAuthenticated remains true
    expect(result.current.isAuthenticated).toBe(true);

    // LD7: Verify that user contains the updated data
    expect(result.current.user).toEqual({ id: '1', email: 'test@example.com', role: 'pharma' });

    // LD8: Verify that loading becomes false
    expect(result.current.loading).toBe(false);

    // LD9: Verify that error remains null
    expect(result.current.error).toBeNull();
  });

  it('should clear authentication errors', () => {
    // LD1: Render the useAuth hook with an error state
    const { result } = renderAuthHook();
    act(() => {
      result.current.login({ email: 'test@example.com', password: 'wrongpassword' } as LoginCredentials);
    });
    waitFor(() => expect(result.current.error).not.toBeNull());

    // LD2: Call the clearError function
    act(() => {
      result.current.clearError();
    });

    // LD3: Verify that error becomes null
    expect(result.current.error).toBeNull();
  });

  it('should check permissions correctly', async () => {
    // LD1: Render the useAuth hook with authenticated state for different user roles
    const { result: pharmaResult } = renderAuthHook();
    act(() => {
      pharmaResult.current.login({ email: 'pharma@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(pharmaResult.current.isAuthenticated).toBe(true));

    const { result: croResult } = renderAuthHook();
    act(() => {
      croResult.current.login({ email: 'cro@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(croResult.current.isAuthenticated).toBe(true));

    const { result: adminResult } = renderAuthHook();
    act(() => {
      adminResult.current.login({ email: 'admin@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(adminResult.current.isAuthenticated).toBe(true));

    const { result: unauthenticatedResult } = renderAuthHook();

    // LD2: Test hasPermission function with various permissions
    const pharmaPermission = 'molecules:view';
    const croPermission = 'submissions:quote';
    const adminPermission = 'users:edit';

    // LD3: Verify that admin users have all permissions
    expect(adminResult.current.hasPermission(pharmaPermission)).toBe(true);
    expect(adminResult.current.hasPermission(croPermission)).toBe(true);
    expect(adminResult.current.hasPermission(adminPermission)).toBe(true);

    // LD4: Verify that pharma users have only pharma permissions
    expect(pharmaResult.current.hasPermission(pharmaPermission)).toBe(true);
    expect(pharmaResult.current.hasPermission(croPermission)).toBe(false);
    expect(pharmaResult.current.hasPermission(adminPermission)).toBe(false);

    // LD5: Verify that CRO users have only CRO permissions
    expect(croResult.current.hasPermission(pharmaPermission)).toBe(false);
    expect(croResult.current.hasPermission(croPermission)).toBe(true);
    expect(croResult.current.hasPermission(adminPermission)).toBe(false);

    // LD6: Verify that unauthenticated users have no permissions
    expect(unauthenticatedResult.current.hasPermission(pharmaPermission)).toBe(false);
    expect(unauthenticatedResult.current.hasPermission(croPermission)).toBe(false);
    expect(unauthenticatedResult.current.hasPermission(adminPermission)).toBe(false);
  });

  it('should check roles correctly', async () => {
    // LD1: Render the useAuth hook with authenticated state for different user roles
    const { result: pharmaResult } = renderAuthHook();
    act(() => {
      pharmaResult.current.login({ email: 'pharma@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(pharmaResult.current.isAuthenticated).toBe(true));

    const { result: croResult } = renderAuthHook();
    act(() => {
      croResult.current.login({ email: 'cro@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(croResult.current.isAuthenticated).toBe(true));

    const { result: adminResult } = renderAuthHook();
    act(() => {
      adminResult.current.login({ email: 'admin@example.com', password: 'password' } as LoginCredentials);
    });
    await waitFor(() => expect(adminResult.current.isAuthenticated).toBe(true));

    const { result: unauthenticatedResult } = renderAuthHook();

    // LD2: Test hasRole function with various roles
    // LD3: Verify that users with a specific role return true for that role
    expect(pharmaResult.current.hasRole(UserRole.PHARMA)).toBe(true);
    expect(croResult.current.hasRole(UserRole.CRO)).toBe(true);
    expect(adminResult.current.hasRole(UserRole.ADMIN)).toBe(true);

    // LD4: Verify that users with a different role return false
    expect(pharmaResult.current.hasRole(UserRole.CRO)).toBe(false);
    expect(croResult.current.hasRole(UserRole.PHARMA)).toBe(false);
    expect(adminResult.current.hasRole(UserRole.PHARMA)).toBe(false);

    // LD5: Verify that unauthenticated users return false for all roles
    expect(unauthenticatedResult.current.hasRole(UserRole.PHARMA)).toBe(false);
    expect(unauthenticatedResult.current.hasRole(UserRole.CRO)).toBe(false);
    expect(unauthenticatedResult.current.hasRole(UserRole.ADMIN)).toBe(false);
  });
});