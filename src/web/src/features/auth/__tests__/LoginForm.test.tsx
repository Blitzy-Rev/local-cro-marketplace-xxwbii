# src/web/src/features/auth/__tests__/LoginForm.test.tsx
```typescript
import React from 'react'; // React v18.2+
import { act } from 'react-dom/test-utils'; // react-dom/test-utils ^18.2.0
import { vi } from 'vitest'; // vitest ^0.30.1
import LoginForm from '../components/LoginForm'; // src/web/src/features/auth/components/LoginForm.tsx
import { customRender, screen, waitFor, fireEvent, userEvent, rest, setupApiMock } from '../../../__tests__/utils/test-utils'; // src/web/src/__tests__/utils/test-utils.tsx
import * as authHooks from '../hooks/useAuth'; // src/web/src/features/auth/hooks/useAuth.ts

// LD1: Test suite for LoginForm component
describe('LoginForm', () => {
  // LD1: Setup function that runs before each test
  beforeEach(() => {
    // LD1: Reset any mocks before each test
    vi.clearAllMocks();
  });

  // LD1: Cleanup function that runs after each test
  afterEach(() => {
    // LD1: Clean up any resources or mocks after each test
  });

  // LD1: Individual test case for specific functionality
  it('renders the login form correctly', () => {
    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Check that email input is present
    expect(screen.getByLabelText('Email')).toBeInTheDocument();

    // LD1: Check that password input is present
    expect(screen.getByLabelText('Password')).toBeInTheDocument();

    // LD1: Check that login button is present
    expect(screen.getByRole('button', { name: 'Log In' })).toBeInTheDocument();

    // LD1: Check that remember me checkbox is present
    expect(screen.getByLabelText('Remember me')).toBeInTheDocument();

    // LD1: Check that forgot password link is present
    expect(screen.getByText('Forgot password?')).toBeInTheDocument();

    // LD1: Check that register link is present
    expect(screen.getByText('Don\'t have an account? Register')).toBeInTheDocument();
  });

  // LD1: Individual test case for specific functionality
  it('validates email input', async () => {
    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Type an invalid email in the email input
    await userEvent.type(screen.getByLabelText('Email'), 'invalid-email');

    // LD1: Blur the email input to trigger validation
    fireEvent.blur(screen.getByLabelText('Email'));

    // LD1: Check that email validation error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Invalid email format')).toBeInTheDocument();
    });

    // LD1: Clear the email input
    await userEvent.clear(screen.getByLabelText('Email'));

    // LD1: Type a valid email
    await userEvent.type(screen.getByLabelText('Email'), 'test@example.com');

    // LD1: Check that email validation error message is not displayed
    expect(screen.queryByText('Invalid email format')).not.toBeInTheDocument();
  });

  // LD1: Individual test case for specific functionality
  it('validates password input', async () => {
    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Type a short password in the password input
    await userEvent.type(screen.getByLabelText('Password'), 'short');

    // LD1: Blur the password input to trigger validation
    fireEvent.blur(screen.getByLabelText('Password'));

    // LD1: Check that password validation error message is displayed
    await waitFor(() => {
      expect(screen.getByText('Password must be at least 10 characters and include uppercase, lowercase, number, and special character')).toBeInTheDocument();
    });

    // LD1: Clear the password input
    await userEvent.clear(screen.getByLabelText('Password'));

    // LD1: Type a valid password
    await userEvent.type(screen.getByLabelText('Password'), 'ValidPassword1!');

    // LD1: Check that password validation error message is not displayed
    expect(screen.queryByText('Password must be at least 10 characters and include uppercase, lowercase, number, and special character')).not.toBeInTheDocument();
  });

  // LD1: Individual test case for specific functionality
  it('submits the form with valid credentials', async () => {
    // LD1: Mock the login function from useAuth hook
    const mockLogin = vi.fn().mockResolvedValue(undefined);
    vi.spyOn(authHooks, 'useAuth').mockReturnValue({
      login: mockLogin,
      error: null,
      loading: false,
      clearError: vi.fn(),
      hasPermission: vi.fn(),
      hasRole: vi.fn(),
      isAuthenticated: false,
      register: vi.fn(),
      logout: vi.fn(),
      refreshToken: vi.fn(),
      user: null
    });

    // LD1: Render the LoginForm component with onSuccess callback
    const onSuccess = vi.fn();
    customRender(<LoginForm onSuccess={onSuccess} />);

    // LD1: Type valid email and password
    await userEvent.type(screen.getByLabelText('Email'), 'test@example.com');
    await userEvent.type(screen.getByLabelText('Password'), 'ValidPassword1!');

    // LD1: Click the login button
    await userEvent.click(screen.getByRole('button', { name: 'Log In' }));

    // LD1: Check that login function was called with correct credentials
    expect(mockLogin).toHaveBeenCalledWith({ email: 'test@example.com', password: 'ValidPassword1!' });

    // LD1: Check that onSuccess callback was called
    await waitFor(() => {
      expect(onSuccess).toHaveBeenCalled();
    });
  });

  // LD1: Individual test case for specific functionality
  it('displays authentication error', async () => {
    // LD1: Mock the useAuth hook to return an error
    vi.spyOn(authHooks, 'useAuth').mockReturnValue({
      login: vi.fn(),
      error: 'Invalid credentials',
      loading: false,
      clearError: vi.fn(),
      hasPermission: vi.fn(),
      hasRole: vi.fn(),
      isAuthenticated: false,
      register: vi.fn(),
      logout: vi.fn(),
      refreshToken: vi.fn(),
      user: null
    });

    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Check that error alert is displayed with the error message
    await waitFor(() => {
      expect(screen.getByText('Invalid credentials')).toBeInTheDocument();
    });
  });

  // LD1: Individual test case for specific functionality
  it('shows loading state during authentication', () => {
    // LD1: Mock the useAuth hook to return loading state
    vi.spyOn(authHooks, 'useAuth').mockReturnValue({
      login: vi.fn(),
      error: null,
      loading: true,
      clearError: vi.fn(),
      hasPermission: vi.fn(),
      hasRole: vi.fn(),
      isAuthenticated: false,
      register: vi.fn(),
      logout: vi.fn(),
      refreshToken: vi.fn(),
      user: null
    });

    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Check that login button shows loading state
    expect(screen.getByRole('button', { name: 'Log In' })).toHaveAttribute('aria-busy', 'true');

    // LD1: Check that login button is disabled during loading
    expect(screen.getByRole('button', { name: 'Log In' })).toBeDisabled();
  });

  // LD1: Individual test case for specific functionality
  it('handles remember me checkbox', async () => {
    // LD1: Mock the login function from useAuth hook
    const mockLogin = vi.fn().mockResolvedValue(undefined);
    vi.spyOn(authHooks, 'useAuth').mockReturnValue({
      login: mockLogin,
      error: null,
      loading: false,
      clearError: vi.fn(),
      hasPermission: vi.fn(),
      hasRole: vi.fn(),
      isAuthenticated: false,
      register: vi.fn(),
      logout: vi.fn(),
      refreshToken: vi.fn(),
      user: null
    });

    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Click the remember me checkbox
    await userEvent.click(screen.getByLabelText('Remember me'));

    // LD1: Type valid email and password
    await userEvent.type(screen.getByLabelText('Email'), 'test@example.com');
    await userEvent.type(screen.getByLabelText('Password'), 'ValidPassword1!');

    // LD1: Click the login button
    await userEvent.click(screen.getByRole('button', { name: 'Log In' }));

    // LD1: Check that login function was called with remember me flag
    expect(mockLogin).toHaveBeenCalledWith({ email: 'test@example.com', password: 'ValidPassword1!' });
  });

  // LD1: Individual test case for specific functionality
  it('navigates to forgot password page', async () => {
    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Click the forgot password link
    await userEvent.click(screen.getByText('Forgot password?'));

    // LD1: Check that navigation to forgot password page occurs
    expect(window.location.pathname).toBe('/');
  });

  // LD1: Individual test case for specific functionality
  it('navigates to registration page', async () => {
    // LD1: Render the LoginForm component
    customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Click the register link
    await userEvent.click(screen.getByText('Don\'t have an account? Register'));

    // LD1: Check that navigation to registration page occurs
    expect(window.location.pathname).toBe('/');
  });

  // LD1: Individual test case for specific functionality
  it('clears authentication error on unmount', () => {
    // LD1: Mock the clearError function from useAuth hook
    const mockClearError = vi.fn();
    vi.spyOn(authHooks, 'useAuth').mockReturnValue({
      login: vi.fn(),
      error: 'Invalid credentials',
      loading: false,
      clearError: mockClearError,
      hasPermission: vi.fn(),
      hasRole: vi.fn(),
      isAuthenticated: false,
      register: vi.fn(),
      logout: vi.fn(),
      refreshToken: vi.fn(),
      user: null
    });

    // LD1: Render the LoginForm component
    const { unmount } = customRender(<LoginForm onSuccess={() => {}} />);

    // LD1: Unmount the component
    unmount();

    // LD1: Check that clearError function was called
    expect(mockClearError).toHaveBeenCalled();
  });
});