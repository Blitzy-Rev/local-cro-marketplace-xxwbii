# Molecular Data Management and CRO Integration Platform - Frontend

Frontend component of the Molecular Data Management and CRO Integration Platform, providing an intuitive user interface for molecular data management, organization, and experimental workflow submission. This React-based application is designed for local deployment with no external dependencies, enabling small to mid-cap pharmaceutical companies to streamline their molecular data management and CRO interactions.

## Features

- Interactive molecular data visualization and management
- Drag-and-drop interface for molecule organization and library creation
- Real-time filtering and sorting of molecular data based on properties
- CSV upload with flexible column mapping
- Experiment creation and submission to CROs
- Result visualization and analysis
- Role-specific interfaces for Pharma users, CRO users, and Administrators
- Responsive design for various screen sizes

## Technology Stack

- **Framework:** React 18.2+
- **Language:** TypeScript 4.9+
- **State Management:** Redux Toolkit 1.9+
- **UI Components:** Material-UI 5.13+
- **Data Fetching:** React Query 4.28+
- **Data Visualization:** D3.js 7.8+, Chart.js 4.3+
- **Form Handling:** Formik 2.2.9, Yup 1.1.1
- **Drag and Drop:** React DnD 16.0.1
- **File Upload:** React Dropzone 14.2.3
- **Routing:** React Router 6.11.1
- **Build Tool:** Vite 4.3.5
- **Testing:** Jest 29.5.0, React Testing Library 14.0.0

## Project Structure

```
src/web/
├── src/                  # Source code
│   ├── api/              # API client functions
│   ├── components/       # Shared UI components
│   │   ├── common/       # Generic UI components
│   │   ├── molecular/    # Molecule-specific components
│   │   └── data-visualization/ # Charts and visualizations
│   ├── features/         # Feature modules
│   │   ├── auth/         # Authentication
│   │   ├── dashboard/    # Dashboard
│   │   ├── molecules/    # Molecule management
│   │   ├── libraries/    # Library management
│   │   ├── experiments/  # Experiment management
│   │   ├── submissions/  # CRO submissions
│   │   ├── results/      # Experimental results
│   │   ├── cro-interface/# CRO user interface
│   │   ├── admin/        # Admin interface
│   │   └── communications/ # User communications
│   ├── hooks/            # Custom React hooks
│   ├── layouts/          # Page layouts
│   ├── store/            # Redux store configuration
│   ├── theme/            # UI theme configuration
│   ├── types/            # TypeScript type definitions
│   ├── utils/            # Utility functions
│   ├── __tests__/        # Test setup and mocks
│   ├── App.tsx           # Main application component
│   ├── index.tsx         # Application entry point
│   └── routes.tsx        # Routing configuration
├── public/               # Static assets
├── docker/               # Docker configuration
├── .env.development      # Development environment variables
├── .env.production       # Production environment variables
├── package.json          # Dependencies and scripts
├── tsconfig.json         # TypeScript configuration
├── vite.config.ts        # Vite configuration
└── README.md             # This file
```

## Setup and Installation

### Prerequisites

- Node.js 16.0+
- npm 8.0+
- Git

### Local Development Setup

1. Clone the repository:
   ```bash
   git clone <repository-url>
   cd <repository-directory>
   ```

2. Navigate to the frontend directory:
   ```bash
   cd src/web
   ```

3. Install dependencies:
   ```bash
   npm install
   ```

4. Create a `.env` file based on `.env.development`:
   ```bash
   cp .env.development .env
   # Edit .env with your configuration
   ```

5. Start the development server:
   ```bash
   npm run dev
   ```

6. Access the application at http://localhost:5173

### Running with Docker

To run the complete application with Docker Compose:

```bash
cd infrastructure
docker-compose up -d
```

The frontend will be available at http://localhost

## Available Scripts

The following npm scripts are available:

- `npm run dev` - Start the development server
- `npm run build` - Build the production-ready application
- `npm run preview` - Preview the production build locally
- `npm run test` - Run tests
- `npm run test:watch` - Run tests in watch mode
- `npm run test:coverage` - Run tests with coverage report
- `npm run lint` - Lint the code
- `npm run lint:fix` - Lint and fix the code
- `npm run format` - Format the code with Prettier
- `npm run typecheck` - Check TypeScript types

## Key Features Implementation

### Molecule Visualization

Molecule structures are rendered using the `MoleculeViewer` component, which converts SMILES strings to visual representations:

```tsx
// src/components/molecular/MoleculeViewer.tsx
import { useEffect, useRef } from 'react';
import { renderMolecule } from '../../utils/molecularUtils';

interface MoleculeViewerProps {
  smiles: string;
  width?: number;
  height?: number;
}

export const MoleculeViewer: React.FC<MoleculeViewerProps> = ({
  smiles,
  width = 200,
  height = 200,
}) => {
  const containerRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (containerRef.current && smiles) {
      renderMolecule(smiles, containerRef.current, { width, height });
    }
  }, [smiles, width, height]);

  return <div ref={containerRef} data-testid="molecule-viewer" />;
};
```

### Drag and Drop Library Management

The drag-and-drop functionality for organizing molecules into libraries is implemented using React DnD:

```tsx
// src/features/libraries/components/DragDropLibrary.tsx
import { useDrop } from 'react-dnd';
import { MoleculeDragItem } from './MoleculeDragItem';
import { useLibraries } from '../hooks/useLibraries';

export const DragDropLibrary = ({ libraryId }) => {
  const { addMoleculeToLibrary } = useLibraries();
  
  const [{ isOver }, drop] = useDrop(() => ({
    accept: 'molecule',
    drop: (item: { id: number }) => {
      addMoleculeToLibrary(libraryId, item.id);
    },
    collect: (monitor) => ({
      isOver: !!monitor.isOver(),
    }),
  }));

  return (
    <div 
      ref={drop} 
      className={`library-dropzone ${isOver ? 'is-over' : ''}`}
    >
      {/* Library content */}
    </div>
  );
};
```

### Real-time Filtering

Molecule filtering is implemented with debounced inputs for performance:

```tsx
// src/features/molecules/components/MoleculeFilter.tsx
import { useState } from 'react';
import { useDebounce } from '../../../hooks/useDebounce';

export const MoleculeFilter = ({ onFilterChange }) => {
  const [filters, setFilters] = useState({});
  const debouncedFilters = useDebounce(filters, 300);

  // Effect to apply filters
  useEffect(() => {
    onFilterChange(debouncedFilters);
  }, [debouncedFilters, onFilterChange]);

  // Filter handling logic
  // ...
};
```

## State Management

The application uses Redux Toolkit for global state management with the following slices:

- `auth` - Authentication state (current user, token, permissions)
- `molecules` - Molecular data and filtering state
- `libraries` - Library management state
- `experiments` - Experiment configuration and tracking
- `submissions` - CRO submission management
- `results` - Experimental result data
- `ui` - UI state (active view, sidebar, modals, notifications)

Example Redux slice:

```typescript
// src/store/molecules/moleculesSlice.ts
import { createSlice, createAsyncThunk } from '@reduxjs/toolkit';
import { getMolecules } from '../../api/molecules';

export const fetchMolecules = createAsyncThunk(
  'molecules/fetchMolecules',
  async (filters, { rejectWithValue }) => {
    try {
      const response = await getMolecules(filters);
      return response;
    } catch (error) {
      return rejectWithValue(error.response.data);
    }
  }
);

const moleculesSlice = createSlice({
  name: 'molecules',
  initialState: {
    molecules: [],
    total: 0,
    loading: false,
    error: null,
    filters: {},
    selectedMolecules: [],
  },
  reducers: {
    // Reducers
  },
  extraReducers: (builder) => {
    // Async thunk handling
  },
});

export const { actions, reducer } = moleculesSlice;
```

For data fetching, the application uses React Query to handle caching, loading states, and error handling.

## API Integration

The frontend communicates with the backend API using Axios. API client functions are organized by domain:

```typescript
// src/api/client.ts
import axios from 'axios';

export const client = axios.create({
  baseURL: import.meta.env.VITE_API_URL || '/api',
  headers: {
    'Content-Type': 'application/json',
  },
});

// Add request interceptor for authentication
client.interceptors.request.use((config) => {
  const token = localStorage.getItem('token');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  return config;
});

// Add response interceptor for error handling
client.interceptors.response.use(
  (response) => response,
  (error) => {
    // Handle authentication errors
    if (error.response && error.response.status === 401) {
      // Redirect to login or refresh token
    }
    return Promise.reject(error);
  }
);
```

Example API module:

```typescript
// src/api/molecules.ts
import { client } from './client';
import { Molecule, MoleculeFilter } from '../types/molecule';

export const getMolecules = async (filters: MoleculeFilter) => {
  const response = await client.get<{ molecules: Molecule[], total: number }>('/v1/molecules', { params: filters });
  return response.data;
};

export const getMoleculeById = async (id: number) => {
  const response = await client.get<Molecule>(`/v1/molecules/${id}`);
  return response.data;
};

export const createMolecule = async (molecule: Partial<Molecule>) => {
  const response = await client.post<Molecule>('/v1/molecules', molecule);
  return response.data;
};

// Additional API functions
```

## Testing

The frontend uses Jest and React Testing Library for testing. Tests are organized alongside the components they test.

```typescript
// src/features/molecules/__tests__/MoleculeTable.test.tsx
import { render, screen } from '@testing-library/react';
import { MoleculeTable } from '../components/MoleculeTable';
import { mockMolecules } from '../../../__tests__/mocks/data';

describe('MoleculeTable', () => {
  it('renders molecule table with data', () => {
    render(<MoleculeTable molecules={mockMolecules} loading={false} />);
    expect(screen.getByText(mockMolecules[0].smiles)).toBeInTheDocument();
  });

  it('shows loading state', () => {
    render(<MoleculeTable molecules={[]} loading={true} />);
    expect(screen.getByTestId('loading-indicator')).toBeInTheDocument();
  });

  it('shows empty state when no molecules', () => {
    render(<MoleculeTable molecules={[]} loading={false} />);
    expect(screen.getByText('No molecules found')).toBeInTheDocument();
  });
});
```

API requests are mocked using Mock Service Worker (MSW):

```typescript
// src/__tests__/mocks/handlers.ts
import { rest } from 'msw';
import { mockMolecules } from './data';

export const handlers = [
  rest.get('/api/v1/molecules', (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json({
        molecules: mockMolecules,
        total: mockMolecules.length,
      })
    );
  }),
  // Additional mock handlers
];
```

## Styling and Theming

The application uses Material-UI for styling with a custom theme:

```typescript
// src/theme/index.ts
import { createTheme } from '@mui/material/styles';
import { palette } from './palette';
import { typography } from './typography';
import { shadows } from './shadows';
import { overrides } from './overrides';

export const theme = createTheme({
  palette,
  typography,
  shadows,
  components: overrides,
});
```

The theme includes a responsive design system:

```typescript
// src/theme/responsive.ts
export const breakpoints = {
  values: {
    xs: 0,
    sm: 600,
    md: 960,
    lg: 1280,
    xl: 1920,
  },
};

export const getResponsiveValue = (value, breakpoint) => {
  // Responsive value logic
};
```

## Deployment

The frontend can be deployed in several ways:

### Production Build

```bash
npm run build
```

This creates a production-ready build in the `dist` directory that can be served by any static file server.

### Docker Deployment

The frontend is containerized using Docker:

```dockerfile
# src/web/docker/Dockerfile
FROM node:18-alpine as build

WORKDIR /app
COPY package*.json ./
RUN npm ci
COPY . .
RUN npm run build

FROM nginx:1.24-alpine
COPY --from=build /app/dist /usr/share/nginx/html
COPY docker/nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
```

The complete application can be deployed using Docker Compose from the `infrastructure` directory:

```bash
cd infrastructure
docker-compose up -d
```

## Contributing

We welcome contributions to the frontend! Please follow these guidelines:

1. Fork the repository and create a feature branch
2. Install dependencies: `npm install`
3. Make your changes following our coding standards
4. Write tests for your changes
5. Run tests: `npm test`
6. Run linting: `npm run lint`
7. Submit a pull request

### Coding Standards

- Use TypeScript for all new code
- Follow the existing code style (enforced by ESLint and Prettier)
- Write comprehensive tests for new features
- Document complex logic with comments
- Use functional components with hooks
- Follow the feature-based organization structure

## License

This project is licensed under the terms specified in the LICENSE file.