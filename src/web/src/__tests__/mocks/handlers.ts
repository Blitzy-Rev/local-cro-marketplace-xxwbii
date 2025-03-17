import { rest, ResponseComposition, RestContext, RestRequest } from 'msw'; // msw ^1.2.1

const API_BASE_URL = '/api/v1';

// Helper functions for creating standardized responses
function createSuccessResponse(data: any, message?: string) {
  return {
    success: true,
    data,
    message: message || 'Operation successful',
  };
}

function createErrorResponse(message: string, statusCode: number = 400, details?: any) {
  return {
    success: false,
    error: message,
    statusCode,
    details,
  };
}

// Mock Data
// Mock users for authentication tests
const mockUsers = [
  { id: '1', email: 'pharma@example.com', role: 'pharma' },
  { id: '2', email: 'cro@example.com', role: 'cro' },
  { id: '3', email: 'admin@example.com', role: 'admin' }
];

// Mock molecules for testing molecule-related functionality
const mockMolecules = [
  {
    id: '1',
    smiles: 'CCO',
    properties: {
      molecular_weight: 46.07,
      logp: -0.14,
      solubility: 3.45,
    },
    flag_status: 'important',
    created_by: '1',
    created_at: '2023-06-01T10:00:00Z'
  },
  {
    id: '2',
    smiles: 'CCCCO',
    properties: {
      molecular_weight: 74.12,
      logp: 0.88,
      solubility: 2.17,
    },
    flag_status: null,
    created_by: '1',
    created_at: '2023-06-02T10:00:00Z'
  },
  {
    id: '3',
    smiles: 'c1ccccc1',
    properties: {
      molecular_weight: 78.11,
      logp: 1.9,
      solubility: 0.12,
    },
    flag_status: null,
    created_by: '1',
    created_at: '2023-06-03T10:00:00Z'
  }
];

// Define request handlers
const handlers = [
  // Authentication endpoints
  rest.post(`${API_BASE_URL}/auth/login`, (req, res, ctx) => {
    const { email, password } = req.body as { email: string; password: string };
    
    // Basic validation
    if (!email || !password) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Email and password are required'))
      );
    }
    
    // Check credentials (for testing purposes)
    const user = mockUsers.find(u => u.email === email);
    if (!user || password !== 'password') { // Simple password check for testing
      return res(
        ctx.status(401),
        ctx.json(createErrorResponse('Invalid credentials', 401))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        user: {
          id: user.id,
          email: user.email,
          role: user.role,
        },
        tokens: {
          access_token: 'mock-access-token',
          refresh_token: 'mock-refresh-token',
          expires_in: 900, // 15 minutes
        }
      }))
    );
  }),
  
  rest.post(`${API_BASE_URL}/auth/register`, (req, res, ctx) => {
    const { email, password, role } = req.body as { email: string; password: string; role: string };
    
    // Basic validation
    if (!email || !password) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Email and password are required'))
      );
    }
    
    // Check if email already exists
    if (mockUsers.some(u => u.email === email)) {
      return res(
        ctx.status(409),
        ctx.json(createErrorResponse('Email already registered', 409))
      );
    }
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(
        { email, role },
        'Registration successful. Please verify your email.'
      ))
    );
  }),
  
  rest.post(`${API_BASE_URL}/auth/refresh-token`, (req, res, ctx) => {
    const { refresh_token } = req.body as { refresh_token: string };
    
    if (!refresh_token) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Refresh token is required'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        access_token: 'new-mock-access-token',
        refresh_token: 'new-mock-refresh-token',
        expires_in: 900, // 15 minutes
      }))
    );
  }),
  
  rest.post(`${API_BASE_URL}/auth/logout`, (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, 'Logout successful'))
    );
  }),
  
  rest.post(`${API_BASE_URL}/auth/verify-email`, (req, res, ctx) => {
    const { token } = req.body as { token: string };
    
    if (!token) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Verification token is required'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, 'Email verification successful'))
    );
  }),
  
  // Molecule endpoints
  rest.get(`${API_BASE_URL}/molecules`, (req, res, ctx) => {
    const page = parseInt(req.url.searchParams.get('page') || '1');
    const limit = parseInt(req.url.searchParams.get('limit') || '10');
    const sortBy = req.url.searchParams.get('sortBy') || 'created_at';
    const sortOrder = req.url.searchParams.get('sortOrder') || 'desc';
    
    // Filter logic would go here in a real implementation
    // For now, just return paginated molecules
    const startIndex = (page - 1) * limit;
    const endIndex = Math.min(startIndex + limit, mockMolecules.length);
    const paginatedMolecules = mockMolecules.slice(startIndex, endIndex);
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        molecules: paginatedMolecules,
        pagination: {
          page,
          limit,
          total: mockMolecules.length,
          pages: Math.ceil(mockMolecules.length / limit),
        }
      }))
    );
  }),
  
  rest.get(`${API_BASE_URL}/molecules/:id`, (req, res, ctx) => {
    const { id } = req.params;
    
    const molecule = mockMolecules.find(m => m.id === id);
    if (!molecule) {
      return res(
        ctx.status(404),
        ctx.json(createErrorResponse('Molecule not found', 404))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(molecule))
    );
  }),
  
  rest.post(`${API_BASE_URL}/molecules`, (req, res, ctx) => {
    const moleculeData = req.body as {
      smiles: string;
      properties: Record<string, any>;
    };
    
    if (!moleculeData.smiles) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('SMILES is required'))
      );
    }
    
    // Create a new molecule with mock data
    const newMolecule = {
      id: (mockMolecules.length + 1).toString(),
      smiles: moleculeData.smiles,
      properties: moleculeData.properties || {},
      flag_status: null,
      created_by: '1',
      created_at: new Date().toISOString(),
    };
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(newMolecule, 'Molecule created successfully'))
    );
  }),
  
  rest.put(`${API_BASE_URL}/molecules/:id`, (req, res, ctx) => {
    const { id } = req.params;
    const updateData = req.body as {
      properties?: Record<string, any>;
      flag_status?: string | null;
    };
    
    const moleculeIndex = mockMolecules.findIndex(m => m.id === id);
    if (moleculeIndex === -1) {
      return res(
        ctx.status(404),
        ctx.json(createErrorResponse('Molecule not found', 404))
      );
    }
    
    const updatedMolecule = {
      ...mockMolecules[moleculeIndex],
      ...updateData,
    };
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(updatedMolecule, 'Molecule updated successfully'))
    );
  }),
  
  rest.delete(`${API_BASE_URL}/molecules/:id`, (req, res, ctx) => {
    const { id } = req.params;
    
    const moleculeIndex = mockMolecules.findIndex(m => m.id === id);
    if (moleculeIndex === -1) {
      return res(
        ctx.status(404),
        ctx.json(createErrorResponse('Molecule not found', 404))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, 'Molecule deleted successfully'))
    );
  }),
  
  rest.get(`${API_BASE_URL}/molecules/properties`, (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse([
        'molecular_weight',
        'logp',
        'solubility',
        'activity',
        'pka',
        'polar_surface_area',
        'h_bond_donors',
        'h_bond_acceptors',
        'rotatable_bonds',
      ]))
    );
  }),
  
  // CSV upload endpoints
  rest.post(`${API_BASE_URL}/csv/upload`, (req, res, ctx) => {
    // For testing, we'll just return a success response
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        file_id: 'mock-file-id',
        headers: ['SMILES', 'MW', 'LogP', 'Activity', 'Solubility', 'Notes'],
        row_count: 1245,
      }))
    );
  }),
  
  rest.post(`${API_BASE_URL}/csv/map`, (req, res, ctx) => {
    const { file_id, mapping } = req.body as {
      file_id: string;
      mapping: Record<string, string>;
    };
    
    if (!file_id || !mapping) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('File ID and mapping are required'))
      );
    }
    
    // Check if SMILES is mapped
    if (!Object.values(mapping).includes('smiles')) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('SMILES column must be mapped'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        job_id: 'mock-job-id',
        status: 'PROCESSING',
      }))
    );
  }),
  
  rest.get(`${API_BASE_URL}/csv/status/:jobId`, (req, res, ctx) => {
    const { jobId } = req.params;
    
    // Simulate different statuses based on jobId
    let status, progress;
    
    if (jobId === 'complete-job') {
      status = 'COMPLETED';
      progress = 100;
    } else if (jobId === 'error-job') {
      status = 'ERROR';
      progress = 50;
    } else {
      // Default to in progress
      status = 'PROCESSING';
      progress = 75;
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        job_id: jobId,
        status,
        progress,
        summary: status === 'COMPLETED' ? {
          total_rows: 1245,
          processed_rows: 1245,
          successful_rows: 1240,
          failed_rows: 5,
          molecules_created: 1240,
        } : null,
        error: status === 'ERROR' ? 'Failed to process CSV file' : null,
      }))
    );
  }),
  
  rest.get(`${API_BASE_URL}/csv/properties`, (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse([
        { name: 'smiles', label: 'SMILES', required: true },
        { name: 'molecular_weight', label: 'Molecular Weight', unit: 'g/mol' },
        { name: 'logp', label: 'LogP' },
        { name: 'solubility', label: 'Solubility', unit: 'mg/mL' },
        { name: 'activity', label: 'Activity', custom: true },
        { name: 'notes', label: 'Notes', custom: true },
      ]))
    );
  }),
  
  // Library endpoints
  rest.get(`${API_BASE_URL}/libraries`, (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse([
        {
          id: '1',
          name: 'High Activity',
          description: 'Molecules with high activity',
          created_by: '1',
          created_at: '2023-04-15T10:00:00Z',
          molecule_count: 24,
        },
        {
          id: '2',
          name: 'Alcohols',
          description: 'Molecules containing alcohol groups',
          created_by: '1',
          created_at: '2023-05-01T10:00:00Z',
          molecule_count: 12,
        },
      ]))
    );
  }),
  
  rest.post(`${API_BASE_URL}/libraries`, (req, res, ctx) => {
    const { name, description } = req.body as {
      name: string;
      description?: string;
    };
    
    if (!name) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Library name is required'))
      );
    }
    
    const newLibrary = {
      id: '3',
      name,
      description: description || '',
      created_by: '1',
      created_at: new Date().toISOString(),
      molecule_count: 0,
    };
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(newLibrary, 'Library created successfully'))
    );
  }),
  
  rest.get(`${API_BASE_URL}/libraries/:id`, (req, res, ctx) => {
    const { id } = req.params;
    
    const library = {
      id,
      name: id === '1' ? 'High Activity' : 'Alcohols',
      description: id === '1' ? 'Molecules with high activity' : 'Molecules containing alcohol groups',
      created_by: '1',
      created_at: '2023-05-01T10:00:00Z',
      molecule_count: id === '1' ? 24 : 12,
      molecules: mockMolecules.slice(0, 2),
    };
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(library))
    );
  }),
  
  rest.post(`${API_BASE_URL}/libraries/:id/molecules`, (req, res, ctx) => {
    const { id } = req.params;
    const { molecule_ids } = req.body as { molecule_ids: string[] };
    
    if (!molecule_ids || !molecule_ids.length) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Molecule IDs are required'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, `${molecule_ids.length} molecules added to library`))
    );
  }),
  
  // Experiment endpoints
  rest.get(`${API_BASE_URL}/experiments`, (req, res, ctx) => {
    const status = req.url.searchParams.get('status');
    
    const experiments = [
      {
        id: '1',
        name: 'Binding Study',
        type: 'Binding Assay',
        status: 'QUOTE_RECEIVED',
        created_by: '1',
        created_at: '2023-06-01T10:00:00Z',
        parameters: {
          target: 'Protein A',
          concentration: '10 μM',
          temperature: '25 °C',
        },
        molecule_count: 3,
      },
      {
        id: '2',
        name: 'ADME Screening',
        type: 'ADME Panel',
        status: 'IN_PROGRESS',
        created_by: '1',
        created_at: '2023-06-01T10:00:00Z',
        parameters: {
          concentration: '5 μM',
        },
        molecule_count: 5,
      },
    ];
    
    let filteredExperiments = [...experiments];
    if (status) {
      filteredExperiments = filteredExperiments.filter(e => e.status === status);
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(filteredExperiments))
    );
  }),
  
  rest.post(`${API_BASE_URL}/experiments`, (req, res, ctx) => {
    const { name, type, parameters } = req.body as {
      name: string;
      type: string;
      parameters?: Record<string, any>;
    };
    
    if (!name || !type) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Experiment name and type are required'))
      );
    }
    
    const newExperiment = {
      id: '3',
      name,
      type,
      status: 'DRAFT',
      created_by: '1',
      created_at: new Date().toISOString(),
      parameters: parameters || {},
      molecule_count: 0,
    };
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(newExperiment, 'Experiment created successfully'))
    );
  }),
  
  rest.post(`${API_BASE_URL}/experiments/:id/molecules`, (req, res, ctx) => {
    const { id } = req.params;
    const { molecule_ids } = req.body as { molecule_ids: string[] };
    
    if (!molecule_ids || !molecule_ids.length) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Molecule IDs are required'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, `${molecule_ids.length} molecules added to experiment`))
    );
  }),
  
  // Submission endpoints
  rest.get(`${API_BASE_URL}/submissions`, (req, res, ctx) => {
    const status = req.url.searchParams.get('status');
    
    const submissions = [
      {
        id: '1',
        experiment_id: '1',
        cro_id: '2',
        status: 'QUOTE_RECEIVED',
        submitted_at: '2023-06-01T12:00:00Z',
        quote: {
          price: 1250.0,
          currency: 'USD',
          estimated_days: 7,
        },
      },
      {
        id: '2',
        experiment_id: '2',
        cro_id: '2',
        status: 'IN_PROGRESS',
        submitted_at: '2023-06-01T14:00:00Z',
        approved_at: '2023-06-02T10:00:00Z',
      },
    ];
    
    let filteredSubmissions = [...submissions];
    if (status) {
      filteredSubmissions = filteredSubmissions.filter(s => s.status === status);
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(filteredSubmissions))
    );
  }),
  
  rest.post(`${API_BASE_URL}/submissions`, (req, res, ctx) => {
    const { experiment_id, cro_id } = req.body as {
      experiment_id: string;
      cro_id: string;
    };
    
    if (!experiment_id || !cro_id) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Experiment ID and CRO ID are required'))
      );
    }
    
    const newSubmission = {
      id: '3',
      experiment_id,
      cro_id,
      status: 'SUBMITTED',
      submitted_at: new Date().toISOString(),
    };
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(newSubmission, 'Submission created successfully'))
    );
  }),
  
  rest.post(`${API_BASE_URL}/submissions/:id/quote`, (req, res, ctx) => {
    const { id } = req.params;
    const { price, estimated_days } = req.body as {
      price: number;
      estimated_days: number;
    };
    
    if (!price || !estimated_days) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Price and estimated days are required'))
      );
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(null, 'Quote provided successfully'))
    );
  }),
  
  rest.post(`${API_BASE_URL}/submissions/:id/approve`, (req, res, ctx) => {
    const { id } = req.params;
    const { approved } = req.body as { approved: boolean };
    
    if (approved === undefined) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Approval status is required'))
      );
    }
    
    const newStatus = approved ? 'IN_PROGRESS' : 'QUOTE_REJECTED';
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(
        { status: newStatus },
        approved ? 'Quote approved successfully' : 'Quote rejected'
      ))
    );
  }),
  
  // Results endpoints
  rest.get(`${API_BASE_URL}/results`, (req, res, ctx) => {
    const submission_id = req.url.searchParams.get('submission_id');
    
    const results = [
      {
        id: '1',
        submission_id: '1',
        status: 'NEW',
        uploaded_at: '2023-05-20T10:00:00Z',
        summary: {
          average_binding: 65.45,
          successful_molecules: 3,
        },
        files: [
          {
            id: '1',
            name: 'binding_results.xlsx',
            type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
            size: 2457600,
            url: '/mock-files/binding_results.xlsx',
          },
          {
            id: '2',
            name: 'methodology.pdf',
            type: 'application/pdf',
            size: 1258291,
            url: '/mock-files/methodology.pdf',
          },
        ],
        molecule_results: [
          {
            molecule_id: '1',
            binding: 85.2,
            ic50: 12.3,
          },
          {
            molecule_id: '2',
            binding: 45.7,
            ic50: 78.5,
          },
        ],
      },
    ];
    
    let filteredResults = [...results];
    if (submission_id) {
      filteredResults = filteredResults.filter(r => r.submission_id === submission_id);
    }
    
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse(filteredResults))
    );
  }),
  
  rest.post(`${API_BASE_URL}/results`, (req, res, ctx) => {
    const { submission_id, files, structured_data } = req.body as {
      submission_id: string;
      files?: any[];
      structured_data?: Record<string, any>;
    };
    
    if (!submission_id) {
      return res(
        ctx.status(400),
        ctx.json(createErrorResponse('Submission ID is required'))
      );
    }
    
    return res(
      ctx.status(201),
      ctx.json(createSuccessResponse(
        { result_id: '2' },
        'Results uploaded successfully'
      ))
    );
  }),
  
  // Health check endpoint
  rest.get(`${API_BASE_URL}/health`, (req, res, ctx) => {
    return res(
      ctx.status(200),
      ctx.json(createSuccessResponse({
        status: 'UP',
        version: '1.0.0',
        timestamp: new Date().toISOString(),
      }))
    );
  }),
];

export default handlers;