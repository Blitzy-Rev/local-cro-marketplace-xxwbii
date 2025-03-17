# API Documentation - Molecular Data Management and CRO Integration Platform

## Introduction

This document provides comprehensive information about the REST API for the Molecular Data Management and CRO Integration Platform. The API enables programmatic access to all platform features including user management, molecular data manipulation, experiment submission, and result retrieval.

## API Conventions

### Base URL

All API endpoints are relative to the base URL: 

```
http://localhost/api
```

For production deployments, replace `localhost` with your server domain.

### HTTP Methods

The API uses standard HTTP methods:

- `GET`: Retrieve resources
- `POST`: Create resources
- `PUT`: Update resources
- `DELETE`: Remove resources

### Request Format

All request bodies should be in JSON format with the `Content-Type` header set to `application/json`.

### Response Format

All responses are in JSON format with the following standard structure:

```json
{
  "status": "success",
  "data": { ... }
}
```

For errors:

```json
{
  "status": "error",
  "error": {
    "code": "ERROR_CODE",
    "message": "Human-readable error message",
    "details": { ... }
  }
}
```

### HTTP Status Codes

| Code | Description |
|------|-------------|
| 200 | OK - Request succeeded |
| 201 | Created - Resource created successfully |
| 400 | Bad Request - Invalid request format or parameters |
| 401 | Unauthorized - Authentication required or failed |
| 403 | Forbidden - Authenticated user lacks permission |
| 404 | Not Found - Resource not found |
| 422 | Unprocessable Entity - Validation errors |
| 429 | Too Many Requests - Rate limit exceeded |
| 500 | Internal Server Error - Server error |

## Authentication

The API uses JWT (JSON Web Token) based authentication.

### Obtaining Tokens

To authenticate, obtain a JWT token by making a POST request to the `/auth/login` endpoint:

```
POST /api/auth/login
```

Request body:

```json
{
  "email": "user@example.com",
  "password": "your_password"
}
```

Response:

```json
{
  "status": "success",
  "data": {
    "token": "eyJhbGciOiJSUzI1...",
    "refresh_token": "eyJhbGciOiJSUzI1...",
    "user": {
      "id": 123,
      "email": "user@example.com",
      "role": "pharma"
    }
  }
}
```

### Using Tokens

For all authenticated requests, include the JWT token in the `Authorization` header:

```
Authorization: Bearer eyJhbGciOiJSUzI1...
```

### Token Refresh

When your access token expires, use the refresh token to obtain a new token:

```
POST /api/auth/refresh
```

Request body:

```json
{
  "refresh_token": "eyJhbGciOiJSUzI1..."
}
```

Response:

```json
{
  "status": "success",
  "data": {
    "token": "eyJhbGciOiJSUzI1...",
    "refresh_token": "eyJhbGciOiJSUzI1..."
  }
}
```

## Rate Limiting

API requests are rate-limited to protect the service. Limits vary by endpoint:

| Endpoint Category | Rate Limit |
|-------------------|------------|
| Authentication | 10 requests per minute |
| General API | 100 requests per minute |
| CSV Processing | 5 requests per minute |

Rate limit information is included in response headers:

- `X-RateLimit-Limit`: Maximum requests per window
- `X-RateLimit-Remaining`: Remaining requests in current window
- `X-RateLimit-Reset`: Seconds until window reset

When rate limits are exceeded, the API returns a 429 status code with a Retry-After header.

## Versioning

The API uses URL path versioning:

```
/api/v1/resource
```

For minor versions, use the `Accept` header:

```
Accept: application/json;version=1.2
```

## API Endpoints

### User Management

#### Register User

```
POST /api/auth/register
```

Create a new user account.

**Request Body:**

```json
{
  "email": "user@example.com",
  "password": "secure_password",
  "role": "pharma"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "User registered successfully. Please check your email for verification.",
    "user_id": 123
  }
}
```

#### Verify Email

```
GET /api/auth/verify?token=verification_token
```

Verify a user's email address.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Email verified successfully."
  }
}
```

#### Login

```
POST /api/auth/login
```

Authenticate a user and receive access token.

**Request Body:**

```json
{
  "email": "user@example.com",
  "password": "secure_password"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "token": "eyJhbGciOiJSUzI1...",
    "refresh_token": "eyJhbGciOiJSUzI1...",
    "user": {
      "id": 123,
      "email": "user@example.com",
      "role": "pharma"
    }
  }
}
```

#### Refresh Token

```
POST /api/auth/refresh
```

Obtain a new access token using a refresh token.

**Request Body:**

```json
{
  "refresh_token": "eyJhbGciOiJSUzI1..."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "token": "eyJhbGciOiJSUzI1...",
    "refresh_token": "eyJhbGciOiJSUzI1..."
  }
}
```

#### Forgot Password

```
POST /api/auth/forgot-password
```

Request a password reset email.

**Request Body:**

```json
{
  "email": "user@example.com"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "If the email exists in our system, a password reset link has been sent."
  }
}
```

#### Reset Password

```
POST /api/auth/reset-password
```

Reset password using token from email.

**Request Body:**

```json
{
  "token": "reset_token_from_email",
  "password": "new_secure_password"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Password has been reset successfully."
  }
}
```

#### Get Current User

```
GET /api/auth/me
```

Get the current authenticated user's details.

**Response:**

```json
{
  "status": "success",
  "data": {
    "user": {
      "id": 123,
      "email": "user@example.com",
      "role": "pharma",
      "created_at": "2023-01-15T12:00:00Z",
      "last_login": "2023-06-01T09:30:00Z"
    }
  }
}
```

#### Update User Profile

```
PUT /api/auth/me
```

Update the current user's profile information.

**Request Body:**

```json
{
  "name": "John Doe",
  "organization": "PharmaCorp"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "user": {
      "id": 123,
      "email": "user@example.com",
      "name": "John Doe",
      "organization": "PharmaCorp",
      "role": "pharma"
    }
  }
}
```

#### Change Password

```
PUT /api/auth/password
```

Change the current user's password.

**Request Body:**

```json
{
  "current_password": "current_secure_password",
  "new_password": "new_secure_password"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Password changed successfully."
  }
}
```

### Molecule Management

#### List Molecules

```
GET /api/molecules
```

Retrieve a paginated list of molecules.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 50, max: 100) |
| sort | string | Field to sort by (e.g., `created_at`) |
| order | string | Sort order (`asc` or `desc`) |
| filter | object | Filter criteria (e.g., `{"property_name": {"min": 0, "max": 100}}`) |

**Response:**

```json
{
  "status": "success",
  "data": {
    "molecules": [
      {
        "id": 1,
        "smiles": "CCO",
        "properties": {
          "molecular_weight": 46.07,
          "logp": -0.14,
          "solubility": 3.45
        },
        "created_at": "2023-05-10T14:30:00Z",
        "created_by": 123,
        "flag_status": null
      },
      // ... more molecules
    ],
    "pagination": {
      "total": 1245,
      "pages": 25,
      "page": 1,
      "limit": 50,
      "prev": null,
      "next": 2
    }
  }
}
```

#### Get Molecule

```
GET /api/molecules/{id}
```

Retrieve a specific molecule by ID.

**Response:**

```json
{
  "status": "success",
  "data": {
    "molecule": {
      "id": 1,
      "smiles": "CCO",
      "properties": {
        "molecular_weight": 46.07,
        "logp": -0.14,
        "solubility": 3.45,
        "custom_property": "value"
      },
      "created_at": "2023-05-10T14:30:00Z",
      "created_by": 123,
      "flag_status": null,
      "libraries": [
        {
          "id": 5,
          "name": "High Activity"
        }
      ],
      "experiments": [
        {
          "id": 10,
          "name": "Binding Assay",
          "status": "completed"
        }
      ]
    }
  }
}
```

#### Create Molecule

```
POST /api/molecules
```

Create a new molecule.

**Request Body:**

```json
{
  "smiles": "CCO",
  "properties": {
    "custom_property": "value"
  }
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "molecule": {
      "id": 1,
      "smiles": "CCO",
      "properties": {
        "molecular_weight": 46.07,
        "logp": -0.14,
        "solubility": 3.45,
        "custom_property": "value"
      },
      "created_at": "2023-06-05T14:30:00Z",
      "created_by": 123,
      "flag_status": null
    }
  }
}
```

#### Update Molecule

```
PUT /api/molecules/{id}
```

Update an existing molecule.

**Request Body:**

```json
{
  "properties": {
    "custom_property": "new_value"
  },
  "flag_status": "important"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "molecule": {
      "id": 1,
      "smiles": "CCO",
      "properties": {
        "molecular_weight": 46.07,
        "logp": -0.14,
        "solubility": 3.45,
        "custom_property": "new_value"
      },
      "created_at": "2023-05-10T14:30:00Z",
      "created_by": 123,
      "flag_status": "important"
    }
  }
}
```

#### Delete Molecule

```
DELETE /api/molecules/{id}
```

Delete a molecule.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Molecule deleted successfully."
  }
}
```

#### Flag Molecule

```
PUT /api/molecules/{id}/flag
```

Flag a molecule for quick reference.

**Request Body:**

```json
{
  "flag_status": "important"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "molecule": {
      "id": 1,
      "flag_status": "important"
    }
  }
}
```

#### Get Molecule Image

```
GET /api/molecules/{id}/image
```

Retrieve a molecular structure image.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| format | string | Image format (`svg` or `png`, default: `svg`) |
| width | integer | Image width in pixels (default: 300) |
| height | integer | Image height in pixels (default: 200) |

**Response:**

Binary image data with appropriate Content-Type header (`image/svg+xml` or `image/png`).

### CSV Management

#### Upload CSV

```
POST /api/csv/upload
```

Upload a CSV file containing molecule data.

**Request:**

Multipart form-data with a `file` field containing the CSV file.

**Response:**

```json
{
  "status": "success",
  "data": {
    "file_id": "a1b2c3d4",
    "filename": "molecules_batch_12.csv",
    "size": 2359045,
    "headers": ["SMILES", "MW", "LogP", "Activity", "Solubility", "Notes"]
  }
}
```

#### Map CSV Headers

```
POST /api/csv/map
```

Map CSV headers to system properties.

**Request Body:**

```json
{
  "file_id": "a1b2c3d4",
  "mapping": {
    "SMILES": "smiles",
    "MW": "molecular_weight",
    "LogP": "logp",
    "Activity": "custom:activity",
    "Solubility": "solubility",
    "Notes": "custom:notes"
  }
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "job_id": "job123",
    "message": "CSV processing started."
  }
}
```

#### Check Import Status

```
GET /api/csv/status/{job_id}
```

Check the status of a CSV import job.

**Response:**

```json
{
  "status": "success",
  "data": {
    "job_id": "job123",
    "status": "processing", // "pending", "processing", "completed", "failed"
    "progress": 65, // Percentage complete
    "processed_rows": 6500,
    "total_rows": 10000,
    "message": "Processing CSV data."
  }
}
```

#### Get Import Results

```
GET /api/csv/results/{job_id}
```

Get the results of a completed CSV import job.

**Response:**

```json
{
  "status": "success",
  "data": {
    "job_id": "job123",
    "status": "completed",
    "processed_rows": 10000,
    "total_rows": 10000,
    "successful_imports": 9850,
    "failed_imports": 150,
    "errors": [
      {
        "row": 45,
        "message": "Invalid SMILES string."
      },
      // More errors...
    ],
    "molecules": [
      {
        "id": 1001,
        "smiles": "CCO"
      },
      // Limited list of imported molecules
    ]
  }
}
```

### Library Management

#### List Libraries

```
GET /api/libraries
```

Retrieve a list of the user's molecule libraries.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |

**Response:**

```json
{
  "status": "success",
  "data": {
    "libraries": [
      {
        "id": 1,
        "name": "High Activity",
        "description": "Molecules with high binding activity",
        "molecule_count": 24,
        "created_at": "2023-04-15T10:20:30Z",
        "created_by": 123
      },
      // More libraries...
    ],
    "pagination": {
      "total": 5,
      "pages": 1,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": null
    }
  }
}
```

#### Create Library

```
POST /api/libraries
```

Create a new molecule library.

**Request Body:**

```json
{
  "name": "Alcohols",
  "description": "Collection of alcohol compounds"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "library": {
      "id": 2,
      "name": "Alcohols",
      "description": "Collection of alcohol compounds",
      "molecule_count": 0,
      "created_at": "2023-06-05T15:10:20Z",
      "created_by": 123
    }
  }
}
```

#### Get Library

```
GET /api/libraries/{id}
```

Retrieve a specific library.

**Response:**

```json
{
  "status": "success",
  "data": {
    "library": {
      "id": 1,
      "name": "High Activity",
      "description": "Molecules with high binding activity",
      "molecule_count": 24,
      "created_at": "2023-04-15T10:20:30Z",
      "created_by": 123
    }
  }
}
```

#### Update Library

```
PUT /api/libraries/{id}
```

Update a library's information.

**Request Body:**

```json
{
  "name": "High Activity Compounds",
  "description": "Updated description for high activity molecules"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "library": {
      "id": 1,
      "name": "High Activity Compounds",
      "description": "Updated description for high activity molecules",
      "molecule_count": 24,
      "created_at": "2023-04-15T10:20:30Z",
      "created_by": 123
    }
  }
}
```

#### Delete Library

```
DELETE /api/libraries/{id}
```

Delete a library.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Library deleted successfully."
  }
}
```

#### Get Library Molecules

```
GET /api/libraries/{id}/molecules
```

Retrieve molecules in a library.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 50, max: 100) |
| sort | string | Field to sort by (e.g., `created_at`) |
| order | string | Sort order (`asc` or `desc`) |

**Response:**

```json
{
  "status": "success",
  "data": {
    "library": {
      "id": 1,
      "name": "High Activity"
    },
    "molecules": [
      {
        "id": 1,
        "smiles": "CCO",
        "properties": {
          "molecular_weight": 46.07,
          "logp": -0.14
        }
      },
      // More molecules...
    ],
    "pagination": {
      "total": 24,
      "pages": 1,
      "page": 1,
      "limit": 50,
      "prev": null,
      "next": null
    }
  }
}
```

#### Add Molecules to Library

```
POST /api/libraries/{id}/molecules
```

Add molecules to a library.

**Request Body:**

```json
{
  "molecule_ids": [1, 2, 3, 4, 5]
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "5 molecules added to library.",
    "library": {
      "id": 1,
      "name": "High Activity",
      "molecule_count": 29
    }
  }
}
```

#### Remove Molecules from Library

```
DELETE /api/libraries/{id}/molecules
```

Remove molecules from a library.

**Request Body:**

```json
{
  "molecule_ids": [4, 5]
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "2 molecules removed from library.",
    "library": {
      "id": 1,
      "name": "High Activity",
      "molecule_count": 27
    }
  }
}
```

#### Export Library

```
GET /api/libraries/{id}/export
```

Export library molecules to a CSV file.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| format | string | Export format (`csv` or `sdf`, default: `csv`) |
| properties | string | Comma-separated list of properties to include |

**Response:**

Binary file data with appropriate Content-Type header (`text/csv` or `chemical/x-mdl-sdfile`).

### Experiment Management

#### List Experiment Types

```
GET /api/experiment-types
```

Retrieve available experiment types.

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiment_types": [
      {
        "id": 1,
        "name": "Binding Assay",
        "description": "Measures binding affinity to target protein",
        "category": "Target Binding",
        "parameters": [
          {
            "name": "target",
            "type": "string",
            "description": "Target protein",
            "required": true,
            "options": ["Protein A", "Protein B", "Protein C"]
          },
          {
            "name": "concentration",
            "type": "number",
            "description": "Concentration in μM",
            "required": true
          }
        ]
      },
      // More experiment types...
    ]
  }
}
```

#### List Experiments

```
GET /api/experiments
```

Retrieve a list of experiments.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |
| status | string | Filter by status |
| type_id | integer | Filter by experiment type |

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiments": [
      {
        "id": 1,
        "name": "Binding Study",
        "type": {
          "id": 1,
          "name": "Binding Assay"
        },
        "status": "completed",
        "molecule_count": 3,
        "created_at": "2023-05-10T09:15:00Z",
        "created_by": 123
      },
      // More experiments...
    ],
    "pagination": {
      "total": 12,
      "pages": 1,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": null
    }
  }
}
```

#### Create Experiment

```
POST /api/experiments
```

Create a new experiment.

**Request Body:**

```json
{
  "name": "ADME Screening",
  "type_id": 2,
  "parameters": {
    "concentration": 10,
    "temperature": 25
  },
  "molecule_ids": [1, 2, 3]
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiment": {
      "id": 2,
      "name": "ADME Screening",
      "type": {
        "id": 2,
        "name": "ADME Panel"
      },
      "status": "draft",
      "parameters": {
        "concentration": 10,
        "temperature": 25
      },
      "molecule_count": 3,
      "created_at": "2023-06-05T16:20:30Z",
      "created_by": 123
    }
  }
}
```

#### Get Experiment

```
GET /api/experiments/{id}
```

Retrieve a specific experiment.

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiment": {
      "id": 1,
      "name": "Binding Study",
      "type": {
        "id": 1,
        "name": "Binding Assay"
      },
      "status": "completed",
      "parameters": {
        "target": "Protein A",
        "concentration": 10,
        "temperature": 25
      },
      "molecule_count": 3,
      "created_at": "2023-05-10T09:15:00Z",
      "created_by": 123,
      "submission": {
        "id": 1,
        "status": "completed",
        "cro_id": 456
      },
      "results": {
        "id": 1,
        "uploaded_at": "2023-05-20T14:30:00Z"
      }
    }
  }
}
```

#### Update Experiment

```
PUT /api/experiments/{id}
```

Update an experiment.

**Request Body:**

```json
{
  "name": "Binding Study - Updated",
  "parameters": {
    "temperature": 30
  }
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiment": {
      "id": 1,
      "name": "Binding Study - Updated",
      "type": {
        "id": 1,
        "name": "Binding Assay"
      },
      "status": "draft",
      "parameters": {
        "target": "Protein A",
        "concentration": 10,
        "temperature": 30
      },
      "molecule_count": 3,
      "created_at": "2023-05-10T09:15:00Z",
      "updated_at": "2023-06-05T16:30:40Z",
      "created_by": 123
    }
  }
}
```

#### Delete Experiment

```
DELETE /api/experiments/{id}
```

Delete an experiment.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Experiment deleted successfully."
  }
}
```

#### Get Experiment Molecules

```
GET /api/experiments/{id}/molecules
```

Retrieve molecules in an experiment.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 50, max: 100) |

**Response:**

```json
{
  "status": "success",
  "data": {
    "experiment": {
      "id": 1,
      "name": "Binding Study"
    },
    "molecules": [
      {
        "id": 1,
        "smiles": "CCO",
        "properties": {
          "molecular_weight": 46.07,
          "logp": -0.14
        }
      },
      // More molecules...
    ],
    "pagination": {
      "total": 3,
      "pages": 1,
      "page": 1,
      "limit": 50,
      "prev": null,
      "next": null
    }
  }
}
```

#### Add Molecules to Experiment

```
POST /api/experiments/{id}/molecules
```

Add molecules to an experiment.

**Request Body:**

```json
{
  "molecule_ids": [4, 5]
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "2 molecules added to experiment.",
    "experiment": {
      "id": 1,
      "name": "Binding Study",
      "molecule_count": 5
    }
  }
}
```

#### Remove Molecules from Experiment

```
DELETE /api/experiments/{id}/molecules
```

Remove molecules from an experiment.

**Request Body:**

```json
{
  "molecule_ids": [4, 5]
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "2 molecules removed from experiment.",
    "experiment": {
      "id": 1,
      "name": "Binding Study",
      "molecule_count": 3
    }
  }
}
```

#### Submit Experiment to Queue

```
POST /api/experiments/{id}/queue
```

Add an experiment to the queue for CRO submission.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Experiment added to queue.",
    "experiment": {
      "id": 1,
      "name": "Binding Study",
      "status": "queued"
    }
  }
}
```

### CRO Submission Management

#### List CROs

```
GET /api/cros
```

Retrieve a list of available CROs.

**Response:**

```json
{
  "status": "success",
  "data": {
    "cros": [
      {
        "id": 456,
        "name": "BioAssay Labs",
        "supported_experiments": [1, 2, 3]
      },
      {
        "id": 457,
        "name": "ChemTest Inc",
        "supported_experiments": [2, 4, 5]
      }
    ]
  }
}
```

#### List Submissions

```
GET /api/submissions
```

Retrieve a list of submissions.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |
| status | string | Filter by status |
| cro_id | integer | Filter by CRO ID |

**Response:**

```json
{
  "status": "success",
  "data": {
    "submissions": [
      {
        "id": 1,
        "experiment": {
          "id": 1,
          "name": "Binding Study"
        },
        "cro": {
          "id": 456,
          "name": "BioAssay Labs"
        },
        "status": "completed",
        "submitted_at": "2023-05-12T10:30:00Z",
        "updated_at": "2023-05-20T14:30:00Z"
      },
      // More submissions...
    ],
    "pagination": {
      "total": 5,
      "pages": 1,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": null
    }
  }
}
```

#### Create Submission

```
POST /api/submissions
```

Create a new CRO submission.

**Request Body:**

```json
{
  "experiment_id": 2,
  "cro_id": 456,
  "details": {
    "notes": "Please process with high priority",
    "contact_name": "John Doe"
  }
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "submission": {
      "id": 2,
      "experiment": {
        "id": 2,
        "name": "ADME Screening"
      },
      "cro": {
        "id": 456,
        "name": "BioAssay Labs"
      },
      "status": "submitted",
      "details": {
        "notes": "Please process with high priority",
        "contact_name": "John Doe"
      },
      "submitted_at": "2023-06-05T16:45:00Z",
      "created_by": 123
    }
  }
}
```

#### Get Submission

```
GET /api/submissions/{id}
```

Retrieve a specific submission.

**Response:**

```json
{
  "status": "success",
  "data": {
    "submission": {
      "id": 1,
      "experiment": {
        "id": 1,
        "name": "Binding Study",
        "type": {
          "id": 1,
          "name": "Binding Assay"
        },
        "parameters": {
          "target": "Protein A",
          "concentration": 10,
          "temperature": 25
        },
        "molecule_count": 3
      },
      "cro": {
        "id": 456,
        "name": "BioAssay Labs"
      },
      "status": "completed",
      "details": {
        "notes": "Process as standard priority"
      },
      "pricing": {
        "amount": 1250.00,
        "currency": "USD",
        "turnaround_days": 7
      },
      "submitted_at": "2023-05-12T10:30:00Z",
      "updated_at": "2023-05-20T14:30:00Z",
      "created_by": 123
    }
  }
}
```

#### Update Submission

```
PUT /api/submissions/{id}
```

Update a submission.

**Request Body:**

```json
{
  "details": {
    "notes": "Updated notes for this submission"
  }
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "submission": {
      "id": 1,
      "status": "submitted",
      "details": {
        "notes": "Updated notes for this submission"
      },
      "updated_at": "2023-06-05T17:00:00Z"
    }
  }
}
```

#### Cancel Submission

```
POST /api/submissions/{id}/cancel
```

Cancel a submission.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Submission cancelled successfully.",
    "submission": {
      "id": 1,
      "status": "cancelled",
      "updated_at": "2023-06-05T17:10:00Z"
    }
  }
}
```

#### Provide Quote (CRO)

```
POST /api/submissions/{id}/quote
```

Provide a price quote for a submission (CRO role only).

**Request Body:**

```json
{
  "pricing": {
    "amount": 1250.00,
    "currency": "USD",
    "turnaround_days": 7
  },
  "notes": "Quote includes all requested analyses."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Quote provided successfully.",
    "submission": {
      "id": 1,
      "status": "quote_provided",
      "pricing": {
        "amount": 1250.00,
        "currency": "USD",
        "turnaround_days": 7
      },
      "updated_at": "2023-06-05T17:20:00Z"
    }
  }
}
```

#### Approve Quote (Pharma)

```
POST /api/submissions/{id}/approve
```

Approve a quote and authorize work to begin (Pharma role only).

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Quote approved successfully.",
    "submission": {
      "id": 1,
      "status": "in_progress",
      "updated_at": "2023-06-05T17:30:00Z"
    }
  }
}
```

#### Reject Quote (Pharma)

```
POST /api/submissions/{id}/reject
```

Reject a quote (Pharma role only).

**Request Body:**

```json
{
  "reason": "Price is too high, please provide a more competitive quote."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Quote rejected successfully.",
    "submission": {
      "id": 1,
      "status": "quote_rejected",
      "updated_at": "2023-06-05T17:40:00Z"
    }
  }
}
```

#### Update Submission Status (CRO)

```
POST /api/submissions/{id}/status
```

Update the status of a submission (CRO role only).

**Request Body:**

```json
{
  "status": "in_progress",
  "progress": 75,
  "notes": "75% complete, on track for delivery."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Status updated successfully.",
    "submission": {
      "id": 1,
      "status": "in_progress",
      "progress": 75,
      "updated_at": "2023-06-05T17:50:00Z"
    }
  }
}
```

### Result Management

#### List Results

```
GET /api/results
```

Retrieve a list of experiment results.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |
| status | string | Filter by status |
| experiment_id | integer | Filter by experiment ID |

**Response:**

```json
{
  "status": "success",
  "data": {
    "results": [
      {
        "id": 1,
        "submission": {
          "id": 1,
          "experiment": {
            "id": 1,
            "name": "Binding Study"
          }
        },
        "status": "available",
        "uploaded_at": "2023-05-20T14:30:00Z",
        "uploaded_by": 456
      },
      // More results...
    ],
    "pagination": {
      "total": 3,
      "pages": 1,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": null
    }
  }
}
```

#### Get Result

```
GET /api/results/{id}
```

Retrieve a specific result.

**Response:**

```json
{
  "status": "success",
  "data": {
    "result": {
      "id": 1,
      "submission": {
        "id": 1,
        "experiment": {
          "id": 1,
          "name": "Binding Study",
          "type": {
            "id": 1,
            "name": "Binding Assay"
          }
        },
        "cro": {
          "id": 456,
          "name": "BioAssay Labs"
        }
      },
      "status": "available",
      "files": [
        {
          "id": 1,
          "filename": "binding_results.xlsx",
          "file_type": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
          "file_size": 2457600,
          "upload_date": "2023-05-20T14:30:00Z"
        },
        {
          "id": 2,
          "filename": "methodology.pdf",
          "file_type": "application/pdf",
          "file_size": 1258291,
          "upload_date": "2023-05-20T14:30:00Z"
        }
      ],
      "molecule_data": [
        {
          "molecule_id": 1,
          "molecule_smiles": "CCO",
          "data": {
            "binding": 85.2,
            "ic50": 12.3
          }
        },
        {
          "molecule_id": 2,
          "molecule_smiles": "CCCCO",
          "data": {
            "binding": 45.7,
            "ic50": 78.5
          }
        },
        {
          "molecule_id": 3,
          "molecule_smiles": "c1ccccc1",
          "data": {
            "binding": 92.1,
            "ic50": 5.6
          }
        }
      ],
      "notes": "All experiments performed according to protocol.",
      "uploaded_at": "2023-05-20T14:30:00Z",
      "uploaded_by": 456
    }
  }
}
```

#### Upload Result (CRO)

```
POST /api/results
```

Upload experiment results (CRO role only).

**Request:**

Multipart form-data with:
- `submission_id`: ID of the submission
- `notes`: Notes about the results
- `files[]`: Result files
- `molecule_data`: JSON string containing structured result data

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Results uploaded successfully.",
    "result": {
      "id": 2,
      "submission_id": 2,
      "status": "available",
      "uploaded_at": "2023-06-05T18:00:00Z",
      "uploaded_by": 456
    }
  }
}
```

#### Update Result Status

```
PUT /api/results/{id}/status
```

Update the status of a result.

**Request Body:**

```json
{
  "status": "reviewed",
  "notes": "Results reviewed and accepted."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Result status updated successfully.",
    "result": {
      "id": 1,
      "status": "reviewed",
      "updated_at": "2023-06-05T18:10:00Z"
    }
  }
}
```

#### Download Result File

```
GET /api/results/{id}/files/{file_id}
```

Download a specific result file.

**Response:**

Binary file data with appropriate Content-Type header.

#### Request Additional Data

```
POST /api/results/{id}/request-data
```

Request additional data for a result.

**Request Body:**

```json
{
  "request": "Please provide additional data on the binding mechanism."
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Additional data request sent successfully.",
    "result": {
      "id": 1,
      "status": "data_requested",
      "updated_at": "2023-06-05T18:20:00Z"
    }
  }
}
```

### Notification Management

#### List Notifications

```
GET /api/notifications
```

Retrieve user notifications.

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |
| read | boolean | Filter by read status |

**Response:**

```json
{
  "status": "success",
  "data": {
    "notifications": [
      {
        "id": 1,
        "type": "quote_received",
        "message": "Quote received for submission #1",
        "read_status": false,
        "data": {
          "submission_id": 1
        },
        "created_at": "2023-06-01T10:30:00Z"
      },
      // More notifications...
    ],
    "pagination": {
      "total": 15,
      "pages": 1,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": null
    }
  }
}
```

#### Mark Notification as Read

```
PUT /api/notifications/{id}/read
```

Mark a notification as read.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Notification marked as read.",
    "notification": {
      "id": 1,
      "read_status": true,
      "read_at": "2023-06-05T18:30:00Z"
    }
  }
}
```

#### Mark All Notifications as Read

```
PUT /api/notifications/read-all
```

Mark all notifications as read.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "All notifications marked as read.",
    "count": 10
  }
}
```

#### Delete Notification

```
DELETE /api/notifications/{id}
```

Delete a notification.

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Notification deleted successfully."
  }
}
```

### Admin Endpoints

#### List Users (Admin)

```
GET /api/admin/users
```

Retrieve a list of users (Admin role only).

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 20, max: 50) |
| role | string | Filter by role |
| status | string | Filter by status |

**Response:**

```json
{
  "status": "success",
  "data": {
    "users": [
      {
        "id": 123,
        "email": "user@example.com",
        "role": "pharma",
        "status": "active",
        "created_at": "2023-01-15T12:00:00Z",
        "last_login": "2023-06-01T09:30:00Z"
      },
      // More users...
    ],
    "pagination": {
      "total": 45,
      "pages": 3,
      "page": 1,
      "limit": 20,
      "prev": null,
      "next": 2
    }
  }
}
```

#### Get User (Admin)

```
GET /api/admin/users/{id}
```

Retrieve a specific user (Admin role only).

**Response:**

```json
{
  "status": "success",
  "data": {
    "user": {
      "id": 123,
      "email": "user@example.com",
      "role": "pharma",
      "status": "active",
      "name": "John Doe",
      "organization": "PharmaCorp",
      "created_at": "2023-01-15T12:00:00Z",
      "last_login": "2023-06-01T09:30:00Z",
      "login_count": 45,
      "failed_login_attempts": 2,
      "molecules_count": 156,
      "experiments_count": 12
    }
  }
}
```

#### Update User (Admin)

```
PUT /api/admin/users/{id}
```

Update a user (Admin role only).

**Request Body:**

```json
{
  "role": "admin",
  "status": "active"
}
```

**Response:**

```json
{
  "status": "success",
  "data": {
    "user": {
      "id": 123,
      "email": "user@example.com",
      "role": "admin",
      "status": "active",
      "updated_at": "2023-06-05T18:40:00Z"
    }
  }
}
```

#### Delete User (Admin)

```
DELETE /api/admin/users/{id}
```

Delete a user (Admin role only).

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "User deleted successfully."
  }
}
```

#### Reset User Password (Admin)

```
POST /api/admin/users/{id}/reset-password
```

Reset a user's password (Admin role only).

**Response:**

```json
{
  "status": "success",
  "data": {
    "message": "Password reset email sent to user."
  }
}
```

#### System Statistics (Admin)

```
GET /api/admin/stats
```

Retrieve system statistics (Admin role only).

**Response:**

```json
{
  "status": "success",
  "data": {
    "users": {
      "total": 45,
      "pharma": 35,
      "cro": 8,
      "admin": 2,
      "active_now": 12
    },
    "molecules": {
      "total": 10542,
      "last_24h": 156
    },
    "experiments": {
      "total": 245,
      "draft": 45,
      "queued": 15,
      "submitted": 30,
      "in_progress": 20,
      "completed": 135
    },
    "system": {
      "database_size": "2.4 GB",
      "storage_used": "15.7 GB",
      "storage_capacity": "50 GB",
      "cpu_usage": 32,
      "memory_usage": 45
    }
  }
}
```

#### System Logs (Admin)

```
GET /api/admin/logs
```

Retrieve system logs (Admin role only).

**Query Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| page | integer | Page number (default: 1) |
| limit | integer | Items per page (default: 50, max: 100) |
| level | string | Filter by log level |
| start_date | string | Filter by start date (ISO format) |
| end_date | string | Filter by end date (ISO format) |

**Response:**

```json
{
  "status": "success",
  "data": {
    "logs": [
      {
        "timestamp": "2023-06-05T18:32:15Z",
        "level": "INFO",
        "message": "User login: user@example.com",
        "service": "auth",
        "additional_data": {
          "user_id": 123,
          "ip": "192.168.1.1"
        }
      },
      // More logs...
    ],
    "pagination": {
      "total": 1245,
      "pages": 25,
      "page": 1,
      "limit": 50,
      "prev": null,
      "next": 2
    }
  }
}
```

## Client Code Examples

### Authentication Example

**JavaScript (Fetch API):**

```javascript
// Login and get token
async function login(email, password) {
  try {
    const response = await fetch('http://localhost/api/auth/login', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ email, password }),
    });
    
    const data = await response.json();
    
    if (data.status === 'success') {
      // Store token for future requests
      localStorage.setItem('token', data.data.token);
      localStorage.setItem('refreshToken', data.data.refresh_token);
      return data.data.user;
    } else {
      throw new Error(data.error.message);
    }
  } catch (error) {
    console.error('Login failed:', error);
    throw error;
  }
}

// Make authenticated request
async function fetchMolecules(page = 1, limit = 50) {
  try {
    const token = localStorage.getItem('token');
    
    if (!token) {
      throw new Error('Not authenticated');
    }
    
    const response = await fetch(`http://localhost/api/molecules?page=${page}&limit=${limit}`, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${token}`,
      },
    });
    
    // Handle 401 (token expired)
    if (response.status === 401) {
      const refreshed = await refreshToken();
      if (refreshed) {
        return fetchMolecules(page, limit);
      } else {
        throw new Error('Session expired');
      }
    }
    
    const data = await response.json();
    
    if (data.status === 'success') {
      return data.data.molecules;
    } else {
      throw new Error(data.error.message);
    }
  } catch (error) {
    console.error('Fetch failed:', error);
    throw error;
  }
}

// Refresh token
async function refreshToken() {
  try {
    const refreshToken = localStorage.getItem('refreshToken');
    
    if (!refreshToken) {
      return false;
    }
    
    const response = await fetch('http://localhost/api/auth/refresh', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ refresh_token: refreshToken }),
    });
    
    const data = await response.json();
    
    if (data.status === 'success') {
      localStorage.setItem('token', data.data.token);
      localStorage.setItem('refreshToken', data.data.refresh_token);
      return true;
    } else {
      localStorage.removeItem('token');
      localStorage.removeItem('refreshToken');
      return false;
    }
  } catch (error) {
    console.error('Token refresh failed:', error);
    localStorage.removeItem('token');
    localStorage.removeItem('refreshToken');
    return false;
  }
}
```

**Python (Requests):**

```python
import requests

class MolecularAPI:
    def __init__(self, base_url="http://localhost/api"):
        self.base_url = base_url
        self.token = None
        self.refresh_token = None
        
    def login(self, email, password):
        """Authenticate and get token"""
        url = f"{self.base_url}/auth/login"
        response = requests.post(url, json={
            "email": email,
            "password": password
        })
        
        data = response.json()
        
        if data["status"] == "success":
            self.token = data["data"]["token"]
            self.refresh_token = data["data"]["refresh_token"]
            return data["data"]["user"]
        else:
            raise Exception(data["error"]["message"])
            
    def refresh(self):
        """Refresh the access token"""
        if not self.refresh_token:
            return False
            
        url = f"{self.base_url}/auth/refresh"
        response = requests.post(url, json={
            "refresh_token": self.refresh_token
        })
        
        data = response.json()
        
        if data["status"] == "success":
            self.token = data["data"]["token"]
            self.refresh_token = data["data"]["refresh_token"]
            return True
        else:
            self.token = None
            self.refresh_token = None
            return False
            
    def get_molecules(self, page=1, limit=50):
        """Get molecules with authentication"""
        if not self.token:
            raise Exception("Not authenticated")
            
        url = f"{self.base_url}/molecules"
        headers = {"Authorization": f"Bearer {self.token}"}
        params = {"page": page, "limit": limit}
        
        response = requests.get(url, headers=headers, params=params)
        
        # Handle expired token
        if response.status_code == 401:
            refreshed = self.refresh()
            if refreshed:
                return self.get_molecules(page, limit)
            else:
                raise Exception("Session expired")
                
        data = response.json()
        
        if data["status"] == "success":
            return data["data"]["molecules"]
        else:
            raise Exception(data["error"]["message"])
```

### CSV Upload Example

**JavaScript:**

```javascript
async function uploadCSV(file) {
  try {
    const token = localStorage.getItem('token');
    
    if (!token) {
      throw new Error('Not authenticated');
    }
    
    // Step 1: Upload the file
    const formData = new FormData();
    formData.append('file', file);
    
    const uploadResponse = await fetch('http://localhost/api/csv/upload', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
      },
      body: formData,
    });
    
    const uploadData = await uploadResponse.json();
    
    if (uploadData.status !== 'success') {
      throw new Error(uploadData.error.message);
    }
    
    const fileId = uploadData.data.file_id;
    const headers = uploadData.data.headers;
    
    // Step 2: Map the headers
    const mapping = {};
    headers.forEach(header => {
      // Simple auto-mapping
      if (header.toLowerCase() === 'smiles') {
        mapping[header] = 'smiles';
      } else if (header.toLowerCase() === 'mw' || header.toLowerCase() === 'molecular weight') {
        mapping[header] = 'molecular_weight';
      } else if (header.toLowerCase() === 'logp') {
        mapping[header] = 'logp';
      } else {
        mapping[header] = `custom:${header.toLowerCase()}`;
      }
    });
    
    const mapResponse = await fetch('http://localhost/api/csv/map', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        file_id: fileId,
        mapping: mapping,
      }),
    });
    
    const mapData = await mapResponse.json();
    
    if (mapData.status !== 'success') {
      throw new Error(mapData.error.message);
    }
    
    const jobId = mapData.data.job_id;
    
    // Step 3: Poll for job status
    return pollJobStatus(jobId);
  } catch (error) {
    console.error('CSV upload failed:', error);
    throw error;
  }
}

async function pollJobStatus(jobId) {
  try {
    const token = localStorage.getItem('token');
    
    if (!token) {
      throw new Error('Not authenticated');
    }
    
    const response = await fetch(`http://localhost/api/csv/status/${jobId}`, {
      method: 'GET',
      headers: {
        'Authorization': `Bearer ${token}`,
      },
    });
    
    const data = await response.json();
    
    if (data.status !== 'success') {
      throw new Error(data.error.message);
    }
    
    const jobStatus = data.data;
    
    if (jobStatus.status === 'completed') {
      return jobStatus;
    } else if (jobStatus.status === 'failed') {
      throw new Error(`Job failed: ${jobStatus.message}`);
    } else {
      // Still processing, wait and try again
      await new Promise(resolve => setTimeout(resolve, 2000));
      return pollJobStatus(jobId);
    }
  } catch (error) {
    console.error('Job status check failed:', error);
    throw error;
  }
}
```

**Python:**

```python
import time
import requests

def upload_csv(api, file_path):
    """Upload a CSV file and process it"""
    if not api.token:
        raise Exception("Not authenticated")
        
    # Step 1: Upload the file
    url = f"{api.base_url}/csv/upload"
    headers = {"Authorization": f"Bearer {api.token}"}
    
    with open(file_path, 'rb') as file:
        files = {'file': file}
        upload_response = requests.post(url, headers=headers, files=files)
        
    upload_data = upload_response.json()
    
    if upload_data["status"] != "success":
        raise Exception(upload_data["error"]["message"])
        
    file_id = upload_data["data"]["file_id"]
    headers_list = upload_data["data"]["headers"]
    
    # Step 2: Map the headers
    mapping = {}
    for header in headers_list:
        # Simple auto-mapping
        if header.lower() == "smiles":
            mapping[header] = "smiles"
        elif header.lower() in ["mw", "molecular weight"]:
            mapping[header] = "molecular_weight"
        elif header.lower() == "logp":
            mapping[header] = "logp"
        else:
            mapping[header] = f"custom:{header.lower()}"
            
    map_url = f"{api.base_url}/csv/map"
    map_response = requests.post(
        map_url,
        headers={
            "Authorization": f"Bearer {api.token}",
            "Content-Type": "application/json"
        },
        json={
            "file_id": file_id,
            "mapping": mapping
        }
    )
    
    map_data = map_response.json()
    
    if map_data["status"] != "success":
        raise Exception(map_data["error"]["message"])
        
    job_id = map_data["data"]["job_id"]
    
    # Step 3: Poll for job status
    return poll_job_status(api, job_id)
    
def poll_job_status(api, job_id):
    """Poll for the status of a CSV import job"""
    if not api.token:
        raise Exception("Not authenticated")
        
    url = f"{api.base_url}/csv/status/{job_id}"
    headers = {"Authorization": f"Bearer {api.token}"}
    
    while True:
        response = requests.get(url, headers=headers)
        data = response.json()
        
        if data["status"] != "success":
            raise Exception(data["error"]["message"])
            
        job_status = data["data"]
        
        if job_status["status"] == "completed":
            return job_status
        elif job_status["status"] == "failed":
            raise Exception(f"Job failed: {job_status['message']}")
        else:
            # Still processing, wait and try again
            time.sleep(2)
```

### Experiment Submission Example

**JavaScript:**

```javascript
async function createAndSubmitExperiment(name, typeId, parameters, moleculeIds, croId) {
  try {
    const token = localStorage.getItem('token');
    
    if (!token) {
      throw new Error('Not authenticated');
    }
    
    // Step 1: Create the experiment
    const createResponse = await fetch('http://localhost/api/experiments', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        name,
        type_id: typeId,
        parameters,
        molecule_ids: moleculeIds
      }),
    });
    
    const createData = await createResponse.json();
    
    if (createData.status !== 'success') {
      throw new Error(createData.error.message);
    }
    
    const experimentId = createData.data.experiment.id;
    
    // Step 2: Add to queue
    const queueResponse = await fetch(`http://localhost/api/experiments/${experimentId}/queue`, {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json',
      },
    });
    
    const queueData = await queueResponse.json();
    
    if (queueData.status !== 'success') {
      throw new Error(queueData.error.message);
    }
    
    // Step 3: Submit to CRO
    const submitResponse = await fetch('http://localhost/api/submissions', {
      method: 'POST',
      headers: {
        'Authorization': `Bearer ${token}`,
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({
        experiment_id: experimentId,
        cro_id: croId,
        details: {
          notes: 'Please process as standard priority'
        }
      }),
    });
    
    const submitData = await submitResponse.json();
    
    if (submitData.status !== 'success') {
      throw new Error(submitData.error.message);
    }
    
    return submitData.data.submission;
  } catch (error) {
    console.error('Experiment submission failed:', error);
    throw error;
  }
}
```

**Python:**

```python
def create_and_submit_experiment(api, name, type_id, parameters, molecule_ids, cro_id):
    """Create and submit an experiment to a CRO"""
    if not api.token:
        raise Exception("Not authenticated")
        
    # Step 1: Create the experiment
    create_url = f"{api.base_url}/experiments"
    headers = {
        "Authorization": f"Bearer {api.token}",
        "Content-Type": "application/json"
    }
    
    create_payload = {
        "name": name,
        "type_id": type_id,
        "parameters": parameters,
        "molecule_ids": molecule_ids
    }
    
    create_response = requests.post(create_url, headers=headers, json=create_payload)
    create_data = create_response.json()
    
    if create_data["status"] != "success":
        raise Exception(create_data["error"]["message"])
        
    experiment_id = create_data["data"]["experiment"]["id"]
    
    # Step 2: Add to queue
    queue_url = f"{api.base_url}/experiments/{experiment_id}/queue"
    queue_response = requests.post(queue_url, headers=headers)
    queue_data = queue_response.json()
    
    if queue_data["status"] != "success":
        raise Exception(queue_data["error"]["message"])
        
    # Step 3: Submit to CRO
    submit_url = f"{api.base_url}/submissions"
    submit_payload = {
        "experiment_id": experiment_id,
        "cro_id": cro_id,
        "details": {
            "notes": "Please process as standard priority"
        }
    }
    
    submit_response = requests.post(submit_url, headers=headers, json=submit_payload)
    submit_data = submit_response.json()
    
    if submit_data["status"] != "success":
        raise Exception(submit_data["error"]["message"])
        
    return submit_data["data"]["submission"]
```