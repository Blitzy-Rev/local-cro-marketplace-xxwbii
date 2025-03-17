-- PostgreSQL initialization script for the Molecular Data Management and CRO Integration Platform
-- This script creates the database schema, tables, indexes, and initial data

-- ================================================================
-- EXTENSIONS
-- ================================================================
CREATE EXTENSION IF NOT EXISTS "uuid-ossp"; -- For UUID generation
CREATE EXTENSION IF NOT EXISTS pg_trgm;     -- For text search capabilities
CREATE EXTENSION IF NOT EXISTS btree_gin;   -- For GIN indexing

-- ================================================================
-- ENUM TYPES
-- ================================================================
-- User roles
CREATE TYPE user_role AS ENUM ('PHARMA', 'CRO', 'ADMIN');

-- User account status
CREATE TYPE user_status AS ENUM ('PENDING', 'ACTIVE', 'INACTIVE', 'LOCKED');

-- Experiment status
CREATE TYPE experiment_status AS ENUM (
    'DRAFT',           -- Initial creation, not yet queued
    'QUEUED',          -- Added to queue, ready for submission
    'SUBMITTED',       -- Submitted to CRO
    'REJECTED',        -- Rejected by CRO
    'QUOTE_PENDING',   -- CRO accepted, waiting for quote approval
    'QUOTE_REJECTED',  -- Quote rejected by Pharma
    'IN_PROGRESS',     -- Approved and in progress at CRO
    'RESULTS_PENDING', -- CRO completed work
    'RESULTS_AVAILABLE', -- CRO uploaded results
    'RESULTS_REJECTED',  -- Results rejected by Pharma
    'COMPLETED',       -- Results reviewed and accepted
    'CANCELLED'        -- Cancelled at any stage
);

-- Submission status
CREATE TYPE submission_status AS ENUM (
    'PENDING',         -- Initial submission
    'REJECTED',        -- Rejected by CRO
    'QUOTE_PROVIDED',  -- CRO provided pricing
    'QUOTE_REJECTED',  -- Quote rejected by pharma
    'APPROVED',        -- Quote approved, work can begin
    'IN_PROGRESS',     -- Work started
    'COMPLETED',       -- Work completed
    'CANCELLED'        -- Cancelled by either party
);

-- Result status
CREATE TYPE result_status AS ENUM (
    'PENDING',         -- Results not yet provided
    'UPLOADED',        -- Results uploaded by CRO
    'APPROVED',        -- Results approved by pharma
    'REJECTED'         -- Results rejected by pharma
);

-- Notification types
CREATE TYPE notification_type AS ENUM (
    'EXPERIMENT_STATUS_CHANGE',
    'SUBMISSION_CREATED',
    'QUOTE_PROVIDED',
    'RESULTS_UPLOADED',
    'SYSTEM_ALERT',
    'USER_MENTION'
);

-- ================================================================
-- TABLES
-- ================================================================

-- Users table
CREATE TABLE users (
    id SERIAL PRIMARY KEY,
    email VARCHAR(255) UNIQUE NOT NULL,
    password_hash VARCHAR(255) NOT NULL,
    role user_role NOT NULL,
    status user_status NOT NULL DEFAULT 'PENDING',
    is_active BOOLEAN NOT NULL DEFAULT FALSE,
    email_verified BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP,
    last_login TIMESTAMP,
    password_history TEXT[] -- Store previous password hashes
);

-- Molecules table
CREATE TABLE molecules (
    id SERIAL PRIMARY KEY,
    smiles VARCHAR(4000) UNIQUE NOT NULL, -- SMILES string representation
    created_by INTEGER NOT NULL REFERENCES users(id),
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP,
    flag_status VARCHAR(50) -- Optional flag for prioritization
);

-- Molecule properties table
CREATE TABLE molecule_properties (
    molecule_id INTEGER NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
    property_name VARCHAR(100) NOT NULL,
    property_value NUMERIC,
    property_unit VARCHAR(50),
    is_calculated BOOLEAN NOT NULL DEFAULT FALSE,
    PRIMARY KEY (molecule_id, property_name)
);

-- Libraries table
CREATE TABLE libraries (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    description TEXT,
    created_by INTEGER NOT NULL REFERENCES users(id),
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP
);

-- Library molecules junction table
CREATE TABLE library_molecules (
    library_id INTEGER NOT NULL REFERENCES libraries(id) ON DELETE CASCADE,
    molecule_id INTEGER NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
    added_at TIMESTAMP NOT NULL DEFAULT NOW(),
    PRIMARY KEY (library_id, molecule_id)
);

-- Experiment types table
CREATE TABLE experiment_types (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL UNIQUE,
    description TEXT,
    category VARCHAR(100)
);

-- Experiments table
CREATE TABLE experiments (
    id SERIAL PRIMARY KEY,
    name VARCHAR(255) NOT NULL,
    type_id INTEGER NOT NULL REFERENCES experiment_types(id),
    status experiment_status NOT NULL DEFAULT 'DRAFT',
    created_by INTEGER NOT NULL REFERENCES users(id),
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP,
    description TEXT
);

-- Experiment parameters table
CREATE TABLE experiment_parameters (
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    parameter_name VARCHAR(100) NOT NULL,
    parameter_value TEXT,
    PRIMARY KEY (experiment_id, parameter_name)
);

-- Experiment molecules junction table
CREATE TABLE experiment_molecules (
    experiment_id INTEGER NOT NULL REFERENCES experiments(id) ON DELETE CASCADE,
    molecule_id INTEGER NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
    added_at TIMESTAMP NOT NULL DEFAULT NOW(),
    PRIMARY KEY (experiment_id, molecule_id)
);

-- Submissions table
CREATE TABLE submissions (
    id SERIAL PRIMARY KEY,
    experiment_id INTEGER NOT NULL REFERENCES experiments(id),
    cro_id INTEGER NOT NULL REFERENCES users(id), -- CRO user assigned to this submission
    status submission_status NOT NULL DEFAULT 'PENDING',
    submitted_at TIMESTAMP NOT NULL DEFAULT NOW(),
    updated_at TIMESTAMP
);

-- Submission details table
CREATE TABLE submission_details (
    submission_id INTEGER NOT NULL REFERENCES submissions(id) ON DELETE CASCADE,
    detail_name VARCHAR(100) NOT NULL,
    detail_value TEXT,
    PRIMARY KEY (submission_id, detail_name)
);

-- Results table
CREATE TABLE results (
    id SERIAL PRIMARY KEY,
    submission_id INTEGER NOT NULL REFERENCES submissions(id),
    status result_status NOT NULL DEFAULT 'PENDING',
    uploaded_at TIMESTAMP,
    approved_at TIMESTAMP
);

-- Result files table
CREATE TABLE result_files (
    id SERIAL PRIMARY KEY,
    result_id INTEGER NOT NULL REFERENCES results(id) ON DELETE CASCADE,
    file_name VARCHAR(255) NOT NULL,
    file_path VARCHAR(1000) NOT NULL,
    file_size INTEGER,
    file_type VARCHAR(50),
    uploaded_at TIMESTAMP NOT NULL DEFAULT NOW()
);

-- Result data table
CREATE TABLE result_data (
    result_id INTEGER NOT NULL REFERENCES results(id) ON DELETE CASCADE,
    molecule_id INTEGER NOT NULL REFERENCES molecules(id) ON DELETE CASCADE,
    data_name VARCHAR(100) NOT NULL,
    data_value NUMERIC,
    data_unit VARCHAR(50),
    PRIMARY KEY (result_id, molecule_id, data_name)
);

-- Notifications table
CREATE TABLE notifications (
    id SERIAL PRIMARY KEY,
    user_id INTEGER NOT NULL REFERENCES users(id) ON DELETE CASCADE,
    type notification_type NOT NULL,
    message TEXT NOT NULL,
    read_status BOOLEAN NOT NULL DEFAULT FALSE,
    created_at TIMESTAMP NOT NULL DEFAULT NOW(),
    read_at TIMESTAMP
);

-- ================================================================
-- INDEXES
-- ================================================================

-- User indexes
CREATE INDEX idx_users_email ON users(email);

-- Molecule indexes
CREATE INDEX idx_molecules_smiles ON molecules(smiles);
CREATE INDEX idx_molecules_created_by ON molecules(created_by);

-- Molecule properties indexes
CREATE INDEX idx_mol_prop_name_value ON molecule_properties(molecule_id, property_name, property_value);
CREATE INDEX idx_property_range ON molecule_properties(property_name, property_value);

-- Library indexes
CREATE INDEX idx_libraries_user ON libraries(created_by);

-- Experiment indexes
CREATE INDEX idx_experiments_status ON experiments(status);
CREATE INDEX idx_experiments_user ON experiments(created_by);

-- Submission indexes
CREATE INDEX idx_submissions_cro ON submissions(cro_id);
CREATE INDEX idx_submissions_status ON submissions(status);

-- Notification indexes
CREATE INDEX idx_notifications_user ON notifications(user_id);
CREATE INDEX idx_notifications_unread ON notifications(user_id, read_status);

-- ================================================================
-- TRIGGERS AND FUNCTIONS
-- ================================================================

-- Function to automatically update updated_at timestamp
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
   NEW.updated_at = NOW();
   RETURN NEW;
END;
$$ LANGUAGE plpgsql;

-- Add the trigger to relevant tables
CREATE TRIGGER update_users_updated_at
BEFORE UPDATE ON users
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_molecules_updated_at
BEFORE UPDATE ON molecules
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_libraries_updated_at
BEFORE UPDATE ON libraries
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_experiments_updated_at
BEFORE UPDATE ON experiments
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

CREATE TRIGGER update_submissions_updated_at
BEFORE UPDATE ON submissions
FOR EACH ROW
EXECUTE FUNCTION update_updated_at_column();

-- ================================================================
-- INITIAL DATA
-- ================================================================

-- Create admin user with password 'admin'
-- Note: In production, this should be changed immediately
-- The hash used here is for the password 'admin' using bcrypt
INSERT INTO users (
    email,
    password_hash,
    role,
    status,
    is_active,
    email_verified
) VALUES (
    'admin@example.com',
    '$2b$12$EixZaYVK1fsbw1ZfbX3OXePaWxn96p36WQoeG6Lruj3vjPGga31lW', -- 'admin'
    'ADMIN',
    'ACTIVE',
    TRUE,
    TRUE
);

-- Insert experiment types
INSERT INTO experiment_types (name, description, category)
VALUES 
    ('Binding Assay', 'Measures binding affinity between molecules and targets', 'Binding'),
    ('ADME Panel', 'Absorption, Distribution, Metabolism, and Excretion assessment', 'ADME'),
    ('Toxicity Assay', 'Evaluates potential toxic effects of compounds', 'Toxicity'),
    ('Solubility Test', 'Measures compound solubility in various solvents', 'Physicochemical');