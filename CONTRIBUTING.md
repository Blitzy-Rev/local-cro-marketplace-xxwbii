# Contributing to the Molecular Data Management and CRO Integration Platform

Thank you for your interest in contributing to the Molecular Data Management and CRO Integration Platform! This document provides guidelines and instructions for contributing to the project.

## Code of Conduct

By participating in this project, you agree to abide by our Code of Conduct. We expect all contributors to be respectful, inclusive, and considerate of others. We are committed to providing a welcoming and inspiring community for all.

## Getting Started

### Prerequisites

Before you begin, ensure you have the following installed:

- **Git**: Version control system
- **Docker**: Version 23.0 or higher
- **Docker Compose**: Version 2.17 or higher
- **Python**: Version 3.10 or higher (for backend development without Docker)
- **Poetry**: Version 1.4 or higher (for Python dependency management)
- **Node.js**: Version 16.0 or higher (for frontend development without Docker)
- **npm**: Version 8.0 or higher (for JavaScript package management)

### Setting Up the Development Environment

1. **Fork the repository** to your GitHub account

2. **Clone your fork** to your local machine:
   ```bash
   git clone https://github.com/your-username/molecular-platform.git
   cd molecular-platform
   ```

3. **Set up the development environment** using Docker Compose:
   ```bash
   cd infrastructure
   cp .env.example .env
   # Edit .env with your preferred configuration
   docker-compose -f docker-compose.yml -f docker-compose.dev.yml up -d
   ```

   Alternatively, you can set up the backend and frontend separately:

   **Backend Setup:**
   ```bash
   cd src/backend
   poetry install
   cp .env.example .env
   # Edit .env with your local configuration
   poetry run uvicorn main:app --reload
   ```

   **Frontend Setup:**
   ```bash
   cd src/web
   npm install
   cp .env.development .env
   # Edit .env with your local configuration
   npm run dev
   ```

4. **Create a feature branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

For more detailed setup instructions, please refer to the [Development Guide](docs/development.md).

## Development Workflow

### Project Structure

The repository is organized into the following main directories:

```
/
├── docs/                  # Documentation
├── infrastructure/        # Docker and deployment configuration
├── scripts/              # Utility scripts for deployment and maintenance
├── src/                  # Source code
│   ├── backend/          # Backend API and services
│   └── web/              # Frontend application
└── .github/              # GitHub workflows and templates
```

### Backend Development

The backend uses FastAPI with Python and follows these key principles:

1. **Create or modify API endpoints** in `src/backend/app/api/api_v1/endpoints/`
2. **Implement business logic** in `src/backend/app/services/`
3. **Define data models** in `src/backend/app/models/` and schemas in `src/backend/app/schemas/`
4. **Add database operations** in `src/backend/app/crud/`
5. **Run tests** to verify your changes:
   ```bash
   # Using Poetry
   cd src/backend
   poetry run pytest
   
   # Using Docker
   docker-compose exec backend pytest
   ```

### Frontend Development

The frontend uses React with TypeScript and follows these key principles:

1. **Create or modify components** in the appropriate feature directory in `src/web/src/features/`
2. **Update Redux state** as needed in `src/web/src/store/`
3. **Add API client functions** for new endpoints in `src/web/src/api/`
4. **Run tests** to verify your changes:
   ```bash
   # Using npm
   cd src/web
   npm test
   
   # Using Docker
   docker-compose exec frontend npm test
   ```

### Database Changes

When making changes to the database schema:

1. **Modify the models** in `src/backend/app/models/`
2. **Create a migration** using Alembic:
   ```bash
   # Using Poetry
   cd src/backend
   poetry run alembic revision --autogenerate -m "description of changes"
   
   # Using Docker
   docker-compose exec backend alembic revision --autogenerate -m "description of changes"
   ```
3. **Apply the migration**:
   ```bash
   # Using Poetry
   poetry run alembic upgrade head
   
   # Using Docker
   docker-compose exec backend alembic upgrade head
   ```

## Coding Standards

### Python Coding Standards

1. **Code Formatting**:
   - Use Black for code formatting
   - Use isort for import sorting
   - Maximum line length: 88 characters (Black default)

2. **Naming Conventions**:
   - Classes: CamelCase (e.g., `MoleculeProcessor`)
   - Functions and variables: snake_case (e.g., `calculate_properties`)
   - Constants: UPPER_CASE (e.g., `MAX_MOLECULES`)
   - Private methods/variables: leading underscore (e.g., `_validate_input`)

3. **Documentation**:
   - Use docstrings for all modules, classes, and functions
   - Follow Google-style docstring format
   - Include type hints for function parameters and return values

4. **Testing**:
   - Write unit tests for all functions and classes
   - Maintain minimum code coverage of 85%
   - Follow the Arrange-Act-Assert pattern

### TypeScript/JavaScript Coding Standards

1. **Code Formatting**:
   - Use Prettier for code formatting
   - Use ESLint for linting
   - Maximum line length: 100 characters

2. **Naming Conventions**:
   - Components: PascalCase (e.g., `MoleculeCard`)
   - Functions and variables: camelCase (e.g., `calculateProperties`)
   - Constants: UPPER_CASE (e.g., `MAX_MOLECULES`)
   - Interfaces/Types: PascalCase with descriptive names (e.g., `MoleculeProperties`)

3. **React Components**:
   - Use functional components with hooks
   - Define prop types using TypeScript interfaces
   - Keep components focused on a single responsibility

4. **Testing**:
   - Write tests for all components and hooks
   - Maintain minimum code coverage of 80%
   - Test behavior, not implementation

### Git Commit Standards

1. **Branch Naming**:
   - Feature branches: `feature/short-description`
   - Bug fixes: `fix/issue-description`
   - Documentation: `docs/topic-description`
   - Performance improvements: `perf/description`
   - Refactoring: `refactor/description`

2. **Commit Messages**:
   - Follow the Conventional Commits specification
   - Format: `type(scope): description`
   - Types: feat, fix, docs, style, refactor, perf, test, chore
   - Keep descriptions concise and descriptive

   Examples:
   ```
   feat(molecules): add property calculation for LogP
   fix(auth): resolve token refresh issue
   docs(readme): update development setup instructions
   refactor(api): simplify error handling
   ```

## Pull Request Process

1. **Ensure your code meets the coding standards** and passes all tests

2. **Update documentation** if necessary, including:
   - Code comments and docstrings
   - README updates if you've changed functionality
   - API documentation for new endpoints

3. **Run the test suite** to ensure all tests pass:
   ```bash
   # Backend tests
   cd src/backend
   poetry run pytest
   
   # Frontend tests
   cd src/web
   npm test
   ```

4. **Check code coverage** to ensure it meets the minimum requirements:
   ```bash
   # Backend coverage
   cd src/backend
   poetry run pytest --cov=app tests/
   
   # Frontend coverage
   cd src/web
   npm test -- --coverage
   ```

5. **Format your code** according to the project standards:
   ```bash
   # Backend formatting
   cd src/backend
   poetry run black app tests
   poetry run isort app tests
   
   # Frontend formatting
   cd src/web
   npm run format
   npm run lint
   ```

6. **Commit your changes** with descriptive commit messages following the Conventional Commits specification

7. **Push your changes** to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

8. **Create a pull request** from your fork to the main repository

9. **Fill out the pull request template** with all required information

10. **Request reviews** from appropriate team members

11. **Address review comments** promptly and thoroughly

12. **Update your PR** if needed based on feedback

Your pull request will be reviewed by the maintainers, who may request changes or provide feedback. Once approved, your changes will be merged into the main codebase.

## Issue Reporting

### Bug Reports

If you find a bug in the project, please report it by creating an issue on GitHub. When reporting bugs, please include:

1. **A clear and descriptive title**
2. **Steps to reproduce the issue**
3. **Expected behavior**
4. **Actual behavior**
5. **Environment details** (OS, browser, etc.)
6. **Screenshots or logs** if applicable

### Feature Requests

If you have an idea for a new feature or enhancement, please create an issue using the feature request template. When submitting feature requests, please include:

1. **A clear and descriptive title**
2. **A detailed description of the feature**
3. **The problem this feature would solve**
4. **How this feature would benefit users**
5. **Any alternative solutions you've considered**

### Security Issues

If you discover a security vulnerability, please do NOT open an issue. Instead, please email [security@example.com] with details about the vulnerability. We take security issues very seriously and will address them promptly.

## Code Review Guidelines

When reviewing code, please follow these guidelines:

1. **Be respectful and constructive**:
   - Focus on the code, not the person
   - Explain why changes are suggested
   - Provide examples when possible

2. **Check for**:
   - Functionality: Does the code work as intended?
   - Code quality: Is the code clean, readable, and maintainable?
   - Performance: Are there any performance concerns?
   - Security: Are there any security vulnerabilities?
   - Tests: Are there appropriate tests?
   - Documentation: Is the code properly documented?

3. **Use GitHub's review features**:
   - Comment on specific lines
   - Suggest changes
   - Approve or request changes

4. **Be timely** with reviews to avoid blocking progress

5. **Follow up** on addressed comments

## License

By contributing to this project, you agree that your contributions will be licensed under the same [MIT License](LICENSE) that covers the project.

## Questions and Support

If you have questions about contributing or need help with the development environment, please:

1. Check the [documentation](docs/)
2. Look for existing issues that might address your question
3. Create a new issue with the question tag if you can't find an answer

Thank you for contributing to the Molecular Data Management and CRO Integration Platform!