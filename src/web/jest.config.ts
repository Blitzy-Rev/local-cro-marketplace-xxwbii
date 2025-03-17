import type { JestConfigWithTsJest } from 'ts-jest'; // v29.1.0

/**
 * Jest configuration for the frontend application of the Molecular Data Management and CRO Integration Platform.
 * This configuration defines testing environment settings, file transformations, module mappings, and coverage reporting.
 */
const config: JestConfigWithTsJest = {
  // Use ts-jest preset for TypeScript support
  preset: 'ts-jest',
  
  // Use jsdom to simulate browser environment for React component testing
  testEnvironment: 'jsdom',
  
  // Setup files to run after Jest is initialized
  setupFilesAfterEnv: ['<rootDir>/src/__tests__/setup.ts'],
  
  // Configure transformers for TypeScript files
  transform: {
    '^.+\\.tsx?$': 'ts-jest',
  },
  
  // Module name mappers for handling various file types and path aliases
  moduleNameMapper: {
    // Handle CSS imports (mock them as they're not needed for tests)
    '\\.(css|less|scss|sass)$': 'identity-obj-proxy',
    
    // Handle image and other file imports
    '\\.(jpg|jpeg|png|gif|eot|otf|webp|svg|ttf|woff|woff2|mp4|webm|wav|mp3|m4a|aac|oga)$':
      '<rootDir>/src/__tests__/mocks/fileMock.ts',
    
    // Support @ alias for src directory
    '^@/(.*)$': '<rootDir>/src/$1',
  },
  
  // Regex pattern for identifying test files
  testRegex: '(/__tests__/.*|(\\.|/)(test|spec))\\.(jsx?|tsx?)$',
  
  // File extensions Jest should recognize
  moduleFileExtensions: ['ts', 'tsx', 'js', 'jsx', 'json', 'node'],
  
  // Files to include in coverage reports
  collectCoverageFrom: [
    'src/**/*.{ts,tsx}',
    '!src/**/*.d.ts',
    '!src/__tests__/**/*',
    '!src/types/**/*',
    '!src/vite-env.d.ts',
    '!src/index.tsx',
    '!src/App.tsx',
  ],
  
  // Minimum coverage thresholds to enforce
  coverageThreshold: {
    global: {
      branches: 80,
      functions: 85,
      lines: 85,
      statements: 85,
    },
  },
  
  // Patterns for paths to ignore during testing
  testPathIgnorePatterns: ['/node_modules/', '/dist/'],
  
  // Plugins for watch mode to improve developer experience
  watchPlugins: [
    'jest-watch-typeahead/filename',
    'jest-watch-typeahead/testname',
  ],
  
  // Global configuration settings for ts-jest
  globals: {
    'ts-jest': {
      tsconfig: 'tsconfig.json',
      isolatedModules: true,
    },
  },
};

export default config;