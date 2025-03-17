/// <reference types="vite/client" />

/**
 * Interface that extends Vite's ImportMetaEnv to add custom environment variables
 */
interface ImportMetaEnv {
  /**
   * Base URL for API requests
   */
  readonly VITE_API_BASE_URL: string;
  
  /**
   * Application title
   */
  readonly VITE_APP_TITLE: string;
  
  /**
   * Application version
   */
  readonly VITE_APP_VERSION: string;
  
  /**
   * Whether to use mock API responses
   */
  readonly VITE_ENABLE_MOCK_API: boolean;
  
  /**
   * Prefix for local storage keys
   */
  readonly VITE_STORAGE_PREFIX: string;
}

/**
 * Interface that extends Vite's ImportMeta to include the env property
 */
interface ImportMeta {
  readonly env: ImportMetaEnv;
}

/**
 * Global variable defined in vite.config.ts for application version
 */
declare const __APP_VERSION__: string;

/**
 * Declaration for SVG imports as React components
 */
declare module '*.svg' {
  import React from 'react';
  const SVGComponent: React.FC<React.SVGProps<SVGSVGElement>>;
  export default SVGComponent;
}