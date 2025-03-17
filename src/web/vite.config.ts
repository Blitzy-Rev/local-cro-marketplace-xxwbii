import { defineConfig } from 'vite'; // vite ^4.3.5
import react from '@vitejs/plugin-react'; // @vitejs/plugin-react ^4.0.0
import svgr from 'vite-plugin-svgr'; // vite-plugin-svgr ^3.2.0
import tsconfigPaths from 'vite-tsconfig-paths'; // vite-tsconfig-paths ^4.2.0
import path from 'path'; // path built-in
import fs from 'fs'; // fs built-in

/**
 * Reads the package.json file to extract the current version
 * @returns The version string from package.json
 */
function getPackageVersion(): string {
  try {
    // First look for package.json in the web directory
    const packageJsonPath = path.resolve(__dirname, './package.json');
    if (fs.existsSync(packageJsonPath)) {
      const packageJson = fs.readFileSync(packageJsonPath, 'utf-8');
      const { version } = JSON.parse(packageJson);
      return version;
    }
    
    // Fallback to the project root package.json
    const rootPackageJsonPath = path.resolve(__dirname, '../../package.json');
    if (fs.existsSync(rootPackageJsonPath)) {
      const packageJson = fs.readFileSync(rootPackageJsonPath, 'utf-8');
      const { version } = JSON.parse(packageJson);
      return version;
    }
    
    return '0.0.0';
  } catch (error) {
    console.error('Error reading package.json:', error);
    return '0.0.0';
  }
}

// https://vitejs.dev/config/
export default defineConfig({
  // Project root directory
  root: __dirname,
  
  // Base public path for assets
  base: '/',
  
  // Directory for static assets
  publicDir: 'public',
  
  // Build configuration
  build: {
    // Output directory (relative to src/web)
    outDir: '../../dist',
    
    // Directory for chunked assets
    assetsDir: 'assets',
    
    // Generate sourcemaps for debugging
    sourcemap: true,
    
    // Minify output for production
    minify: true,
    
    // Clean output directory before build
    emptyOutDir: true,
    
    // Rollup specific options
    rollupOptions: {
      output: {
        // Manual chunk splitting for optimal loading
        manualChunks: {
          // Core React libraries
          vendor: ['react', 'react-dom', 'react-router-dom'],
          
          // UI component libraries
          ui: ['@mui/material', '@mui/icons-material', '@emotion/react', '@emotion/styled'],
          
          // State management
          state: ['@reduxjs/toolkit', 'react-redux'],
          
          // Data fetching
          data: ['react-query', 'axios']
        }
      }
    }
  },
  
  // Plugins
  plugins: [
    // React plugin with fast refresh
    react({ fastRefresh: true }),
    
    // SVG as React components
    svgr(),
    
    // TypeScript path resolution
    tsconfigPaths()
  ],
  
  // Path resolution
  resolve: {
    alias: {
      // Alias for src directory within the web frontend
      '@': path.resolve(__dirname, './src')
    }
  },
  
  // Development server configuration
  server: {
    // Port to run the dev server
    port: 3000,
    
    // Fail if port is already in use
    strictPort: true,
    
    // Listen on all network interfaces
    host: true,
    
    // Proxy configuration for API requests
    proxy: {
      // API endpoints
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        secure: false
      },
      
      // Storage endpoints
      '/storage': {
        target: 'http://localhost:9000',
        changeOrigin: true,
        secure: false
      },
      
      // WebSocket endpoints
      '/ws': {
        target: 'ws://localhost:8000',
        ws: true
      }
    }
  },
  
  // Define global constants
  define: {
    __APP_VERSION__: JSON.stringify(getPackageVersion())
  },
  
  // Dependency optimization
  optimizeDeps: {
    include: ['react', 'react-dom', 'react-router-dom', '@mui/material', 'lodash']
  },
  
  // CSS configuration
  css: {
    devSourcemap: true
  }
});