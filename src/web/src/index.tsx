import React from 'react'; // React v18.2+
import ReactDOM from 'react-dom/client'; // react-dom v18.2+
import App from './App'; // Import the main App component

// LD1: Define the root element where the React app will be mounted
const rootElement = document.getElementById('root');

// LD1: Check if the root element exists to prevent errors
if (!rootElement) {
  console.error('Root element with id "root" not found in the document.');
} else {
  // LD1: Create a React root using createRoot API for React 18
  const root = ReactDOM.createRoot(rootElement);

  // LD1: Render the App component into the root element
  root.render(
    <React.StrictMode>
      <App />
    </React.StrictMode>
  );
}