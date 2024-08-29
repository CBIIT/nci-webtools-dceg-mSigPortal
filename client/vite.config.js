import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react-swc';

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [react()],
  server: {
    port: 3000,
    proxy: {
      '/api': 'http://localhost:8330',
      '/extraction': 'http://localhost:8332',
    },
  },
  define: {
    'process.env': process.env,
  },
});
