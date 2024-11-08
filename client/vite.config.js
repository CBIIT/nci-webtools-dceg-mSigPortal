import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react-swc';

// https://vitejs.dev/config/
const BASE_URL = '/mutational-signatures';
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  return {
    plugins: [react()],
    base: BASE_URL,
    server: {
      port: 3000,
      proxy: {
        [`${BASE_URL}/api`]: {
          target: 'http://localhost:8330',
          rewrite: (path) => path.replace(BASE_URL, ''),
        },
        [`${BASE_URL}/extraction`]: {
          target: 'http://localhost:8332',
          rewrite: (path) => path.replace(BASE_URL, ''),
        },
      },
    },
    define: {
      'process.env': env,
    },
    build: {
      rollupOptions: {
        input: '/index.html', 
      },
    },
    optimizeDeps: {
      include: ['react', 'react-dom'], 
    },
  };
});
