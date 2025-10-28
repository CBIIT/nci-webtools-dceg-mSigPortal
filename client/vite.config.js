import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react-swc';

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  const BASE_URL = env.APP_PATH || '';
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
        [`${BASE_URL}/data`]: {
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
      __APP_VERSION__: JSON.stringify(env.VITE_APP_VERSION || 'dev'),
      __APP_LAST_UPDATE__: JSON.stringify(env.VITE_APP_LAST_UPDATE || 'Unknown'),
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
