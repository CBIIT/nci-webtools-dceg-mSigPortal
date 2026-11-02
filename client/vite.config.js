import { defineConfig, loadEnv } from 'vite';
import react from '@vitejs/plugin-react-swc';
import path from 'path';
import { fileURLToPath } from 'url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));

// https://vitejs.dev/config/
export default defineConfig(({ mode }) => {
  const env = loadEnv(mode, process.cwd(), '');
  const BASE_URL = env.APP_PATH || '';
  return {
    plugins: [
      react(),
      {
        name: 'fix-gz-mime-type',
        configureServer(server) {
          server.middlewares.use((req, res, next) => {
            if (req.url?.endsWith('.gz')) {
              res.setHeader('Content-Type', 'application/gzip');
              res.setHeader('Content-Encoding', 'identity');
            }
            next();
          });
        },
      },
    ],
    resolve: {
      alias: {
        '@': path.resolve(__dirname, 'src'),
      },
    },
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

        [`/api`]: {
          target: 'http://localhost:8330',
          rewrite: (path) => path.replace(BASE_URL, ''),
        },
        [`/data`]: {
          target: 'http://localhost:8330',
          rewrite: (path) => path.replace(BASE_URL, ''),
        },
      },
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
