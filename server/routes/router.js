import express from 'express';
import compression from 'compression';
import webApiRoutes from './web.js';
import publicApiRoutes from './public.js';

const apiRouter = express.Router();

// parse json requests
apiRouter.use(express.json({ limit: '50mb' }));
// compress all responses
apiRouter.use(compression());
// register routes
apiRouter.use('/web', webApiRoutes);
apiRouter.use('/api', publicApiRoutes);

export { apiRouter };
