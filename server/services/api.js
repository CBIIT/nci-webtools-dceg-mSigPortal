const express = require('express');
const compression = require('compression');
const webApiRoutes = require('./routes/analysis');
const publicApiRoutes = require('./routes/public');

const apiRouter = express.Router();

// parse json requests
apiRouter.use(express.json({ limit: '50mb' }));

// compress all responses
apiRouter.use(compression());

// register routes
apiRouter.use('/web', webApiRoutes);
apiRouter.use('/api', publicApiRoutes);

module.exports = { apiRouter };
