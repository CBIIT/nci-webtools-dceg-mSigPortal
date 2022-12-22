const express = require('express');
const compression = require('compression');
// const cors = require('cors');
const extraction = require('./extraction');
const logger = require('../services/logger');

const apiRouter = express.Router();

// parse json requests
apiRouter.use(express.json({ limit: '50mb' }));

// compress all responses
apiRouter.use(compression());

// // enable cors
// apiRouter.use(cors());

// register routes
apiRouter.use('/extraction', extraction);

apiRouter.use((error, req, res, next) => {
  logger.error(error);
  res.status(500).json([error.message]);
});

module.exports = apiRouter;
