import express from 'express';
import Router from 'express-promise-router';
import compression from 'compression';
import webApiRoutes from './web.js';
import publicApiRoutes from './public.js';
import {
  handleValidationErrors,
  logRequests,
  logErrors,
  logFiles,
} from '../services/middleware.js';

export function createWebApi(env) {
  // register middleware
  const router = Router();
  router.use(express.json());
  router.use(compression());
  router.use(logRequests());

  // serve static files under /data
  router.use('/data', express.static(env.DATA_FOLDER));

  // register routes
  router.get('/ping', async (req, res) => res.json(true));
  router.use(webApiRoutes);

  router.use(logErrors());
  return router;
}

export function createPublicApi(env) {
  // register middleware
  const router = Router();
  router.use(express.json());
  router.use(logRequests());

  router.get('/ping', async (req, res) => res.json(true));
  router.use(publicApiRoutes);

  router.use(logErrors());
  return router;
}
