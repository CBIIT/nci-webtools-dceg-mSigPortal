import express from 'express';
import Router from 'express-promise-router';
import compression from 'compression';
import cors from 'cors';
import {xss} from 'express-xss-sanitizer';
import {
  handleValidationErrors,
  logRequests,
  logErrors,
  logFiles,
} from '../services/middleware.js';
import apiSpec from '../apiSpec.json' with { type: 'json' };
import { router as general } from '../services/api/general.js';
import { router as visualizationRoutes } from '../services/api/visualization/visualization.js';
import { router as userVisualizationRoutes } from '../services/api/visualization/userVisualization.js';
import { router as explorationRoutes } from '../services/api/exploration/exploration.js';
import { router as userExplorationRoutes } from '../services/api/exploration/userExploration.js';
import { router as signatureRoutes } from '../services/api/signature/signature.js';
import { router as associationRoutes } from '../services/api/association/association.js';
import { router as etiologyRoutes } from '../services/api/etiology/etiology.js';
import { router as publicationsRoutes } from '../services/api/publications/publications.js';
import { router as patternRoutes } from '../services/api/pattern/pattern.js';
import { router as extractionRoutes } from '../services/api/extraction/extraction.js';
import { router as refittingRoutes } from '../services/api/refitting/refitting.js';

export function createApi(env) {
  // register middleware
  const router = Router();
  router.get('/ping', async (req, res) => res.json(true));
  router.use(express.json());
  router.use(compression());
  router.use(logRequests());
  router.use(cors());
  router.use(xss());
  // serve static files under /data
  router.use('/data', express.static(env.DATA_FOLDER));

  // serve swagger apispec
  router.get('/', (req, res) => {
    const url = process.env.APP_BASE_URL;
    const servers = [{ url }];
    res.json({ ...apiSpec, servers });
  });

  // register routes
  router.use(general);
  router.use(visualizationRoutes);
  router.use(userVisualizationRoutes);
  router.use(signatureRoutes);
  router.use(explorationRoutes);
  router.use(userExplorationRoutes);
  router.use(associationRoutes);
  router.use(etiologyRoutes);
  router.use(publicationsRoutes);
  router.use(patternRoutes);
  router.use(extractionRoutes);
  router.use(refittingRoutes);

  router.use(logErrors());
  return router;
}
