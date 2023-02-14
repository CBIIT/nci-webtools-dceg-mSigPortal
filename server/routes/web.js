import express from 'express';
import logger from '../services/logger.js';
import config from '../config.json' assert { type: 'json' };
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
import { router as extrationRoutes } from '../services/api/extraction/extraction.js';

const router = express.Router();

router.use(
  '/data',
  express.static(config.folders.data, {
    setHeaders: (res, path, stat) => {
      res.set('Cache-Control', 'max-age=0, must-revalidate');
    },
  })
);

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
router.use(extrationRoutes);
router.get('/ping', (req, res) => res.send(true));
router.use((error, req, res, next) => {
  logger.error(error);
  res.status(500).json([error.message]);
});

export default router;
