import { Router } from 'express';
import cors from 'cors';
import config from '../config.json' assert { type: 'json' };
import apiSpec from '../apiSpec.json' assert { type: 'json' };
import { router as visualizationRoutes } from '../services/api/visualization/visualization.js';
import { router as explorationRoutes } from '../services/api/exploration/exploration.js';
import { router as signatureRoutes } from '../services/api/signature/signature.js';
import { router as associationRoutes } from '../services/api/association/association.js';
import { router as etiologyRoutes } from '../services/api/etiology/etiology.js';
import { router as publicationsRoutes } from '../services/api/publications/publications.js';
import { router as patternRoutes } from '../services/api/pattern/pattern.js';

const router = Router();

router.use(cors());
router.get('/', (req, res) => {
  const url = config.email.baseUrl;
  const servers = [{ url }];
  res.json({ ...apiSpec, servers });
});
router.get('/ping', (req, res) => res.send(true));
router.use(visualizationRoutes);
router.use(signatureRoutes);
router.use(explorationRoutes);
router.use(associationRoutes);
router.use(etiologyRoutes);
router.use(publicationsRoutes);
router.use(patternRoutes);

export default router;
