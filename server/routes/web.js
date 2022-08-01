const express = require('express');
const logger = require('../services/logger');
const config = require('../config.json');

const { router: analyisRoutes } = require('../services/analysis/analysis');

const {
  router: visualizationRoutes,
} = require('../services/analysis/visualization/visualization');
const {
  router: explorationRoutes,
} = require('../services/analysis/exploration/exploration');
const {
  router: signatureRoutes,
} = require('../services/analysis/signature/signature');
const {
  router: associationRoutes,
} = require('../services/analysis/association/association');
const {
  router: etiologyRoutes,
} = require('../services/analysis/etiology/etiology');
const {
  router: publicationsRoutes,
} = require('../services/analysis/publications/publications');

const router = express.Router();

router.use(
  '/results',
  express.static(config.results.folder, {
    setHeaders: (res, path, stat) => {
      res.set('Cache-Control', 'max-age=0, must-revalidate');
    },
  })
);

router.use(analyisRoutes);
router.use(visualizationRoutes);
router.use(signatureRoutes);
router.use(explorationRoutes);
router.use(associationRoutes);
router.use(etiologyRoutes);
router.use(publicationsRoutes);

router.use((error, req, res, next) => {
  logger.error(error);
  res.status(500).json(error);
});

module.exports = router;
