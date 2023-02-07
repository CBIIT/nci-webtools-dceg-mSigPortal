const express = require('express');
const logger = require('../services/logger');
const config = require('../config.json');

const { router: analysisRoutes } = require('../services/api/analysis');

const {
  router: visualizationRoutes,
} = require('../services/api/visualization/visualization');
const {
  router: userVisualizationRoutes,
} = require('../services/api/visualization/userVisualization');
const {
  router: explorationRoutes,
} = require('../services/api/exploration/exploration');
const {
  router: userExplorationRoutes,
} = require('../services/api/exploration/userExploration');
const {
  router: signatureRoutes,
} = require('../services/api/signature/signature');
const {
  router: associationRoutes,
} = require('../services/api/association/association');
const { router: etiologyRoutes } = require('../services/api/etiology/etiology');
const {
  router: publicationsRoutes,
} = require('../services/api/publications/publications');
const { router: patternRoutes } = require('../services/api/pattern/pattern');

const router = express.Router();

router.use(
  '/results',
  express.static(config.results.folder, {
    setHeaders: (res, path, stat) => {
      res.set('Cache-Control', 'max-age=0, must-revalidate');
    },
  })
);

router.use(analysisRoutes);
router.use(visualizationRoutes);
router.use(userVisualizationRoutes);
router.use(signatureRoutes);
router.use(explorationRoutes);
router.use(userExplorationRoutes);
router.use(associationRoutes);
router.use(etiologyRoutes);
router.use(publicationsRoutes);
router.use(patternRoutes);

router.get('/ping', (req, res) => res.send(true));

router.use((error, req, res, next) => {
  logger.error(error);
  res.status(500).json([error.message]);
});

module.exports = router;
