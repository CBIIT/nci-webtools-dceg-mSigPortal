const { Router } = require('express');
const cors = require('cors');
const config = require('../config.json');
const apiSpec = require('../apiSpec.json');

const {
  router: visualizationRoutes,
} = require('../services/api/visualization/visualization');
const {
  router: explorationRoutes,
} = require('../services/api/exploration/exploration');
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

module.exports = router;
