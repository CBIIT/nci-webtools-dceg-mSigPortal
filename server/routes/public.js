const { Router } = require('express');
const cors = require('cors');
const config = require('../config.json');
const apiSpec = require('../apiSpec.json');

const visualizationRoutes = require('../services/analysis/visualization/visualization');
const explorationRoutes = require('../services/analysis/exploration/exploration');

const router = Router();

router.use(cors());

router.get('/', (req, res) => {
  const url = config.email.baseUrl;
  const servers = [{ url }];
  res.json({ ...apiSpec, servers });
});

router.get('/ping', (req, res) => res.send(true));

// visualization
router.use(visualizationRoutes);

// exploration
router.use(explorationRoutes);

module.exports = router;
