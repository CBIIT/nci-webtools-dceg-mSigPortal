const { Router } = require('express');
const cors = require('cors');
const config = require('../config.json');
const apiSpec = require('../apiSpec.json');

const {
  visualizationOptions,
  visualizationSamples,
  mutationalProfiles,
  profilerSummary,
  querySeqmatrix,
  querySignature,
} = require('../services/analysis/visualization');

const {
  queryExposure,
  explorationOptions,
  tmb,
} = require('../services/analysis/exploration');

const router = Router();

router.use(cors());

router.get('/', (req, res) => {
  const url = config.email.baseUrl;
  const servers = [{ url }];
  res.json({ ...apiSpec, servers });
});

router.get('/ping', (req, res) => res.send(true));

router.get('/seqmatrix', querySeqmatrix);

router.get('/exposure', queryExposure);

router.get('/signature', querySignature);

// visualization
router.get('/visualizationOptions', visualizationOptions);

router.get('/visualizationSamples', visualizationSamples);

router.get('/mutationalProfiles', mutationalProfiles);

router.get('/profilerSummary', profilerSummary);

// exploration
router.get('/explorationOptions', explorationOptions);

router.get('/tmb', tmb);

module.exports = router;
