const { Router } = require('express');
const cors = require('cors');
const config = require('../../config.json');
const apiSpec = require('../../apiSpec.json');

const {
  visualizationOptions,
  visualizationSamples,
  mutationalProfiles,
  profilerSummary,
  querySeqmatrix,
  querySignature,
} = require('../analysis/visualization');

const {
  queryExposure,
  explorationOptions,
  explorationTmbData,
} = require('../analysis/exploration');

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

router.get('/visualizationOptions', visualizationOptions);

router.get('/visualizationSamples', visualizationSamples);

router.get('/mutationalProfiles', mutationalProfiles);

router.get('/profilerSummary', profilerSummary);

module.exports = router;
