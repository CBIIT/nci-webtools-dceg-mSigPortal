const { Router } = require('express');
const cors = require('cors');
const config = require('../../config.json');
const apiSpec = require('../../apiSpec.json');
const {
  visualizationData,
  querySeqmatrix,
  queryExposure,
  querySignature,
} = require('../apiQuery');

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

router.get('/visualizationData', visualizationData);

module.exports = router;
