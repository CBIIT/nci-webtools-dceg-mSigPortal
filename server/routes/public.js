const { Router } = require('express');
const cors = require('cors');
const config = require('../config.json');
const apiSpec = require('../apiSpec.json');

const {
  querySeqmatrix,
} = require('../services/analysis/visualization/visualization');
const {
  queryExposure,
} = require('../services/analysis/exploration/exploration');
const { querySignature } = require('../services/analysis/signature/signature');
const {
  queryAssociation,
} = require('../services/analysis/association/association');

const router = Router();

router.use(cors());

router.get('/', (req, res) => {
  const url = config.email.baseUrl;
  const servers = [{ url }];
  res.json({ ...apiSpec, servers });
});

router.get('/ping', (req, res) => res.send(true));

router.get('/seqmatrix', querySeqmatrix);
router.get('/signature', querySignature);
router.get('/exposure', queryExposure);
router.get('/association', queryAssociation);

module.exports = router;
