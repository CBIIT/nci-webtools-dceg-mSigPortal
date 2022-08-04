const { Router } = require('express');
const cors = require('cors');
const config = require('../config.json');
const apiSpec = require('../apiSpec.json');

const {
  querySeqmatrix,
} = require('../services/api/visualization/visualization');
const {
  queryExposure,
} = require('../services/api/exploration/exploration');
const { querySignature } = require('../services/api/signature/signature');
const {
  queryAssociation,
} = require('../services/api/association/association');
const {
  queryEtiology,
} = require('../services/api/etiology/etiology');
const {
  queryPublications,
} = require('../services/api/publications/publications');

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
router.use('/etiology', queryEtiology);
router.use('/publications', queryPublications);

module.exports = router;
