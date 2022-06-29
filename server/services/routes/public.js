const { Router } = require('express');
const cors = require('cors');
const config = require('../../config.json');
const apiSpec = require('../../apiSpec.json');
const {
  querySeqmatrix,
  queryExposure,
  querySignature,
} = require('../analysis');

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

module.exports = router;
