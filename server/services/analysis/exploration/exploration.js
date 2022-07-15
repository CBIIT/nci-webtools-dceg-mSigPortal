const { Router } = require('express');
const { v4: uuidv4 } = require('uuid');
const { getExposureData } = require('../../query');
const { calculateTmb, calculateTmbSignature } = require('./tmb');

function alphaNumericSort(array) {
  return array.sort((a, b) => {
    return a.localeCompare(b, undefined, {
      numeric: true,
      sensitivity: 'base',
    });
  });
}
async function queryExposure(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['sample', 'signatureName', 'exposure'];
    const data = await getExposureData(connection, query, columns, limit);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

// query public exploration options for exploration tab
async function explorationOptions(req, res, next) {
  try {
    const { study, strategy, signatureSetName, cancer } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName, cancer };
    const columns = ['study', 'strategy', 'cancer', 'signatureSetName'];
    const data = await getExposureData(connection, query, columns);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function explorationSamples(req, res, next) {
  try {
    const { study, strategy, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName };
    const columns = ['sample', 'signatureName'];
    const data = await getExposureData(connection, query, columns);
    const samples = alphaNumericSort([...new Set(data.map((e) => e.sample))]);
    const signatures = alphaNumericSort([
      ...new Set(data.map((e) => e.signatureName)),
    ]);
    res.json({ samples, signatures });
  } catch (error) {
    next(error);
  }
}

async function tmb(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['cancer', 'sample', 'exposure'];
    const data = await getExposureData(connection, query, columns, limit);
    const tmb = calculateTmb(data, query.study);
    res.json(tmb);
  } catch (error) {
    next(error);
  }
}

async function tmbSignature(req, res, next) {
  try {
    const { limit, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['sample', 'signatureName', 'exposure'];
    const data = await getExposureData(connection, query, columns, limit);
    const tmbSignature = calculateTmbSignature(
      data.filter((e) => e.exposure),
      query.study
    );
    res.json(tmbSignature);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/exposure', queryExposure);
router.get('/explorationOptions', explorationOptions);
router.get('/explorationSamples', explorationSamples);
router.get('/tmb', tmb);
router.get('/tmbSignature', tmbSignature);

module.exports = router;
