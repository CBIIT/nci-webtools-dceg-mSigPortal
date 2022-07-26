const { Router } = require('express');
const { v4: uuidv4 } = require('uuid');
const {
  getExposureData,
  getSeqmatrixData,
  getSignatureData,
} = require('../../query');
const { addBurden } = require('./burden');

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

    const columns = '*';
    const data = await getExposureData(connection, query, columns, limit);
    res.json(addBurden(data));
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
    const signatureNames = alphaNumericSort([
      ...new Set(data.map((e) => e.signatureName)),
    ]);
    res.json({ samples, signatureNames });
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/exposure', queryExposure);
router.get('/explorationOptions', explorationOptions);
router.get('/explorationSamples', explorationSamples);

module.exports = { router, queryExposure };
