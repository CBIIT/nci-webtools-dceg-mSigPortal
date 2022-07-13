const { v4: uuidv4 } = require('uuid');
const { getExposureData } = require('../query');
const { groupBy } = require('lodash');

async function queryExposure(req, res, next) {
  try {
    const { study, strategy, cancer, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = {
      study,
      strategy,
      cancer,
      signatureSetName,
    };
    const columns = ['sample', 'signatureName', 'exposure'];
    const data = await getExposureData(connection, query, columns);
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

async function tmb(req, res, next) {
  try {
    const { study, strategy, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName };
    const columns = ['cancer', 'sample', 'exposure'];
    const data = await getExposureData(connection, query, columns);
    const tmb = calculateTmb(data, study);
    res.json(tmb);
  } catch (error) {
    next(error);
  }
}

function calculateTmb(data, study) {
  // Calculate the number of mutations per megabase for each study
  const genome = { PCAWG: 'GRCh37', TCGA: 'GRCh37' };
  const genomeSize = { GRCh37: 3101976562 / Math.pow(10, 6) };
  const burden = (exposure) => Math.log10(exposure / genomeSize[genome[study]]);

  const groupByCancer = groupBy(data, 'cancer');

  // calculate burden per cancer/sample
  const cancerBurden = Object.entries(groupByCancer)
    .map(([cancer, values]) => {
      const groupBySample = groupBy(values, 'sample');
      const samples = Object.entries(groupBySample).map(([sample, e]) => ({
        sample: sample,
        tmb: burden(e.reduce((sum, e) => e.exposure + sum, 0)),
      }));
      samples.sort((a, b) => a.tmb - b.tmb);

      // find median tmb
      const tmbs = samples.map((e) => e.tmb);
      const medianTmb =
        tmbs.length % 2 == 0
          ? (tmbs[tmbs.length / 2] + tmbs[tmbs.length / 2 - 1]) / 2
          : tmbs[Math.floor(tmbs.length / 2)];

      return { cancer, samples, medianTmb };
    })
    .sort((a, b) => a.medianTmb - b.medianTmb);

  return cancerBurden;
}

module.exports = {
  explorationOptions,
  queryExposure,
  tmb,
};
