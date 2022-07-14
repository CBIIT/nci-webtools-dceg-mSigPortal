const { v4: uuidv4 } = require('uuid');
const { getExposureData } = require('../../query');
const { calculateTmb } = require('./tmb');

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

module.exports = {
  explorationOptions,
  queryExposure,
  tmb,
};
