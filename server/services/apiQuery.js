const { v4: uuidv4 } = require('uuid');
const {
  getAssociationData,
  getSeqmatrixData,
  getExposureData,
  getSignatureData,
} = require('./query');

// query public seqmatrix data for visualization tab
async function visualizationData(req, res, next) {
  try {
    const { study, cancer, strategy } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, cancer };
    const columns = ['profile', 'sample', 'profile', 'matrix'];
    const data = await getSeqmatrixData(connection, query, columns);
    const projectID = uuidv4();

    res.json({ data, projectID });
  } catch (error) {
    next(error);
  }
}

async function querySeqmatrix(req, res, next) {
  try {
    const { type, ...query } = req.query;
    const connection = req.app.locals.connection;

    if (!type) throw "Seqmatrix API: Missing 'type' parameter";

    let columns = [];
    if (type == 'studies') columns = ['study', 'cancer', 'strategy'];
    if (type == 'samples') columns = ['profile', 'sample', 'profile', 'matrix'];
    if (type == 'mutationalProfiles') columns = ['mutationType', 'mutations'];
    if (type == 'profilerSummary') columns = ['sample', 'profile', 'mutations'];

    const data = await getSeqmatrixData(connection, query, columns);
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

// query public exploration data for exploration tab
async function explorationTmbData(req, res, next) {
  try {
    const { study, strategy, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = { study, strategy, signatureSetName };
    const columns = ['cancer', 'sample', 'exposure'];
    const data = await getExposureData(connection, query, columns);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

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

async function querySignature(req, res, next) {
  try {
    const { profile, matrix, signatureSetName } = req.query;
    const connection = req.app.locals.connection;

    const query = { profile, matrix, signatureSetName };
    const columns = !signatureSetName
      ? ['signatureSetName']
      : [
          'strandInfo',
          'strand',
          'signatureName',
          'mutationType',
          'contribution',
        ];
    const data = await getSignatureData(connection, query, columns);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

module.exports = {
  visualizationData,
  explorationOptions,
  explorationTmbData,
  querySeqmatrix,
  queryExposure,
  querySignature,
};
