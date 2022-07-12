const { v4: uuidv4 } = require('uuid');
const { getSeqmatrixData, getSignatureData } = require('../query');

// query public seqmatrix data for visualization tab
async function visualizationOptions(req, res, next) {
  try {
    const connection = req.app.locals.connection;

    const query = {};
    const columns = ['study', 'cancer', 'strategy'];
    const data = await getSeqmatrixData(connection, query, columns);
    const projectID = uuidv4();

    res.json({ data, projectID });
  } catch (error) {
    next(error);
  }
}

async function visualizationSamples(req, res, next) {
  try {
    // const { study, cancer, strategy } = req.query;
    const connection = req.app.locals.connection;

    // const query = { study, strategy, cancer };
    const columns = ['profile', 'sample', 'profile', 'matrix'];
    const data = await getSeqmatrixData(connection, req.query, columns);
    const projectID = uuidv4();

    res.json({ data, projectID });
  } catch (error) {
    next(error);
  }
}

async function mutationalProfiles(req, res, next) {
  try {
    // const { study,strategy,cancer, sample, profile, matrix } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['mutationType', 'mutations'];
    const data = await getSeqmatrixData(connection, req.query, columns);
    res.json(data);
  } catch (error) {
    next(error);
  }
}

async function profilerSummary(req, res, next) {
  try {
    // const { study, cancer, strategy } = req.query;
    const connection = req.app.locals.connection;

    const columns = ['sample', 'profile', 'mutations'];
    const data = await getSeqmatrixData(connection, req.query, columns);
    res.json(data);
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
  visualizationOptions,
  visualizationSamples,
  mutationalProfiles,
  profilerSummary,
  querySeqmatrix,
  querySignature,
};
