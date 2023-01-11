const { Router } = require('express');
const { isArguments } = require('lodash');
const {
  getSignatureData,
  getSignatureOptions,
  getSignatureSummary,
} = require('../../query');

async function querySignature(req, res, next) {
  try {
    const { limit, offset, scalarValue, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const rowMode = 'object';
    const distinct = true;
    const data = await getSignatureData(
      connection,
      query,
      columns,
      limit,
      offset,
      rowMode,
      distinct
    );
    const signatureName = req.query.signatureName;
    if (scalarValue) {
      const scalarArray = scalarValue.split(';');
      const signatureNameArray = signatureName.split(';');
      console.log(scalarArray);
      for (var i = 0; i < signatureNameArray.length; i++) {
        for (var j = 0; j < data.length; j++) {
          if (data[j].signatureName === signatureNameArray[i]) {
            data[j].contribution =
              data[j].contribution * parseInt(scalarArray[i]);
          }
        }
      }
    }
    res.json(data);
  } catch (error) {
    next(error);
  }
}

// query public exploration options for exploration tab
async function signatureOptions(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getSignatureOptions(
      connection,
      query,
      columns,
      limit,
      offset
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

// query signature summary data for catalog RS tab
async function signatureSummary(req, res, next) {
  try {
    const { limit, offset, ...query } = req.query;
    const connection = req.app.locals.connection;

    const columns = '*';
    const data = await getSignatureSummary(
      connection,
      query,
      columns,
      limit,
      offset
    );
    res.json(data);
  } catch (error) {
    next(error);
  }
}

const router = Router();

router.get('/mutational_signature', querySignature);
router.get('/mutational_signature_options', signatureOptions);
router.get('/mutational_signature_summary', signatureSummary);

module.exports = { router, querySignature };
