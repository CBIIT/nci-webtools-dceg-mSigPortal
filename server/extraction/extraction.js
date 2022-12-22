const { Router } = require('express');
const { randomUUID } = require('crypto');
const { validate } = require('uuid');
const formidable = require('formidable');
const fs = require('fs-extra');
const path = require('path');
const { unparse } = require('papaparse');
const AWS = require('aws-sdk');
const tar = require('tar');
const { groupBy } = require('lodash');
const logger = require('../services/logger');
const config = require('../config.json');
const { getSignatureData } = require('../services/query');

function upload(req, res, next) {
  const id = randomUUID();
  const form = formidable({
    uploadDir: path.join(config.results.folder, id),
    multiples: true,
  });

  // create upload directory
  fs.mkdirSync(form.uploadDir);

  form
    .on('fileBegin', (field, file) => {
      uploadPath = path.join(form.uploadDir, file.name);
      file.path = uploadPath;
    })
    .on('error', (error) => {
      logger.error('error at /upload');
      next(error);
    })
    .on('end', () => {
      res.json({ id });
    });

  form.parse(req);
}

async function submit(req, res, next) {
  try {
    const { args, signatureQuery, id, email } = req.body;
    if (!validate(id)) throw Error('Invalid ID');

    const workPath = path.resolve(config.results.folder, id);

    // query signature data
    const connection = req.app.locals.connection;
    const columns = ['signatureName', 'mutationType', 'contribution'];
    const limit = false;
    const signatureData = await getSignatureData(
      connection,
      signatureQuery,
      columns,
      limit
    );

    // transform data into format accepted by SigProfiler
    const groupByType = groupBy(signatureData, (e) => e.mutationType);
    const transposeSignature = Object.values(groupByType).map((signatures) =>
      signatures.reduce((obj, e) => {
        return {
          Type: e.mutationType,
          [e.signatureName]: e.contribution,
          ...obj,
        };
      }, {})
    );

    // write data to tsv file
    const signatureTSV = unparse(transposeSignature, {
      delimiter: '\t',
    });
    const signatureFilePath = path.join(workPath, 'signature.tsv');
    fs.writeFileSync(signatureFilePath, signatureTSV);

    // modify and include parameters
    const transformArgs = {
      ...args,
      input_data: path.join(workPath, args.input_data),
      output: path.join(workPath, 'output'),
      signature_database: signatureFilePath,
    };
    logger.debug(transformArgs);
    const cliArgs = Object.entries(transformArgs)
      .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
      .join(' ');

    const { execa } = await import('execa');
    const { stdout, stderr } = await execa('python3', [
      'services/python/mSigPortal_Profiler_Extraction.py',
      cliArgs,
    ]);
    res.json({ stdout, stderr });
  } catch (error) {
    logger.error('error at /submit');
    next(error);
  }
}
const router = Router();

router.get('/ping', (req, res) => res.send(true));
router.post('/upload', upload);
router.post('/submit', submit);

module.exports = router;
