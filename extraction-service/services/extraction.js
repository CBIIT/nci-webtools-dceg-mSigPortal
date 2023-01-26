import fs from 'fs-extra';
import path from 'path';
import { stringify } from 'csv-stringify';
// const tar = require('tar');
import { groupBy } from 'lodash-es';
import { getSignatureData } from './query.js';
import { execa } from 'execa';

export async function extraction(params, app, logger2, env = process.env) {
  const id = params.id;
  const logger = app.locals.logger;
  try {
    const { args, signatureQuery, id, email } = params;

    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    // query signature data
    const connection = app.locals.connection;
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
    const signatureFilePath = path.join(outputFolder, 'signature.tsv');
    stringify(
      transposeSignature,
      {
        header: true,
        delimiter: '\t',
      },
      (error, output) => fs.writeFileSync(signatureFilePath, output)
    );

    // modify and include parameters
    const transformArgs = {
      ...args,
      input_data: path.join(inputFolder, args.input_data),
      output: path.join(outputFolder),
      signature_database: signatureFilePath,
    };
    logger.debug(transformArgs);
    const cliArgs = Object.entries(transformArgs)
      .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
      .join(' ');

    const { all } = await execa(
      'python3',
      ['services/python/mSigPortal-SigProfilerExtractor.py', cliArgs],
      { all: true, shell: true }
    );
    logger.debug(all);
    return { id };
  } catch (error) {
    logger.debug(error);
    logger.error('error at /submit');
  }
}
