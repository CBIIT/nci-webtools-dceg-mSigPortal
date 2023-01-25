import fs from 'fs-extra';
import path from 'path';
import { unparse } from 'papaparse';
// const tar = require('tar');
import { groupBy } from 'lodash';
import { getSignatureData } from '../../server/services/query';

export async function extraction(params, app, logger2, env = process.env) {
  const id = params.id;
  const logger = app.locals.logger;
  try {
    const { args, signatureQuery, id, email } = params;
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
