import { readdir, unlinkSync, writeFileSync } from 'fs';
import path from 'path';
import { stringify } from 'csv-stringify';
// const tar = require('tar');
import { groupBy } from 'lodash-es';
import { getSignatureData } from './query.js';
import { execa } from 'execa';
import validator from 'validator';
import mapValues from 'lodash/mapValues.js';
import { readJson, writeJson, mkdirs } from './utils.js';
import { sendNotification } from './notifications.js';

export async function extraction(
  params,
  logger,
  dbConnection,
  env = process.env
) {
  const { args, signatureQuery, id, email } = params;
  const paths = await getPaths(params, env);
  const submittedTime = new Date();
  logger.info(paths);

  try {
    if (!id) throw new Error('Missing id');
    if (!validator.isUUID(id)) throw new Error('Invalid id');

    await mkdirs([paths.inputFolder, paths.outputFolder]);
    await writeJson(paths.paramsFile, params);
    await writeJson(paths.statusFile, { id, status: 'IN_PROGRESS' });
    await writeJson(
      paths.manifestFile,
      mapValues(paths, (value) => path.parse(value).base)
    );

    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    // query signature data
    const connection = dbConnection;
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
      (error, output) => writeFileSync(signatureFilePath, output)
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

    logger.info(`[${id}] Run extraction`);
    const { all } = await execa(
      'python3',
      ['services/python/mSigPortal-SigProfilerExtractor.py', cliArgs],
      { all: true, shell: true }
    );
    logger.debug(all);

    // write success status
    const status = { id, status: 'COMPLETED' };
    await writeJson(paths.statusFile, status);

    // send success notification if email was provided
    if (params.email) {
      //delete input files
      readdir(paths.inputFolder, (err, files) => {
        if (err) {
          console.log(err);
        }

        files.forEach((file) => {
          const fileDir = path.join(paths.inputFolder, file);

          if (file !== 'params.json') {
            unlinkSync(fileDir);
          }
        });
      });

      await sendNotification(
        params.email,
        `Extraction Complete - ${params.jobName}`,
        'templates/user-success-email.html',
        {
          jobName: params.jobName,
          submittedAt: submittedTime.toISOString(),
          executionTime:
            (new Date().getTime() - submittedTime.getTime()) / 1000,
          resultsUrl: `${env.APP_BASE_URL}/#/extraction/${id}`,
        }
      );
    }
    return { id };
  } catch (error) {
    // send error notification if email was provided
    console.log(error);
    logger.error(error);
    const status = { id, status: 'FAILED', error: { ...error } };
    await writeJson(paths.statusFile, status);
    //delete input files
    readdir(paths.inputFolder, (err, files) => {
      if (err) {
        console.log(err);
      }

      files.forEach((file) => {
        const fileDir = path.join(paths.inputFolder, file);

        if (file !== 'params.json') {
          unlinkSync(fileDir);
        }
      });
    });

    if (params.email) {
      await sendNotification(
        params.email,
        `Analysis Failed - ${params.jobName}`,
        'templates/user-failure-email.html',
        {
          jobName: params.jobName,
          submittedAt: submittedTime.toISOString(),
          executionTime:
            (new Date().getTime() - submittedTime.getTime()) / 1000,
          error: formatObject(error),
        }
      );

      return false;
    }
  }
}

/**
 * Returns a list of all file paths for a set of parameters
 * @param {any} params
 * @param {any} env
 * @returns {any} paths
 */
export async function getPaths(params, env = process.env) {
  const { id } = params;
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFile = path.resolve(inputFolder, 'params.json');
  const statusFile = path.resolve(outputFolder, 'status.json');
  const manifestFile = path.resolve(outputFolder, 'manifest.json');
  const databaseFile = path.resolve(outputFolder, 'results.db');

  return {
    inputFolder,
    outputFolder,
    paramsFile,
    statusFile,
    manifestFile,
    databaseFile,
  };
}

export async function waitUntilComplete(
  id,
  env = process.env,
  checkInterval = 1000
) {
  const start = Date.now();
  const statusFilePath = path.resolve(env.OUTPUT_FOLDER, id, 'status.json');
  const isComplete = ({ status }) => ['COMPLETED', 'FAILED'].includes(status);
  let statusFile = await readJson(statusFilePath);

  if (!statusFile) {
    throw new Error('No status file found');
  }

  while (!isComplete(statusFile)) {
    statusFile = await readJson(statusFilePath);
    await setTimeout(checkInterval);
  }

  return Date.now() - start;
}
