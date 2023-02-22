import { readdir, unlinkSync, writeFileSync, createReadStream } from 'fs';
import path from 'path';
import { stringify } from 'csv-stringify';
// const tar = require('tar');
import { groupBy } from 'lodash-es';
import { getSignatureData } from './query.js';
import { execa } from 'execa';
import validator from 'validator';
import mapValues from 'lodash/mapValues.js';
import { randomUUID } from 'crypto';
import { readJson, writeJson, mkdirs } from './utils.js';
import { sendNotification } from './notifications.js';
import { formatObject } from './logger.js';
import axios from 'axios';
import FormData from 'form-data';
import { uploadDirectory } from './s3.js';

async function uploadWorkingDirectory(inputFolder, outputFolder, id, env) {
  // upload input folder
  await uploadDirectory(
    inputFolder,
    path.join(env.INPUT_KEY_PREFIX, id),
    env.DATA_BUCKET,
    { region: env.AWS_DEFAULT_REGION }
  );

  // upload output folder
  await uploadDirectory(
    outputFolder,
    path.join(env.OUTPUT_KEY_PREFIX, id),
    env.DATA_BUCKET,
    { region: env.AWS_DEFAULT_REGION }
  );
}

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

    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    await mkdirs([paths.inputFolder, paths.outputFolder]);
    await writeJson(paths.paramsFile, params);
    await writeJson(paths.statusFile, { id, status: 'IN_PROGRESS' });
    await writeJson(
      paths.manifestFile,
      mapValues(paths, (value) => path.parse(value).base)
    );
    await uploadWorkingDirectory(inputFolder, outputFolder, id, env);

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

    let denovoId, decomposedId;
    // run exploration calculation on denovo and decomposed solutions
    try {
      logger.info(`[${id}] Run Denovo Exploration`);
      const denovoFormData = new FormData();
      denovoFormData.append('matrixFile', createReadStream(paths.matrixFile));
      denovoFormData.append(
        'exposureFile',
        createReadStream(paths.denovoExposureInput)
      );
      denovoFormData.append(
        'signatureFile',
        createReadStream(paths.denovoSignatureInput)
      );

      const denovoUpload = await axios.post(
        `${env.API_BASE_URL}/web/upload/${randomUUID()}`,
        denovoFormData,
        { headers: denovoFormData.getHeaders() }
      );

      const denovoExploration = await axios.post(
        `${env.API_BASE_URL}/web/submitExploration/${denovoUpload.data.id}`,
        {
          matrixFile: path.parse(paths.matrixFile).base,
          exposureFile: path.parse(paths.denovoExposureInput).base,
          signatureFile: path.parse(paths.denovoSignatureInput).base,
        }
      );

      denovoId = denovoExploration.data;
    } catch (error) {
      logger.error('Denovo Exploration Error');
      console.log(error);
      throw error.data;
    }

    try {
      logger.info(`[${id}] Run Decomposed Exploration`);
      const decomposedFormData = new FormData();
      decomposedFormData.append(
        'matrixFile',
        createReadStream(paths.matrixFile)
      );
      decomposedFormData.append(
        'exposureFile',
        createReadStream(paths.decomposedExposureInput)
      );
      decomposedFormData.append(
        'signatureFile',
        createReadStream(paths.decomposedSignatureInput)
      );

      const decomposedUpload = await axios.post(
        `${env.API_BASE_URL}/web/upload/${randomUUID()}`,
        decomposedFormData,
        { headers: decomposedFormData.getHeaders() }
      );

      const decomposedExploration = await axios.post(
        `${env.API_BASE_URL}/web/submitExploration/${decomposedUpload.data.id}`,
        {
          matrixFile: path.parse(paths.matrixFile).base,
          exposureFile: path.parse(paths.decomposedExposureInput).base,
          signatureFile: path.parse(paths.decomposedSignatureInput).base,
        }
      );

      decomposedId = decomposedExploration.data;
    } catch (error) {
      logger.error('Decomposed Exploration Error');
      throw error.data;
    }

    // add exploration ids to manifest
    await writeJson(paths.manifestFile, {
      ...mapValues(paths, (value) => path.parse(value).base),
      denovoId,
      decomposedId,
    });

    // write success status
    const status = { id, status: 'COMPLETED' };
    await writeJson(paths.statusFile, status);

    await uploadWorkingDirectory(inputFolder, outputFolder, id, env);

    // // upload denovo output
    // await uploadDirectory(
    //   path.resolve(env.OUTPUT_FOLDER, denovoId),
    //   path.join(env.OUTPUT_KEY_PREFIX, denovoId),
    //   env.DATA_BUCKET,
    //   { region: env.AWS_DEFAULT_REGION }
    // );

    // // upload decmoposed output
    // await uploadDirectory(
    //   path.resolve(env.OUTPUT_FOLDER, decomposedId),
    //   path.join(env.OUTPUT_KEY_PREFIX, decomposedId),
    //   env.DATA_BUCKET,
    //   { region: env.AWS_DEFAULT_REGION }
    // );

    // send success notification if email was provided
    if (params.email) {
      logger.info(`[${id}] Sending success notificaiton`);
      //delete input files
      // readdir(paths.inputFolder, (err, files) => {
      //   if (err) {
      //     console.log(err);
      //   }

      //   files.forEach((file) => {
      //     const fileDir = path.join(paths.inputFolder, file);

      //     if (file !== 'params.json') {
      //       unlinkSync(fileDir);
      //     }
      //   });
      // });
      await sendNotification(
        params.email,
        `Extraction Complete - ${params.jobName}`,
        'templates/user-success-email.html',
        {
          jobName: params.jobName,
          submittedAt: submittedTime.toISOString(),
          executionTime:
            (new Date().getTime() - submittedTime.getTime()) / 1000,
          resultsUrl: path.join(env.APP_BASE_URL, '#', 'extraction', id),
        }
      );
    }
    return { id };
  } catch (error) {
    // send error notification if email was provided
    logger.error(`[${id}] Sending error notificaiton`);
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

    await uploadWorkingDirectory(
      paths.inputFolder,
      paths.outputFolder,
      id,
      env
    );

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
  const { id, args } = params;
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFile = path.resolve(inputFolder, 'params.json');
  const statusFile = path.resolve(outputFolder, 'status.json');
  const manifestFile = path.resolve(outputFolder, 'manifest.json');
  const databaseFile = path.resolve(outputFolder, 'results.db');

  // map files to be used as input for exploration module
  const solutionsFolder = path.resolve(
    outputFolder,
    args.context_type,
    'Suggested_Solution'
  );
  const denovoFolder = path.resolve(
    solutionsFolder,
    `${args.context_type}_De-Novo_Solution`
  );
  const decomposedFolder = path.resolve(
    solutionsFolder,
    `COSMIC_${args.context_type}_Decomposed_Solution`
  );
  // SigProfilerExtraction log
  const extractionLog = path.resolve(outputFolder, 'JOB_METADATA.txt');

  // matrix file - input for extraction and exploration
  const matrixFile = path.resolve(inputFolder, args.input_data);

  // files for denovo exploration input
  const denovoExposureInput = path.resolve(
    denovoFolder,
    'Activities',
    `${args.context_type}_De-Novo_Activities_refit.txt`
  );
  const denovoSignatureInput = path.resolve(
    denovoFolder,
    'Signatures',
    `${args.context_type}_De-Novo_Signatures.txt`
  );

  // files for decomposed exploration input
  const decomposedExposureInput = path.resolve(
    decomposedFolder,
    'Activities',
    `COSMIC_${args.context_type}_Activities.txt`
  );
  const decomposedSignatureInput = path.resolve(outputFolder, 'signature.tsv');

  return {
    inputFolder,
    outputFolder,
    paramsFile,
    statusFile,
    manifestFile,
    databaseFile,
    extractionLog,
    matrixFile,
    denovoExposureInput,
    denovoSignatureInput,
    decomposedExposureInput,
    decomposedSignatureInput,
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
