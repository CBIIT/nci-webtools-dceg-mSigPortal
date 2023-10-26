import { readdirSync, unlinkSync, writeFileSync, createReadStream } from 'fs';
import path from 'path';
import { stringify } from 'csv-stringify';
import { groupBy } from 'lodash-es';
import { getSignatureData, getSeqmatrixData } from './query.js';
import { execa } from 'execa';
import validator from 'validator';
import mapValues from 'lodash/mapValues.js';
import { randomUUID } from 'crypto';
import Papa from 'papaparse';
import knex from 'knex';
import { readJson, writeJson, mkdirs, getFiles } from './utils.js';
import { sendNotification } from './notifications.js';
import { formatObject } from './logger.js';
import axios from 'axios';
import FormData from 'form-data';
import { uploadDirectory } from './s3.js';
import { importUserSession } from './importSignatures.js';

function parseCSV(filepath) {
  const file = createReadStream(filepath);
  return new Promise((resolve, reject) => {
    Papa.parse(file, {
      header: true,
      transformHeader: (e) => e.trim(),
      transform: (e) => e.trim(),
      complete(results, file) {
        resolve(results.data);
      },
      error(err, file) {
        reject(err);
      },
    });
  });
}

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
  const { args, signatureQuery, seqmatrixQuery, id, email } = params;
  let paths = await getPaths(params, env);
  const submittedTime = new Date(
    (await readJson(paths.statusFile)).submittedAt
  );

  try {
    logger.info(id);
    logger.info(params);
    if (!id) throw new Error('Missing id');
    if (!validator.isUUID(id)) throw new Error('Invalid id');

    const inputFolder = path.resolve(env.INPUT_FOLDER, id);

    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    await mkdirs([paths.inputFolder, paths.outputFolder]);
    await writeJson(paths.paramsFile, params);
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      id,
      status: 'IN_PROGRESS',
    });
    await writeJson(
      paths.manifestFile,
      mapValues(paths, (value) => path.parse(value).base)
    );

    // await uploadWorkingDirectory(inputFolder, outputFolder, id, env);

    const connection = dbConnection;
    const limit = false;

    // query signature data
    const signatureColumns = ['signatureName', 'mutationType', 'contribution'];

    const signatureData = await getSignatureData(
      connection,
      signatureQuery,
      signatureColumns,
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

    ////////////////// query seqmatrix data ////////////////////////
    let seqmatrixData;
    let seqmatrixFilePath;
    let seqmatrixFileName;
    let tsvString;
    if (params.form.source === 'public') {
      const seqmatrixColumns = [
        'sample',
        'profile',
        'matrix',
        'mutationType',
        'mutations',
      ];
      seqmatrixData = await getSeqmatrixData(
        connection,
        seqmatrixQuery,
        seqmatrixColumns,
        limit
      );

      // transform data into format accepted by SigProfiler
      // Extract unique mutation types and samples
      const mutationTypes = [
        ...new Set(seqmatrixData.map((d) => d.mutationType)),
      ];
      const samples = [...new Set(seqmatrixData.map((d) => d.sample))];

      // Initialize the result object with mutation types as keys
      const result = mutationTypes.reduce((acc, mt) => {
        acc[mt] = {};
        return acc;
      }, {});

      // Fill in the result object with mutation counts for each sample
      samples.forEach((s) => {
        mutationTypes.forEach((mt) => {
          const count = seqmatrixData.reduce((acc, d) => {
            if (d.sample === s && d.mutationType === mt) {
              return acc + d.mutations;
            } else {
              return acc;
            }
          }, 0);
          result[mt][s] = count;
        });
      });

      // Write the result to a TSV file using csv-stringify
      const headers = ['MutationType', ...samples];
      const rows = Object.entries(result).map(([mt, counts]) => [
        mt,
        ...samples.map((s) => counts[s] || 0),
      ]);
      const tsvData = rows.map((row) => {
        return {
          MutationType: row[0],
          ...Object.fromEntries(
            samples.map((sample, i) => [sample, row[i + 1] || 0])
          ),
        };
      });

      seqmatrixFilePath = path.join(outputFolder, 'seqmatrix.tsv');
      tsvString = await new Promise((resolve, reject) => {
        stringify(
          tsvData,
          { delimiter: '\t', header: true },
          (error, tsvString) => {
            if (error) {
              reject(error);
            } else {
              resolve(tsvString);
            }
          }
        );
      });
      // Write the TSV string to the file using writeFileSync
      writeFileSync(seqmatrixFilePath, tsvString);
      logger.info('Result written to signature.tsv');

      seqmatrixFileName = path.basename(seqmatrixFilePath);

      // Write the TSV data to the input_data file
      // const inputFilePath = path.join(inputFolder, args.input_data);
      // writeFileSync(inputFilePath, tsvString);
      // console.log('Data written to ExtractionData.all');
    }

    // modify and include parameters
    const transformArgs = {
      ...args,
      // input_data: path.join(inputFolder, args.input_data),
      input_data:
        params.form.source === 'public'
          ? seqmatrixFilePath
          : path.join(inputFolder, args.input_data),
      output: path.join(outputFolder),
      signature_database: signatureFilePath,
      gpu: env.EXTRACTION_WORKER_TYPE === 'batch' ? 'True' : 'False',
    };

    const cliArgs = Object.entries(transformArgs)
      .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
      .join(' ');

    logger.info(`[${id}] Run extraction`);
    await execa(
      'python3',
      ['services/python/mSigPortal-SigProfilerExtractor.py', cliArgs],
      { shell: true, all: true }
    )
      .pipeStdout(process.stdout)
      .pipeStderr(process.stderr)
      .pipeAll(path.resolve(outputFolder, 'extraction_log.txt'));

    logger.info('Finished Extraction');

    // update paths
    paths = { ...paths, ...getResultsPaths(params, env) };

    // import signatures data to database
    const decomposedSignatures = await parseCSV(paths.decomposedSignatureFile);
    const denovoSignatures = await parseCSV(paths.denovoSignatureInput);
    function signatureMapping(e) {
      const { MutationType, ...signatures } = e;
      return Object.entries(signatures).map(([signatureName, mutations]) => ({
        signatureName,
        MutationType,
        mutations,
      }));
    }
    const transformSignatures = [
      ...decomposedSignatures.map(signatureMapping).flat(),
      ...denovoSignatures.map(signatureMapping).flat(),
    ];

    const localDb = knex({
      client: 'better-sqlite3',
      connection: {
        filename: path.join(outputFolder, `local.sqlite3`),
      },
      useNullAsDefault: true,
    });

    await importUserSession(localDb, { signature: transformSignatures });

    // parse signatureMap csv to JSON
    const signatureMap = await parseCSV(paths.signatureMapFile);
    await writeJson(paths.signatureMapJson, signatureMap);

    // run exploration calculation on denovo and decomposed solutions
    let denovoId, decomposedId;
    try {
      logger.info(`[${id}] Run Denovo Exploration`);
      const denovoFormData = new FormData();
      if (params.form.source === 'public') {
        denovoFormData.append(
          'matrixFile',
          createReadStream(seqmatrixFilePath)
        );
      } else {
        denovoFormData.append('matrixFile', createReadStream(paths.matrixFile));
      }

      denovoFormData.append(
        'exposureFile',
        createReadStream(paths.denovoExposureInput)
      );
      denovoFormData.append(
        'signatureFile',
        createReadStream(paths.denovoSignatureInput)
      );

      const denovoUpload = await axios.post(
        `${env.API_BASE_URL}/api/upload/${randomUUID()}`,
        denovoFormData,
        { headers: denovoFormData.getHeaders() }
      );

      const denovoExploration = await axios.post(
        `${env.API_BASE_URL}/api/submitExploration/${denovoUpload.data.id}`,
        {
          matrixFile:
            params.form.source === 'public'
              ? path.parse(seqmatrixFilePath).base
              : path.parse(paths.matrixFile).base,
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
      if (params.form.source === 'public') {
        decomposedFormData.append(
          'matrixFile',
          createReadStream(seqmatrixFilePath)
        );
      } else {
        decomposedFormData.append(
          'matrixFile',
          createReadStream(paths.matrixFile)
        );
      }

      decomposedFormData.append(
        'exposureFile',
        createReadStream(paths.decomposedExposureInput)
      );
      decomposedFormData.append(
        'signatureFile',
        createReadStream(paths.decomposedSignatureInput)
      );

      const decomposedUpload = await axios.post(
        `${env.API_BASE_URL}/api/upload/${randomUUID()}`,
        decomposedFormData,
        { headers: decomposedFormData.getHeaders() }
      );

      const decomposedExploration = await axios.post(
        `${env.API_BASE_URL}/api/submitExploration/${decomposedUpload.data.id}`,
        {
          matrixFile:
            params.form.source === 'public'
              ? path.parse(seqmatrixFilePath).base
              : path.parse(paths.matrixFile).base,
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
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'COMPLETED',
      stopped: new Date(),
    });

    // await uploadWorkingDirectory(inputFolder, outputFolder, id, env);

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
    logger.debug(
      `Execution Time: ${
        (new Date().getTime() - submittedTime.getTime()) / 1000
      }`
    );

    // send success notification if email was provided
    if (params.email) {
      logger.info(`[${id}] Sending success notification`);
      await sendNotification(
        params.email,
        `mSigPortal - Extraction Complete - ${params.jobName}`,
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
    logger.error(`[${id}] Sending error notification`);
    logger.error(error);
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'FAILED',
      error: { ...error },
      stopped: new Date(),
    });

    // await uploadWorkingDirectory(
    //   paths.inputFolder,
    //   paths.outputFolder,
    //   id,
    //   env
    // );
    logger.debug(
      `Execution Time: ${
        (new Date().getTime() - submittedTime.getTime()) / 1000
      }`
    );
    if (params.email) {
      await sendNotification(
        params.email,
        `mSigPortal - Extraction Analysis Failed - ${params.jobName}`,
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
  } finally {
    // delete input files
    for await (const file of getFiles(paths.inputFolder)) {
      if (path.basename(file) !== 'params.json') {
        unlinkSync(file);
      }
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
  let inputFolder;
  if (params.form.source === 'user' || params.form.source === 'public') {
    inputFolder = path.resolve(env.INPUT_FOLDER, id);
  } else {
    inputFolder = '';
  }

  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFile = path.resolve(inputFolder, 'params.json');
  const statusFile = path.resolve(outputFolder, 'status.json');
  const manifestFile = path.resolve(outputFolder, 'manifest.json');
  const databaseFile = path.resolve(outputFolder, 'results.db');

  // SigProfilerExtraction log
  const extractionLog = path.resolve(outputFolder, 'JOB_METADATA.txt');

  // matrix file - input for extraction and exploration
  const matrixFile =
    params.form.source === 'public'
      ? ''
      : path.resolve(inputFolder, args.input_data);

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
    decomposedSignatureInput,
  };
}

/**
 * Dynamically generate paths to results files/folders
 * @param {any} params
 * @param {any} env
 * @returns {any} paths
 */
function getResultsPaths(params, env) {
  const { id } = params;
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

  const getDirectories = (source) =>
    readdirSync(source, { withFileTypes: true })
      .filter((dirent) => dirent.isDirectory())
      .map((dirent) => dirent.name);
  // assumes the only folder at /outputFolder is the name of the context type
  const [contextType] = getDirectories(outputFolder);

  // map files to be used as input for exploration module
  const solutionsFolder = path.resolve(
    outputFolder,
    contextType,
    'Suggested_Solution'
  );

  const denovoFolder = path.resolve(
    solutionsFolder,
    `${contextType}_De-Novo_Solution`
  );

  const decomposedFolder = path.resolve(
    solutionsFolder,
    `COSMIC_${contextType}_Decomposed_Solution`
  );

  // files for denovo exploration input
  const denovoExposureInput = path.resolve(
    denovoFolder,
    'Activities',
    `${contextType}_De-Novo_Activities_refit.txt`
  );

  const denovoSignatureInput = path.resolve(
    denovoFolder,
    'Signatures',
    `${contextType}_De-Novo_Signatures.txt`
  );

  // files for decomposed exploration input
  const decomposedExposureInput = path.resolve(
    decomposedFolder,
    'Activities',
    `COSMIC_${contextType}_Activities.txt`
  );

  // signature map file
  const signatureMapFile = path.resolve(
    decomposedFolder,
    `De_Novo_map_to_COSMIC_${contextType}.csv`
  );

  // signature map file
  const signatureMapJson = path.resolve(
    decomposedFolder,
    `De_Novo_map_to_COSMIC_${contextType}.json`
  );

  const decomposedSignatureFile = path.resolve(
    decomposedFolder,
    'Signatures',
    `COSMIC_${contextType}_Signatures.txt`
  );

  return {
    denovoExposureInput,
    denovoSignatureInput,
    decomposedExposureInput,
    signatureMapFile,
    signatureMapJson,
    decomposedSignatureFile,
    contextType,
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
