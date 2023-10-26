import fs from 'fs';
import path from 'path';
import rWrapper from 'r-wrapper';
import { execa } from 'execa';
import isUUID from 'validator/lib/isUUID.js';
import mapValues from 'lodash/mapValues.js';
import { parseCSV, parseTSV } from '../general.js';
import { schema } from './userSchema.js';
import { readJson, writeJson, mkdirs, getFiles } from '../../utils.js';
import { sqliteImport } from '../../sqlite.js';
import { sendNotification } from '../../notifications.js';
import { formatObject } from '../../logger.js';
const r = rWrapper.async;

// transform all matrix files into a single json object
export async function getMatrices(matrixFolder) {
  let files = [];
  const fileRegex = /\.(all|region)$/;
  for await (const f of getFiles(matrixFolder)) {
    if (fileRegex.test(f)) files = [...files, f];
  }
  const profileRegex = /\.([A-Z]+)\d+\.(all|region)$/;
  const matrixRegex = /\.[A-Z]+(\d+)\.(all|region)$/;

  return (
    await Promise.all(
      files.map(async (file) => {
        const profile = file.match(profileRegex)[1];
        const matrix = file.match(matrixRegex)[1];
        const { data } = await parseTSV(file);
        return data
          .map((e) => {
            const { MutationType, ...samples } = e;
            return Object.entries(samples).map(([sample, mutations]) => {
              const splitSample = sample.split('@');
              return {
                sample: splitSample[0],
                filter: splitSample[1] || '',
                profile,
                matrix,
                mutationType: MutationType,
                mutations,
              };
            });
          })
          .flat();
      })
    )
  ).flat();
}

export async function profilerExtraction(
  params,
  logger,
  dbConnection,
  env = process.env
) {
  const { args, email } = params;
  const id = args.Project_ID;
  const paths = await getPaths(args, env);
  const submittedTime = new Date(
    (await readJson(paths.statusFile)).submittedAt
  );

  logger.info(paths);

  try {
    if (!id) throw new Error('Missing id');
    if (!isUUID(id)) throw new Error('Invalid id');

    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
    const profilerExtractionOutput = path.resolve(
      outputFolder,
      'profilerExtraction'
    );

    await mkdirs([paths.inputFolder, paths.outputFolder]);
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      id,
      status: 'IN_PROGRESS',
    });
    await writeJson(
      paths.manifestFile,
      mapValues(paths, (value) => path.parse(value).base)
    );

    // modify and include parameters
    const transformArgs = {
      ...args,
      Input_Path: path.resolve(inputFolder, args.Input_Path),
      Output_Dir: profilerExtractionOutput,
      ...(args?.Bed && {
        Bed: path.resolve(inputFolder, args.Bed),
      }),
    };
    const cliArgs = Object.entries(transformArgs)
      .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
      .join(' ');

    logger.info(`[${id}] Run Profiler Extraction`);
    const script = await execa(
      'python3',
      ['services/python/mSigPortal_Profiler_Extraction.py', cliArgs],
      { shell: true, all: true }
    )
      .pipeStdout(process.stdout)
      .pipeStderr(process.stderr)
      .pipeAll(path.resolve(outputFolder, 'profiler_extraction_log.txt'));

    // parse all matrix files and to json file
    const matrixFiles = path.resolve(profilerExtractionOutput, 'output');
    if (!fs.existsSync(matrixFiles)) {
      const { stdout, stderr } = script;
      if (stdout) logger.error(stdout);
      if (stderr) logger.error(stderr);

      const message = stdout.split('\n').filter((e) => e.includes('Error'));
      if (message.length) throw new Error(message[0]);
      else throw new Error('Error occurred during profiler extraction.');
    }
    const seqmatrix = await getMatrices(matrixFiles);

    // parse cluster data if option was used
    let cluster = [];
    if (args?.Cluster == 'True') {
      const clusterFile = path.resolve(
        profilerExtractionOutput,
        `Cluster/Result/${id}_clustered_class_All.txt`
      );
      const clusterData = fs.existsSync(clusterFile)
        ? await parseCSV(clusterFile)
        : [];
      // filter cluster objects for relevant values and rename keys if needed
      if (clusterData.length) {
        cluster = clusterData.map((e) => {
          const { project, samples, ID, __parsed_extra, ...rest } = e;
          return { sample: samples, geneId: ID, ...rest };
        });
      }
    }

    logger.info(`[${id}] Create sqlite database`);
    const importStatus = await sqliteImport(
      dbConnection,
      { seqmatrix, cluster },
      schema
    );
    if (!importStatus)
      throw new Error('Failed to import matrix data into database');

    // write success status
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'COMPLETED',
    });

    logger.debug(
      `Execution Time: ${
        (new Date().getTime() - submittedTime.getTime()) / 1000
      }`
    );

    // send success notification if email was provided
    if (email) {
      logger.info(`[${id}] Sending success notification`);
      await sendNotification(
        params.email,
        `mSigPortal - Visualization Complete - ${params.jobName}`,
        'templates/user-success-email.html',
        {
          jobName: params.jobName,
          submittedAt: submittedTime.toISOString(),
          executionTime:
            (new Date().getTime() - submittedTime.getTime()) / 1000,
          resultsUrl: `${env.APP_BASE_URL}/#/visualization/${id}`,
        }
      );
    }
    return { id };
  } catch (error) {
    logger.error(error);
    // write failed status
    await writeJson(paths.statusFile, {
      ...(await readJson(paths.statusFile)),
      status: 'FAILED',
      error: error.message,
    });

    if (email) {
      logger.info(`[${id}] Sending error notification`);
      await sendNotification(
        params.email,
        `mSigPortal - Visualization Analysis Failed - ${params.jobName}`,
        'templates/user-failure-email.html',
        {
          jobName: params.jobName,
          submittedAt: submittedTime.toISOString(),
          executionTime:
            (new Date().getTime() - submittedTime.getTime()) / 1000,
          error: formatObject(error),
        }
      );
    }
    return { error };
  } finally {
    // delete input files
    for await (const file of getFiles(paths.inputFolder)) {
      if (path.basename(file) !== 'params.json') {
        fs.unlinkSync(file);
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
export async function getPaths(args, env = process.env) {
  const id = args.Project_ID;
  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFile = path.resolve(inputFolder, 'params.json');
  const statusFile = path.resolve(outputFolder, 'status.json');
  const manifestFile = path.resolve(outputFolder, 'manifest.json');
  // const databaseFile = path.resolve(outputFolder, 'results.db');

  return {
    inputFolder,
    outputFolder,
    paramsFile,
    statusFile,
    manifestFile,
    // databaseFile,
  };
}
