import { Router } from 'express';
import { randomUUID } from 'crypto';
import { validate } from 'uuid';
import fs from 'fs-extra';
import path from 'path';
import glob from 'glob';
import archiver from 'archiver';
import rWrapper from 'r-wrapper';
import AWS from 'aws-sdk';
import tar from 'tar';
import config from '../../../config.json' assert { type: 'json' };
import logger from '../../logger.js';
import { parseCSV, importUserSession } from '../general.js';
import { schema } from './userSchema.js';
import { getExposureData, getSignatureData } from '../../query.js';
import { createCacheMiddleware } from '../../cache.js';
import { execa } from 'execa';
const r = rWrapper.async;

// config info for R functions
const rConfig = {
  s3Data: config.data.s3,
  bucket: config.data.bucket,
  localData: path.resolve(config.data.localData),
  wd: path.resolve(config.results.folder),
};

async function* getFiles(dir) {
  const dirents = await fs.promises.readdir(dir, { withFileTypes: true });
  for (const dirent of dirents) {
    const res = path.resolve(dir, dirent.name);
    if (dirent.isDirectory()) {
      yield* getFiles(res);
    } else {
      yield res;
    }
  }
}

// transform all matrix files into a single json object
async function getMatrices(matrixFolder) {
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
        const data = await parseCSV(file);
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

function getRelativePath(paths, id = '') {
  let newPaths = {};
  const resultsPath = path.resolve(config.results.folder, id);
  const dataPath = path.resolve(config.data.localData);

  Object.keys(paths).map((key) => {
    const fullPath = path.resolve(paths[key]);

    if (fullPath.includes(resultsPath))
      newPaths[key] = fullPath.replace(resultsPath, '');
    else if (fullPath.includes(dataPath))
      newPaths[key] = fullPath.replace(dataPath + '/', '');
  });
  return newPaths;
}

async function profilerExtraction(paramsObj) {
  return new Promise(async (resolve, reject) => {
    const outputFolder = paramsObj.Output_Dir;

    const cliArgs = Object.entries(paramsObj)
      .reduce((params, [key, value]) => [...params, `--${key} ${value}`], [])
      .join(' ');

    const command = [
      'services/python/mSigPortal_Profiler_Extraction.py',
      cliArgs,
    ];
    try {
      const { all: scriptOutput } = await execa('python3', command, {
        all: true,
        shell: true,
      });
      try {
        // parse matrix files and transform into single json file
        const matrixFiles = path.resolve(outputFolder, 'output');
        const seqmatrix = await getMatrices(matrixFiles);
        const matricesFile = path.join(outputFolder, 'matrices.json');
        fs.writeFileSync(matricesFile, JSON.stringify(seqmatrix));

        // parse cluster data if option was used
        let cluster = [];
        if (paramsObj?.Cluster == 'True') {
          const clusterFile = path.join(
            outputFolder,
            `Cluster/Result/${paramsObj.Project_ID}_clustered_class_All.txt`
          );
          const clusterData = fs.existsSync(clusterFile)
            ? await parseCSV(clusterFile)
            : [];
          // filter cluster objects for relevent values and rename keys if needed
          if (clusterData.length) {
            cluster = clusterData.map((e) => {
              const { project, samples, ID, __parsed_extra, ...rest } = e;
              return { sample: samples, geneId: ID, ...rest };
            });
          }
        }
        resolve({
          scriptOutput,
          seqmatrix,
          ...(cluster && { cluster }),
        });
      } catch (error) {
        logger.error('Error parsing matrix data. ');
        reject(error);
      }
    } catch (error) {
      reject(error);
    }
  });
}
// retrieves parsed data files - modify paths if needed
async function getResultsFiles(resultsPath, id = '') {
  const statisticsPath = path.join(resultsPath, 'Statistics.txt');
  const matrixPath = path.join(resultsPath, 'matrix_files_list.txt');
  const downloadsPath = path.join(resultsPath, 'output');
  let matrixList = [];
  let statistics = '';
  let downloads = [];
  if (fs.existsSync(matrixPath)) matrixList = await parseCSV(matrixPath);
  if (fs.existsSync(statisticsPath))
    statistics = fs.readFileSync(statisticsPath, 'utf8');
  if (fs.existsSync(downloadsPath)) {
    downloads = fs
      .readdirSync(downloadsPath)
      .filter((file) => file.endsWith('.zip'))
      .map((file) => file);
  }
  if (fs.existsSync(path.join(downloadsPath, 'vcf_files_zip'))) {
    downloads = [
      ...downloads,
      ...fs
        .readdirSync(path.join(downloadsPath, 'vcf_files_zip'))
        .filter((file) => file.endsWith('.zip'))
        .map((file) => path.join('vcf_files_zip', file)),
    ];
  }
  // convert to relative paths
  matrixList.forEach(
    (plot) => (plot.Path = getRelativePath({ Path: plot.Path }, id).Path)
  );
  return {
    statistics: statistics,
    matrixList: matrixList,
    downloads: downloads,
  };
}

async function submitVisualization(req, res, next) {
  res.setTimeout(15 * 60 * 1000, () => {
    res.status(504).send('request timed out');
  });

  try {
    const userId = req.body.Project_ID;
    // resolve file paths
    const params = {
      ...req.body,
      Input_Path: path.resolve(
        config.folders.input,
        userId,
        req.body.Input_Path
      ),
      Output_Dir: path.resolve(config.folders.output, userId),
      ...(req.body?.Bed && {
        Bed: path.resolve(config.folders.input, userId, req.body.Bed),
      }),
    };
    fs.mkdirSync(params.Output_Dir, { recursive: true });

    const { scriptOutput, seqmatrix, cluster } = await profilerExtraction(
      params
    );
    const outputFolder = path.resolve(params.Output_Dir);
    // fs.mkdirSync(resultsPath, { recursive: true });

    // import data into user session table
    const connection = req.app.locals.sqlite(userId, 'local');
    const importStatus = await importUserSession(
      connection,
      { seqmatrix, cluster },
      schema
    );
    if (!importStatus)
      next(new Error('Failed to import matrix data into database'));
    if (fs.existsSync(path.join(outputFolder, 'svg_files_list.txt'))) {
      res.json({
        scriptOutput,
        id: userId,
        ...(await getResultsFiles(outputFolder, userId)),
      });
    } else {
      logger.error(
        '/profilerExtraction: An error occured while extracting profiles'
      );
      res.status(500).json({ scriptOutput });
    }
  } catch (error) {
    logger.error('/profilerExtraction: Caught error while extracting profiles');
    logger.error(error);
    res.status(500).json(error);
  }
}

async function getResults(req, res, next) {
  logger.info(`/getResults: Retrieving Results for ${req.body.id}`);
  const userResults = path.resolve(
    config.results.folder,
    req.body.id,
    'results'
  );
  if (fs.existsSync(path.join(userResults, 'svg_files_list.txt'))) {
    res.json(await getResultsFiles(userResults, req.body.id));
  } else {
    logger.info('/getResults: Results not found');
    res.status(500).json('Results not found');
  }
}
async function wrapper(fn, args) {
  return await JSON.parse(await r('services/R/visualizeWrapper.R', fn, args));
}
// Visualization Calculation functions
async function visualizationWrapper(req, res, next) {
  const { fn, args, id } = req.body;
  logger.debug(`/visualizationWrapper: %o`, req.body);
  // create directory for results if needed
  const savePath = id ? path.join(id, 'results', fn, '/') : null;
  if (id) fs.mkdirSync(path.join(rConfig.wd, savePath), { recursive: true });
  try {
    const { scriptOutput, ...rest } = await wrapper('wrapper', {
      fn,
      args,
      config: {
        ...rConfig,
        savePath,
      },
    });
    // generate an id if getPublicData is called
    res.json({
      id: fn == 'getPublicData' ? randomUUID() : id,
      scriptOutput,
      ...rest,
    });
  } catch (err) {
    logger.error(`/visualizationWrapper: An error occured with fn: ${fn}`);
    next(err);
  }
}

async function getSignaturesUser(req, res, next) {
  logger.info('/getSignaturesUser: Parsing File');
  try {
    const file = path.resolve(req.body.path);
    logger.debug(file);
    if (file.indexOf(path.resolve(config.results.folder)) == 0) {
      const data = await parseCSV(file);
      res.json(data);
    } else {
      logger.info('traversal error');
      res.status(500).end('Not found');
    }
  } catch (err) {
    logger.info('/getSignaturesUser: An error occured');
    logger.error(err);
    res.status(500).json(err.message);
  }
}

function visualizationDownload(req, res, next) {
  logger.info(
    `/visualization/download: id:${req.query.id} file:${req.query.file}`
  );
  const file = path.resolve(
    path.join(
      config.results.folder,
      req.query.id,
      'results/output',
      req.query.file
    )
  );
  if (file.indexOf(path.resolve(config.results.folder)) == 0) {
    res.download(file);
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
}

// Generate public data files for download
async function visualizationDownloadPublic(req, res, next) {
  logger.info(`/visualization/downloadPublic`);
  const { id, ...args } = req.body.args;
  const { study, cancerType, experimentalStrategy } = args;
  const savePath = path.join(
    config.results.folder,
    id,
    `results/msigportal-${study}-${cancerType}-${experimentalStrategy}`
  );
  try {
    await r('services/R/visualizeWrapper.R', 'downloadPublicData', {
      args,
      config: { ...rConfig, savePath },
    });
    const file = path.resolve(
      path.join(
        config.results.folder,
        id,
        `results/msigportal-${study}-${cancerType}-${experimentalStrategy}.tar.gz`
      )
    );
    if (file.indexOf(path.resolve(config.results.folder)) == 0) {
      res.download(file);
    } else {
      logger.error('visualizationDownloadPublic failed');
      res.status(500).end('Not found');
    }
  } catch (err) {
    logger.error(err);
    res.status(500).json(err.message);
    next(err);
  }
}

async function submitQueue(req, res, next) {
  const id = req.body.args.id[1];
  const sqs = new AWS.SQS();
  try {
    // upload archived project directory
    await new AWS.S3()
      .upload({
        Body: tar
          .c({ sync: true, gzip: true, C: config.results.folder }, [id])
          .read(),
        Bucket: config.queue.bucket,
        Key: `${config.queue.inputKeyPrefix}${id}/${id}.tgz`,
      })
      .promise();
    const { QueueUrl } = await sqs
      .getQueueUrl({ QueueName: config.queue.url })
      .promise();
    await sqs
      .sendMessage({
        QueueUrl: QueueUrl,
        MessageDeduplicationId: id,
        MessageGroupId: id,
        MessageBody: JSON.stringify({
          ...req.body,
          timestamp: new Date().toLocaleString('en-US', {
            timeZone: 'America/New_York',
          }),
        }),
      })
      .promise();
    logger.info('Queue submitted ID: ' + id);
    res.json({ id });
  } catch (err) {
    logger.info('Queue failed to submit ID: ' + id);
    next(err);
  }
}

async function getQueueResults(req, res, next) {
  try {
    const s3 = new AWS.S3();
    const { id } = req.params;
    logger.info(`Fetch Queue Result: ${id}`);
    // validate id format
    if (!validate(id)) next(new Error(`Invalid request`));
    // ensure output directory exists
    const resultsFolder = path.resolve(config.results.folder, id);
    await fs.promises.mkdir(resultsFolder, { recursive: true });
    // find objects which use the specified id as the prefix
    const objects = await s3
      .listObjectsV2({
        Bucket: config.queue.bucket,
        Prefix: `${config.queue.outputKeyPrefix}${id}/`,
      })
      .promise();
    // download results
    for (let { Key } of objects.Contents) {
      const filename = path.basename(Key);
      const filepath = path.resolve(resultsFolder, filename);
      // download results if they do not exist
      if (!fs.existsSync(filepath)) {
        logger.info(`Downloading result: ${Key}`);
        const object = await s3
          .getObject({
            Bucket: config.queue.bucket,
            Key,
          })
          .promise();
        await fs.promises.writeFile(filepath, object.Body);
        // extract and delete archive
        if (path.extname(filename) == '.tgz') {
          await new Promise((resolve, reject) => {
            fs.createReadStream(filepath)
              .on('end', () =>
                fs.unlink(filepath, (err) => {
                  if (err) {
                    reject(err);
                  } else {
                    resolve();
                  }
                })
              )
              .pipe(tar.x({ strip: 1, C: resultsFolder }));
          });
        }
      }
    }
    let paramsPath = path.resolve(resultsFolder, `params.json`);
    if (fs.existsSync(paramsPath)) {
      const data = JSON.parse(String(await fs.promises.readFile(paramsPath)));
      res.json(data);
    } else {
      next(new Error(`Params not found`));
    }
  } catch (error) {
    next(error);
  }
}

async function getVisExample(req, res, next) {
  try {
    const { example } = req.params;
    logger.info(`Fetching example: ${example}`);
    // check exists
    const examplePath = path.resolve(
      config.data.examples,
      'visualization',
      `${example}.tgz`
    );
    if (fs.existsSync(examplePath)) {
      // copy example to results with unique id
      const id = randomUUID();
      const resultsPath = path.resolve(config.results.folder, id);
      await fs.promises.mkdir(resultsPath, { recursive: true });
      // await fs.copy(examplePath, resultsPath);
      await new Promise((resolve, reject) => {
        fs.createReadStream(examplePath)
          .on('end', () => resolve())
          .on('error', (err) => reject(err))
          .pipe(tar.x({ strip: 1, C: resultsPath }));
      });
      const paramsPath = path.join(resultsPath, `params.json`);
      // rename file paths with new ID
      let params = JSON.parse(String(await fs.promises.readFile(paramsPath)));
      const oldID = params.visualization.main.id;
      await replace({
        files: paramsPath,
        from: new RegExp(oldID, 'g'),
        to: id,
      });
      params = JSON.parse(String(await fs.promises.readFile(paramsPath)));
      const svgPath = path.join(resultsPath, 'results', 'svg_files_list.txt');
      if (fs.existsSync(svgPath)) {
        const matrixPath = path.join(
          resultsPath,
          'results',
          'matrix_files_list.txt'
        );
        await replace({
          files: [svgPath, matrixPath],
          from: new RegExp(oldID, 'g'),
          to: id,
        });
      }
      // rename files with new ID
      glob(path.join(resultsPath, 'results/**'), (error, files) => {
        files.forEach((file) => {
          if (file.includes(oldID)) {
            fs.renameSync(file, file.replace(oldID, id));
          }
        });
      });
      res.json({ id, state: params.visualization });
    } else {
      throw `Invalid example`;
    }
  } catch (error) {
    next(error);
  }
}

async function downloadWorkspace(req, res, next) {
  logger.info(`/visualization/downloadWorkspace`);
  const { state, id } = req.body;
  const session = path.resolve(config.results.folder, id);
  const archive = archiver('zip', {
    zlib: { level: 6 }, // Sets the compression level.
  });
  if (fs.existsSync(session)) {
    try {
      await fs.promises.writeFile(
        path.join(session, 'params.json'),
        JSON.stringify(state)
      );
      archive.on('finish', function (error) {
        return res.end();
      });
      archive
        .directory(session, false)
        .on('error', function (err) {
          throw err;
        })
        .pipe(res);
      archive.finalize();
    } catch (err) {
      next(err);
    }
  } else {
    logger.info('traversal error');
    res.status(500).end('Not found');
  }
}

async function getPublicTreeLeafData(req, res, next) {
  try {
    const { connection } = req.app.locals;
    const { study, strategy, cancer, signatureSetName, profileMatrix } =
      req.body;
    const exposureData = await getExposureData(
      connection,
      { study, strategy },
      '*',
      1e8
    );
    const signatureData = await getSignatureData(
      connection,
      { strategy, signatureSetName },
      '*',
      1e8
    );
    let params = {};
    for (let key of ['study', 'strategy', 'cancer']) {
      if (req.body[key]) {
        params[key] = req.body[key];
      }
    }
    let seqmatrixData = await connection
      .select('*', connection.raw('concat(profile, matrix) as "profileMatrix"'))
      .from('seqmatrix')
      .where(params)
      .andWhere(connection.raw('concat(profile, matrix)'), 'in', profileMatrix);
    // console.log('exposureData', exposureData[0], exposureData?.length);
    // console.log('seqmatrixData', seqmatrixData[0], seqmatrixData?.length);
    // console.log('signatureData', signatureData[0], signatureData?.length);
    if (
      !exposureData?.length ||
      !seqmatrixData?.length ||
      !signatureData?.length
    ) {
      throw new Error('No data found');
    }
    const args = { exposureData, seqmatrixData, signatureData };
    const results = await wrapper('wrapper', { fn: 'getTreeLeaf', args });
    res.json(results);
  } catch (error) {
    console.log(error);
    next(error);
  }
}
const router = Router();
router.post('/submitVisualization', submitVisualization);
router.post('/getResults', getResults);
router.post('/getSignaturesUser', getSignaturesUser);
router.post('/visualizationWrapper', visualizationWrapper);
router.get('/visualization/download', visualizationDownload);
router.post('/visualization/downloadPublic', visualizationDownloadPublic);
router.post('/queue', submitQueue);
router.get('/getQueueResults/:id', getQueueResults);
router.get('/getVisExample/:example', getVisExample);
router.post('/downloadWorkspace', downloadWorkspace);
router.post(
  '/treeLeaf',
  createCacheMiddleware((req) =>
    [
      'treeLeaf',
      req.body.study,
      req.body.strategy,
      req.body.cancer,
      req.body.signatureSetName,
      req.body.profileMatrix,
    ].join(':')
  ),
  getPublicTreeLeafData
);

export {
  router,
  getFiles,
  getMatrices,
  getRelativePath,
  profilerExtraction,
  getResultsFiles,
  wrapper,
  visualizationWrapper,
  visualizationDownload,
  visualizationDownloadPublic,
  submitQueue,
  getQueueResults,
  getVisExample,
  downloadWorkspace,
};
