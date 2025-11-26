import { Router } from 'express';
import { randomUUID } from 'crypto';
import fs from 'fs';
import path from 'path';
import archiver from 'archiver';
import rWrapper from 'r-wrapper';
import { parseCSV } from '../general.js';
import { getExposureData, getSignatureData } from '../../query.js';
import { createCacheMiddleware } from '../../cache.js';
import isUUID from 'validator/lib/isUUID.js';
import { mkdirs, writeJson } from '../../utils.js';
import { getWorker } from '../../workers.js';

const r = rWrapper.async;
const env = process.env;

export async function submit(req, res, next) {
  const { id } = req.params;
  if (!isUUID(id)) res.status(500).json('Invalid ID');

  const inputFolder = path.resolve(env.INPUT_FOLDER, id);
  const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);
  const paramsFilePath = path.resolve(inputFolder, 'params.json');
  const statusFilePath = path.resolve(outputFolder, 'status.json');
  await mkdirs([inputFolder, outputFolder]);

  const status = {
    id,
    status: 'SUBMITTED',
    submittedAt: new Date(),
  };

  // use fargate worker if email is provided, otherwise use local worker
  const type =
    env.NODE_ENV === 'development' || !req.body?.email ? 'local' : 'fargate';
  const worker = getWorker(type);

  await writeJson(paramsFilePath, req.body);
  await writeJson(statusFilePath, status);

  worker(id, req.app, 'visualization', env);
  res.json(status);
}

async function getJobStatus(id) {
  if (!isUUID(id)) return `${id} is not a valid ID`;
  try {
    const inputFolder = path.resolve(env.INPUT_FOLDER, id);
    const outputFolder = path.resolve(env.OUTPUT_FOLDER, id);

    const paramsFilePath = path.resolve(inputFolder, 'params.json');
    const statusFilePath = path.resolve(outputFolder, 'status.json');
    const manifestFilePath = path.resolve(outputFolder, 'manifest.json');
    const params = await readJson(paramsFilePath);
    const status = await readJson(statusFilePath);
    const manifest = await readJson(manifestFilePath);

    return { params, status, manifest };
  } catch (error) {
    return error.message;
  }
}

export async function refresh(req, res, next) {
  const { logger } = req.app.locals;
  try {
    const id = req.params.id;
    const data = await getJobStatus(id);
    res.json(data);
  } catch (error) {
    logger.error('/refreshExtraction Error');
    next(error);
  }
}

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
  const resultsPath = path.resolve(env.OUTPUT_FOLDER, id);

  Object.keys(paths).map((key) => {
    const fullPath = path.resolve(paths[key]);

    if (fullPath.includes(resultsPath))
      newPaths[key] = fullPath.replace(resultsPath, '');
  });
  return newPaths;
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

async function getResults(req, res, next) {
  const { logger } = req.app.locals;
  logger.info(`/getResults: Retrieving Results for ${req.body.id}`);
  const userResults = path.resolve(env.OUTPUT_FOLDER, req.body.id, 'results');
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
  const { logger } = req.app.locals;
  const { fn, args, id = randomUUID() } = req.body;

  // config info for R functions
  const rConfig = {
    prefix: env.DATA_BUCKET_PREFIX,
    bucket: env.DATA_BUCKET,
    wd: path.resolve(env.DATA_FOLDER),
  };

  // create directory for results if needed
  const savePath = path.join('output', id, 'results', fn, '/');
  await mkdirs([path.resolve(rConfig.wd, savePath)]);

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
    logger.error(`/visualizationWrapper: An error occurred with fn: ${fn}`);
    next(err);
  }
}

function visualizationDownload(req, res, next) {
  // const { logger } = req.app.locals;
  // logger.info(
  //   `/visualization/download: id:${req.query.id} file:${req.query.file}`
  // );
  // const file = path.resolve(
  //   path.join(
  //     config.results.folder,
  //     req.query.id,
  //     'results/output',
  //     req.query.file
  //   )
  // );
  // if (file.indexOf(path.resolve(config.results.folder)) == 0) {
  //   res.download(file);
  // } else {
  //   logger.info('traversal error');
  res.status(500).end('Not found');
  // }
}

// Generate public data files for download
async function visualizationDownloadPublic(req, res, next) {
  // logger.info(`/visualization/downloadPublic`);
  // const { id, ...args } = req.body.args;
  // const { study, cancerType, experimentalStrategy } = args;
  // const savePath = path.join(
  //   config.results.folder,
  //   id,
  //   `results/msigportal-${study}-${cancerType}-${experimentalStrategy}`
  // );
  // // config info for R functions
  // const rConfig = {
  //   prefix: env.DATA_BUCKET_PREFIX,
  //   bucket: env.DATA_BUCKET,
  //   wd: path.resolve(env.OUTPUT_FOLDER),
  // };
  // try {
  //   await r('services/R/visualizeWrapper.R', 'downloadPublicData', {
  //     args,
  //     config: { ...rConfig, savePath },
  //   });
  //   const file = path.resolve(
  //     path.join(
  //       config.results.folder,
  //       id,
  //       `results/msigportal-${study}-${cancerType}-${experimentalStrategy}.tar.gz`
  //     )
  //   );
  //   if (file.indexOf(path.resolve(config.results.folder)) == 0) {
  //     res.download(file);
  //   } else {
  //     logger.error('visualizationDownloadPublic failed');
  res.status(500).end('Not found');
  //   }
  // } catch (err) {
  //   logger.error(err);
  //   res.status(500).json(err.message);
  //   next(err);
  // }
}

async function downloadWorkspace(req, res, next) {
  const { logger } = req.app.locals;
  logger.info(`/visualization/downloadWorkspace`);
  const { state, id } = req.body;
  const session = path.resolve(env.OUTPUT_FOLDER, id);
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
    let { study, strategy, cancer, signatureSetName, profile, matrix } =
      req.body;
    const pickNotNull = (obj, nulls = [null, undefined, '']) => 
      Object.fromEntries(Object.entries(obj).filter(([_, v]) => !nulls.includes(v)));

    // determine signature set names from profile and matrix (e.g. SBS96 -> *_SBS96)
    const signatureSetNames = await connection
      .select('signatureSetName')
      .from('exposure')
      .where({study})
      .andWhere('signatureSetName', 'like', `%${profile + matrix}`);

    const signatureSetNameValues = signatureSetNames.map(s => s.signatureSetName);
    if (!signatureSetName || !signatureSetNameValues.includes(signatureSetName)) {
      signatureSetName = signatureSetNameValues[0];
    }

    const exposureData = await getExposureData(
      connection,
      pickNotNull({ study, strategy, signatureSetName }),
      '*',
      1e8
    );
    
    const signatureData = await getSignatureData(
      connection,
      pickNotNull({ strategy, signatureSetName }),
      '*',
      1e8
    );

    let seqmatrixData = await connection
      .select('*', connection.raw('concat(profile, matrix) as "profileMatrix"'))
      .from('seqmatrix')
      .where(pickNotNull({ study, strategy, cancer, profile, matrix }))

    const estimatedMutations = connection
      .select('e.cancer', 'e.sample', 'mutationType', connection.raw('SUM(contribution * exposure) as estimated'))
      .from('signature as s')
      .innerJoin('exposure as e', function() {
          this.on('s.strategy', '=', 'e.strategy')
              .andOn('s.signatureSetName', '=', 'e.signatureSetName')
              .andOn('s.signatureName', '=', 'e.signatureName')
      })
      .where(pickNotNull({
          'e.study': study,
          'e.strategy': strategy,
          'e.cancer': cancer,
          's.profile': profile,
          's.matrix': matrix,
          's.signatureSetName': signatureSetName
      }))
      .andWhere('exposure', '>', 0)
      .groupBy('e.cancer', 'e.sample', 'mutationType')
  
    const mutations = connection('seqmatrix as s')
      .select(
        's.study', 's.strategy', 's.cancer', 's.sample', 's.profile', 's.matrix', 's.mutationType', 
          connection.raw('cast(s.mutations as bigint) as mutations'),
          connection.raw(`cast(floor(a.estimated) as bigint) as estimated`),
        // connection.raw(`1.0 * s.mutations * a.estimated as mutations_estimated`),
        // connection.raw(`1.0 * s.mutations * s.mutations as mutations2`), // we need to materialize these columns types to avoid underflow/overflow when calculating cosine similarity
        // connection.raw(`1.0 * a.estimated * a.estimated as estimated2`),
      )
      .join('estimated_mutations as a', function() {
          this.on('s.mutationType', '=', 'a.mutationType')
              .andOn('s.sample', '=', 'a.sample')
              .andOn('s.cancer', '=', 'a.cancer')
      })
      .where(pickNotNull({
          's.study': study,
          's.strategy': strategy,
          's.cancer': cancer,
          's.profile': profile,
          's.matrix': matrix
      }))
  
    const cosineSimilarityData = await connection
      .with('estimated_mutations', estimatedMutations)
      .with('mutations', mutations)
      .select(
          'b.cancer', 
          'b.sample', 
          connection.raw('SUM(b.mutations * b.estimated) as numerator'),
          connection.raw('SUM(b.mutations * b.mutations) as denominator1'),
          connection.raw('SUM(b.estimated * b.estimated) as denominator2')
      )
      .from('mutations as b')
      .groupBy('b.cancer', 'b.sample');
    
    // todo: execute this in the database (requires that we resolve underflow/overflow issues)
    const cosineSimilarityMap = cosineSimilarityData.reduce((acc, curr) => ({
      ...acc,
      [curr.sample]: +curr.numerator / (Math.sqrt(+curr.denominator1) * Math.sqrt(curr.denominator2))
    }), {});

    if (
      !exposureData?.length ||
      !seqmatrixData?.length ||
      !signatureData?.length
    ) {
      console.log('exposureData', exposureData[0], exposureData?.length);
      console.log('seqmatrixData', seqmatrixData[0], seqmatrixData?.length);
      console.log('signatureData', signatureData[0], signatureData?.length);
      throw new Error('No data found');
    }
    const args = { exposureData, seqmatrixData, signatureData };
    const results = await wrapper('wrapper', { fn: 'getTreeLeaf', args });
    results.output.params = { study, strategy, cancer, signatureSetName, profile, matrix };
    for (let record of results.output?.attributes || []) {
      record.Cosine_similarity = cosineSimilarityMap[record.Sample] || 0;
    }
    res.json(results);
  } catch (error) {
    console.log(error);
    next(error);
  }
}

const router = Router();
router.post('/submitVisualization/:id?', submit);
router.post('/getResults', getResults);
router.post('/visualizationWrapper', visualizationWrapper);
router.get('/visualization/download', visualizationDownload);
router.post('/visualization/downloadPublic', visualizationDownloadPublic);
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
      req.body.profile,
      req.body.matrix,
    ].join(':')
  ),
  getPublicTreeLeafData
);

export {
  router,
  getFiles,
  getMatrices,
  getRelativePath,
  getResultsFiles,
  wrapper,
  visualizationWrapper,
  visualizationDownload,
  visualizationDownloadPublic,
  downloadWorkspace,
};
