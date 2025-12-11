import Router from 'express-promise-router';
import path from 'path';
import fs from 'fs-extra';
import { v4 as uuid, validate as validateUUID } from 'uuid';
import { execa } from 'execa';
import { body, validationResult } from 'express-validator';
import { getWorker } from '../../workers.js';
import { mkdirs, writeJson, readJson } from '../../utils.js';

export const router = Router();
const env = process.env;

// Validation middleware
const validateRefittingInput = [
  body('genome')
    .optional()
    .isIn(['hg19', 'hg38'])
    .withMessage('Genome must be hg19 or hg38'),
  body('matchOnOncotree')
    .optional()
    .isBoolean()
    .withMessage('matchOnOncotree must be boolean'),
  body('outputFilename')
    .optional()
    .isString()
    .withMessage('outputFilename must be a string'),
];

/**
 * POST /submitRefitting/:id
 * Submit a refitting job
 */
router.post(
  '/submitRefitting/:id',
  validateRefittingInput,
  async (req, res) => {
    const logger = req.app.locals.logger;
    const id = req.params.id;

    if (!validateUUID(id)) {
      return res.status(400).json({
        success: false,
        error: 'Invalid job ID format. Must be a valid UUID.',
      });
    }

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

    const type = env.WORKER_TYPE || 'local';
    const worker = getWorker(type);

    await writeJson(paramsFilePath, req.body);
    await writeJson(statusFilePath, status);

    worker(id, req.app, 'refitting', env);
    res.json(status);
  }
);

/**
 * GET /refreshRefitting/:id
 * Check the status of a refitting job
 */
router.get('/refreshRefitting/:id', async (req, res) => {
  const logger = req.app.locals.logger;
  const { id: jobId } = req.params;

  if (!validateUUID(jobId)) {
    return res.status(400).json({
      success: false,
      error: 'Invalid job ID format. Must be a valid UUID.',
    });
  }

  try {
    const data = await getJobStatus(jobId);
    if (typeof data === 'string') {
      return res.status(404).json({
        success: false,
        error: data,
      });
    }
    res.json(data);
  } catch (error) {
    logger.error(`Error checking status for job ${jobId}:`, error);
    res.status(500).json({
      success: false,
      error: 'Internal server error',
      message: error.message,
    });
  }
});

/**
 * GET /refitting/run/:id
 * Execute the refitting job (called by worker)
 */
router.get('/refitting/run/:id', async (req, res) => {
  const logger = req.app.locals.logger;
  const { id: jobId } = req.params;

  // Validate that jobId is a valid UUID
  if (!validateUUID(jobId)) {
    return res.status(400).json({
      success: false,
      error: 'Invalid job ID format. Must be a valid UUID.',
    });
  }

  try {
    const inputPath = path.join(env.INPUT_FOLDER || './data/input', jobId);
    const outputPath = path.join(env.OUTPUT_FOLDER || './data/output', jobId);

    const paramsFile = path.join(inputPath, 'params.json');
    if (!fs.existsSync(paramsFile)) {
      return res.status(404).json({
        success: false,
        error: `Job not found: ${jobId}`,
      });
    }

    const params = await readJson(paramsFile);

    // Get file paths from input directory
    const files = fs.readdirSync(inputPath).filter((f) => f !== 'params.json');

    // Find the files based on common file extensions
    const mafFile = files.find(
      (f) =>
        f.toLowerCase().includes('maf') ||
        f.endsWith('.maf') ||
        f.endsWith('.txt')
    );
    const genomicFile = files.find((f) => f.toLowerCase().includes('genomic'));
    const clinicalFile = files.find((f) =>
      f.toLowerCase().includes('clinical')
    );

    if (!mafFile || !genomicFile || !clinicalFile) {
      return res.status(400).json({
        success: false,
        error: 'Missing required input files',
        foundFiles: files,
      });
    }

    const mafFilePath = path.join(inputPath, mafFile);
    const genomicFilePath = path.join(inputPath, genomicFile);
    const clinicalFilePath = path.join(inputPath, clinicalFile);

    // Start the refitting process asynchronously
    startRefittingJob({
      jobId,
      mafFilePath,
      genomicFilePath,
      clinicalFilePath,
      inputPath,
      outputPath,
      params,
      logger,
    });

    res.json({
      success: true,
      jobId: jobId,
      message: 'Refitting job started',
    });
  } catch (error) {
    logger.error(`Error running refitting job ${jobId}:`, error);
    res.status(500).json({
      success: false,
      error: 'Internal server error',
      message: error.message,
    });
  }
});

/**
 * Helper function to get job status, params, and manifest
 */
async function getJobStatus(id) {
  if (!validateUUID(id)) {
    return `${id} is not a valid ID`;
  }

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

/**
 * POST /refreshRefittingMulti
 * Check the status of multiple refitting jobs
 */
router.post('/refreshRefittingMulti', async (req, res) => {
  const logger = req.app.locals.logger;
  try {
    const ids = req.body;
    const statuses = await Promise.all(ids.map(getJobStatus));
    res.json(statuses);
  } catch (error) {
    logger.error('/refreshRefittingMulti Error');
    res.status(500).json({
      success: false,
      error: 'Internal server error',
      message: error.message,
    });
  }
});

/**
 * Start a refitting job asynchronously
 */
async function startRefittingJob({
  jobId,
  mafFilePath,
  genomicFilePath,
  clinicalFilePath,
  inputPath,
  outputPath,
  params,
  logger,
}) {
  const statusFile = path.join(outputPath, 'status.json');

  try {
    // Update status to processing
    await writeJson(statusFile, {
      id: jobId,
      status: 'IN_PROGRESS',
      submittedAt: new Date().toISOString(),
      params,
    });

    // Prepare arguments for the refitting service
    const refittingServicePath = path.join(
      process.cwd(),
      '../refitting-service'
    );

    const nodeArgs = [
      'app.js',
      jobId,
      '--mafFile',
      mafFilePath,
      '--genomicFile',
      genomicFilePath,
      '--clinicalFile',
      clinicalFilePath,
      '--outputPath',
      outputPath,
      '--genome',
      params.genome,
      '--matchOnOncotree',
      params.matchOnOncotree.toString(),
      '--outputFilename',
      params.outputFilename,
      '--jobName',
      params.jobName || jobId,
      '--email',
      params.email || '',
    ];

    logger.info(`Starting refitting process for job ${jobId}`);

    // Execute the refitting service
    const { stdout, stderr } = await execa('node', nodeArgs, {
      cwd: refittingServicePath,
      timeout: 30 * 60 * 1000, // 30 minutes timeout
    });

    console.log(`[${jobId}] Refitting service completed!`);
    console.log(`[${jobId}] STDOUT:`, stdout);
    if (stderr) {
      console.log(`[${jobId}] STDERR:`, stderr);
    }

    logger.info(`Refitting job ${jobId} completed successfully`);

    // Create manifest.json similar to extraction service
    const manifestData = {
      jobId,
      completedAt: new Date().toISOString(),
      files: [params.outputFilename],
    };
    const manifestFile = path.join(outputPath, 'manifest.json');
    await writeJson(manifestFile, manifestData);

    // Update status to completed
    const currentStatus = await readJson(statusFile);
    await writeJson(statusFile, {
      id: jobId,
      status: 'COMPLETED',
      submittedAt: currentStatus.submittedAt,
      stopped: new Date().toISOString(),
      params,
      outputFilename: params.outputFilename,
      stdout: stdout,
      stderr: stderr,
    });
  } catch (error) {
    logger.error(`Refitting job ${jobId} failed:`, error);

    // Update status to failed
    const currentStatus = await readJson(statusFile).catch(() => ({}));
    await writeJson(statusFile, {
      id: jobId,
      status: 'FAILED',
      submittedAt: currentStatus.submittedAt || new Date().toISOString(),
      stopped: new Date().toISOString(),
      params,
      error: error.message,
      stderr: error.stderr,
    });
  }
}
