import Router from 'express-promise-router';
import multer from 'multer';
import path from 'path';
import fs from 'fs-extra';
import { v4 as uuid, validate as validateUUID } from 'uuid';
import { execa } from 'execa';
import { body, validationResult } from 'express-validator';

export const router = Router();

// Configure multer for file uploads
const upload = multer({
  dest: './data/input/', // Temporary location, will be moved to proper job folder
  limits: {
    fileSize: 100 * 1024 * 1024, // 100MB limit
  },
  fileFilter: function (req, file, cb) {
    // Accept only text files and csv files
    const allowedTypes = ['.txt', '.csv', '.tsv', '.maf'];
    const ext = path.extname(file.originalname).toLowerCase();
    if (allowedTypes.includes(ext)) {
      cb(null, true);
    } else {
      cb(new Error(`Invalid file type: ${ext}. Allowed types: ${allowedTypes.join(', ')}`));
    }
  }
});

// Validation middleware
const validateRefittingInput = [
  body('genome').optional().isIn(['hg19', 'hg38']).withMessage('Genome must be hg19 or hg38'),
  body('matchOnOncotree').optional().isBoolean().withMessage('matchOnOncotree must be boolean'),
  body('outputFilename').optional().isString().withMessage('outputFilename must be a string'),
];

/**
 * POST /submitRefitting/:id
 * Submit a refitting job with 3 required files
 */
router.post('/submitRefitting/:id', 
  upload.fields([
    { name: 'mafFile', maxCount: 1 },
    { name: 'genomicFile', maxCount: 1 },
    { name: 'clinicalFile', maxCount: 1 }
  ]),
  validateRefittingInput,
  async (req, res) => {
    const logger = req.app.locals.logger;
    
    console.log(`=== [${new Date().toISOString()}] REFITTING REQUEST RECEIVED ===`);
    console.log('Request params:', req.params);
    console.log('Request body:', req.body);
    console.log('Request files:', req.files ? Object.keys(req.files) : 'No files');
    
    try {
      // Check validation errors
      const errors = validationResult(req);
      if (!errors.isEmpty()) {
        console.log('Validation errors:', errors.array());
        return res.status(400).json({
          success: false,
          error: 'Validation failed',
          details: errors.array()
        });
      }

      // Check that all required files are uploaded
      const files = req.files;
      if (!files || !files.mafFile || !files.genomicFile || !files.clinicalFile) {
        console.log('Missing files. Received files:', files ? Object.keys(files) : 'none');
        return res.status(400).json({
          success: false,
          error: 'Missing required files. Please upload mafFile, genomicFile, and clinicalFile.'
        });
      }

      const jobId = req.params.id || uuid();
      console.log('Processing job ID:', jobId);
      
      // Validate that jobId is a valid UUID
      if (!validateUUID(jobId)) {
        console.log('Invalid UUID format:', jobId);
        return res.status(400).json({
          success: false,
          error: 'Invalid job ID format. Must be a valid UUID.'
        });
      }

      const inputPath = path.join(process.env.INPUT_FOLDER || './data/input', jobId);
      const outputPath = path.join(process.env.OUTPUT_FOLDER || './data/output', jobId);
      
      console.log('Input path:', inputPath);
      console.log('Output path:', outputPath);
      
      // Ensure both input and output directories exist
      fs.ensureDirSync(inputPath);
      fs.ensureDirSync(outputPath);
      
      // Move files from temp directory to job directory
      const mafFilePath = path.join(inputPath, `mafFile_${files.mafFile[0].originalname}`);
      const genomicFilePath = path.join(inputPath, `genomicFile_${files.genomicFile[0].originalname}`);
      const clinicalFilePath = path.join(inputPath, `clinicalFile_${files.clinicalFile[0].originalname}`);
      
      console.log('Moving files to:');
      console.log('  MAF:', mafFilePath);
      console.log('  Genomic:', genomicFilePath);
      console.log('  Clinical:', clinicalFilePath);
      
      await fs.move(files.mafFile[0].path, mafFilePath);
      await fs.move(files.genomicFile[0].path, genomicFilePath);
      await fs.move(files.clinicalFile[0].path, clinicalFilePath);
      
      console.log('Files moved successfully');

      // Extract parameters with defaults
      const params = {
        jobId,
        genome: req.body.genome || 'hg19',
        matchOnOncotree: req.body.matchOnOncotree === 'true' || false,
        outputFilename: req.body.outputFilename || 'H_Burden_est.csv'
      };

      // Write params.json to input folder (similar to extraction service)
      const paramsFilePath = path.join(inputPath, 'params.json');
      await fs.writeJson(paramsFilePath, params);

      // Write initial status.json to output folder
      const statusFilePath = path.join(outputPath, 'status.json');
      const initialStatus = {
        id: jobId,
        status: 'SUBMITTED',
        startTime: new Date().toISOString()
      };
      await fs.writeJson(statusFilePath, initialStatus);

      logger.info(`Starting refitting job ${jobId} with files:`, {
        mafFile: {
          originalname: files.mafFile[0].originalname,
          path: files.mafFile[0].path,
          size: files.mafFile[0].size
        },
        genomicFile: {
          originalname: files.genomicFile[0].originalname,
          path: files.genomicFile[0].path,
          size: files.genomicFile[0].size
        },
        clinicalFile: {
          originalname: files.clinicalFile[0].originalname,
          path: files.clinicalFile[0].path,
          size: files.clinicalFile[0].size
        },
        inputPath,
        outputPath,
        params
      });

      console.log('Starting refitting job asynchronously...');
      // Start the refitting process asynchronously
      startRefittingJob({
        jobId,
        mafFilePath,
        genomicFilePath,
        clinicalFilePath,
        inputPath,
        outputPath,
        params,
        logger
      });

      console.log('Responding with success...');
      // Return job ID immediately
      const response = {
        success: true,
        jobId: jobId,
        message: 'Refitting job started successfully',
        status: 'processing'
      };
      console.log('Response:', response);
      res.json(response);

    } catch (error) {
      logger.error('Error in refitting endpoint:', error);
      res.status(500).json({
        success: false,
        error: 'Internal server error',
        message: error.message
      });
    }
  }
);

/**
 * GET /refreshRefitting/:id
 * Check the status of a refitting job
 */
router.get('/refreshRefitting/:id', async (req, res) => {
  const logger = req.app.locals.logger;
  const { id: jobId } = req.params;

  // Validate that jobId is a valid UUID
  if (!validateUUID(jobId)) {
    return res.status(400).json({
      success: false,
      error: 'Invalid job ID format. Must be a valid UUID.'
    });
  }

  try {
    const jobStatus = await getJobStatus(jobId);
    if (typeof jobStatus === 'string') {
      return res.status(404).json({
        success: false,
        error: jobStatus
      });
    }
    res.json({
      success: true,
      ...jobStatus.status
    });
  } catch (error) {
    logger.error(`Error checking status for job ${jobId}:`, error);
    res.status(500).json({
      success: false,
      error: 'Internal server error',
      message: error.message
    });
  }
});

/**
 * Helper function to get job status, params
 */
async function getJobStatus(id) {
  if (!validateUUID(id)) {
    return `${id} is not a valid ID`;
  }

  try {
    const inputPath = path.join(process.env.INPUT_FOLDER || './data/input', id);
    const outputPath = path.join(process.env.OUTPUT_FOLDER || './data/output', id);
    
    const paramsFile = path.join(inputPath, 'params.json');
    const statusFile = path.join(outputPath, 'status.json');
    
    // Check if files exist
    if (!fs.existsSync(statusFile)) {
      return `Job not found: ${id}`;
    }

    const params = fs.existsSync(paramsFile) ? await fs.readJson(paramsFile) : {};
    const status = await fs.readJson(statusFile);
    
    // If job is completed, check for output file
    if (status.status === 'COMPLETED') {
      const outputFile = path.join(outputPath, status.outputFilename || 'H_Burden_est.csv');
      if (fs.existsSync(outputFile)) {
        status.downloadUrl = `/refitting/download/${id}/${status.outputFilename || 'H_Burden_est.csv'}`;
      }
    }

    return { params, status };
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
      message: error.message
    });
  }
});

/**
 * Start a refitting job asynchronously
 */
async function startRefittingJob({ jobId, mafFilePath, genomicFilePath, clinicalFilePath, inputPath, outputPath, params, logger }) {
  const statusFile = path.join(outputPath, 'status.json');
  
  console.log(`=== [${new Date().toISOString()}] STARTING REFITTING JOB ${jobId} ===`);
  console.log('Input parameters:', { jobId, mafFilePath, genomicFilePath, clinicalFilePath, inputPath, outputPath, params });
  
  try {
    console.log(`[${jobId}] Updating status to PROCESSING...`);
    // Update status to processing
    await fs.writeJson(statusFile, {
      id: jobId,
      status: 'PROCESSING',
      startTime: new Date().toISOString(),
      params
    });
    console.log(`[${jobId}] Status updated to PROCESSING`);

    // Prepare arguments for the refitting service
    const refittingServicePath = path.join(process.cwd(), '../refitting-service');
    console.log(`[${jobId}] Refitting service path: ${refittingServicePath}`);
    
    const nodeArgs = [
      'app.js',
      jobId,
      '--mafFile', mafFilePath,
      '--genomicFile', genomicFilePath,
      '--clinicalFile', clinicalFilePath,
      '--outputPath', outputPath,
      '--genome', params.genome,
      '--matchOnOncotree', params.matchOnOncotree.toString(),
      '--outputFilename', params.outputFilename
    ];

    console.log(`[${jobId}] Node.js arguments:`, nodeArgs);
    console.log(`[${jobId}] Working directory: ${refittingServicePath}`);
    console.log(`[${jobId}] About to execute: node ${nodeArgs.join(' ')}`);
    
    logger.info(`Starting refitting process for job ${jobId}`);

    // Execute the refitting service
    console.log(`[${jobId}] Executing refitting service...`);
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
      files: [params.outputFilename]
    };
    const manifestFile = path.join(outputPath, 'manifest.json');
    await fs.writeJson(manifestFile, manifestData);

    // Update status to completed
    await fs.writeJson(statusFile, {
      id: jobId,
      status: 'COMPLETED',
      startTime: (await fs.readJson(statusFile)).startTime,
      endTime: new Date().toISOString(),
      params,
      outputFilename: params.outputFilename,
      stdout: stdout,
      stderr: stderr
    });

  } catch (error) {
    console.log(`=== [${new Date().toISOString()}] ERROR in refitting job ${jobId} ===`);
    console.error(`[${jobId}] Error details:`, error);
    console.error(`[${jobId}] Error message:`, error.message);
    console.error(`[${jobId}] Error stack:`, error.stack);
    if (error.stdout) console.log(`[${jobId}] Error STDOUT:`, error.stdout);
    if (error.stderr) console.log(`[${jobId}] Error STDERR:`, error.stderr);
    
    logger.error(`Refitting job ${jobId} failed:`, error);

    // Update status to failed
    const currentStatus = await fs.readJson(statusFile).catch(() => ({}));
    await fs.writeJson(statusFile, {
      id: jobId,
      status: 'FAILED',
      startTime: currentStatus.startTime || new Date().toISOString(),
      endTime: new Date().toISOString(),
      params,
      error: error.message,
      stderr: error.stderr
    });
    
    console.log(`[${jobId}] Status updated to FAILED`);
  }
}