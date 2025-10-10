import Router from 'express-promise-router';
import multer from 'multer';
import path from 'path';
import fs from 'fs-extra';
import { v4 as uuid } from 'uuid';
import { execa } from 'execa';
import { body, validationResult } from 'express-validator';

export const router = Router();

// Configure multer for file uploads
const storage = multer.diskStorage({
  destination: function (req, file, cb) {
    // Create temp directory for this job
    const jobId = req.jobId || uuid();
    req.jobId = jobId;
    const uploadPath = path.join(process.env.UPLOAD_FOLDER || './uploads', jobId);
    fs.ensureDirSync(uploadPath);
    cb(null, uploadPath);
  },
  filename: function (req, file, cb) {
    // Use original filename with field name prefix for identification
    cb(null, `${file.fieldname}_${file.originalname}`);
  }
});

const upload = multer({
  storage: storage,
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
 * POST /refitting/sbs
 * Submit a refitting job with 3 required files
 */
router.post('/refitting/sbs', 
  upload.fields([
    { name: 'mafFile', maxCount: 1 },
    { name: 'genomicFile', maxCount: 1 },
    { name: 'clinicalFile', maxCount: 1 }
  ]),
  validateRefittingInput,
  async (req, res) => {
    const logger = req.app.locals.logger;
    
    try {
      // Check validation errors
      const errors = validationResult(req);
      if (!errors.isEmpty()) {
        return res.status(400).json({
          success: false,
          error: 'Validation failed',
          details: errors.array()
        });
      }

      // Check that all required files are uploaded
      const files = req.files;
      if (!files || !files.mafFile || !files.genomicFile || !files.clinicalFile) {
        return res.status(400).json({
          success: false,
          error: 'Missing required files. Please upload mafFile, genomicFile, and clinicalFile.'
        });
      }

      const jobId = req.jobId || uuid();
      const uploadPath = path.join(process.env.UPLOAD_FOLDER || './uploads', jobId);
      
      // Extract file paths
      const mafFilePath = files.mafFile[0].path;
      const genomicFilePath = files.genomicFile[0].path;
      const clinicalFilePath = files.clinicalFile[0].path;

      // Extract parameters with defaults
      const params = {
        genome: req.body.genome || 'hg19',
        matchOnOncotree: req.body.matchOnOncotree === 'true' || false,
        outputFilename: req.body.outputFilename || 'H_Burden_est.csv'
      };

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
        uploadPath,
        params
      });

      // Create output directory
      const outputPath = path.join(uploadPath, 'output');
      fs.ensureDirSync(outputPath);

      // Start the refitting process asynchronously
      startRefittingJob({
        jobId,
        mafFilePath,
        genomicFilePath,
        clinicalFilePath,
        outputPath,
        params,
        logger
      });

      // Return job ID immediately
      res.json({
        success: true,
        jobId: jobId,
        message: 'Refitting job started successfully',
        status: 'processing'
      });

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
 * GET /refitting/status/:jobId
 * Check the status of a refitting job
 */
router.get('/refitting/status/:jobId', async (req, res) => {
  const logger = req.app.locals.logger;
  const { jobId } = req.params;

  try {
    const uploadPath = path.join(process.env.UPLOAD_FOLDER || './uploads', jobId);
    const outputPath = path.join(uploadPath, 'output');
    const statusFile = path.join(uploadPath, 'status.json');
    
    // Check if status file exists
    if (!fs.existsSync(statusFile)) {
      return res.status(404).json({
        success: false,
        error: 'Job not found'
      });
    }

    const status = await fs.readJson(statusFile);
    
    // If job is completed, check for output file
    if (status.status === 'completed') {
      const outputFile = path.join(outputPath, status.outputFilename || 'H_Burden_est.csv');
      if (fs.existsSync(outputFile)) {
        status.downloadUrl = `/refitting/download/${jobId}/${status.outputFilename || 'H_Burden_est.csv'}`;
      }
    }

    res.json({
      success: true,
      ...status
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
 * GET /refitting/download/:jobId/:filename
 * Download the results of a completed refitting job
 */
router.get('/refitting/download/:jobId/:filename', async (req, res) => {
  const logger = req.app.locals.logger;
  const { jobId, filename } = req.params;

  try {
    const uploadPath = path.join(process.env.UPLOAD_FOLDER || './uploads', jobId);
    const outputPath = path.join(uploadPath, 'output');
    const filePath = path.join(outputPath, filename);

    // Security check: ensure the file is within the job's output directory
    const resolvedPath = path.resolve(filePath);
    const resolvedOutputPath = path.resolve(outputPath);
    if (!resolvedPath.startsWith(resolvedOutputPath)) {
      return res.status(403).json({
        success: false,
        error: 'Access denied'
      });
    }

    if (!fs.existsSync(filePath)) {
      return res.status(404).json({
        success: false,
        error: 'File not found'
      });
    }

    res.download(filePath, filename);

  } catch (error) {
    logger.error(`Error downloading file for job ${jobId}:`, error);
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
async function startRefittingJob({ jobId, mafFilePath, genomicFilePath, clinicalFilePath, outputPath, params, logger }) {
  const uploadPath = path.dirname(mafFilePath);
  const statusFile = path.join(uploadPath, 'status.json');
  
  try {
    // Update status to processing
    await fs.writeJson(statusFile, {
      jobId,
      status: 'processing',
      startTime: new Date().toISOString(),
      params
    });

    // Prepare arguments for the refitting service
    const refittingServicePath = path.join(process.cwd(), '../refitting-service');
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

    logger.info(`Starting refitting process for job ${jobId}`);

    // Execute the refitting service
    const { stdout, stderr } = await execa('node', nodeArgs, {
      cwd: refittingServicePath,
      timeout: 30 * 60 * 1000, // 30 minutes timeout
    });

    logger.info(`Refitting job ${jobId} completed successfully`);

    // Update status to completed
    await fs.writeJson(statusFile, {
      jobId,
      status: 'completed',
      startTime: (await fs.readJson(statusFile)).startTime,
      endTime: new Date().toISOString(),
      params,
      outputFilename: params.outputFilename,
      stdout: stdout,
      stderr: stderr
    });

  } catch (error) {
    logger.error(`Refitting job ${jobId} failed:`, error);

    // Update status to failed
    await fs.writeJson(statusFile, {
      jobId,
      status: 'failed',
      startTime: (await fs.readJson(statusFile)).startTime,
      endTime: new Date().toISOString(),
      params,
      error: error.message,
      stderr: error.stderr
    });
  }
}