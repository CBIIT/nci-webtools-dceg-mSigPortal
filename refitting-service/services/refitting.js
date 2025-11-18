import rWrapper from "r-wrapper";
import path from 'path';
import fs from 'fs';
import { readJson, writeJson } from './utils.js';
import { sendNotification } from './notifications.js';

const r = rWrapper.async;

/**
 * Main refitting function that processes a refitting job
 * @param {object} params - Job parameters including file paths, genome, etc.
 * @param {object} logger - Logger instance
 * @param {object} env - Environment variables
 * @returns {Promise<object>} Result of the refitting operation
 */
export async function refitting(params, logger, env = process.env) {
  const start = new Date();
  const { id, mafFile, genomicFile, clinicalFile, outputPath, genome, outputFilename, matchOnOncotree, jobName, email, profileType } = params;
  
  // Get paths for status file
  const statusFilePath = path.join(outputPath, 'status.json');
  
  try {
    logger.info(`Starting refitting job ${id}`);
    
    // Determine profile type - require either profileType or params.profileType
    const profile = (profileType || params.profileType);
    
    if (!profile) {
      throw new Error(`Missing required parameter: profileType. Must be either 'SBS' or 'DBS'`);
    }
    
    const normalizedProfile = profile.toUpperCase();
    
    // Validate profile type
    if (!['SBS', 'DBS'].includes(normalizedProfile)) {
      throw new Error(`Invalid profile type: ${normalizedProfile}. Must be either 'SBS' or 'DBS'`);
    }
    
    // Update status to IN_PROGRESS
    const currentStatus = await readJson(statusFilePath).catch(() => ({}));
    await writeJson(statusFilePath, {
      ...currentStatus,
      id,
      status: 'IN_PROGRESS',
      submittedAt: currentStatus.submittedAt || new Date().toISOString(),
    });
    
    // Validate required parameters (profileType already validated above)
    const required = ['mafFile', 'genomicFile', 'clinicalFile', 'outputPath'];
    const missing = required.filter(param => !params[param]);
    if (missing.length > 0) {
      throw new Error(`Missing required parameters: ${missing.join(', ')}`);
    }

    // Validate files exist
    const fileParams = ['mafFile', 'genomicFile', 'clinicalFile'];
    for (const param of fileParams) {
      if (!fs.existsSync(params[param])) {
        throw new Error(`File not found: ${params[param]}`);
      }
    }

    // Ensure output directory exists
    if (!fs.existsSync(outputPath)) {
      fs.mkdirSync(outputPath, { recursive: true });
    }

    // Set up refitting files directory (reference files) - use DATA_FOLDER env var with refitting subfolder
    const dataFolder = env.DATA_FOLDER;
    if (!dataFolder) {
      throw new Error('DATA_FOLDER environment variable must be set');
    }
    const commonFilesDir = path.join(dataFolder, 'refitting');
    
    // Determine which R function to call based on profile type
    const rFunctionName = normalizedProfile === 'DBS' ? 'run_dbs_refitting' : 'run_sbs_refitting';
    const defaultOutputFilename = normalizedProfile === 'DBS' ? 'H_Burden_est_DBS.csv' : 'H_Burden_est.csv';
    
    // Call the appropriate R refitting function
    logger.info(`Running ${normalizedProfile} refitting for job ${id}`);
    const result = await r("./refitting.R", rFunctionName, {
      maf_file: mafFile,
      genomic_file: genomicFile,
      clinical_file: clinicalFile,
      output_dir: outputPath,
      common_files_dir: commonFilesDir,
      genome: genome || 'hg19',
      save_csv: true,
      out_file: outputFilename || defaultOutputFilename,
      match_on_oncotree: matchOnOncotree === 'true' || matchOnOncotree === true
    });
    
    const end = new Date();
    const duration = (end - start) / 1000;
    logger.info(`Job ${id} completed successfully in ${duration}s - ${result?.H_Burden?.length || 0} entries processed`);
    
    // Update status to COMPLETED
    const statusData = await readJson(statusFilePath).catch(() => ({}));
    await writeJson(statusFilePath, {
      ...statusData,
      id,
      status: 'COMPLETED',
      stopped: new Date().toISOString(),
      params,
      result: {
        entriesProcessed: result?.H_Burden?.length || 0
      }
    });

    // Send success notification if email was provided
    if (params.email) {
      try {
        const submittedTime = new Date(statusData.submittedAt || start);
        const executionTime = (new Date().getTime() - submittedTime.getTime()) / 1000;
        
        const emailData = {
          jobName: jobName || id,
          submittedAt: submittedTime.toISOString(),
          executionTime: executionTime,
          resultsUrl: `${env.APP_BASE_URL || 'http://localhost:8330'}/#/refitting/${id}`,
        };
        
        await sendNotification(
          email,
          `mSigPortal - Refitting Complete - ${jobName || id}`,
          'templates/user-success-email.html',
          emailData,
          env
        );
        logger.info(`[${id}] Success notification sent to ${email}`);
      } catch (emailError) {
        logger.error(`[${id}] Failed to send success notification: ${emailError.message}`);
      }
    }
    
    return result;
    
  } catch (error) {
    const done = new Date();
    const duration = (done - start) / 1000;
    logger.error(`Refitting job ${id} failed after ${duration} seconds`);
    logger.error(error);
    
    // Update status to FAILED
    const statusData = await readJson(statusFilePath).catch(() => ({}));
    await writeJson(statusFilePath, {
      ...statusData,
      id,
      status: 'FAILED',
      stopped: new Date().toISOString(),
      params,
      error: error.message || String(error),
    });

    // Send failure notification if email was provided
    if (params.email) {
      try {
        const submittedTime = new Date(statusData.submittedAt || start);
        const executionTime = (new Date().getTime() - submittedTime.getTime()) / 1000;
        
        const emailData = {
          jobName: jobName || id,
          submittedAt: submittedTime.toISOString(),
          executionTime: executionTime,
          error: error.message || JSON.stringify(error, null, 2),
        };
        
        await sendNotification(
          email,
          `mSigPortal - Refitting Analysis Failed - ${jobName || id}`,
          'templates/user-failure-email.html',
          emailData,
          env
        );
        logger.info(`[${id}] Failure notification sent to ${email}`);
      } catch (emailError) {
        logger.error(`[${id}] Failed to send failure notification: ${emailError.message}`);
      }
    }
    
    throw error;
  }
}
