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
    logger.info(`Starting refitting job ${id} with params:`, params);
    
    // Determine profile type (default to SBS for backward compatibility)
    const profile = (profileType || 'SBS').toUpperCase();
    logger.info(`Profile type: ${profile}`);
    
    // Validate profile type
    if (!['SBS', 'DBS'].includes(profile)) {
      throw new Error(`Invalid profile type: ${profile}. Must be either 'SBS' or 'DBS'`);
    }
    
    // Update status to IN_PROGRESS
    const currentStatus = await readJson(statusFilePath).catch(() => ({}));
    await writeJson(statusFilePath, {
      ...currentStatus,
      id,
      status: 'IN_PROGRESS',
      submittedAt: currentStatus.submittedAt || new Date().toISOString(),
    });
    
    // Validate required parameters
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

    // Set up common files directory (reference files)
    const commonFilesDir = path.join(process.cwd(), 'data');
    
    // Determine which R function to call based on profile type
    const rFunctionName = profile === 'DBS' ? 'run_dbs_refitting' : 'run_sbs_refitting';
    const defaultOutputFilename = profile === 'DBS' ? 'H_Burden_est_DBS.csv' : 'H_Burden_est.csv';
    
    // Call the appropriate R refitting function
    logger.info(`Calling R refitting script for job ${id} with function: ${rFunctionName}`);
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
    
    logger.info(`R ${profile} refitting completed successfully for job ${id}`);
    logger.info(`Results summary: ${result?.H_Burden?.length || 0} entries processed`);
    
    const end = new Date();
    const duration = (end - start) / 1000;
    logger.info(`Refitting job ${id} completed in ${duration} seconds`);
    
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
    logger.info(`[${id}] Checking email notification - email param: ${email}, params.email: ${params.email}`);
    if (params.email) {
      logger.info(`[${id}] Email is provided, preparing to send success notification`);
      
      try {
        const submittedTime = new Date(statusData.submittedAt || start);
        const executionTime = (new Date().getTime() - submittedTime.getTime()) / 1000;
        
        const emailData = {
          jobName: jobName || id,
          submittedAt: submittedTime.toISOString(),
          executionTime: executionTime,
          resultsUrl: `${env.APP_BASE_URL || 'http://localhost:8330'}/#/refitting/${id}`,
        };
        
        logger.info(`[${id}] Email notification data:`, emailData);
        
        await sendNotification(
          email,
          `mSigPortal - Refitting Complete - ${jobName || id}`,
          'templates/user-success-email.html',
          emailData,
          env
        );
        logger.info(`[${id}] Success notification sent successfully to ${email}`);
      } catch (emailError) {
        logger.error(`[${id}] Failed to send success notification:`, emailError);
        logger.error(`[${id}] Email error details:`, emailError.message);
        logger.error(`[${id}] Email error stack:`, emailError.stack);
      }
    } else {
      logger.info(`[${id}] No email provided, skipping notification. params:`, { email, 'params.email': params.email });
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
    logger.info(`[${id}] Checking failure email notification - email param: ${email}, params.email: ${params.email}`);
    if (params.email) {
      logger.info(`[${id}] Email is provided, preparing to send failure notification`);
      
      try {
        const submittedTime = new Date(statusData.submittedAt || start);
        const executionTime = (new Date().getTime() - submittedTime.getTime()) / 1000;
        
        const emailData = {
          jobName: jobName || id,
          submittedAt: submittedTime.toISOString(),
          executionTime: executionTime,
          error: error.message || JSON.stringify(error, null, 2),
        };
        
        logger.info(`[${id}] Failure email notification data:`, emailData);
        
        await sendNotification(
          email,
          `mSigPortal - Refitting Analysis Failed - ${jobName || id}`,
          'templates/user-failure-email.html',
          emailData,
          env
        );
        logger.info(`[${id}] Failure notification sent successfully to ${email}`);
      } catch (emailError) {
        logger.error(`[${id}] Failed to send failure notification:`, emailError);
        logger.error(`[${id}] Email error details:`, emailError.message);
        logger.error(`[${id}] Email error stack:`, emailError.stack);
      }
    } else {
      logger.info(`[${id}] No email provided, skipping failure notification. params:`, { email, 'params.email': params.email });
    }
    
    throw error;
  }
}
