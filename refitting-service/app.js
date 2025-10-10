import { createLogger } from "./logger.js";
import rWrapper from "r-wrapper";
import path from 'path';
import fs from 'fs';
const r = rWrapper.async;

// Parse command line arguments
const args = process.argv.slice(2);
const jobId = args[0];

if (!jobId) {
  console.error("No jobId specified");
  process.exit(1);
}

// Parse additional arguments
const argMap = {};
for (let i = 1; i < args.length; i += 2) {
  if (args[i].startsWith('--') && i + 1 < args.length) {
    argMap[args[i].substring(2)] = args[i + 1];
  }
}

const logger = createLogger(`Job ${jobId}`, process.env?.LOG_LEVEL || "debug");

try {
  await refitting(jobId, argMap);
  process.exit(0);
} catch (e) {
  logger.error(e);
  process.exit(1);
}

async function refitting(jobId, params) {
  const start = new Date();
  try {
    logger.info(`Starting refitting job ${jobId} with params:`, params);
    
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
    if (!fs.existsSync(params.outputPath)) {
      fs.mkdirSync(params.outputPath, { recursive: true });
    }

    // Set up common files directory (reference files)
    const commonFilesDir = path.join(process.cwd(), 'data');
    
    // Call the R refitting function
    const result = await r("./refitting.R", "run_sbs_refitting", {
      maf_file: params.mafFile,
      genomic_file: params.genomicFile,
      clinical_file: params.clinicalFile,
      output_dir: params.outputPath,
      common_files_dir: commonFilesDir,
      genome: params.genome || 'hg19',
      save_csv: true,
      out_file: params.outputFilename || 'H_Burden_est.csv',
      match_on_oncotree: params.matchOnOncotree === 'true'
    });
    
    logger.info(`R refitting completed successfully for job ${jobId}`);
    logger.info(`Results summary: ${result?.H_Burden?.length || 0} entries processed`);
    
    const end = new Date();
    const duration = (end - start) / 1000;
    logger.info(`Refitting job ${jobId} completed in ${duration} seconds`);
    
    return result;
  } catch (e) {
    const done = new Date();
    const duration = (done - start) / 1000;
    logger.info(`Refitting job ${jobId} failed after ${duration} seconds`);
    logger.error(e);
    throw e;
  }
}
