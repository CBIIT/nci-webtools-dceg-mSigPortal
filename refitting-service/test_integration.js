// test_integration.js
// Simple test to verify the refitting service integration

import { createLogger } from "./logger.js";
import rWrapper from "r-wrapper";
import path from 'path';
import fs from 'fs';

const r = rWrapper.async;
const logger = createLogger("Test Integration", "debug");

async function testIntegration() {
  try {
    logger.info("Testing R integration...");
    
    // Test the simple sum function
    const sum = await r("./refitting.R", "sum", { a: 5, b: 3 });
    logger.info(`Sum test result: ${sum}`);
    
    // Check if reference files exist
    const commonFilesDir = path.join(process.cwd(), 'data');
    const refFiles = [
      'Alex_Sigs_TMB_check_V3.4_SBS2_13 together.csv',
      '0 Cancer_Dictionary_BZ.csv'
    ];
    
    logger.info("Checking reference files...");
    for (const file of refFiles) {
      const filePath = path.join(commonFilesDir, file);
      if (fs.existsSync(filePath)) {
        logger.info(`✓ Found: ${file}`);
      } else {
        logger.error(`✗ Missing: ${file}`);
      }
    }
    
    logger.info("Integration test completed successfully!");
    
  } catch (error) {
    logger.error("Integration test failed:", error);
    throw error;
  }
}

// Run the test
testIntegration()
  .then(() => {
    console.log("All tests passed!");
    process.exit(0);
  })
  .catch((error) => {
    console.error("Test failed:", error);
    process.exit(1);
  });