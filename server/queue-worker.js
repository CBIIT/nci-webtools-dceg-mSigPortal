import fs from 'fs';
import path from 'path';
import AWS from 'aws-sdk';
import nodemailer from 'nodemailer';
import rWrapper from 'r-wrapper';
import tar from 'tar';
import config from './config.json' assert { type: 'json' };
import logger from './services/logger.js';
import { parseCSV } from './services/api/analysis.js';
import {
  getRelativePath,
  profilerExtraction,
} from './services/api/visualization/userVisualization.js';
import {
  defaultProfile,
  defaultMatrix,
  defaultFilter,
} from './services/utils.js';
const r = rWrapper.async;

(async function main() {
  // update aws configuration if all keys are supplied, otherwise
  // fall back to default credentials/IAM role
  if (config.aws) {
    AWS.config.update(config.aws);
  }
  // create required folders
  for (let folder of [config.logs.folder, config.results.folder]) {
    fs.mkdirSync(folder, { recursive: true });
  }
  receiveMessage();
})();
/**
 * Reads a template, substituting {tokens} with data values
 * @param {string} filepath
 * @param {object} data
 */
async function readTemplate(filePath, data) {
  const template = await fs.promises.readFile(path.resolve(filePath));
  // replace {tokens} with data values or removes them if not found
  return String(template).replace(
    /{[^{}]+}/g,
    (key) => data[key.replace(/[{}]+/g, '')] || ''
  );
}
/**
 * Writes the contents of a stream to a file and resolves once complete
 * @param {*} readStream
 * @param {*} filePath
 */
function streamToFile(readStream, filePath) {
  return new Promise((resolve, reject) => {
    const file = fs.createWriteStream(filePath);
    const stream = readStream.pipe(file);
    stream.on('error', (error) => reject(error));
    stream.on('close', (_) => resolve());
  });
}
/**
 * Downloads work files from s3 for calculation
 * @param {string} id
 * @param {string} savePath
 */
async function downloadS3(id, savePath) {
  const s3 = new AWS.S3();
  const objects = await s3
    .listObjectsV2({
      Bucket: config.queue.bucket,
      Prefix: `${config.queue.inputKeyPrefix}${id}/`,
    })
    .promise();
  // download work files
  for (let { Key } of objects.Contents) {
    const filename = path.basename(Key);
    const filepath = path.resolve(savePath, filename);
    logger.info(`Downloading: ${Key}`);
    const object = await s3
      .getObject({
        Bucket: config.queue.bucket,
        Key,
      })
      .promise();
    await fs.promises.writeFile(filepath, object.Body);
    // extract and delete archive
    if (path.extname(filename) == '.tgz') {
      fs.createReadStream(filepath)
        .pipe(tar.x({ strip: 1, C: savePath }))
        .once('finish', () =>
          fs.unlink(filepath, (e) => {
            if (e) logger.error(e);
          })
        );
    }
  }
}
/**
 * Processes a message and sends emails when finished
 * @param {object} params
 */
async function processMessage(params) {
  const { args, state: visualizationStore, timestamp } = params;
  const id = args.id[1];
  const s3 = new AWS.S3();
  const mailer = nodemailer.createTransport(config.email.smtp);
  try {
    // Setup folders
    const directory = path.resolve(config.results.folder, id);
    await fs.promises.mkdir(directory, { recursive: true });
    const rConfig = {
      s3Data: config.data.s3,
      bucket: config.data.bucket,
      localData: path.resolve(config.data.localData),
      wd: path.resolve(config.results.folder),
    };
    // python extraction
    const start = new Date().getTime();
    await downloadS3(id, directory);
    const { stdout, stderr, projectPath } = await profilerExtraction(args);
    // logger.debug('stdout:' + stdout);
    // logger.debug('stderr:' + stderr);
    const matrixPath = path.join(directory, 'results/matrix_files_list.txt');
    const svgPath = path.join(directory, 'results/svg_files_list.txt');
    if (!fs.existsSync(svgPath)) {
      logger.error(stdout);
      logger.error(stderr);
      let error = new Error(
        'profilerExtraction: An Error Occured While Extracting Profiles'
      );
      error.stdout = stdout;
      error.stderr = stderr;
      throw error;
    }
    if (!fs.existsSync(matrixPath))
      throw `matrix file lists does not exist at ${matrixPath}`;
    if (!fs.existsSync(svgPath))
      throw `svg file list does not exist at ${svgPath}`;
    let matrixList = await parseCSV(matrixPath);
    let svgList = await parseCSV(svgPath);
    svgList = svgList.map(({ Profile_Type, Matrix_Size, Path, ...e }) => ({
      ...e,
      profileType: Profile_Type,
      matrixSize: Matrix_Size,
      Path: getRelativePath({ Path: Path }, id).Path,
    }));
    matrixList = matrixList.map(
      ({ Profile_Type, Matrix_Size, Path, ...e }) => ({
        ...e,
        profileType: Profile_Type,
        matrixSize: Matrix_Size,
        Path: getRelativePath({ Path: Path }, id).Path,
      })
    );
    let newState = {
      ...visualizationStore,
      state: {
        ...visualizationStore.state,
        submitted: true,
        id,
      },
    };
    // get dropdown options
    // Mutational Profiles
    const nameOptions = [
      ...new Set(svgList.map(({ Sample_Name }) => Sample_Name)),
    ];
    const selectName = nameOptions[0];
    const filteredPlots = svgList.filter(
      (row) => row.Sample_Name == selectName
    );
    const filteredProfileOptions = [
      ...new Set(filteredPlots.map(({ profileType }) => profileType)),
    ];
    const profile = defaultProfile(filteredProfileOptions);
    const filteredMatrixOptions = [
      ...new Set(
        filteredPlots
          .filter((row) => row.profileType == profile)
          .map(({ matrixSize }) => matrixSize)
      ),
    ].sort((a, b) => a - b);
    const matrix = defaultMatrix(profile, filteredMatrixOptions);
    const filteredFilterOptions = [
      ...new Set(
        filteredPlots
          .filter((row) => row.matrixSize == matrix)
          .map(({ Filter }) => Filter)
      ),
    ];
    const filteredMatrixList = [
      ...new Set(
        matrixList
          .filter((row) => row.profileType == profile)
          .map(({ matrixSize }) => matrixSize)
      ),
    ];
    newState = {
      ...newState,
      state: { ...newState.state, profileOptions: filteredFilterOptions },
    };
    // Cosine Similarity - Profile Comparison - PCA - Kataegis
    const sampleNameOptions = [
      ...new Set(
        svgList.map((row) => {
          if (row.Filter != 'NA') return `${row.Sample_Name}@${row.Filter}`;
          else return row.Sample_Name;
        })
      ),
    ];
    const profileOptions = [...new Set(svgList.map((row) => row.profileType))];
    const selectProfile = defaultProfile(profileOptions);
    const selectMatrix = defaultMatrix(selectProfile, filteredMatrixOptions);
    const refSignatureSetOptions = await r(
      'services/R/visualizeWrapper.R',
      'getReferenceSignatureSets',
      {
        args: { profileType: selectProfile },
        config,
      }
    );
    newState = {
      ...newState,
      cosineSimilarity: {
        ...newState.cosineSimilarity,
        withinProfileType: selectProfile,
        refProfileType: selectProfile,
        refSignatureSet: refSignatureSetOptions[0],
        refSignatureSetOptions: refSignatureSetOptions,
        withinMatrixSize: selectMatrix,
        withinMatrixOptions: filteredMatrixList,
        userProfileType: selectProfile,
        userMatrixSize: selectMatrix,
        userMatrixOptions: filteredMatrixOptions,
      },
      profileComparison: {
        ...newState.profileComparison,
        withinProfileType: selectProfile,
        withinSampleName1: sampleNameOptions[0],
        withinSampleName2: sampleNameOptions[1],
        sampleOptions: sampleNameOptions,
        refProfileType: selectProfile,
        refSampleName: sampleNameOptions[0],
        refSignatureSet: refSignatureSetOptions[0],
        refSignatureSetOptions: refSignatureSetOptions,
        userProfileType: selectProfile,
        userMatrixSize: selectMatrix,
        userMatrixOptions: filteredMatrixOptions,
        userSampleName: sampleNameOptions[0],
      },
      pca: {
        ...newState.pca,
        profileType: selectProfile,
        signatureSet: refSignatureSetOptions[0],
        signatureSetOptions: refSignatureSetOptions,
        userProfileType: selectProfile,
        userMatrixSize: selectMatrix,
        userMatrixOptions: filteredMatrixOptions,
      },
      kataegis: {
        ...newState.kataegis,
        sample: sampleNameOptions[0],
        sampleOptions: sampleNameOptions,
      },
    };
    // R calculations
    // Profiler Summary
    try {
      const profilerSummaryPath = path.join(id, 'results/profilerSummary/');
      await fs.promises.mkdir(path.join(config.wd, profilerSummaryPath), {
        recursive: true,
      });
      const profilerSummary = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'profilerSummary',
          args: {
            matrixList: JSON.stringify(matrixList),
          },
          config: {
            ...rConfig,
            savePath: profilerSummaryPath,
          },
        }
      );
      const { stdout, output } = JSON.parse(profilerSummary);
      if (output.plotPath) {
        newState = {
          ...newState,
          profilerSummary: {
            ...newState.profilerSummary,
            plotPath: output.plotPath,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // cosinse similarity within
    try {
      const csWithinPath = path.join(id, 'results/cosineSimilarityWithin/');
      await fs.promises.mkdir(path.join(config.wd, csWithinPath), {
        recursive: true,
      });
      const csWithin = await r('services/R/visualizeWrapper.R', 'wrapper', {
        fn: 'cosineSimilarityWithin',
        args: {
          matrixFile: matrixList.filter(
            (row) =>
              row.profileType == selectProfile && row.matrixSize == selectMatrix
          )[0].Path,
        },
        config: {
          ...rConfig,
          savePath: csWithinPath,
        },
      });
      let { output, stdout } = JSON.parse(csWithin);
      if (output.plotPath) {
        newState = {
          ...newState,
          cosineSimilarity: {
            ...newState.cosineSimilarity,
            withinPlotPath: output.plotPath,
            withinTxtPath: output.txtPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // cosinse similarity ref
    try {
      const csRefPath = path.join(id, 'results/cosineSimilarityRefSig/');
      await fs.promises.mkdir(path.join(config.wd, csRefPath), {
        recursive: true,
      });
      const csRef = await r('services/R/visualizeWrapper.R', 'wrapper', {
        fn: 'cosineSimilarityRefSig',
        args: {
          profileType: selectProfile,
          signatureSet: refSignatureSetOptions[0],
          matrixFile: matrixList.filter(
            (row) =>
              row.profileType == selectProfile &&
              row.matrixSize == defaultMatrix(selectProfile, ['96', '78', '83'])
          )[0].Path,
        },
        config: {
          ...rConfig,
          savePath: csRefPath,
        },
      });
      let { output, stdout } = JSON.parse(csRef);
      if (output.plotPath) {
        newState = {
          ...newState,
          cosineSimilarity: {
            ...newState.cosineSimilarity,
            refPlotPath: output.plotPath,
            refTxtPath: output.txtPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // cosinse similarity public
    try {
      const csPublicPath = path.join(id, 'results/cosineSimilarityPublic/');
      await fs.promises.mkdir(path.join(config.wd, csPublicPath), {
        recursive: true,
      });
      const csPublic = await r('services/R/visualizeWrapper.R', 'wrapper', {
        fn: 'cosineSimilarityPublic',
        args: {
          matrixFile: matrixList.filter(
            (row) =>
              row.profileType == selectProfile && row.matrixSize == selectMatrix
          )[0].Path,
          study: 'PCAWG',
          cancerType: 'Lung-AdenoCA',
          profileName: selectProfile + selectMatrix,
        },
        config: {
          ...rConfig,
          savePath: csPublicPath,
        },
      });
      let { output, stdout } = JSON.parse(csPublic);
      if (output.plotPath) {
        newState = {
          ...newState,
          cosineSimilarity: {
            ...newState.cosineSimilarity,
            pubPlotPath: output.plotPath,
            pubTxtPath: output.txtPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // mutational pattern
    try {
      const mpPath = path.join(id, 'results/mutationalPattern/');
      await fs.promises.mkdir(path.join(config.wd, mpPath), {
        recursive: true,
      });
      const mutationalPattern = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'mutationalPattern',
          args: {
            matrixFile: matrixList.filter(
              (row) => row.profileType == 'SBS' && row.matrixSize == '96'
            )[0].Path,
            proportion: parseFloat(0.8),
            pattern: 'NCG>NTG',
          },
          config: {
            ...rConfig,
            savePath: mpPath,
          },
        }
      );
      let { output, stdout } = JSON.parse(mutationalPattern);
      if (output.plotPath) {
        newState = {
          ...newState,
          mutationalPattern: {
            ...newState.mutationalPattern,
            ...output,
            debugR: stdout,
          },
        };
      } else {
        logger.dbeug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // profile comparison within
    try {
      const pcWithinPath = path.join(id, 'results/profileComparisonWithin/');
      await fs.promises.mkdir(path.join(config.wd, pcWithinPath), {
        recursive: true,
      });
      const profileComparisonWithin = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'profileComparisonWithin',
          args: {
            profileType: selectProfile,
            sampleName1: sampleNameOptions[0],
            sampleName2: sampleNameOptions[1],
            matrixFile: matrixList.filter(
              (row) =>
                row.profileType == selectProfile &&
                row.matrixSize ==
                  defaultMatrix(selectProfile, ['96', '78', '83'])
            )[0].Path,
          },
          config: {
            ...rConfig,
            savePath: pcWithinPath,
          },
        }
      );
      let { output, stdout } = JSON.parse(profileComparisonWithin);
      if (output.plotPath) {
        newState = {
          ...newState,
          profileComparison: {
            ...newState.profileComparison,
            withinPlotPath: output.plotPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // profile comparison ref
    try {
      const pcRefPath = path.join(id, 'results/profileComparisonRefSig/');
      await fs.promises.mkdir(path.join(config.wd, pcRefPath), {
        recursive: true,
      });
      const refCompare = await r(
        'services/R/visualizeWrapper.R',
        'getSignaturesR',
        {
          args: {
            profileType: selectProfile,
            signatureSetName: refSignatureSetOptions[0],
          },
          config,
        }
      );
      const profileComparisonRefSig = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'profileComparisonRefSig',
          args: {
            profileType: selectProfile,
            sampleName: sampleNameOptions[0],
            signatureSet: refSignatureSetOptions[0],
            compare: refCompare[0],
            matrixFile: matrixList.filter(
              (row) =>
                row.profileType == selectProfile &&
                row.matrixSize ==
                  defaultMatrix(selectProfile, ['96', '78', '83'])
            )[0].Path,
          },
          config: {
            ...rConfig,
            savePath: pcRefPath,
          },
        }
      );
      let { output, stdout } = JSON.parse(profileComparisonRefSig);
      if (output.plotPath) {
        newState = {
          ...newState,
          profileComparison: {
            ...newState.profileComparison,
            refSignatures: refCompare,
            refCompare: refCompare[0],
            refPlotPath: output.plotPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // profile comparison public
    try {
      const pcPubPath = path.join(id, 'results/profileComparisonPublic/');
      await fs.promises.mkdir(path.join(config.wd, pcPubPath), {
        recursive: true,
      });
      const profileComparisonPublic = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'profileComparisonPublic',
          args: {
            profileName: selectProfile + selectMatrix,
            matrixFile: matrixList.filter(
              (row) =>
                row.profileType == selectProfile &&
                row.matrixSize == selectMatrix
            )[0].Path,
            userSample: sampleNameOptions[0],
            study: 'PCAWG',
            cancerType: 'Lung-AdenoCA',
            publicSample: 'SP53073',
          },
          config: {
            ...rConfig,
            savePath: pcPubPath,
          },
        }
      );
      let { output, stdout } = JSON.parse(profileComparisonPublic);
      if (output.plotPath) {
        newState = {
          ...newState,
          profileComparison: {
            ...newState.profileComparison,
            pubPlotPath: output.plotPath,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // pca
    try {
      const pcaPath = path.join(id, 'results/pca/');
      await fs.promises.mkdir(path.join(config.wd, pcaPath), {
        recursive: true,
      });
      const pca = await r('services/R/visualizeWrapper.R', 'wrapper', {
        fn: 'pca',
        args: {
          profileType: selectProfile,
          signatureSet: refSignatureSetOptions[0],
          matrixFile: matrixList.filter(
            (row) =>
              row.profileType == selectProfile &&
              row.matrixSize == defaultMatrix(selectProfile, ['96', '78', '83'])
          )[0].Path,
        },
        config: {
          ...rConfig,
          savePath: pcaPath,
        },
      });
      let { output, stdout } = JSON.parse(pca);
      if (output.pca1) {
        newState = {
          ...newState,
          pca: {
            ...newState.pca,
            ...output,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // pca public
    try {
      const pcaPublicPath = path.join(id, 'results/pcaWithPublic/');
      await fs.promises.mkdir(path.join(config.wd, pcaPublicPath), {
        recursive: true,
      });
      const pcaWithPublic = await r(
        'services/R/visualizeWrapper.R',
        'wrapper',
        {
          fn: 'pcaWithPublic',
          args: {
            matrixFile: matrixList.filter(
              (row) =>
                row.profileType == selectProfile &&
                row.matrixSize == selectMatrix
            )[0].Path,
            study: 'PCAWG',
            cancerType: 'Lung-AdenoCA',
            profileName: selectProfile + selectMatrix,
          },
          config: {
            ...rConfig,
            savePath: pcaPublicPath,
          },
        }
      );
      let { output, stdout } = JSON.parse(pcaWithPublic);
      if (output.pca1) {
        newState = {
          ...newState,
          pca: {
            ...newState.pca,
            pubPca1: output.pca1,
            pubPca2: output.pca2,
            pubPca3: output.pca3,
            pubPca2Data: output.pca2Data,
            pubPca3Data: output.pca3Data,
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    // kataegis
    try {
      const kataegisPath = path.join(id, 'results/kataegis/');
      await fs.promises.mkdir(path.join(config.wd, kataegisPath), {
        recursive: true,
      });
      const kataegis = await r('services/R/visualizeWrapper.R', 'wrapper', {
        fn: 'kataegis',
        args: {
          sample: sampleNameOptions[0],
          highlight: false,
          min: parseInt(5),
          max: parseInt(100),
          chromosome: 'None',
        },
        config: {
          ...rConfig,
          savePath: kataegisPath,
        },
      });
      let { output, stdout } = JSON.parse(kataegis);
      if (output.plotPath) {
        newState = {
          ...newState,
          kataegis: {
            ...newState.kataegis,
            plotPath: output.plotPath,
            txtPath: output.txtPath,
            kataegisData: JSON.parse(output.data),
            debugR: stdout,
          },
        };
      } else {
        logger.debug(stdout);
        logger.error(output);
      }
    } catch (err) {
      logger.error(err);
    }
    const end = new Date().getTime();
    const time = end - start;
    const minutes = Math.floor(time / 60000);
    var seconds = ((time % 60000) / 1000).toFixed(0);
    const runtime = (minutes > 0 ? minutes + ' min ' : '') + seconds + ' secs';
    // upload parameters
    await s3
      .upload({
        Body: JSON.stringify({ ...params, visualization: newState }),
        Bucket: config.queue.bucket,
        Key: `${config.queue.outputKeyPrefix}${id}/params.json`,
      })
      .promise();
    // upload archived project directory
    await s3
      .upload({
        Body: tar
          .c({ sync: true, gzip: true, C: config.results.folder }, [id])
          .read(),
        Bucket: config.queue.bucket,
        Key: `${config.queue.outputKeyPrefix}${id}/${id}.tgz`,
      })
      .promise();
    // specify email template variables
    const templateData = {
      jobName: 'mSigPortal',
      originalTimestamp: timestamp,
      runTime: runtime,
      resultsUrl: `${config.email.baseUrl}/#/visualization/queue/${id}`,
      supportEmail: config.email.adminSupport,
    };
    // send user success email
    logger.info(`Sending user success email`);
    const userEmailResults = await mailer.sendMail({
      from: config.email.adminSupport,
      to: visualizationStore.state.email,
      subject: `mSigPortal Results - ${timestamp} EST`,
      html: await readTemplate(
        __dirname + '/templates/user-success-email.html',
        templateData
      ),
    });
    return true;
  } catch (err) {
    logger.error(err);
    const stdout = err.stdout ? err.stdout.toString() : '';
    const stderr = err.stderr ? err.stderr.toString() : '';
    // template variables
    const templateData = {
      id,
      parameters: JSON.stringify(args, null, 4),
      originalTimestamp: timestamp,
      exception: err.toString(),
      processOutput: !stdout && !stderr ? null : stdout + stderr,
      supportEmail: config.email.adminSupport,
    };
    // send techSupport error email
    logger.error(`Sending techSupport error email`);
    const adminEmailResults = await mailer.sendMail({
      from: config.email.adminSupport,
      to: config.email.techSupport,
      subject: `mSigPortal Error: ${id} - ${timestamp} EST`,
      html: await readTemplate(
        __dirname + '/templates/admin-failure-email.html',
        templateData
      ),
    });
    // send user error email
    if (visualizationStore.state.email) {
      logger.error(`Sending user error email`);
      const userEmailResults = await mailer.sendMail({
        from: config.email.adminSupport,
        to: visualizationStore.state.email,
        subject: 'mSigPortal Error',
        html: await readTemplate(
          __dirname + '/templates/user-failure-email.html',
          templateData
        ),
      });
    }
    return false;
  }
}
/**
 * Receives messages from the queue at regular intervals,
 * specified by config.pollInterval
 */
async function receiveMessage() {
  const sqs = new AWS.SQS();
  try {
    // to simplify running multiple workers in parallel,
    // fetch one message at a time
    const { QueueUrl } = await sqs
      .getQueueUrl({ QueueName: config.queue.url })
      .promise();
    const data = await sqs
      .receiveMessage({
        QueueUrl: QueueUrl,
        MaxNumberOfMessages: 1,
        VisibilityTimeout: config.queue.visibilityTimeout,
        WaitTimeSeconds: 20,
      })
      .promise();
    if (data.Messages && data.Messages.length > 0) {
      const message = data.Messages[0];
      const params = JSON.parse(message.Body);
      logger.info(`Received Message ${params.args.id[1]}`);
      // logger.debug(message.Body);
      // while processing is not complete, update the message's visibilityTimeout
      const intervalId = setInterval((_) => {
        logger.debug('refreshing visibility timeout');
        sqs
          .changeMessageVisibility({
            QueueUrl: QueueUrl,
            ReceiptHandle: message.ReceiptHandle,
            VisibilityTimeout: config.queue.visibilityTimeout,
          })
          .send();
      }, 1000 * (config.queue.visibilityTimeout - 5));
      // processMessage should return a boolean status indicating success or failure
      const status = await processMessage(params);
      clearInterval(intervalId);
      // if message was not processed successfully, send it to the
      // error queue (add metadata in future if needed)
      //   if (!status && config.queue.errorUrl) {
      //     // generate new unique id for error message
      //     const id = crypto.randomBytes(16).toString('hex');
      //     await sqs
      //       .sendMessage({
      //         QueueUrl: config.queue.errorUrl,
      //         MessageDeduplicationId: id,
      //         MessageGroupId: id,
      //         MessageBody: JSON.stringify(params),
      //       })
      //       .promise();
      //   }
      // remove original message from queue once processed
      await sqs
        .deleteMessage({
          QueueUrl: QueueUrl,
          ReceiptHandle: message.ReceiptHandle,
        })
        .promise();
    }
  } catch (e) {
    // catch exceptions related to sqs
    logger.error(e);
  } finally {
    // schedule receiving next message
    setTimeout(receiveMessage, 1000 * (config.queue.pollInterval || 60));
  }
}
