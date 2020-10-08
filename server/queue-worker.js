const fs = require('fs');
const path = require('path');
const AWS = require('aws-sdk');
const nodemailer = require('nodemailer');
const r = require('r-wrapper').async;
const tar = require('tar');
const config = require('./config.json');
const logger = require('./logger');
const { profilerExtraction } = require('./controllers');

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
 * Processes a message and sends emails when finished
 * @param {object} params
 */
async function processMessage({ args, email, timestamp }) {
  const id = args.projectID[1];
  const s3 = new AWS.S3();
  const mailer = nodemailer.createTransport(config.email.smtp);

  try {
    // get calculation results
    const directory = path.resolve(config.results.folder, id);
    await fs.promises.mkdir(directory, { recursive: true });

    const start = new Date().getTime();
    const { stdout, stderr, projectPath } = await profilerExtraction(args);
    const end = new Date().getTime();

    const time = end - start;
    const minutes = Math.floor(time / 60000);
    var seconds = ((time % 60000) / 1000).toFixed(0);

    const runtime = (minutes > 0 ? minutes + ' min ' : '') + seconds + ' secs';

    // upload archived project directory
    await s3
      .upload({
        Body: tar.c({ gzip: true, C: config.results.folder }, [id]),
        Bucket: config.s3.bucket,
        Key: `${config.aws.key}${id}.tgz`,
      })
      .promise();

    // specify email template variables
    const templateData = {
      jobName: 'mSigPortal',
      originalTimestamp: timestamp,
      runTime: runtime,
      resultsUrl: `${config.email.baseUrl}/mutational-signatures/${id}`,
    };

    // send user success email
    logger.info(`Sending user success email`);
    const userEmailResults = await mailer.sendMail({
      from: config.email.sender,
      to: email,
      subject: 'mSigPortal Results - ' + timestamp,
      html: await readTemplate(
        __dirname + '/templates/user-success-email.html',
        templateData
      ),
    });

    return true;
  } catch (e) {
    logger.error(e);

    // template variables
    const templateData = {
      id: id,
      parameters: JSON.stringify(args, null, 4),
      originalTimestamp: timestamp,
      exception: e.toString(),
      processOutput: e.stdout ? e.stdout.toString() : null,
      supportEmail: config.email.admin,
    };

    // send admin error email
    logger.info(`Sending admin error email`);
    const adminEmailResults = await mailer.sendMail({
      from: config.email.sender,
      to: config.email.admin,
      subject: `mSigPortal Error: ${id}`, // searchable calculation error subject
      html: await readTemplate(
        __dirname + '/templates/admin-failure-email.html',
        templateData
      ),
    });

    // send user error email
    if (email) {
      logger.info(`Sending user error email`);
      const userEmailResults = await mailer.sendMail({
        from: config.email.sender,
        to: email,
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

      logger.info(`Received Message : ${message.Body}`);

      // while processing is not complete, update the message's visibilityTimeout
      const intervalId = setInterval(
        (_) =>
          sqs
            .changeMessageVisibility({
              QueueUrl: QueueUrl,
              ReceiptHandle: message.ReceiptHandle,
              VisibilityTimeout: config.queue.visibilityTimeout,
            })
            .send(),
        1000 * (config.queue.visibilityTimeout - 1)
      );

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