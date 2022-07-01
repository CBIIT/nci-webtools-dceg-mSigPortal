import {
  DeleteMessageCommand,
  GetQueueUrlCommand,
  ReceiveMessageCommand,
} from "@aws-sdk/client-sqs";

export async function processMessages({
  sqs,
  queueName,
  visibilityTimeout = 60,
  waitTime = 20,
  pollInterval = 60,
  messageHandler,
  messageParser = JSON.parse,
  errorHandler = console.error,
}) {
  const args = {
    sqs,
    queueName,
    visibilityTimeout,
    waitTime,
    pollInterval,
    messageHandler,
    messageParser,
    errorHandler,
  };

  try {
    await processMessage(args);
  } catch (exception) {
    errorHandler(exception);
  } finally {
    setTimeout(() => processMessages(args), pollInterval * 1000);
  }
}

export async function processMessage({
  sqs,
  queueName,
  visibilityTimeout = 60,
  waitTime = 20,
  messageHandler,
  messageParser = JSON.parse,
  errorHandler = console.error,
}) {
  try {
    const { QueueUrl: queueUrl } = await sqs.send(
      new GetQueueUrlCommand({
        QueueName: queueName,
      })
    );

    // to simplify running multiple workers in parallel,
    // fetch one message at a time
    const data = await sqs.send(
      new ReceiveMessageCommand({
        QueueUrl: queueUrl,
        MaxNumberOfMessages: 1,
        VisibilityTimeout: visibilityTimeout,
        WaitTimeSeconds: waitTime,
      })
    );

    if (data.Messages && data.Messages.length > 0) {
      const message = data.Messages[0];
      const messageBody = messageParser(message.Body);

      try {
        // delete message immediately since processing usually takes longer than the message's validity period
        await sqs.send(
          new DeleteMessageCommand({
            QueueUrl: queueUrl,
            ReceiptHandle: message.ReceiptHandle,
          })
        );
        await messageHandler(messageBody);
      } catch (e) {
        await errorHandler(e, message);
      }
    }
  } catch (e) {
    await errorHandler(e);
  }
}
