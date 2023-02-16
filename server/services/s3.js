import {
  S3Client,
  GetObjectCommand,
  PutObjectCommand,
  ListObjectsCommand,
} from '@aws-sdk/client-s3';
import { createWriteStream, createReadStream, readdirSync, statSync } from 'fs';
import path from 'path';

export async function getObject(
  destination,
  key,
  bucket,
  config = { region: 'us-east-1' }
) {
  const s3 = new S3Client(config);
  const params = { Bucket: bucket, Key: key };
  const { Body } = await s3.send(new GetObjectCommand(params));
  return writeStream(destination, Body);
}

export function putObject(file, key, bucket, config = { region: 'us-east-1' }) {
  const s3 = new S3Client(config);
  const params = {
    Bucket: bucket,
    Key: key,
    Body: createReadStream(file),
  };
  return s3.send(new PutObjectCommand(params));
}

export function listObjects(prefix, bucket, config = { region: 'us-east-1' }) {
  const s3 = new S3Client(config);
  return s3.send(new ListObjectsCommand({ Prefix: prefix, Bucket: bucket }));
}

function writeStream(destination, stream) {
  return new Promise((resolve, reject) => {
    const writer = createWriteStream(destination);
    writer.on('finish', () => resolve(destination));
    writer.on('error', reject);
    stream.pipe(writer);
  });
}

export function getDirectory(destination, key, bucket, config = {}) {
  return new Promise(async (resolve, reject) => {
    try {
      const { Contents: files } = await listObjects(key, bucket, config);
      if (!files)
        reject(
          new Error(`${path.join(bucket, key)} is empty or does not exist`)
        );
      const response = await Promise.all(
        files.map((e) => {
          return getObject(
            path.join(destination, e.Key.replace(key, '')),
            e.Key,
            bucket,
            config
          );
        })
      );
      resolve(response);
    } catch (error) {
      reject(error);
    }
  });
}

export function putDirectory(directory, key, bucket, config = {}) {
  return new Promise(async (resolve, reject) => {
    try {
      const files = readdirSync(directory);
      if (!files) reject(new Error(`${directory} is empty or does not exist`));

      const response = await Promise.all(
        files.map((file) => {
          const filePath = path.resolve(directory, file);
          if (statSync(filePath).isFile()) {
            return putObject(filePath, path.join(key, file), bucket, config);
          }
        })
      );
      resolve(response);
    } catch (error) {
      reject(error);
    }
  });
}
