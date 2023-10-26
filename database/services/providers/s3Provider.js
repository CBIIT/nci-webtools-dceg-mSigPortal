import {
  GetObjectAttributesCommand,
  GetObjectCommand,
  paginateListObjectsV2,
} from "@aws-sdk/client-s3";

export class S3Provider {
  constructor(client, basePath) {
    this.type = "s3";
    this.client = client;
    this.basePath = basePath;
    this.listFiles = this.listFiles.bind(this);
    this.readFile = this.readFile.bind(this);
    this.readFileMetadata = this.readFileMetadata.bind(this);
    this.parseS3Path = this.parseS3Path.bind(this);
    this.isAbsoluteS3Path = this.isAbsoluteS3Path.bind(this);
  }

  async listFiles(path) {
    const targetPath = this.isAbsoluteS3Path(path) ? path : this.basePath + path;
    const { bucket, key } = this.parseS3Path(targetPath);
    let keys = [];

    const paginator = paginateListObjectsV2({ client: this.client }, { Bucket: bucket, Prefix: key });

    for await (const page of paginator) {
      const items = (page.Contents ?? []).map((content) => `s3://${bucket}/${content.Key}`);
      keys = keys.concat(items);
    }

    return keys.filter((key) => !key.endsWith("/"));
  }

  async readFile(path) {
    const targetPath = this.isAbsoluteS3Path(path) ? path : this.basePath + path;
    const { bucket, key } = this.parseS3Path(targetPath);
    const s3Response = await this.client.send(
      new GetObjectCommand({
        Bucket: bucket,
        Key: key,
      })
    );
    return s3Response.Body;
  }

  async readFileMetadata(path) {
    const targetPath = this.isAbsoluteS3Path(path) ? path : this.basePath + path;
    const { bucket, key } = this.parseS3Path(targetPath);
    const s3Response = await this.client.send(
      new GetObjectAttributesCommand({
        Bucket: bucket,
        Key: key,
      })
    );
    return s3Response;
  }

  parseS3Path(path) {
    const delimiter = "/";
    const parts = path.replace("s3://", "").split(delimiter);
    const bucket = parts.shift();
    const key = parts.join(delimiter);
    return { bucket, key };
  }

  isAbsoluteS3Path(path) {
    return /^s3:\/\//i.test(path);
  }
}
