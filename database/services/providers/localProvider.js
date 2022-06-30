import { resolve } from "path";
import { createReadStream } from "fs";
import { stat, readdir } from "fs/promises";

export class LocalProvider {
  constructor(basePath) {
    this.type = "local";
    this.basePath = basePath || "";
    this.listFiles = this.listFiles.bind(this);
    this.readFile = this.readFile.bind(this);
    this.readFileMetadata = this.readFileMetadata.bind(this);
  }

  /**
   * Recursively lists all files within a directory.
   * @param {*} path
   * @returns
   */
  async listFiles(path) {
    const targetPath = resolve(this.basePath, path);
    const metadata = await stat(targetPath);
    if (!metadata.isDirectory()) {
      return [];
    }

    let filePaths = [];
    for (const filePath of await readdir(targetPath)) {
      const absolutePath = resolve(targetPath, filePath);
      const metadata = await stat(absolutePath);
      if (metadata.isFile()) {
        filePaths.push(absolutePath);
      } else if (metadata.isDirectory()) {
        const subFilePaths = await this.listFiles(filePath);
        filePaths = filePaths.concat(subFilePaths);
      }
    }
    return filePaths;
  }

  async readFile(path) {
    const targetPath = resolve(this.basePath, path);
    return createReadStream(targetPath);
  }

  async readFileMetadata(path) {
    const targetPath = resolve(this.basePath, path);
    return await stat(targetPath);
  }
}
