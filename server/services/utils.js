import { fileURLToPath } from 'url';
import { existsSync } from 'fs';
import { mkdir, writeFile, readFile, copyFile, readdir } from 'fs/promises';
import path from 'path';
import { promisify } from 'util';
import { execFile } from 'child_process';
import template from 'lodash/template.js';
import { pickBy } from 'lodash-es';
import { validate as validateUuid } from 'uuid';

// promisified executeFile
export const execFileAsync = promisify(execFile);

/**
 * Validates that an identifier is safe to use in a filesystem path.
 * Accepts UUIDs or a strict allowlist of filename-safe characters.
 * @param {unknown} id
 * @returns {boolean}
 */
export function isValidId(id) {
  if (typeof id !== 'string' || id.length === 0 || id.length > 255) return false;
  return validateUuid(id) || /^[A-Za-z0-9_-]+$/.test(id);
}

/**
 * Reduces a user-supplied name to a safe basename, rejecting traversal.
 * @param {unknown} name
 * @returns {string} A safe filename
 * @throws {Error} If the name is missing or resolves to a traversal segment
 */
export function sanitizeFilename(name) {
  if (typeof name !== 'string' || name.includes('\0')) {
    throw new Error('Invalid filename');
  }
  const base = path.basename(name);
  if (!base || base === '.' || base === '..') {
    throw new Error('Invalid filename');
  }
  return base;
}

/**
 * Resolves path segments against a base directory and guarantees the result
 * stays within that base, preventing path traversal.
 * @param {string} base Base directory
 * @param {...string} segments Path segments to join
 * @returns {string} The contained, resolved path
 * @throws {Error} If any segment contains a null byte or escapes the base
 */
export function resolveWithin(base, ...segments) {
  const resolvedBase = path.resolve(base);
  for (const segment of segments) {
    if (typeof segment !== 'string' || segment.includes('\0')) {
      throw new Error('Invalid path segment');
    }
  }
  const target = path.resolve(resolvedBase, ...segments);
  if (target !== resolvedBase && !target.startsWith(resolvedBase + path.sep)) {
    throw new Error('Resolved path escapes base directory');
  }
  return target;
}

/**
 * Checks if the current module is the main module.
 * @param {ImportMeta} importMeta
 * @param {NodeJS.ProcessEnv} env
 * @returns
 */
export function isMainModule(importMeta, env = process.env) {
  const mainModulePath = env.pm_exec_path || process.argv[1];
  const currentModulePath = fileURLToPath(importMeta.url);
  return mainModulePath === currentModulePath;
}

/**
 * Creates directories if they don't exist.
 * @param {string[]} dirs
 * @returns
 */
export async function mkdirs(dirs) {
  return await Promise.all(dirs.map((dir) => mkdir(dir, { recursive: true })));
}

/**
 * Writes json to a file.
 * @param {string} filepath
 * @param {any} data
 * @returns {Promise<void>} fulfilled when the file is written
 */
export async function writeJson(filepath, data) {
  if (typeof filepath !== 'string' || filepath.includes('\0')) {
    throw new Error('Invalid filepath');
  }
  return await writeFile(filepath, JSON.stringify(data), 'utf-8');
}

/**
 * Reads json from a file.
 * @param {string} filepath
 * @returns {any} data
 */
export async function readJson(filepath) {
  try {
    if (typeof filepath !== 'string' || filepath.includes('\0')) {
      throw new Error('Invalid filepath');
    }
    const data = await readFile(filepath, 'utf8');
    return JSON.parse(data);
  } catch (e) {
    return null;
  }
}

export async function renderTemplate(filepath, data) {
  if (typeof filepath !== 'string' || filepath.includes('\0')) {
    throw new Error('Invalid filepath');
  }
  const templateContents = await readFile(filepath, 'utf8');
  return template(templateContents)(data);
}

/**
 * Selects the first file which exists from the given list of files.
 * @param {string[]} filePaths
 * @returns string
 */
export function coalesceFilePaths(filePaths) {
  for (const filePath of filePaths) {
    if (typeof filePath !== 'string' || filePath.includes('\0')) continue;
    if (existsSync(filePath)) {
      return filePath;
    }
  }
}

/**
 * Removes the file extension from the given file path.
 * @param {string} filePath
 * @returns {string} The file path without the extension.
 */
export function stripExtension(filePath) {
  if (!filePath) return null;
  const { dir, name } = path.parse(filePath);
  return path.join(dir, name);
}

/**
 * Copies files from the given source directory to the given destination directory.
 * @param {string} source Source directory path
 * @param {string} destination Destination directory path
 * @param {boolean} overwrite Overwrite existing files
 */
export async function copyFiles(source, destination, overwrite = false) {
  const sourceFiles = await readdir(source, { withFileTypes: true });
  for (const file of sourceFiles.filter((f) => f.isFile())) {
    const sourceFilePath = path.resolve(source, file.name);
    const destinationFilePath = path.resolve(destination, file.name);
    if (overwrite || !existsSync(destinationFilePath)) {
      await copyFile(sourceFilePath, destinationFilePath);
    }
  }
}

export function pickNonNullValues(object) {
  return pickBy(object, (v) => v !== null);
}

//
/**
 * async generator for retrieving paths for all files under a given directory
 * @param {string} filePath
 * @returns {AsyncGenerator} async generator to consume
 * consume the generator like so:
 * for await (const f of getFiles(filePath)) {
    if (f) ...
  }
 */
export async function* getFiles(filePath) {
  const dirents = await readdir(filePath, { withFileTypes: true });
  for (const dirent of dirents) {
    const res = path.resolve(filePath, dirent.name);
    if (dirent.isDirectory()) {
      yield* getFiles(res);
    } else {
      yield res;
    }
  }
}
