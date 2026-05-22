// ── Field Cache ───────────────────────────────────────────────────────
// Disk cache for extracted field JSON. Cache key = sha256(filename:field)[0:16].json
// Invalidation via mtime comparison: cache is stale if source HDF5 is newer.

import { readFile, writeFile, stat, mkdir } from 'node:fs/promises';
import { join } from 'node:path';
import { createHash } from 'node:crypto';

function cacheKey(filename, field) {
  const hash = createHash('sha256').update(`${filename}:${field}`).digest('hex');
  return hash.slice(0, 16) + '.json';
}

function cacheDir(dataDir) {
  return join(dataDir, '.cache');
}

/**
 * Try to read a valid cache entry.
 * @returns {string|null}  Raw JSON string if valid, null if miss/stale.
 */
export async function getCached(dataDir, filename, field) {
  const dir = cacheDir(dataDir);
  const cachePath = join(dir, cacheKey(filename, field));
  const sourcePath = join(dataDir, filename);

  try {
    const [cacheStat, sourceStat] = await Promise.all([stat(cachePath), stat(sourcePath)]);
    if (cacheStat.mtimeMs >= sourceStat.mtimeMs) {
      return await readFile(cachePath, 'utf8');
    }
  } catch (_) {
    // cache file doesn't exist or source missing → miss
  }
  return null;
}

/**
 * Write extracted field JSON to cache.
 * @param {string} dataDir
 * @param {string} filename
 * @param {string} field
 * @param {string} json  Serialised JSON string.
 */
export async function writeCache(dataDir, filename, field, json) {
  const dir = cacheDir(dataDir);
  try {
    await mkdir(dir, { recursive: true });
  } catch (_) { /* directory exists */ }
  const cachePath = join(dir, cacheKey(filename, field));
  await writeFile(cachePath, json, 'utf8');
}
