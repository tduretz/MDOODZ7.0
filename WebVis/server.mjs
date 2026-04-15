import { createServer } from 'node:http';
import { readFile, readdir, stat } from 'node:fs/promises';
import { createGzip } from 'node:zlib';
import { pipeline } from 'node:stream/promises';
import { Readable } from 'node:stream';
import { join, extname, resolve, normalize } from 'node:path';
import { fileURLToPath } from 'node:url';

import { readMetadata, extractField, datasetExists, readAllMetadata } from './hdf5-reader.mjs';
import { getField, listAvailableFields, getFieldLabel, getFieldUnit } from './fields.mjs';
import { getCached, writeCache } from './field-cache.mjs';

// ── Configuration ──────────────────────────────────────────────────────
const DATA_DIR = resolve(process.argv[2] || '.');
const PORT     = parseInt(process.env.PORT || '3000', 10);
const __dirname = fileURLToPath(new URL('.', import.meta.url));
const PUBLIC    = join(__dirname, 'public');

// ── MIME types ─────────────────────────────────────────────────────────
const MIME = {
  '.html': 'text/html; charset=utf-8',
  '.css':  'text/css; charset=utf-8',
  '.mjs':  'application/javascript; charset=utf-8',
  '.js':   'application/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.png':  'image/png',
  '.svg':  'image/svg+xml',
  '.ico':  'image/x-icon',
};

// ── Validation ─────────────────────────────────────────────────────────
const FILENAME_RE = /^Output\d+\.gzip\.h5$/;

function isPathSafe(p) {
  return !p.includes('..') && !p.includes('\0');
}

// ── Helpers ────────────────────────────────────────────────────────────
function jsonResponse(res, statusCode, data) {
  const body = JSON.stringify(data);
  res.writeHead(statusCode, { 'Content-Type': 'application/json; charset=utf-8' });
  res.end(body);
}

async function gzipJsonResponse(req, res, statusCode, data) {
  const body = JSON.stringify(data);
  const acceptGzip = (req.headers['accept-encoding'] || '').includes('gzip');
  res.setHeader('Content-Type', 'application/json; charset=utf-8');
  res.setHeader('Access-Control-Allow-Origin', '*');
  if (acceptGzip) {
    res.writeHead(statusCode, { 'Content-Encoding': 'gzip' });
    await pipeline(Readable.from(body), createGzip(), res);
  } else {
    res.writeHead(statusCode);
    res.end(body);
  }
}

// ── Route handlers ─────────────────────────────────────────────────────
let cachedFileList = null; // in-memory cache for file-list metadata

async function handleApiFiles(req, res) {
  if (!cachedFileList) {
    const entries = await readdir(DATA_DIR);
    const filenames = entries.filter(f => FILENAME_RE.test(f)).sort();
    cachedFileList = await readAllMetadata(DATA_DIR, filenames);
  }
  await gzipJsonResponse(req, res, 200, { files: cachedFileList });
}

async function handleApiFields(req, res, filename) {
  if (!FILENAME_RE.test(filename)) {
    jsonResponse(res, 400, { error: 'Invalid filename' });
    return;
  }
  const filePath = join(DATA_DIR, filename);
  try {
    await stat(filePath);
  } catch {
    jsonResponse(res, 404, { error: 'File not found' });
    return;
  }
  const params = await readMetadata(filePath);
  const fieldNames = await listAvailableFields(filePath, datasetExists);
  const fields = fieldNames.map(name => ({
    name,
    label: getFieldLabel(name),
    formattedUnit: getFieldUnit(name),
  }));
  await gzipJsonResponse(req, res, 200, { fields, params });
}

async function handleApiFieldData(req, res, filename, fieldName) {
  if (!FILENAME_RE.test(filename)) {
    jsonResponse(res, 400, { error: 'Invalid filename' });
    return;
  }
  const fieldDef = getField(fieldName);
  if (!fieldDef) {
    jsonResponse(res, 400, { error: `Unknown field: ${fieldName}` });
    return;
  }
  const filePath = join(DATA_DIR, filename);
  try {
    await stat(filePath);
  } catch {
    jsonResponse(res, 404, { error: 'File not found' });
    return;
  }

  // Try disk cache first
  const cached = await getCached(DATA_DIR, filename, fieldName);
  if (cached) {
    const acceptGzip = (req.headers['accept-encoding'] || '').includes('gzip');
    res.setHeader('Content-Type', 'application/json; charset=utf-8');
    res.setHeader('Access-Control-Allow-Origin', '*');
    if (acceptGzip) {
      res.writeHead(200, { 'Content-Encoding': 'gzip' });
      await pipeline(Readable.from(cached), createGzip(), res);
    } else {
      res.writeHead(200);
      res.end(cached);
    }
    return;
  }

  const data = await extractField(filePath, fieldDef);
  const json = JSON.stringify(data);
  // Write cache asynchronously — don't block response
  writeCache(DATA_DIR, filename, fieldName, json).catch(() => {});
  await gzipJsonResponse(req, res, 200, data);
}

async function handleApiParams(req, res, filename) {
  if (!FILENAME_RE.test(filename)) {
    jsonResponse(res, 400, { error: 'Invalid filename' });
    return;
  }
  // Try to serve from in-memory metadata cache first
  if (cachedFileList) {
    const entry = cachedFileList.find(f => f.name === filename);
    if (entry) {
      await gzipJsonResponse(req, res, 200, {
        step: entry.step,
        time: entry.time,
        Nx: entry.nx ?? null,
        Nz: entry.nz ?? null,
        dt: entry.dt ?? null,
      });
      return;
    }
  }
  // Fallback: read from HDF5 file
  const filePath = join(DATA_DIR, filename);
  try {
    await stat(filePath);
  } catch {
    jsonResponse(res, 404, { error: 'File not found' });
    return;
  }
  const params = await readMetadata(filePath);
  await gzipJsonResponse(req, res, 200, {
    step: parseInt(filename.match(/\d+/)?.[0] || '0', 10),
    time: params.time,
    Nx: params.nx,
    Nz: params.nz,
    dt: params.dt,
  });
}

async function handleStatic(req, res, urlPath) {
  let relPath = urlPath === '/' ? '/index.html' : urlPath;
  if (!isPathSafe(relPath)) {
    res.writeHead(400);
    res.end('Bad request');
    return;
  }
  const absPath = normalize(join(PUBLIC, relPath));
  // Ensure resolved path is still under PUBLIC
  if (!absPath.startsWith(PUBLIC)) {
    res.writeHead(400);
    res.end('Bad request');
    return;
  }
  try {
    const content = await readFile(absPath);
    const ext = extname(absPath);
    const mime = MIME[ext] || 'application/octet-stream';
    res.writeHead(200, { 'Content-Type': mime });
    res.end(content);
  } catch {
    res.writeHead(404);
    res.end('Not found');
  }
}

// ── Router ─────────────────────────────────────────────────────────────
const server = createServer(async (req, res) => {
  const url = new URL(req.url, `http://${req.headers.host}`);
  const path = decodeURIComponent(url.pathname);

  if (!isPathSafe(path)) {
    res.writeHead(400);
    res.end('Bad request');
    return;
  }

  try {
    if (path === '/api/files') {
      await handleApiFiles(req, res);
    } else if (path.startsWith('/api/params/')) {
      const filename = path.slice('/api/params/'.length);
      await handleApiParams(req, res, filename);
    } else if (path.startsWith('/api/fields/')) {
      const filename = path.slice('/api/fields/'.length);
      await handleApiFields(req, res, filename);
    } else if (path.startsWith('/api/field-data/')) {
      const rest = path.slice('/api/field-data/'.length);
      const slashIdx = rest.indexOf('/');
      if (slashIdx === -1) {
        jsonResponse(res, 400, { error: 'Missing field name in URL' });
        return;
      }
      const filename = rest.slice(0, slashIdx);
      const fieldName = rest.slice(slashIdx + 1);
      await handleApiFieldData(req, res, filename, fieldName);
    } else {
      await handleStatic(req, res, path);
    }
  } catch (err) {
    console.error('Server error:', err);
    if (!res.headersSent) {
      jsonResponse(res, 500, { error: 'Internal server error' });
    }
  }
});

server.listen(PORT, () => {
  console.log(`MDOODZ WebVis server`);
  console.log(`  Data dir: ${DATA_DIR}`);
  console.log(`  URL:      http://localhost:${PORT}`);
});
