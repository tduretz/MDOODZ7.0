import { createServer } from 'node:http';
import { readFile, readdir, stat } from 'node:fs/promises';
import { createGzip } from 'node:zlib';
import { pipeline } from 'node:stream/promises';
import { Readable } from 'node:stream';
import { join, extname, resolve, normalize } from 'node:path';
import { fileURLToPath } from 'node:url';

import { readMetadata, extractField, datasetExists, readAllMetadata, extractOverlayData } from './hdf5-reader.mjs';
import { getField, listAvailableFields, getFieldLabel, getFieldUnit } from './fields.mjs';
import { getCached, writeCache } from './field-cache.mjs';

// ── Configuration ──────────────────────────────────────────────────────
const __dirname = fileURLToPath(new URL('.', import.meta.url));
const ROOT_DIR  = resolve(process.argv[2] || join(__dirname, '..'));
const PORT      = parseInt(process.env.PORT || '3000', 10);
const PUBLIC    = join(__dirname, 'public');

// ── Dataset (directory) scanning ───────────────────────────────────────
// Recursively scans ROOT_DIR for directories containing Output*.gzip.h5.
// Dataset names are paths relative to ROOT_DIR.
const SKIP_DIRS = new Set(['node_modules', '.git', '_deps', 'CMakeFiles', '.cache', '__pycache__']);
let datasetList = [];   // [{name, path, fileCount}]
let activeDataset = null; // current dataset path (absolute)

async function scanDatasets() {
  datasetList = [];

  async function walk(dirPath, relPath) {
    let entries;
    try { entries = await readdir(dirPath, { withFileTypes: true }); } catch { return; }

    const h5Files = entries.filter(e => !e.isDirectory() && FILENAME_RE.test(e.name));
    if (h5Files.length > 0) {
      datasetList.push({ name: relPath || '.', path: dirPath, fileCount: h5Files.length });
    }

    for (const entry of entries) {
      if (!entry.isDirectory() || SKIP_DIRS.has(entry.name)) continue;
      const childPath = join(dirPath, entry.name);
      const childRel = relPath ? `${relPath}/${entry.name}` : entry.name;
      await walk(childPath, childRel);
    }
  }

  await walk(ROOT_DIR, '');

  // Sort: '.' first, then alphabetically
  datasetList.sort((a, b) => {
    if (a.name === '.') return -1;
    if (b.name === '.') return 1;
    return a.name.localeCompare(b.name);
  });
  // Default to first dataset
  if (datasetList.length > 0 && !activeDataset) {
    activeDataset = datasetList[0].path;
  }
  // If active dataset was removed, reset to first
  if (activeDataset && !datasetList.find(d => d.path === activeDataset)) {
    activeDataset = datasetList.length > 0 ? datasetList[0].path : null;
  }
}

function getDataDir() {
  return activeDataset || ROOT_DIR;
}

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
let cachedDataDir  = null; // which dataset the cache belongs to

async function handleApiFiles(req, res) {
  const dataDir = getDataDir();
  if (!cachedFileList || cachedDataDir !== dataDir) {
    const entries = await readdir(dataDir);
    const filenames = entries.filter(f => FILENAME_RE.test(f)).sort();
    cachedFileList = await readAllMetadata(dataDir, filenames);
    cachedDataDir  = dataDir;
  }
  await gzipJsonResponse(req, res, 200, { files: cachedFileList });
}

async function handleApiFields(req, res, filename) {
  if (!FILENAME_RE.test(filename)) {
    jsonResponse(res, 400, { error: 'Invalid filename' });
    return;
  }
  const filePath = join(getDataDir(), filename);
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
  const filePath = join(getDataDir(), filename);
  try {
    await stat(filePath);
  } catch {
    jsonResponse(res, 404, { error: 'File not found' });
    return;
  }

  // Try disk cache first
  const dataDir = getDataDir();
  const cached = await getCached(dataDir, filename, fieldName);
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
  writeCache(getDataDir(), filename, fieldName, json).catch(() => {});
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
  const filePath = join(getDataDir(), filename);
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

const VALID_OVERLAY_LAYERS = new Set([
  'phases', 'temperature', 'velocity', 'topo', 'director', 'sigma1', 'edot1', 'melt',
]);

async function handleApiOverlayData(req, res, filename, searchParams) {
  if (!FILENAME_RE.test(filename)) {
    jsonResponse(res, 400, { error: 'Invalid filename' });
    return;
  }
  const filePath = join(getDataDir(), filename);
  try { await stat(filePath); } catch {
    jsonResponse(res, 404, { error: 'File not found' });
    return;
  }
  const raw = searchParams.get('layers') || '';
  const layers = raw.split(',').map(s => s.trim()).filter(s => VALID_OVERLAY_LAYERS.has(s));
  if (layers.length === 0) {
    jsonResponse(res, 400, { error: 'No valid layers requested' });
    return;
  }
  const data = await extractOverlayData(filePath, layers);
  await gzipJsonResponse(req, res, 200, data);
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
    if (path === '/api/datasets') {
      await gzipJsonResponse(req, res, 200, {
        datasets: datasetList.map(d => ({ name: d.name, fileCount: d.fileCount })),
        active: datasetList.find(d => d.path === activeDataset)?.name || null,
      });
    } else if (path === '/api/dataset' && req.method === 'POST') {
      // Read JSON body
      const chunks = [];
      for await (const chunk of req) chunks.push(chunk);
      const body = JSON.parse(Buffer.concat(chunks).toString());
      const target = datasetList.find(d => d.name === body.name);
      if (!target) {
        jsonResponse(res, 404, { error: 'Dataset not found' });
        return;
      }
      activeDataset = target.path;
      cachedFileList = null;  // invalidate file-list cache
      cachedDataDir = null;
      jsonResponse(res, 200, { active: target.name });
    } else if (path === '/api/rescan' && req.method === 'POST') {
      await scanDatasets();
      cachedFileList = null;
      cachedDataDir = null;
      await gzipJsonResponse(req, res, 200, {
        datasets: datasetList.map(d => ({ name: d.name, fileCount: d.fileCount })),
        active: datasetList.find(d => d.path === activeDataset)?.name || null,
      });
    } else if (path === '/api/files') {
      await handleApiFiles(req, res);
    } else if (path.startsWith('/api/params/')) {
      const filename = path.slice('/api/params/'.length);
      await handleApiParams(req, res, filename);
    } else if (path.startsWith('/api/overlay-data/')) {
      const filename = path.slice('/api/overlay-data/'.length);
      await handleApiOverlayData(req, res, filename, url.searchParams);
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

server.listen(PORT, async () => {
  await scanDatasets();
  console.log(`MDOODZ WebVis server`);
  console.log(`  Root dir:  ${ROOT_DIR}`);
  console.log(`  Datasets:  ${datasetList.map(d => d.name).join(', ') || '(none)'}`);
  console.log(`  Active:    ${datasetList.find(d => d.path === activeDataset)?.name || '(none)'}`);
  console.log(`  URL:       http://localhost:${PORT}`);
});
