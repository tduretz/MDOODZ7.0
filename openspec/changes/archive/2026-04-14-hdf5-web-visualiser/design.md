## Context

MDOODZ produces `Output*.gzip.h5` files (8–24 MB each, gzip-compressed datasets) containing staggered-grid fields at cell centres, vertices, Vx-nodes, and Vz-nodes. Currently, visualisation requires Julia/Makie or Python/h5py — both heavy to install and non-interactive. The goal is a zero-compilation, npm-only web viewer that reads these files directly and renders fields in the browser canvas.

The staggered grid is the main complexity: different fields live on different grids (centres = `(Nx-1)×(Nz-1)`, vertices = `Nx×Nz`, Vx = `Nx×(Nz+1)`, Vz = `(Nx+1)×Nz`). Derived fields like stress/strain-rate invariants require averaging vertex values to centres. Air cells (phase = -1) must be masked as NaN.

## Goals / Non-Goals

**Goals:**

- Browse a directory of HDF5 output files in the browser
- Select any field from the file and render it as a colour-mapped 2D image
- Support all grid locations (centres, vertices, Vx, Vz, VizGrid phase maps)
- Compute derived fields (τ_II, ε̇_II) server-side
- Sub-second response for a 600×400 grid
- Install with `npm install` only — no C compiler, no system HDF5, no build step for the frontend
- Clean MVC frontend extensible for future overlays

**Non-Goals:**

- 3D visualisation
- Particle file reading (Particles*.gzip.h5) — future extension
- Vector field overlays (velocity quivers, principal stress) — future extension
- Multi-panel comparison layouts — future extension
- Video/animation export — future extension
- Editing or writing HDF5 files
- Authentication or multi-user access

## Decisions

### 1. HDF5 reading: h5wasm (WebAssembly) in-process

**Choice**: Use `h5wasm` (npm package, zero native deps) to read HDF5 files directly in the Node.js server process.

**Alternatives considered**:
- **C++ CLI extractor** (spawn child process per request): adds CMake build step, process spawn latency (~50ms per request), JSON serialisation overhead. More complex deployment.
- **hdf5-io / hdf5lib** (native Node bindings): require `node-gyp`, system HDF5 headers, platform-specific compilation. Fragile on macOS/Linux differences.
- **Python subprocess**: adds Python dependency — exactly what we're trying to avoid.

**Rationale**: h5wasm compiles HDF5 C library to WASM. It reads gzip-compressed datasets natively, returns typed arrays (`Float32Array`, `Float64Array`). No native compilation. Node.js ≥16 with `NODERAWFS=1` gives direct filesystem access. The API (`file.get("/Centers/P").value`) maps cleanly to our use case. ~6k weekly npm downloads, maintained by NIST.

### 2. Server: plain Node.js `http` module, no framework

**Choice**: Single `server.mjs` using `node:http`, `node:fs`, `node:path`. No Express, no Koa.

**Alternatives considered**:
- **Express**: unnecessary abstraction for 3-4 routes. Adds a dependency tree.
- **Fastify**: same concern — overkill for a file browser API.

**Rationale**: The server has exactly 4 routes (`/api/files`, `/api/fields/:file`, `/api/field-data/:file/:field`, and static file serving). Plain `http` handles this trivially. Zero dependencies beyond h5wasm.

### 3. Field extraction: server-side reshape + derive, send JSON to browser

**Choice**: Server reads HDF5, reshapes flat arrays to 2D, computes derived fields, sends JSON `{ xCoords, zCoords, values: [[...], ...], nx, nz, unit, log }` to the browser.

**Alternatives considered**:
- **Browser-side h5wasm** (send raw .h5 file to browser, read in WASM there): transfers 8-24 MB per request vs ~1 MB for a single field. Also requires Emscripten FS mounting in browser — finicky.
- **Binary transfer** (ArrayBuffer): faster than JSON but harder to debug, version, and extend. JSON overhead is acceptable for 600×400 grids (~2 MB for a float field with coordinates).

**Rationale**: Server-side extraction keeps the browser thin. A single field at 600×400 = 240k floats ≈ 1.5 MB JSON (gzip to ~300 KB with HTTP compression). Response time dominated by HDF5 read (~50ms) not JSON serialisation.

### 4. Field registry: declarative field definitions

**Choice**: Define a `fields.mjs` registry mapping logical field names to HDF5 paths, grid types, scaling, units, and optional derivation functions.

```
{
  "Pressure":       { path: "/Centers/P",    grid: "center", unit: "Pa" },
  "Viscosity":      { path: "/Centers/eta_n", grid: "center", unit: "Pa·s", log: true },
  "Stress II":      { derive: "stress_invariant", grid: "center", unit: "Pa" },
  "Strain rate II": { derive: "strainrate_invariant", grid: "center", unit: "1/s", log: true },
  "Temperature":    { path: "/Centers/T",    grid: "center", unit: "K", offset: -273.15 },
  "Phases":         { path: "/VizGrid/compo", grid: "vizgrid", discrete: true },
  ...
}
```

**Rationale**: Adding new fields means adding one entry — no code changes. Derived fields use named derivation functions that combine multiple datasets. This is the extension point for future fields.

### 5. Grid handling: coordinate arrays from HDF5 + dimension table

**Choice**: Read coordinate arrays (`xc_coord`, `zc_coord`, `xg_coord`, `zg_coord`, `zvx_coord`, `xvz_coord`) from `/Model/` group. Use `Params[3]` (Nx) and `Params[4]` (Nz) to derive grid dimensions per grid type.

Dimension lookup:
| Grid type | rows | cols | x-coords | z-coords |
|-----------|------|------|----------|----------|
| `center`  | Nx-1 | Nz-1 | `xc_coord` | `zc_coord` |
| `vertex`  | Nx   | Nz   | `xg_coord` | `zg_coord` |
| `vx`      | Nx   | Nz+1 | `xg_coord` | `zvx_coord` |
| `vz`      | Nx+1 | Nz   | `xvz_coord` | `zg_coord` |
| `vizgrid` | ncx  | ncz  | from `VizGrid/xviz` | from `VizGrid/zviz` |

**Rationale**: This is the canonical mapping from the MDOODZ staggered grid. Reading coordinate arrays from HDF5 (already in SI/metres) avoids reconstruction errors.

### 6. Derived field computation (invariants)

**Choice**: Compute stress and strain-rate second invariants server-side using typed array arithmetic:

```
τ_II = sqrt(0.5 * (τ_xx² + τ_zz² + (τ_xx + τ_zz)²) + τ_xz_c²)
```

where `τ_xz_c` is the average of 4 surrounding vertex values interpolated to each cell centre.

**Rationale**: This matches the Julia/Python reference implementations. The vertex-to-centre averaging pattern (2×2 block average) is the standard MDOODZ staggered-grid interpolation.

### 7. Frontend architecture: ES6 MVC with observer pattern

**Choice**: Three layers, all plain ES6 modules (`.mjs` files served as-is, no bundler):

- **Model** (`model.mjs`): holds app state (current file, field, data, colour range). Emits events via `EventTarget`/`CustomEvent`.
- **Views** (`views/*.mjs`): each view is a class that owns a DOM element, subscribes to model events, re-renders on change. Views: `FieldCanvas`, `ControlPanel`, `ColourBar`, `FileList`, `StatusBar`.
- **Controller** (`controller.mjs`): wires UI events (dropdown change, slider input) to model mutations and API fetches. Injects model into views.

**Alternatives considered**:
- **Lit / Preact**: lightweight frameworks but still a dependency + possible build step.
- **Vanilla DOM manipulation without structure**: leads to spaghetti; MVC gives extension points.

**Rationale**: Observer pattern (`model.addEventListener('field-loaded', ...)`) gives reactive updates without a framework. Each view is self-contained — adding a new overlay view later means adding one file and subscribing to model events.

### 8. Canvas rendering: ImageData pixel-fill with colour map LUT

**Choice**: Render fields by writing pixels directly to a `<canvas>` via `ImageData`. Build a 256-entry colour map lookup table (LUT), map normalised data values to LUT indices, write RGBA bytes.

**Alternatives considered**:
- **WebGL**: faster for huge grids but complex setup for simple 2D heatmaps. Overkill for 600×400.
- **SVG rectangles**: extremely slow for >10k cells.
- **Server-side image generation** (sharp/canvas): adds native dependency, moves rendering off the client.

**Rationale**: `ImageData` pixel-fill for 240k cells takes <10ms. Colour maps (viridis, turbo, inferno, etc.) stored as 256×3 arrays. The LUT approach is cache-friendly and trivially supports log scaling (take log10 before normalisation).

### 9. Colour maps: bundled as JSON arrays

**Choice**: Ship 5-6 common scientific colour maps (viridis, turbo, inferno, plasma, coolwarm, discrete phase palette) as JSON arrays of `[r, g, b]` triplets (256 entries each, ~5 KB total). Bundled in `public/colormaps/`.

**Rationale**: No external dependency. Same maps used by gnuplot/matplotlib. Easy to add more.

### 10. Directory layout

```
WebVis/
├── package.json          # { "type": "module", dependencies: { "h5wasm": "^0.10" } }
├── server.mjs            # Node.js HTTP server + API routes
├── fields.mjs            # Field registry (name → HDF5 path, grid type, unit)
├── hdf5-reader.mjs       # h5wasm wrapper: open file, extract field, derive invariants
├── public/               # Static files served to browser
│   ├── index.html        # Single page, loads ES6 modules
│   ├── style.css         # Minimal CSS (grid layout, dark/light theme)
│   ├── app.mjs           # Entry point: creates controller, model, views
│   ├── model.mjs         # App state + EventTarget
│   ├── controller.mjs    # Orchestration: UI events → model → API
│   ├── views/
│   │   ├── FieldCanvas.mjs
│   │   ├── ColourBar.mjs
│   │   ├── ControlPanel.mjs
│   │   ├── FileList.mjs
│   │   └── StatusBar.mjs
│   └── colormaps/
│       └── maps.json     # Colour map LUTs
└── README.md
```

### 11. API endpoints

| Method | Path | Response |
|--------|------|----------|
| GET | `/api/files` | `{ files: ["Output00000.gzip.h5", ...] }` |
| GET | `/api/fields/:file` | `{ fields: ["Pressure", "Temperature", ...], params: { time, nx, nz, dx, dz } }` |
| GET | `/api/field-data/:file/:field` | `{ xCoords: [...], zCoords: [...], values: [[...]], nx, nz, unit, log, min, max }` |
| GET | `/*` | Static files from `public/` |

The `field-data` endpoint does all heavy lifting: opens HDF5, reads datasets, reshapes, derives if needed, masks air, computes min/max, returns JSON. File handles are cached (LRU, max 5 open files) to avoid re-opening for rapid field switching.

### 12. HTTP compression

**Choice**: Enable gzip for JSON responses using `node:zlib` `createGzip()` piped to the response stream when `Accept-Encoding: gzip` is present.

**Rationale**: A 600×400 float field is ~1.5 MB JSON → ~300 KB gzipped. Significant for remote access.

## Risks / Trade-offs

**[h5wasm compatibility]** → h5wasm may not support every HDF5 feature. MDOODZ files use standard float32/float64 datasets with gzip compression — well within h5wasm's support. Mitigation: test against all 29 sample files in `benchmark-results/rifting-600x400-v2-h5/` early.

**[Large grids]** → For grids much larger than 600×400 (e.g., 2000×1000), JSON transfer and canvas rendering may slow down. Mitigation: future work could add server-side downsampling or binary transfer. Current scope targets the common 150×100 to 600×400 range.

**[File locking]** → h5wasm opens files read-only. If MDOODZ is actively writing while the viewer reads, reads may fail or return partial data. Mitigation: the viewer targets post-run analysis, not live monitoring. Document this limitation.

**[No bundler]** → Serving raw ES6 modules means many small HTTP requests on first load (~10 files). Mitigation: total JS is <50 KB; on localhost this is negligible. HTTP/2 would eliminate the concern but is not needed for local use.

**[Typed array precision]** → HDF5 stores floats as float32. JavaScript `Float32Array` preserves this. JSON serialisation converts to double-precision text — slight size inflation. Acceptable for the target grid sizes.
