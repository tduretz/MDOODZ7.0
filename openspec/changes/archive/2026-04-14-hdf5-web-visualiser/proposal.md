## Why

Visualising MDOODZ HDF5 output currently requires either a Julia/Makie environment (heavy, slow to install) or Python/h5py + gnuplot (fragile, no interactivity). Researchers need a fast, dependency-light way to browse simulation output in the browser — select a time step, pick a field, see a colour-mapped image immediately. A lightweight Node.js web app using `h5wasm` (WebAssembly HDF5 reader) would give sub-second renders with minimal setup — no native compilation, no system HDF5 library required.

## What Changes

- Add a **Node.js HTTP server** (zero-framework, plain `http` module + `h5wasm`) that:
  - Reads `Output*.gzip.h5` files directly in-process via h5wasm (WebAssembly HDF5 library) — no native compilation, no child process spawning.
  - Handles staggered-grid reshaping, coordinate mapping, derived field computation (stress/strain-rate invariants, vertex-to-centre averaging), and air-cell masking in JavaScript using typed arrays.
  - Serves a static ES6+ single-page frontend.
  - Exposes REST endpoints to list available HDF5 files, list fields in a file, and return extracted 2D field data as JSON.
- Add an **ES6+ browser frontend** (no build step, no React/Vue) with:
  - MVC architecture: separate Model, View, and Controller modules; reactivity via a listener/observer pattern.
  - A canvas-based colour-mapped field renderer (client-side, no server-side image generation).
  - Controls: file selector, time-step slider, field dropdown, colour-map picker, colour-range controls.
  - Clean, minimal UI that looks professional.
- The system is designed to be **extensible** — future features (overlays, vector fields, multi-panel, topo line, animations) plug into the existing MVC structure.

## Capabilities

### New Capabilities

- `hdf5-reader`: Server-side module using h5wasm (WebAssembly HDF5) to read MDOODZ output files in-process. Handles all staggered-grid layouts (centres, vertices, Vx-nodes, Vz-nodes, VizGrid), coordinate mapping, derived invariant computation (stress/strain-rate second invariants, vertex-to-centre averaging), and air-cell masking. Pure JavaScript with typed arrays — no native compilation.
- `web-server`: Node.js HTTP server that serves the frontend and exposes REST API endpoints (`/api/files`, `/api/fields`, `/api/field-data`). Reads HDF5 directly via the hdf5-reader module.
- `field-renderer`: Browser-side ES6+ canvas renderer that takes 2D field data + coordinates and produces colour-mapped images with proper axis scaling, colour bar, and NaN masking (air cells). Supports multiple colour maps.
- `web-ui`: ES6+ MVC frontend — Model (field data, app state), View (canvas renderer, controls, layout), Controller (user interaction → model updates → view re-renders). Reactivity via observer/listener pattern. No framework, no build step.

### Modified Capabilities

_(none — this is a new standalone tool, no existing specs are affected)_

## Impact

- **New directories**: `WebVis/` (or similar) at the project root containing the server and frontend.
- **Build system**: No CMake changes — the web app is npm-only.
- **Dependencies added**: Node.js (runtime, assumed installed), `h5wasm` (WebAssembly HDF5 reader, zero native bindings). No system HDF5 library needed for the web app.
- **No changes** to MDLIB, SETS, or any existing simulation code.
- **Data files**: uses existing `Output*.gzip.h5` files as-is; no format changes.
