# MDOODZ WebVis

Lightweight web visualiser for MDOODZ HDF5 simulation output.

## Prerequisites

- Node.js ≥ 18

## Installation

```bash
cd WebVis
npm install
```

This installs [h5wasm](https://github.com/usnistgov/h5wasm) — a WebAssembly port of the HDF5 C library (no native dependencies).

## Usage

```bash
node server.mjs /path/to/hdf5/dir
```

Then open http://localhost:3000.

The data directory should contain MDOODZ output files named `Output*.gzip.h5`.

### Options

- **Port**: set via `PORT` env variable (default: 3000)
  ```bash
  PORT=8080 node server.mjs /path/to/data
  ```

## Supported Fields

**47 fields** across all MDOODZ grid types:

| Category | Examples |
|----------|----------|
| Centre fields | Pressure, Temperature, Viscosity, Density, Grain size, Melt fraction, Porosity, Stresses, Strain rates, Accumulated strains, eII partitions, Divergence fields, Cohesion, Friction angle |
| Vertex fields | Shear stress, Shear strain rate, Vertex viscosity/density |
| Velocity fields | Vx, Vz (staggered grids) |
| Derived fields | Stress II (τ\_II), Strain rate II (ε̇\_II) |
| Phase maps | Phases, Phases HR, Phases dual HR |

## Architecture

```
WebVis/
├── server.mjs          # Node.js HTTP server (plain node:http)
├── hdf5-reader.mjs     # h5wasm-based HDF5 reading + derived fields
├── fields.mjs          # Declarative field registry
├── package.json
└── public/
    ├── index.html       # Single-page shell
    ├── app.mjs          # Entry point (wires MVC)
    ├── model.mjs        # EventTarget-based state model
    ├── controller.mjs   # API orchestration
    ├── style.css        # Dark theme, CSS Grid layout
    ├── colormaps/
    │   └── maps.json    # 256-entry LUTs (viridis, turbo, inferno, plasma, coolwarm)
    └── views/
        ├── FieldCanvas.mjs    # Canvas ImageData pixel renderer
        ├── ColourBar.mjs      # Gradient bar + tick labels
        ├── ControlPanel.mjs   # File/field/colourmap selectors
        ├── StatusBar.mjs      # Time, grid dims, field name
        └── LoadingOverlay.mjs # Loading indicator
```

**No build step, no framework, no bundler.** Plain ES6 modules served directly.

## API

| Endpoint | Returns |
|----------|---------|
| `GET /api/files` | `{ files: ["Output00000.gzip.h5", ...] }` |
| `GET /api/fields/:file` | `{ fields: [...], params: { time, nx, nz, dx, dz, dt, ... } }` |
| `GET /api/field-data/:file/:field` | `{ values, xCoords, zCoords, nx, nz, unit, log, min, max }` |

Responses are gzip-compressed when `Accept-Encoding: gzip` is present.
