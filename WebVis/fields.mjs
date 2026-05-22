// ── Field Registry ─────────────────────────────────────────────────────
// Declarative mapping of display names → HDF5 paths, grid types, units,
// rendering hints, and derive functions.
// Each entry has:
//   path      – HDF5 dataset path (null for derived fields)
//   gridType  – 'center' | 'vertex' | 'vx' | 'vz' | 'vizgrid' | 'vizgrid_hr'
//   unit      – display unit string
//   log       – if true, frontend should use log₁₀ colour scale
//   offset    – value added to every element (e.g., -273.15 for K→°C)
//   discrete  – if true, field contains integer categories (skip air mask, use discrete palette)
//   derive    – name of a derive function in hdf5-reader (null for plain datasets)

const REGISTRY = {
  // ── Centre fields (Nx-1 × Nz-1) ──────────────────────────────────────
  'Pressure':           { path: '/Centers/P',           gridType: 'center', unit: 'Pa',    log: false, label: 'P',           formattedUnit: 'Pa'      },
  'Temperature':        { path: '/Centers/T',           gridType: 'center', unit: '°C',    log: false, offset: -273.15, label: 'T',  formattedUnit: '°C'  },
  'Viscosity':          { path: '/Centers/eta_n',       gridType: 'center', unit: 'Pa·s',  log: true,  label: 'η',           formattedUnit: 'Pa·s'    },
  'Density':            { path: '/Centers/rho_n',       gridType: 'center', unit: 'kg/m³', log: false, label: 'ρ',           formattedUnit: 'kg·m⁻³'  },
  'Grain size':         { path: '/Centers/d',           gridType: 'center', unit: 'm',     log: true,  label: 'd',           formattedUnit: 'm'       },
  'Melt fraction':      { path: '/Centers/X',           gridType: 'center', unit: '',      log: false, label: 'ϕₘₑₗₜ',      formattedUnit: ''        },
  'Porosity':           { path: '/Centers/phi',         gridType: 'center', unit: '',      log: false, label: 'φ',           formattedUnit: ''        },
  'Overpressure':       { path: '/Centers/OverS',       gridType: 'center', unit: 'Pa',    log: false, label: 'ΔP',          formattedUnit: 'Pa'      },
  'Stress xx':          { path: '/Centers/sxxd',        gridType: 'center', unit: 'Pa',    log: false, label: 'σ′ₓₓ',       formattedUnit: 'Pa'      },
  'Stress zz':          { path: '/Centers/szzd',        gridType: 'center', unit: 'Pa',    log: false, label: 'σ′zz',        formattedUnit: 'Pa'      },
  'Strain rate xx':     { path: '/Centers/exxd',        gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ₓₓ',        formattedUnit: 's⁻¹'     },
  'Strain rate zz':     { path: '/Centers/ezzd',        gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇zz',        formattedUnit: 's⁻¹'     },
  'Divergence':         { path: '/Centers/divu',        gridType: 'center', unit: '1/s',   log: false, label: '∇·v',         formattedUnit: 's⁻¹'     },
  'Div elastic':        { path: '/Centers/divu_el',     gridType: 'center', unit: '1/s',   log: false, label: '∇·vₑₗ',      formattedUnit: 's⁻¹'     },
  'Div plastic':        { path: '/Centers/divu_pl',     gridType: 'center', unit: '1/s',   log: false, label: '∇·vₚₗ',      formattedUnit: 's⁻¹'     },
  'Div reaction':       { path: '/Centers/divu_r',      gridType: 'center', unit: '1/s',   log: false, label: '∇·vᵣ',       formattedUnit: 's⁻¹'     },
  'Div thermal':        { path: '/Centers/divu_th',     gridType: 'center', unit: '1/s',   log: false, label: '∇·vₜₕ',      formattedUnit: 's⁻¹'     },
  'Strain (total)':     { path: '/Centers/strain',      gridType: 'center', unit: '',      log: false, label: 'εₜₒₜ',       formattedUnit: ''        },
  'Strain (elastic)':   { path: '/Centers/strain_el',   gridType: 'center', unit: '',      log: false, label: 'εₑₗ',        formattedUnit: ''        },
  'Strain (plastic)':   { path: '/Centers/strain_pl',   gridType: 'center', unit: '',      log: false, label: 'εₚₗ',        formattedUnit: ''        },
  'Strain (pwl)':       { path: '/Centers/strain_pwl',  gridType: 'center', unit: '',      log: false, label: 'εₚwₗ',       formattedUnit: ''        },
  'Strain (exp)':       { path: '/Centers/strain_exp',  gridType: 'center', unit: '',      log: false, label: 'εₑₓₚ',       formattedUnit: ''        },
  'Strain (lin)':       { path: '/Centers/strain_lin',  gridType: 'center', unit: '',      log: false, label: 'εₗᵢₙ',       formattedUnit: ''        },
  'Strain (gbs)':       { path: '/Centers/strain_gbs',  gridType: 'center', unit: '',      log: false, label: 'εgbs',        formattedUnit: ''        },
  'Strain (pl vol)':    { path: '/Centers/strain_pl_vol', gridType: 'center', unit: '',    log: false, label: 'εₚₗ ᵥₒₗ',    formattedUnit: ''        },
  'eII elastic':        { path: '/Centers/eII_el',      gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(ₑₗ)',    formattedUnit: 's⁻¹'     },
  'eII plastic':        { path: '/Centers/eII_pl',      gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(ₚₗ)',    formattedUnit: 's⁻¹'     },
  'eII pwl':            { path: '/Centers/eII_pwl',     gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(ₚwₗ)',   formattedUnit: 's⁻¹'     },
  'eII exp':            { path: '/Centers/eII_exp',     gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(ₑₓₚ)',   formattedUnit: 's⁻¹'     },
  'eII lin':            { path: '/Centers/eII_lin',     gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(ₗᵢₙ)',   formattedUnit: 's⁻¹'     },
  'eII gbs':            { path: '/Centers/eII_gbs',     gridType: 'center', unit: '1/s',   log: true,  label: 'ε̇ᵢᵢ(gbs)',   formattedUnit: 's⁻¹'     },
  'Cohesion':           { path: '/Centers/cohesion',    gridType: 'center', unit: 'Pa',    log: false, label: 'C',           formattedUnit: 'Pa'      },
  'Friction angle':     { path: '/Centers/friction',    gridType: 'center', unit: '°',     log: false, label: 'φf',          formattedUnit: '°'       },
  'Anisotropy factor':  { path: '/Centers/ani_fac',     gridType: 'center', unit: '',      log: false, label: 'δₐₙᵢ',       formattedUnit: ''        },
  'Director nx':        { path: '/Centers/nx',          gridType: 'center', unit: '',      log: false, label: 'nₓ',          formattedUnit: ''        },
  'Director nz':        { path: '/Centers/nz',          gridType: 'center', unit: '',      log: false, label: 'nz',          formattedUnit: ''        },

  // ── Vertex fields (Nx × Nz) ──────────────────────────────────────────
  'Shear stress':       { path: '/Vertices/sxz',        gridType: 'vertex', unit: 'Pa',    log: false, label: 'σₓz',         formattedUnit: 'Pa'      },
  'Shear strain rate':  { path: '/Vertices/exz',        gridType: 'vertex', unit: '1/s',   log: true,  label: 'ε̇ₓz',        formattedUnit: 's⁻¹'     },
  'Viscosity (vertex)': { path: '/Vertices/eta_s',      gridType: 'vertex', unit: 'Pa·s',  log: true,  label: 'ηₛ',          formattedUnit: 'Pa·s'    },
  'Density (vertex)':   { path: '/Vertices/rho_s',      gridType: 'vertex', unit: 'kg/m³', log: false, label: 'ρₛ',          formattedUnit: 'kg·m⁻³'  },

  // ── Velocity fields (staggered) ──────────────────────────────────────
  'Vx':                 { path: '/VxNodes/Vx',           gridType: 'vx',     unit: 'm/s',   log: false, label: 'Vₓ',          formattedUnit: 'm·s⁻¹'   },
  'Vz':                 { path: '/VzNodes/Vz',           gridType: 'vz',     unit: 'm/s',   log: false, label: 'Vz',          formattedUnit: 'm·s⁻¹'   },

  // ── Derived fields (computed from multiple datasets) ─────────────────
  'Stress II':          { path: null, gridType: 'center', unit: 'Pa',  log: false, derive: 'stress_invariant',     label: 'τᵢᵢ',   formattedUnit: 'Pa'    },
  'Strain rate II':     { path: null, gridType: 'center', unit: '1/s', log: true,  derive: 'strainrate_invariant', label: 'ε̇ᵢᵢ',   formattedUnit: 's⁻¹'   },

  // ── VizGrid fields ───────────────────────────────────────────────────
  'Phases':             { path: '/VizGrid/compo',        gridType: 'vizgrid',    unit: '', log: false, discrete: true, label: 'Phase',       formattedUnit: '' },
  'Phases HR':          { path: '/VizGrid/compo_hr',     gridType: 'vizgrid_hr', unit: '', log: false, discrete: true, label: 'Phase (HR)',  formattedUnit: '' },
  'Phases dual HR':     { path: '/VizGrid/compo_dual_hr', gridType: 'vizgrid_hr', unit: '', log: false, discrete: true, label: 'Phase (2×HR)', formattedUnit: '' },
};

/** Return the field definition for a given display name, or null. */
export function getField(name) {
  return REGISTRY[name] || null;
}

/** Return all registered field names. */
export function allFieldNames() {
  return Object.keys(REGISTRY);
}

/** Return the scientific display label for a field (Unicode). */
export function getFieldLabel(name) {
  const def = REGISTRY[name];
  return def ? (def.label || name) : name;
}

/** Return the formatted unit string for a field (Unicode). */
export function getFieldUnit(name) {
  const def = REGISTRY[name];
  return def ? (def.formattedUnit || '') : '';
}

/**
 * Return names of fields that actually exist in the given file.
 * @param {string}   filePath
 * @param {Function} datasetExists – async (filePath, hdf5Path) → bool
 */
export async function listAvailableFields(filePath, datasetExists) {
  const available = [];
  const checks = Object.entries(REGISTRY).map(async ([name, def]) => {
    if (def.derive) {
      // derived fields: check that all source datasets exist
      if (def.derive === 'stress_invariant') {
        const ok = await datasetExists(filePath, '/Centers/sxxd')
               && await datasetExists(filePath, '/Centers/szzd')
               && await datasetExists(filePath, '/Vertices/sxz');
        if (ok) available.push(name);
      } else if (def.derive === 'strainrate_invariant') {
        const ok = await datasetExists(filePath, '/Centers/exxd')
               && await datasetExists(filePath, '/Centers/ezzd')
               && await datasetExists(filePath, '/Vertices/exz');
        if (ok) available.push(name);
      }
    } else {
      if (await datasetExists(filePath, def.path)) available.push(name);
    }
  });
  await Promise.all(checks);
  // Sort to keep a stable order matching REGISTRY definition order
  const order = Object.keys(REGISTRY);
  available.sort((a, b) => order.indexOf(a) - order.indexOf(b));
  return available;
}
