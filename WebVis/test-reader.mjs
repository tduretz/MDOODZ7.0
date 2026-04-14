import { readMetadata, extractField, computeMinMax } from './hdf5-reader.mjs';

const TEST_FILE = '../benchmark-results/rifting-600x400-v2-h5/Output00000.gzip.h5';

async function test() {
  console.log('=== readMetadata ===');
  const meta = await readMetadata(TEST_FILE);
  console.log(meta);
  console.assert(meta.nx > 0, 'nx should be > 0');
  console.assert(meta.nz > 0, 'nz should be > 0');

  console.log('\n=== extractField: Pressure (center) ===');
  const P = await extractField(TEST_FILE, { path: '/Centers/P', gridType: 'center', unit: 'Pa' });
  console.log(`shape: ${P.nx} x ${P.nz}, min=${P.min}, max=${P.max}`);
  console.assert(P.nx === meta.nx - 1, `Expected rows=${meta.nx - 1}, got ${P.nx}`);
  console.assert(P.nz === meta.nz - 1, `Expected cols=${meta.nz - 1}, got ${P.nz}`);
  console.assert(P.xCoords.length === meta.nx - 1, 'xCoords length mismatch');
  console.assert(P.zCoords.length === meta.nz - 1, 'zCoords length mismatch');

  console.log('\n=== extractField: sxz (vertex) ===');
  const sxz = await extractField(TEST_FILE, { path: '/Vertices/sxz', gridType: 'vertex', unit: 'Pa' });
  console.log(`shape: ${sxz.nx} x ${sxz.nz}, min=${sxz.min}, max=${sxz.max}`);
  console.assert(sxz.nx === meta.nx, `Expected rows=${meta.nx}, got ${sxz.nx}`);
  console.assert(sxz.nz === meta.nz, `Expected cols=${meta.nz}, got ${sxz.nz}`);

  console.log('\n=== extractField: Vx (vx-nodes) ===');
  const vx = await extractField(TEST_FILE, { path: '/VxNodes/Vx', gridType: 'vx', unit: 'm/s' });
  console.log(`shape: ${vx.nx} x ${vx.nz}, min=${vx.min}, max=${vx.max}`);
  console.assert(vx.nx === meta.nx, `Expected rows=${meta.nx}, got ${vx.nx}`);
  console.assert(vx.nz === meta.nz + 1, `Expected cols=${meta.nz + 1}, got ${vx.nz}`);

  console.log('\n=== extractField: Vz (vz-nodes) ===');
  const vz = await extractField(TEST_FILE, { path: '/VzNodes/Vz', gridType: 'vz', unit: 'm/s' });
  console.log(`shape: ${vz.nx} x ${vz.nz}, min=${vz.min}, max=${vz.max}`);
  console.assert(vz.nx === meta.nx + 1, `Expected rows=${meta.nx + 1}, got ${vz.nx}`);
  console.assert(vz.nz === meta.nz, `Expected cols=${meta.nz}, got ${vz.nz}`);

  console.log('\n=== extractField: Stress II (derived) ===');
  const sII = await extractField(TEST_FILE, { derive: 'stress_invariant', gridType: 'center', unit: 'Pa' });
  console.log(`shape: ${sII.nx} x ${sII.nz}, min=${sII.min}, max=${sII.max}`);
  console.assert(sII.min >= 0, 'Stress II min should be >= 0');

  console.log('\n=== extractField: Strain rate II (derived) ===');
  const eII = await extractField(TEST_FILE, { derive: 'strainrate_invariant', gridType: 'center', unit: '1/s', log: true });
  console.log(`shape: ${eII.nx} x ${eII.nz}, min=${eII.min}, max=${eII.max}`);
  console.assert(eII.min >= 0, 'Strain rate II min should be >= 0');

  console.log('\n=== extractField: Phases (vizgrid) ===');
  const ph = await extractField(TEST_FILE, { path: '/VizGrid/compo', gridType: 'vizgrid', discrete: true });
  console.log(`shape: ${ph.nx} x ${ph.nz}, min=${ph.min}, max=${ph.max}`);

  console.log('\n=== extractField: Temperature (with offset) ===');
  const T = await extractField(TEST_FILE, { path: '/Centers/T', gridType: 'center', unit: '°C', offset: -273.15 });
  console.log(`shape: ${T.nx} x ${T.nz}, min=${T.min.toFixed(1)}°C, max=${T.max.toFixed(1)}°C`);

  console.log('\n✓ All tests passed!');
}

test().catch(e => { console.error(e); process.exit(1); });
