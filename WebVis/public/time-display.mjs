// ── Adaptive Time Display ─────────────────────────────────────────────
// Converts model time (seconds) to human-readable yr / ka / Ma.

const SEC_PER_YR = 3.1558e7;

/**
 * Format model time in seconds to a human-readable string.
 * @param {number} seconds  Model time in seconds.
 * @param {string|null} overrideUnit  Force a specific unit: 'yr', 'ka', 'Ma', or null for auto.
 * @returns {{ value: number, unit: string, formatted: string }}
 */
export function formatTime(seconds, overrideUnit = null) {
  const yr = seconds / SEC_PER_YR;

  let unit;
  if (overrideUnit) {
    unit = overrideUnit;
  } else if (yr < 1000) {
    unit = 'yr';
  } else if (yr < 1e6) {
    unit = 'ka';
  } else {
    unit = 'Ma';
  }

  let value;
  switch (unit) {
    case 'yr': value = yr;       break;
    case 'ka': value = yr / 1e3; break;
    case 'Ma': value = yr / 1e6; break;
    default:   value = yr / 1e6; unit = 'Ma'; break;
  }

  return { value, unit, formatted: `${value.toFixed(2)} ${unit}` };
}
