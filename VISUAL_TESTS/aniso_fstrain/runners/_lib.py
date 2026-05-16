"""Shared helpers for the ani_fstrain visual-tests runners.

The runners share a single C binary (`SETS/AnisoFstrainBox.c` →
`cmake-exec/AnisoFstrainBox/AnisoFstrainBox`) and a single .txt template
(`VISUAL_TESTS/aniso_fstrain/templates/base.txt`). Each runner defines a
sweep grid, builds a per-cell .txt by substituting placeholders, runs
the binary, and captures the HDF5 outputs under
`VISUAL_TESTS/aniso_fstrain/runs/<test>/<cell>/output/`.

This module exposes:
- `ROOT`, `BIN`, `TEMPLATE`, `RUNS_DIR` — common paths
- `render_template(subs, template=TEMPLATE)` — placeholder substitution
- `default_subs(**overrides)` — convenient dict with sensible defaults
- `run_cell(cell_dir, txt, env_threads=8, timeout=600)` — invoke binary
- `run_sweep(cells, max_workers=4, threads_per_cell=2, ...)` — parallel
- `ensure_build()` — verifies the binary exists, with an actionable msg
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

# repo-root resolution: VISUAL_TESTS/aniso_fstrain/runners/_lib.py → ../../..
ROOT = Path(__file__).resolve().parents[3]
EXEC_DIR = ROOT / "cmake-exec/AnisoFstrainBox"
BIN = EXEC_DIR / "AnisoFstrainBox"
TEMPLATE = ROOT / "VISUAL_TESTS/aniso_fstrain/templates/base.txt"
RUNS_DIR = ROOT / "VISUAL_TESTS/aniso_fstrain/runs"


# --------------------------------------------------------------------------- #
# Default substitutions
# --------------------------------------------------------------------------- #

# All placeholders the template understands. Defaults reproduce the
# minimal AnisoFstrainBox.txt (olivine, simple shear, 21x21, Nt=5).
DEFAULTS: dict[str, object] = {
    "WRITER_STEP":       1,
    "SCALE_ETA":         "1e0",
    "SCALE_L":           "1e0",
    "SCALE_V":           "1e0",
    "SCALE_T":           "1",
    "NX":                21,
    "NZ":                21,
    "NT":                5,
    "DT":                "0.5",
    "MECHANICAL":        1,
    "BKG_STRAIN_RATE":   "1.0",
    "BKG_TEMPERATURE":   "300.0",
    "PWLV":              0,
    "ETA0":              "1e0",
    "ANISO_ANGLE":       80,
    "ANISO_FACTOR":      "4.0",
    "ANI_FSTRAIN":       3,
    "ANI_RELAX_EPS_MAX": "-1.0",
}


def render_template(subs: dict, template: Path = TEMPLATE) -> str:
    """Substitute every `__KEY__` placeholder in the template with subs[KEY].

    Missing keys fall back to DEFAULTS; unknown placeholders raise."""
    text = template.read_text()
    merged: dict[str, object] = {**DEFAULTS, **subs}
    for key, value in merged.items():
        placeholder = f"__{key}__"
        text = text.replace(placeholder, str(value))
    # Sanity check — should have no `__FOO__` left over
    import re
    leftovers = re.findall(r"__[A-Z_]+__", text)
    if leftovers:
        raise ValueError(f"Template still has unsubstituted placeholders: {sorted(set(leftovers))}")
    return text


def default_subs(**overrides) -> dict[str, object]:
    """Return a copy of DEFAULTS with the given overrides applied."""
    out = dict(DEFAULTS)
    out.update(overrides)
    return out


# --------------------------------------------------------------------------- #
# Binary invocation
# --------------------------------------------------------------------------- #

def ensure_build() -> None:
    if not BIN.exists():
        msg = (
            f"AnisoFstrainBox binary not found at {BIN}\n"
            f"Build it from {ROOT}/cmake-build with:\n"
            f"  cd {ROOT}/cmake-build && cmake -USET .. && make -j 8 AnisoFstrainBox"
        )
        raise FileNotFoundError(msg)


def run_cell(cell_dir: Path, subs: dict, *,
             threads: int = 2, timeout: int = 600,
             extra_env: dict | None = None) -> dict:
    """Render the per-cell .txt, copy/symlink the binary, run it.

    Returns a dict with rc, wall, log_path, txt_path, output_dir.
    """
    cell_dir.mkdir(parents=True, exist_ok=True)
    txt_path = cell_dir / "input.txt"
    txt_path.write_text(render_template(subs))

    # Symlink the binary so the .gzip outputs go into cell_dir/output/
    cell_bin = cell_dir / BIN.name
    if not cell_bin.exists():
        try:
            os.symlink(BIN, cell_bin)
        except FileExistsError:
            pass

    # Clean previous output
    out = cell_dir / "output"
    if out.exists():
        shutil.rmtree(out)

    log = cell_dir / "stdout.log"
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    if extra_env:
        env.update(extra_env)
    t0 = time.time()
    try:
        rc = subprocess.run(
            [str(cell_bin), str(txt_path)],
            cwd=str(cell_dir), env=env,
            stdout=open(log, "w"), stderr=subprocess.STDOUT,
            timeout=timeout,
        ).returncode
    except subprocess.TimeoutExpired:
        rc = -1
    wall = time.time() - t0
    return {
        "cell_dir":  cell_dir,
        "rc":        rc,
        "wall":      wall,
        "log_path":  log,
        "txt_path":  txt_path,
        "output_dir": out,
    }


@dataclass
class Cell:
    """One sweep cell."""
    name: str
    subs: dict
    meta: dict = field(default_factory=dict)


def run_sweep(test_name: str,
              cells: Iterable[Cell],
              *,
              max_workers: int = 4,
              threads_per_cell: int = 2,
              timeout: int = 600,
              on_result: Callable[[Cell, dict], None] | None = None,
              ) -> list[dict]:
    """Run a list of `Cell`s in parallel under `RUNS_DIR/<test_name>/`."""
    ensure_build()
    test_dir = RUNS_DIR / test_name
    test_dir.mkdir(parents=True, exist_ok=True)
    cells = list(cells)
    print(f"[{test_name}] running {len(cells)} cells "
          f"with {max_workers} workers x OMP={threads_per_cell} threads ...", flush=True)
    results: list[dict] = []
    t_global = time.time()
    with ThreadPoolExecutor(max_workers=max_workers) as ex:
        futs = {
            ex.submit(run_cell, test_dir / c.name, c.subs,
                      threads=threads_per_cell, timeout=timeout): c
            for c in cells
        }
        done = 0
        for fut in as_completed(futs):
            c = futs[fut]
            try:
                r = fut.result()
            except Exception as e:
                r = {"cell_dir": test_dir / c.name, "rc": -2,
                     "wall": 0.0, "exception": repr(e)}
            r["cell"] = c
            results.append(r)
            done += 1
            ok = r.get("rc", -3) == 0
            mark = "OK" if ok else f"FAIL(rc={r.get('rc')})"
            print(f"  [{done}/{len(cells)}] {c.name}: {mark} wall={r.get('wall',0):.1f}s",
                  flush=True)
            if on_result:
                try:
                    on_result(c, r)
                except Exception as e:
                    print(f"    on_result error: {e!r}", flush=True)
    total = time.time() - t_global
    fails = [r for r in results if r.get("rc", -3) != 0]
    print(f"[{test_name}] DONE in {total:.1f}s; {len(results)-len(fails)} ok, {len(fails)} failed",
          flush=True)
    return results
