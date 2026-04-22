#!/usr/bin/env python3
"""
MDOODZ performance sweep: run scenarios at multiple grid sizes and thread counts,
collect perf.csv, and write a Markdown report.

Usage:
  python misc/perf_sweep.py                         # default sweep
  python misc/perf_sweep.py --scenarios BlankenBench RiftingChenin
  python misc/perf_sweep.py --steps 5 --threads 1 4 8
  python misc/perf_sweep.py --out benchmark-results/my-run
  python misc/perf_sweep.py --report-only <dir>     # just generate report from existing CSVs
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# ---- Scenario definitions ----
# src_txt: path to the master config (relative to ROOT)
# grids: list of (nx, nz) to sweep
SCENARIOS = {
    "BlankenBench": {
        "src_txt": "cmake-exec/BlankenBench/BlankenBench.txt",  # no SETS entry
        "exec": "cmake-exec/BlankenBench/BlankenBench",
        "grids": [(41, 41), (81, 81), (161, 161)],
        "extra_patches": [
            (r"^writer\s*=.*", "writer = 0"),
            (r"^writer_step\s*=.*", "writer_step = 1"),
        ],
    },
    "RiftingChenin": {
        "src_txt": "cmake-exec/RiftingChenin/RiftingChenin.txt",
        "exec": "cmake-exec/RiftingChenin/RiftingChenin",
        "grids": [(150, 100), (300, 200)],
        "extra_patches": [
            (r"^writer\s*=.*", "writer = 0"),
        ],
    },
}

PERF_COLS = [
    "step", "wall_s", "time_ma", "rheology_s", "assembly_s", "solve_s",
    "thermal_s", "advection_s", "free_surface_s", "reseeding_s", "melting_s",
    "anisotropy_s", "gse_s", "output_s", "interp_s", "stokes_setup_s",
    "nl_overhead_s", "post_solve_s", "nit", "n_particles", "neq_mom",
    "neq_cont", "peak_rss_mb", "user_cpu_s", "sys_cpu_s",
]

PHASE_COLS = [
    "rheology_s", "assembly_s", "solve_s", "thermal_s", "advection_s",
    "interp_s", "post_solve_s", "nl_overhead_s", "stokes_setup_s",
    "free_surface_s", "reseeding_s", "melting_s", "anisotropy_s", "gse_s",
    "output_s",
]


def patch_txt(src: Path, dst: Path, nx: int, nz: int, steps: int, extra: list):
    text = src.read_text()

    def sub(pattern, replacement, t):
        new_t, n = re.subn(pattern, replacement, t, flags=re.MULTILINE)
        return new_t

    text = sub(r"^Nx\s*=.*", f"Nx      = {nx}", text)
    text = sub(r"^Nz\s*=.*", f"Nz      = {nz}", text)
    text = sub(r"^Nt\s*=.*", f"Nt      = {steps}", text)
    text = sub(r"^t_end\s*=.*", "t_end   = 0", text)
    text = sub(r"^log_timestamp\s*=.*", "log_timestamp = 0", text)
    for pattern, replacement in extra:
        text = sub(pattern, replacement, text)

    dst.write_text(text)


def run_one(exec_path: Path, txt_dir: Path, nx: int, nz: int, threads: int,
            steps: int, src_txt: Path, extra_patches: list):
    """Run one config. Returns summary dict or None on failure."""
    run_txt = txt_dir / exec_path.name / f"{exec_path.name}.txt"
    perf_csv = txt_dir / exec_path.name / "perf.csv"

    patch_txt(src_txt, run_txt, nx, nz, steps, extra_patches)
    if perf_csv.exists():
        perf_csv.unlink()

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)

    t0 = time.perf_counter()
    result = subprocess.run(
        [str(exec_path)],
        cwd=str(exec_path.parent),
        capture_output=True, text=True, env=env,
        timeout=600,
    )
    elapsed = time.perf_counter() - t0

    if result.returncode != 0:
        print(f"  FAILED {nx}x{nz} t{threads}: exit {result.returncode}")
        print(result.stderr[-500:] if result.stderr else "")
        return None

    if not perf_csv.exists():
        print(f"  FAILED {nx}x{nz} t{threads}: no perf.csv written")
        return None

    rows = list(csv.DictReader(perf_csv.read_text().splitlines()))
    if not rows:
        print(f"  FAILED {nx}x{nz} t{threads}: empty perf.csv")
        return None

    def avg(col):
        vals = [float(r[col]) for r in rows if col in r and r[col] not in ("", "nan")]
        return sum(vals) / len(vals) if vals else 0.0

    summary = {
        "scenario": exec_path.name,
        "nx": nx, "nz": nz,
        "grid": f"{nx}x{nz}",
        "threads": threads,
        "steps": len(rows),
        "elapsed_s": round(elapsed, 2),
        "avg_wall_s": round(avg("wall_s"), 4),
        "avg_nit": round(avg("nit"), 2),
        "peak_rss_mb": round(avg("peak_rss_mb"), 1),
    }
    for col in PHASE_COLS:
        summary[f"avg_{col}"] = round(avg(col), 4)

    return summary


def run_sweep(scenarios: list, threads: list, steps: int, out_dir: Path) -> list:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_results = []

    for name in scenarios:
        cfg = SCENARIOS[name]
        exec_path = ROOT / cfg["exec"]
        src_txt = ROOT / cfg["src_txt"]
        txt_dir = ROOT / "cmake-exec"

        if not exec_path.exists():
            print(f"[SKIP] {name}: executable not found at {exec_path}")
            continue
        if not src_txt.exists():
            print(f"[SKIP] {name}: source config not found at {src_txt}")
            continue

        print(f"\n{'='*50}")
        print(f"  Scenario: {name}")
        print(f"{'='*50}")

        for nx, nz in cfg["grids"]:
            for t in threads:
                label = f"{name} {nx}x{nz} {t}t"
                print(f"  Running {label} ({steps} steps)...", end="", flush=True)
                r = run_one(exec_path, txt_dir, nx, nz, t, steps, src_txt,
                            cfg["extra_patches"])
                if r:
                    all_results.append(r)
                    print(f" done  wall={r['avg_wall_s']:.3f}s/step  "
                          f"solve={r['avg_solve_s']:.4f}s  nit={r['avg_nit']:.1f}")

    return all_results


def write_summary_csv(results: list, path: Path):
    if not results:
        return
    cols = list(results[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(results)


def phase_bar(r: dict, width=30) -> str:
    """ASCII bar showing phase fractions of avg_wall_s."""
    wall = r["avg_wall_s"]
    if wall <= 0:
        return ""
    phases = [
        ("interp", r.get("avg_interp_s", 0)),
        ("solve",  r.get("avg_solve_s", 0)),
        ("assem",  r.get("avg_assembly_s", 0)),
        ("rheo",   r.get("avg_rheology_s", 0)),
        ("therm",  r.get("avg_thermal_s", 0)),
        ("adv",    r.get("avg_advection_s", 0)),
        ("post",   r.get("avg_post_solve_s", 0)),
        ("other",  0),
    ]
    accounted = sum(v for _, v in phases[:-1])
    phases[-1] = ("other", max(0.0, wall - accounted))

    chars = "■▪·"
    SYMS = ["I", "S", "A", "R", "T", "V", "P", "O"]
    bar = ""
    for i, (label, val) in enumerate(phases):
        n = max(0, round(val / wall * width))
        bar += SYMS[i] * n
    return bar[:width].ljust(width)


def write_report(results: list, out_path: Path, ts: str):
    if not results:
        out_path.write_text("No results to report.\n")
        return

    by_scenario = defaultdict(list)
    for r in results:
        by_scenario[r["scenario"]].append(r)

    lines = []
    lines += [
        f"# MDOODZ Performance Report",
        f"",
        f"_Generated: {ts}_",
        f"",
        "## Legend",
        "",
        "Phase bar key: `I`=interp, `S`=solve, `A`=assembly, `R`=rheology, "
        "`T`=thermal, `V`=advection, `P`=post-solve, `O`=other",
        "",
    ]

    for scen, rows in by_scenario.items():
        lines += [f"## {scen}", ""]

        # --- phase breakdown table (group by grid, sort by threads) ---
        lines += [
            "### Phase breakdown (avg per time step)",
            "",
            "| Grid | Thr | wall_s | solve_s | assembly_s | interp_s | "
            "post_s | rheol_s | thermal_s | nit | Phase bar (I·S·A·R·T·V·P·O) |",
            "|------|-----|-------:|--------:|-----------:|---------:|"
            "-------:|--------:|----------:|----:|:----------------------------|",
        ]
        for r in sorted(rows, key=lambda x: (x["nx"] * x["nz"], x["threads"])):
            lines.append(
                f"| {r['grid']} | {r['threads']} "
                f"| {r['avg_wall_s']:.3f} "
                f"| {r['avg_solve_s']:.4f} "
                f"| {r['avg_assembly_s']:.4f} "
                f"| {r['avg_interp_s']:.4f} "
                f"| {r['avg_post_solve_s']:.4f} "
                f"| {r['avg_rheology_s']:.4f} "
                f"| {r['avg_thermal_s']:.4f} "
                f"| {r['avg_nit']:.1f} "
                f"| `{phase_bar(r)}` |"
            )
        lines += [""]

        # --- thread scaling (for each grid, show speedup vs 1t) ---
        grids_in_scen = sorted(set((r["nx"], r["nz"]) for r in rows),
                               key=lambda g: g[0] * g[1])
        lines += ["### Thread scaling", ""]
        lines += [
            "| Grid | Thr | wall_s | Speedup vs 1t | Efficiency |",
            "|------|-----|-------:|--------------:|-----------:|",
        ]
        for (nx, nz) in grids_in_scen:
            grid_rows = sorted([r for r in rows if r["nx"] == nx and r["nz"] == nz],
                               key=lambda x: x["threads"])
            t1_wall = next((r["avg_wall_s"] for r in grid_rows if r["threads"] == 1), None)
            for r in grid_rows:
                speedup = f"{t1_wall / r['avg_wall_s']:.2f}×" if t1_wall else "—"
                eff = f"{t1_wall / r['avg_wall_s'] / r['threads'] * 100:.0f}%" if t1_wall else "—"
                lines.append(
                    f"| {r['grid']} | {r['threads']} "
                    f"| {r['avg_wall_s']:.3f} | {speedup} | {eff} |"
                )
        lines += [""]

        # --- solve fraction vs grid size ---
        lines += ["### Stokes solve fraction vs grid size", ""]
        lines += [
            "| Grid | Threads | wall_s | solve_s | solve% | assembly% |",
            "|------|---------|-------:|--------:|-------:|----------:|",
        ]
        for r in sorted(rows, key=lambda x: (x["threads"], x["nx"] * x["nz"])):
            w = r["avg_wall_s"]
            solve_pct = r["avg_solve_s"] / w * 100 if w > 0 else 0
            assem_pct = r["avg_assembly_s"] / w * 100 if w > 0 else 0
            lines.append(
                f"| {r['grid']} | {r['threads']} "
                f"| {r['avg_wall_s']:.3f} "
                f"| {r['avg_solve_s']:.4f} "
                f"| {solve_pct:.1f}% "
                f"| {assem_pct:.1f}% |"
            )
        lines += [""]

    # --- memory ---
    lines += ["## Memory", ""]
    lines += [
        "| Scenario | Grid | Thr | peak_rss_mb |",
        "|----------|------|-----|------------:|",
    ]
    for r in results:
        lines.append(f"| {r['scenario']} | {r['grid']} | {r['threads']} | {r['peak_rss_mb']:.0f} |")
    lines += [""]

    # --- observations ---
    lines += [
        "## Observations",
        "",
        "*(Fill in after reviewing the tables above.)*",
        "",
        "- **Dominant phase:** …",
        "- **Thread scaling:** …",
        "- **Stokes solve fraction:** …",
        "- **Next bottleneck after current optimizations:** …",
        "",
    ]

    out_path.write_text("\n".join(lines) + "\n")
    print(f"\nReport written to: {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scenarios", nargs="+", default=list(SCENARIOS.keys()),
                    help="Scenarios to run (default: all)")
    ap.add_argument("--threads", nargs="+", type=int, default=[1, 2, 4, 8],
                    help="Thread counts")
    ap.add_argument("--steps", type=int, default=5,
                    help="Time steps per run")
    ap.add_argument("--out", default=None,
                    help="Output directory (default: benchmark-results/<timestamp>)")
    ap.add_argument("--report-only", metavar="DIR",
                    help="Skip runs; generate report from existing summary.csv in DIR")
    args = ap.parse_args()

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    out_dir = Path(args.out) if args.out else (
        ROOT / "benchmark-results" / datetime.now().strftime("%Y%m%d-%H%M%S"))

    if args.report_only:
        rdir = Path(args.report_only)
        csv_path = rdir / "summary.csv"
        if not csv_path.exists():
            sys.exit(f"No summary.csv in {rdir}")
        results = list(csv.DictReader(csv_path.read_text().splitlines()))
        for r in results:
            for k, v in r.items():
                try:
                    r[k] = float(v) if "." in v else int(v)
                except (ValueError, TypeError):
                    pass
        report_path = rdir / f"REPORT-{datetime.now().strftime('%Y%m%d-%H%M%S')}.md"
        write_report(results, report_path, ts)
        return

    results = run_sweep(args.scenarios, args.threads, args.steps, out_dir)

    csv_path = out_dir / "summary.csv"
    write_summary_csv(results, csv_path)
    print(f"\nSummary CSV: {csv_path}")

    report_path = out_dir / f"REPORT-{datetime.now().strftime('%Y%m%d-%H%M%S')}.md"
    write_report(results, report_path, ts)


if __name__ == "__main__":
    main()
