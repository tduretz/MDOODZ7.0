#!/usr/bin/env python3
"""Aggregate the per-run MDOODZ perf.csv files into a single master CSV.

MDOODZ writes `perf.csv` in its current working directory with one row per
time step and columns covering every subsystem (rheology, assembly, solve,
thermal, advection, free_surface, …). `run_local.sh` now preserves those
files under `benchmarks/ec2/results/perf_detail/<key>/perf.csv`.

This script walks that tree, adds the configuration key as leading columns,
and emits:
  benchmarks/ec2/results/perf_detail.csv   (master merged table)
  benchmarks/ec2/results/perf_detail.md    (per-config subsystem breakdown)

The markdown summary is what the defence consumes when it needs to answer
"where is the time going?" — the wall-time column in results.csv alone
cannot distinguish solve_s from rheology_s from interp_s.

Usage:
  python3 benchmarks/ec2/_aggregate_perf_detail.py
  python3 benchmarks/ec2/_aggregate_perf_detail.py --root <path>
"""
from __future__ import annotations

import argparse
import csv
import pathlib
import re
import sys
from collections import defaultdict

KEY_RE = re.compile(
    r"^(?P<nx>\d+)x(?P<nz>\d+)_lin(?P<lin>-?\d+)_thr(?P<thr>\d+)_rep(?P<rep>.+)$"
)


def parse_key(name: str) -> dict | None:
    m = KEY_RE.match(name)
    if not m:
        return None
    d = m.groupdict()
    return {
        "nx": int(d["nx"]),
        "nz": int(d["nz"]),
        "lin_solver": int(d["lin"]),
        "n_threads": int(d["thr"]),
        "rep": d["rep"],
    }


def gather_perf_rows(perf_detail_dir: pathlib.Path):
    """Yield (key_dict, perf.csv header list, perf.csv data row list)."""
    for child in sorted(perf_detail_dir.iterdir()):
        if not child.is_dir():
            continue
        key = parse_key(child.name)
        if key is None:
            continue
        csv_path = child / "perf.csv"
        if not csv_path.exists():
            continue
        with csv_path.open() as f:
            r = csv.reader(f)
            rows = list(r)
        if len(rows) < 2:
            continue
        header = rows[0]
        for data_row in rows[1:]:
            yield key, header, data_row


def parse_fgmres_stats(mdoodz_log: pathlib.Path) -> dict:
    """Extract FGMRES convergence info from MDOODZ log, if present.

    Returns {final_res, total_iter, restarts, converged}. The GMG-FGMRES
    lines written by StokesAssemblyGMG.c look like:
        GMG-FGMRES: restart N, total_iter K, residual X.XXXe±NN
    with an optional terminating
        GMG-FGMRES converged: iters=K, restarts=N, final_res=X.XXXe±NN
    """
    stats = {"final_res": "", "total_iter": -1, "restarts": -1, "converged": False}
    if not mdoodz_log.exists():
        return stats
    try:
        text = mdoodz_log.read_text()
    except OSError:
        return stats

    # Detect explicit converged line first.
    for line in text.splitlines():
        m = re.search(
            r"GMG-FGMRES converged: iters=(\d+), restarts=(\d+), final_res=([\d\.eE+-]+)",
            line,
        )
        if m:
            stats["total_iter"] = int(m.group(1))
            stats["restarts"] = int(m.group(2))
            stats["final_res"] = m.group(3)
            stats["converged"] = True
            return stats

    # Otherwise, take the last restart line.
    last = None
    for line in text.splitlines():
        m = re.search(
            r"GMG-FGMRES: restart (\d+), total_iter (\d+), residual ([\d\.eE+-]+)",
            line,
        )
        if m:
            last = m
    if last is not None:
        stats["restarts"] = int(last.group(1)) + 1  # 0-indexed → count
        stats["total_iter"] = int(last.group(2))
        stats["final_res"] = last.group(3)
        stats["converged"] = False
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--root",
        default=str(pathlib.Path(__file__).parent / "results"),
        help="results/ directory containing perf_detail/",
    )
    args = ap.parse_args()

    root = pathlib.Path(args.root)
    perf_detail = root / "perf_detail"
    if not perf_detail.exists():
        print(f"{perf_detail} does not exist; nothing to aggregate", file=sys.stderr)
        return 1

    master_path = root / "perf_detail.csv"
    md_path = root / "perf_detail.md"

    # Merged CSV
    merged_rows = []
    master_header: list[str] | None = None
    config_rows = defaultdict(list)  # (nx,nz,lin,thr) -> list of row dicts

    for key, header, row in gather_perf_rows(perf_detail):
        fgmres = parse_fgmres_stats(
            perf_detail / _key_dir_name(key) / "mdoodz.log"
        )
        extra_header = [
            "nx",
            "nz",
            "lin_solver",
            "n_threads",
            "rep",
            "fgmres_total_iter",
            "fgmres_restarts",
            "fgmres_final_res",
            "fgmres_converged",
        ]
        if master_header is None:
            master_header = extra_header + header
        extra_row = [
            key["nx"],
            key["nz"],
            key["lin_solver"],
            key["n_threads"],
            key["rep"],
            fgmres["total_iter"],
            fgmres["restarts"],
            fgmres["final_res"],
            fgmres["converged"],
        ]
        merged_rows.append(extra_row + row)

        # Non-numeric reps ("probe_tol1e-6", "probe_tol1e-11") are kept
        # separate so the subsystem breakdown shows each hand-run probe as
        # its own line instead of median-merging them with the real sweep.
        rep_is_numeric = key["rep"].isdigit()
        rep_suffix = "" if rep_is_numeric else f"_rep:{key['rep']}"
        config_key = (
            key["nx"],
            key["nz"],
            key["lin_solver"],
            key["n_threads"],
            rep_suffix,
        )
        row_dict = dict(zip(header, row))
        row_dict.update({"rep": key["rep"], **fgmres})
        config_rows[config_key].append(row_dict)

    if master_header is None:
        print(f"no perf.csv files found under {perf_detail}", file=sys.stderr)
        return 1

    with master_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(master_header)
        w.writerows(merged_rows)

    # Per-config subsystem breakdown in markdown
    lines = [
        "# MDOODZ perf.csv subsystem breakdown",
        "",
        f"Source: `{perf_detail}` ({len(merged_rows)} rows across "
        f"{len(config_rows)} configurations).",
        "",
        "Each row of MDOODZ's `perf.csv` is one simulated time step; the columns",
        "split the wall-clock into subsystem buckets (rheology update, Jacobian",
        "assembly, Stokes solve, particle interpolation, HDF5 output, …). This",
        "table reports the **step-1 medians** for the single-Picard-iter SolVi",
        "fixtures in the sweep, which is where all the interesting time lives.",
        "",
        "| nx  | nz  | lin_solver | threads | note | wall_s | solve_s | assembly_s | interp_s | output_s | rheology_s | stokes_setup_s | post_solve_s | nit | neq_mom | neq_cont | peak_rss_MB | fgmres_iters | fgmres_restarts | fgmres_final_res | fgmres_converged |",
        "|----:|----:|-----------:|--------:|:-----|-------:|--------:|-----------:|---------:|---------:|-----------:|---------------:|-------------:|----:|--------:|---------:|------------:|-------------:|----------------:|-----------------:|-----------------:|",
    ]

    def _fmt(x, fmt="{:.3f}"):
        try:
            return fmt.format(float(x))
        except (TypeError, ValueError):
            return str(x)

    for (nx, nz, lin, thr, note), rows in sorted(config_rows.items()):
        # Median across reps on step 1 (single-step SolVi fixtures)
        step1 = [r for r in rows if r.get("step") == "1"]
        if not step1:
            step1 = rows
        def med(field, fmt="{:.3f}"):
            vals = []
            for r in step1:
                try:
                    vals.append(float(r.get(field, "")))
                except (TypeError, ValueError):
                    pass
            if not vals:
                return "-"
            vals.sort()
            return fmt.format(vals[len(vals) // 2])

        def majority_bool(field):
            vals = [r.get(field) for r in step1]
            return "yes" if sum(1 for v in vals if v in (True, "True", "true")) > len(vals) / 2 else "no"

        def last_str(field):
            for r in step1[::-1]:
                v = r.get(field)
                if v not in (None, "", "-1"):
                    return str(v)
            return "-"

        lines.append(
            "| "
            + " | ".join(
                [
                    str(nx),
                    str(nz),
                    str(lin),
                    str(thr),
                    note or "-",
                    med("wall_s"),
                    med("solve_s"),
                    med("assembly_s"),
                    med("interp_s"),
                    med("output_s"),
                    med("rheology_s"),
                    med("stokes_setup_s"),
                    med("post_solve_s"),
                    med("nit", "{:.0f}"),
                    med("neq_mom", "{:.0f}"),
                    med("neq_cont", "{:.0f}"),
                    med("peak_rss_mb", "{:.1f}"),
                    last_str("total_iter"),
                    last_str("restarts"),
                    last_str("final_res"),
                    majority_bool("converged") if lin == 3 else "-",
                ]
            )
            + " |"
        )

    md_path.write_text("\n".join(lines) + "\n")
    print(f"wrote {master_path} ({len(merged_rows)} rows)")
    print(f"wrote {md_path}")
    return 0


def _key_dir_name(key: dict) -> str:
    return f"{key['nx']}x{key['nz']}_lin{key['lin_solver']}_thr{key['n_threads']}_rep{key['rep']}"


if __name__ == "__main__":
    sys.exit(main())
