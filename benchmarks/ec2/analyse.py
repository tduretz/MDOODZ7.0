#!/usr/bin/env python3
"""Aggregate results.csv into the defence-ready summary + regression plots.

Read one-row-per-run CSV, group by (sweep, nx, nz, lin_solver, n_threads),
drop rows flagged thermal_throttled, compute median + IQR for wall time
and peak RSS, fit log(RSS) = log(a) + alpha * log(DOF) per solver, write:

    results/summary.md
    results/figs/perf_memory_scaling.png
    results/figs/perf_walltime_comparison.png
    results/figs/perf_strong_scaling.png

Usage:
    python analyse.py results/results.csv [--git-sha SHA] [--out summary.md]

The `--git-sha` filter lets a reviewer re-aggregate a historical subset
(D3 in the perf-harness spec: per-row provenance enables future
re-aggregation without relaunching EC2).
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import sys
from collections import defaultdict
from statistics import median
from typing import Iterable


def _percentile(values: list[float], p: float) -> float:
    if not values:
        return float("nan")
    xs = sorted(values)
    k = (len(xs) - 1) * p
    lo, hi = math.floor(k), math.ceil(k)
    if lo == hi:
        return xs[int(k)]
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def iqr(values: list[float]) -> float:
    return _percentile(values, 0.75) - _percentile(values, 0.25)


def _linfit_loglog(xs: list[float], ys: list[float]) -> tuple[float, float]:
    """Least-squares fit log y = log a + alpha * log x. Returns (alpha, log_a)."""
    if len(xs) < 2:
        return float("nan"), float("nan")
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    n = len(lx)
    mx = sum(lx) / n
    my = sum(ly) / n
    num = sum((lx[i] - mx) * (ly[i] - my) for i in range(n))
    den = sum((lx[i] - mx) ** 2 for i in range(n))
    alpha = num / den if den > 0 else float("nan")
    log_a = my - alpha * mx
    return alpha, log_a


def _linfit_confint(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Return alpha, plus 95% CI half-width computed via Student-t on the slope."""
    if len(xs) < 3:
        a, _ = _linfit_loglog(xs, ys)
        return a, float("nan"), float("nan")
    alpha, log_a = _linfit_loglog(xs, ys)
    n = len(xs)
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    mx = sum(lx) / n
    resid = [ly[i] - (log_a + alpha * lx[i]) for i in range(n)]
    s2 = sum(r * r for r in resid) / (n - 2)
    denom = sum((lx[i] - mx) ** 2 for i in range(n))
    se = math.sqrt(s2 / denom) if denom > 0 else float("nan")
    # Student-t critical value for 95% CI, two-tailed, small n
    t_crit = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447, 7: 2.365,
              8: 2.306, 9: 2.262, 10: 2.228}.get(n - 2, 1.96)
    return alpha, t_crit * se, log_a


def load_rows(path: str, git_sha: str | None) -> list[dict]:
    out: list[dict] = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if git_sha and row.get("git_sha") != git_sha:
                continue
            out.append(row)
    return out


def aggregate(rows: Iterable[dict]):
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        if str(r.get("thermal_throttled", "")).strip().lower() in ("true", "1"):
            continue
        key = (int(r["grid_nx"]), int(r["grid_nz"]),
               int(r["lin_solver"]), int(r["n_threads"]))
        groups[key].append(r)
    summary = {}
    for key, grp in groups.items():
        walls = [float(r["wall_time_s"]) for r in grp]
        rss = [float(r["peak_rss_kb"]) for r in grp]
        summary[key] = {
            "n": len(grp),
            "wall_median": median(walls),
            "wall_iqr": iqr(walls),
            "rss_median": median(rss),
            "rss_iqr": iqr(rss),
            "fgmres_iters": median(int(r.get("n_its_fgmres", 0) or 0) for r in grp),
        }
    return summary


def write_summary(summary: dict, out_path: str, source_csv: str, warnings: list[str]) -> None:
    lines: list[str] = []
    lines.append("# Perf-sweep summary — `add-gmg-stokes-defence`\n")
    lines.append(f"Source: `{source_csv}`\n")
    lines.append("")

    if warnings:
        lines.append("## WARNINGS\n")
        lines.extend(f"- {w}" for w in warnings)
        lines.append("")

    # Memory-scaling alpha fits per solver (single-thread only).
    lines.append("## Memory-scaling exponent\n")
    lines.append("Fit: `log RSS = log a + alpha * log DOF` across memory-scaling grids, n_threads = 1.\n")
    for lin_solver in (0, 3):
        pts = [(nx * nz, summary[(nx, nz, lin_solver, 1)]["rss_median"])
               for (nx, nz, s, t) in summary
               if s == lin_solver and t == 1]
        if len(pts) >= 2:
            alpha, ci, _ = _linfit_confint([p[0] for p in pts], [p[1] for p in pts])
            lines.append(f"- `lin_solver = {lin_solver}`: alpha = {alpha:.3f} +/- {ci:.3f} (95% CI), "
                         f"from {len(pts)} grid points")
        else:
            lines.append(f"- `lin_solver = {lin_solver}`: insufficient points ({len(pts)}) for fit")
    lines.append("")

    # 201 wall-time ratio.
    lines.append("## 201 x 201 wall-time comparison\n")
    cholmod_201 = summary.get((201, 201, 0, 1))
    gmg_201 = summary.get((201, 201, 3, 1))
    if cholmod_201 and gmg_201:
        ratio = cholmod_201["wall_median"] / gmg_201["wall_median"]
        lines.append(f"- CHOLMOD median wall time: {cholmod_201['wall_median']:.3f} s  "
                     f"(IQR {cholmod_201['wall_iqr']:.3f})")
        lines.append(f"- GMG median wall time:     {gmg_201['wall_median']:.3f} s  "
                     f"(IQR {gmg_201['wall_iqr']:.3f})")
        lines.append(f"- Ratio t(CHOLMOD) / t(GMG) = {ratio:.2f}x")
        if ratio < 3.0:
            warnings.append(f"201 wall-time ratio {ratio:.2f}x is below the >= 3x regression bound")
    else:
        lines.append("- 201x201 wall-time data missing; re-run the walltime_comparison sweep.")
    lines.append("")

    # Strong scaling.
    lines.append("## Strong-scaling at 201 x 201, lin_solver = 3\n")
    scaling_pts = sorted([(t, summary[(201, 201, 3, t)]["wall_median"])
                          for (nx, nz, s, t) in summary
                          if s == 3 and nx == 201 and nz == 201])
    if scaling_pts:
        t1 = next((w for (t, w) in scaling_pts if t == 1), None)
        for t, w in scaling_pts:
            eff = 100.0 * (t1 / w) / t if t1 else float("nan")
            lines.append(f"- threads = {t:>3d}: wall = {w:.3f} s, efficiency = {eff:5.1f}%")
    else:
        lines.append("- strong-scaling data missing; re-run the strong_scaling sweep.")
    lines.append("")

    # Per-configuration table (full detail).
    lines.append("## All configurations\n")
    lines.append("| nx  | nz  | lin_solver | n_threads | n | wall_median (s) | wall_iqr | rss_median (kB) | rss_iqr | fgmres_iters |")
    lines.append("|----:|----:|-----------:|----------:|--:|----------------:|---------:|----------------:|--------:|-------------:|")
    for key in sorted(summary):
        nx, nz, s, t = key
        r = summary[key]
        lines.append(f"| {nx} | {nz} | {s} | {t} | {r['n']} | {r['wall_median']:.4f} | {r['wall_iqr']:.4f} | "
                     f"{int(r['rss_median'])} | {int(r['rss_iqr'])} | {int(r['fgmres_iters'])} |")
    lines.append("")

    with open(out_path, "w") as f:
        f.write("\n".join(lines))


def _try_plot(summary: dict, figs_dir: str) -> list[str]:
    """Produce three PNGs if matplotlib is available; return filenames written."""
    try:
        import matplotlib  # type: ignore

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt  # type: ignore
    except ImportError:
        return []
    os.makedirs(figs_dir, exist_ok=True)
    written: list[str] = []

    # 1. Memory scaling.
    fig, ax = plt.subplots(figsize=(5, 4))
    for lin_solver, marker, colour in ((0, "o", "#B03030"), (3, "s", "#2B50A0")):
        pts = sorted([(nx * nz, summary[(nx, nz, lin_solver, 1)]["rss_median"])
                      for (nx, nz, s, t) in summary
                      if s == lin_solver and t == 1])
        if len(pts) >= 2:
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            ax.loglog(xs, ys, marker=marker, color=colour, label=f"lin_solver = {lin_solver}")
            alpha, _, _ = _linfit_confint(xs, ys)
            ax.text(xs[-1], ys[-1] * 1.1, fr"$\alpha = {alpha:.2f}$", color=colour)
    ax.set_xlabel("DOF (nx * nz)")
    ax.set_ylabel("Peak RSS (kB)")
    ax.set_title("Memory scaling")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(figs_dir, "perf_memory_scaling.png")
    fig.savefig(p, dpi=140)
    plt.close(fig)
    written.append(p)

    # 2. Wall-time comparison at 201.
    if (201, 201, 0, 1) in summary and (201, 201, 3, 1) in summary:
        fig, ax = plt.subplots(figsize=(4, 4))
        ws = [summary[(201, 201, 0, 1)]["wall_median"],
              summary[(201, 201, 3, 1)]["wall_median"]]
        ax.bar(["CHOLMOD", "GMG"], ws, color=["#B03030", "#2B50A0"])
        ax.set_ylabel("Wall time (s)")
        ax.set_title("201 x 201, n_threads = 1")
        fig.tight_layout()
        p = os.path.join(figs_dir, "perf_walltime_comparison.png")
        fig.savefig(p, dpi=140)
        plt.close(fig)
        written.append(p)

    # 3. Strong scaling.
    scaling_pts = sorted([(t, summary[(201, 201, 3, t)]["wall_median"])
                          for (nx, nz, s, t) in summary
                          if s == 3 and nx == 201 and nz == 201])
    if len(scaling_pts) >= 2:
        fig, ax = plt.subplots(figsize=(5, 4))
        ts = [p[0] for p in scaling_pts]
        ws = [p[1] for p in scaling_pts]
        eff = [100.0 * (ws[0] / ws[i]) / ts[i] for i in range(len(ts))]
        ax.plot(ts, eff, marker="o", color="#2B50A0")
        ax.axhline(100.0, color="#808080", linestyle="--")
        ax.set_xlabel("n_threads")
        ax.set_ylabel("Efficiency (%)")
        ax.set_title("Strong scaling @ 201 x 201")
        ax.set_xscale("log", base=2)
        fig.tight_layout()
        p = os.path.join(figs_dir, "perf_strong_scaling.png")
        fig.savefig(p, dpi=140)
        plt.close(fig)
        written.append(p)

    return written


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("csv", help="path to results.csv")
    ap.add_argument("--git-sha", default=None, help="filter rows by git_sha")
    ap.add_argument("--out", default=None, help="output summary.md path")
    ap.add_argument("--figs-dir", default=None, help="output figures directory")
    args = ap.parse_args(argv[1:])

    rows = load_rows(args.csv, args.git_sha)
    if not rows:
        print("no rows found (check --git-sha filter)", file=sys.stderr)
        return 1

    summary = aggregate(rows)

    warnings: list[str] = []
    # Throttle-fraction warning per configuration.
    throttled = defaultdict(int)
    total = defaultdict(int)
    for r in rows:
        key = (int(r["grid_nx"]), int(r["grid_nz"]),
               int(r["lin_solver"]), int(r["n_threads"]))
        total[key] += 1
        if str(r.get("thermal_throttled", "")).strip().lower() in ("true", "1"):
            throttled[key] += 1
    for key, t in throttled.items():
        frac = t / total[key]
        if frac > 0.20:
            warnings.append(f"configuration {key} had {frac*100:.0f}% throttled runs "
                            "(> 20% threshold) — consider rescheduling or a different instance type")

    out_path = args.out or os.path.join(os.path.dirname(args.csv) or ".", "summary.md")
    figs_dir = args.figs_dir or os.path.join(os.path.dirname(args.csv) or ".", "figs")
    write_summary(summary, out_path, args.csv, warnings)
    written = _try_plot(summary, figs_dir)
    print(f"wrote {out_path} ({len(summary)} configurations)")
    for w in written:
        print(f"wrote {w}")
    for w in warnings:
        print(f"WARNING: {w}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
