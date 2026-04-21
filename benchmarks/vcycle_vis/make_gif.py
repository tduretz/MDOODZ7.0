#!/usr/bin/env python3
"""V-cycle residual animation for the add-gmg-stokes-defence defence.

Reads an HDF5 snapshot tree produced by `gmg_dump_vcycle = 1`
(MDLIB/MultigridStokes.c §2.3-§2.5) and produces:

  * figs/vcycle_sol_cx_41.gif        — 2 fps animated GIF of the whole
                                       V-cycle, one frame per operator
                                       snapshot (pre-smooth pre, post;
                                       restrict pre, post; coarse-solve
                                       pre, post; prolongate pre, post;
                                       post-smooth pre, post; all levels).
  * figs/vcycle_still_{01..06}.png   — Six key stills for PDF embedding:
                                         01 pre-smooth       (level 0, pre)
                                         02 restrict-1       (level 1, post)
                                         03 coarse-solve     (coarsest, post)
                                         04 prolongate-1     (level 1, post)
                                         05 post-smooth      (level 0, post)
                                         06 converged        (level 0, last)

Each frame has two panels:
  * Left  — residual-magnitude field, symmetric-log colourmap (viridis).
  * Right — radial power spectrum of the residual, log-log axes. The
            spectrum is computed from |FFT2(res_magnitude)| binned
            radially; this is the standard "which frequencies does
            the smoother kill, which does coarse-grid correction
            kill" diagnostic (§10.3 in the defence).

Input dump schema (set by `gmg_vcycle_dump_snapshot` in
MultigridStokes.c):
  /Level/dims       int[4]  = { Nx, Nz, Ncx, Ncz }
  /Level/spacing    dbl[2]  = { dx, dz }
  /Level/level      int     = multigrid level index (0 = fine)
  /Fields/Vx, Vz, P, res_u, res_v, res_p
  /Meta/step        c-str   = {"pre_smooth","restrict","coarse_solve","prolongate","post_smooth"}
  /Meta/phase       c-str   = {"pre","post"}
  /Meta/seq         int     = global traversal order

Usage:
    python3 make_gif.py <path/to/vcycle_dump> [--output-dir FIGS] [--fps 2]
"""

from __future__ import annotations

import argparse
import io
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import SymLogNorm
from PIL import Image


FRAME_RE = re.compile(
    r"^level_(?P<level>\d+)_(?P<step>[a-z_]+?)_(?P<phase>pre|post)_(?P<seq>\d+)\.h5$"
)


@dataclass(frozen=True)
class Frame:
    path: Path
    level: int
    step: str
    phase: str
    seq: int

    @property
    def label(self) -> str:
        return f"L{self.level} {self.step} ({self.phase})  seq={self.seq:03d}"


def discover_frames(dump_dir: Path) -> List[Frame]:
    frames: list[Frame] = []
    for h5_path in sorted(dump_dir.glob("level_*.h5")):
        m = FRAME_RE.match(h5_path.name)
        if not m:
            print(f"warn: skipping non-matching file {h5_path.name}", file=sys.stderr)
            continue
        frames.append(
            Frame(
                path=h5_path,
                level=int(m["level"]),
                step=m["step"],
                phase=m["phase"],
                seq=int(m["seq"]),
            )
        )
    frames.sort(key=lambda f: f.seq)
    return frames


def first_vcycle_slice(frames: List[Frame]) -> List[Frame]:
    """Return the subset of `frames` belonging to the first V-cycle only.

    A V-cycle starts with `level 0 pre_smooth pre` and ends with
    `level 0 post_smooth post`. The dumper arms once per SolveStokesGMG
    call and records every V-cycle that runs inside the FGMRES outer
    loop, so for a typical solve we get 15-30 V-cycles' worth of
    snapshots. For the animation we only want the first traversal so
    the spectrum-collapse story is clean. We scan for the earliest
    (pre_smooth, pre, level=0) → (post_smooth, post, level=0) boundary.
    """
    if not frames:
        return frames
    start_idx: Optional[int] = None
    for i, f in enumerate(frames):
        if f.level == 0 and f.step == "pre_smooth" and f.phase == "pre":
            start_idx = i
            break
    if start_idx is None:
        return frames
    end_idx = len(frames) - 1
    for j in range(start_idx + 1, len(frames)):
        f = frames[j]
        if f.level == 0 and f.step == "post_smooth" and f.phase == "post":
            end_idx = j
            break
    return frames[start_idx : end_idx + 1]


def load_cell_residual(frame: Frame) -> tuple[np.ndarray, float, float]:
    """Interpolate res_u (Vx-nodes), res_v (Vz-nodes), res_p (centres) onto
    a common cell grid and return its magnitude + cell spacing.

    res_u has shape (Nx, Ncz),   res_v has shape (Ncx, Nz),
    res_p has shape (Ncx, Ncz).

    We average res_u onto cell centres along x (Vx-node straddles vertical
    cell faces) and res_v onto cell centres along z, then compute
    sqrt(res_u_c^2 + res_v_c^2 + res_p^2) as a single scalar residual
    magnitude per cell — the standard diagnostic for multigrid plots.
    """
    with h5py.File(frame.path, "r") as f:
        dims = f["Level/dims"][:]
        spacing = f["Level/spacing"][:]
        Nx, Nz, Ncx, Ncz = int(dims[0]), int(dims[1]), int(dims[2]), int(dims[3])
        res_u = np.asarray(f["Fields/res_u"][:]).reshape(Ncz, Nx).T
        res_v = np.asarray(f["Fields/res_v"][:]).reshape(Nz, Ncx).T
        res_p = np.asarray(f["Fields/res_p"][:]).reshape(Ncz, Ncx).T
    # res_u sits on vertical faces: shape (Nx, Ncz). Average the two
    # x-neighbours onto the (Ncx, Ncz) centres.
    res_u_c = 0.5 * (res_u[:-1, :] + res_u[1:, :])
    # res_v sits on horizontal faces: shape (Ncx, Nz). Average the two
    # z-neighbours onto the (Ncx, Ncz) centres.
    res_v_c = 0.5 * (res_v[:, :-1] + res_v[:, 1:])
    mag = np.sqrt(res_u_c**2 + res_v_c**2 + res_p**2)
    return mag, float(spacing[0]), float(spacing[1])


def radial_spectrum(field: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (k_bins, P_k) for the radial power spectrum of `field`.
    Uses 2D FFT; bins radially over integer wavenumber indices."""
    nx, nz = field.shape
    F = np.fft.fftshift(np.fft.fft2(field - field.mean()))
    P = np.abs(F) ** 2
    kx = np.fft.fftshift(np.fft.fftfreq(nx)) * nx
    kz = np.fft.fftshift(np.fft.fftfreq(nz)) * nz
    KX, KZ = np.meshgrid(kx, kz, indexing="ij")
    K = np.sqrt(KX**2 + KZ**2)
    # Integer-binned radial average.
    k_max = int(np.floor(K.max()))
    if k_max < 1:
        return np.array([1.0]), np.array([P.sum()])
    bins = np.arange(k_max + 1)
    radial_sum = np.zeros(k_max + 1)
    radial_cnt = np.zeros(k_max + 1, dtype=int)
    K_idx = np.clip(K.astype(int), 0, k_max)
    for b in range(k_max + 1):
        mask = K_idx == b
        if not np.any(mask):
            continue
        radial_sum[b] = P[mask].sum()
        radial_cnt[b] = int(mask.sum())
    radial_avg = np.where(radial_cnt > 0, radial_sum / np.maximum(radial_cnt, 1), 0)
    # drop k=0 DC and zero-count bins for log-log plotting.
    mask_nz = (radial_avg > 0) & (bins > 0)
    return bins[mask_nz].astype(float), radial_avg[mask_nz]


def render_frame(
    frame: Frame,
    mag: np.ndarray,
    dx: float,
    dz: float,
    vmin: float,
    vmax: float,
    step_idx: int,
    n_steps: int,
) -> Image.Image:
    """Render one 2-panel frame and return it as a PIL Image."""
    fig, (ax_field, ax_spec) = plt.subplots(
        1, 2, figsize=(10, 4), dpi=120, gridspec_kw={"width_ratios": [1.1, 1.0]}
    )
    nx, nz = mag.shape
    # x = j*dx on the centres, z = k*dz on the centres.
    ext = (0.0, nx * dx, 0.0, nz * dz)
    norm = SymLogNorm(
        linthresh=max(vmax, 1e-30) * 1e-6,
        linscale=0.5,
        vmin=vmin,
        vmax=vmax,
    )
    im = ax_field.imshow(
        mag.T,
        origin="lower",
        extent=ext,
        cmap="viridis",
        norm=norm,
        aspect="equal",
    )
    ax_field.set_title(f"residual magnitude — {frame.label}", fontsize=10)
    ax_field.set_xlabel("x")
    ax_field.set_ylabel("z")
    plt.colorbar(im, ax=ax_field, fraction=0.042, pad=0.04)

    k_bins, P_k = radial_spectrum(mag)
    ax_spec.loglog(k_bins, P_k, "-o", color="#4b3b8f", markersize=3)
    ax_spec.set_title("radial power spectrum of residual", fontsize=10)
    ax_spec.set_xlabel("wavenumber  k")
    ax_spec.set_ylabel(r"$\langle|\hat r|^2\rangle_k$")
    ax_spec.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        f"V-cycle step {step_idx+1}/{n_steps}  —  SolCx 41×41  —  "
        f"lin_solver=3, gmg_dump_vcycle=1",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    img = Image.open(buf).convert("RGB").copy()
    buf.close()
    return img


def choose_stills(
    vcycle_frames: List[Frame],
    all_frames: List[Frame],
) -> List[Optional[Frame]]:
    """Select six representative stills:
      01  level 0 pre-smooth  (pre)      — opening residual
      02  level 1 restrict    (post)     — 1st-level restriction
      03  coarsest coarse_solve (post)   — direct solve bottom
      04  level 1 prolongate  (post)     — 1st-level prolongation
      05  level 0 post_smooth (post)     — 1st V-cycle exit residual
      06  level 0 post_smooth (post) LAST — converged tail (final V-cycle)

    Picks 01–05 from the first V-cycle; 06 from the full dump (so the
    "converged" still lives many V-cycles after the opener and the
    spectrum collapse is visible when 01 and 06 sit side-by-side).
    """
    def first_in(src: List[Frame], **kw) -> Optional[Frame]:
        for f in src:
            if all(getattr(f, k) == v for k, v in kw.items()):
                return f
        return None

    def last_in(src: List[Frame], **kw) -> Optional[Frame]:
        match: Optional[Frame] = None
        for f in src:
            if all(getattr(f, k) == v for k, v in kw.items()):
                match = f
        return match

    if not vcycle_frames:
        return [None] * 6
    max_level = max(f.level for f in vcycle_frames)
    return [
        first_in(vcycle_frames, level=0, step="pre_smooth", phase="pre"),
        first_in(vcycle_frames, level=1, step="restrict", phase="post"),
        first_in(vcycle_frames, level=max_level, step="coarse_solve", phase="post"),
        first_in(vcycle_frames, level=1, step="prolongate", phase="post"),
        first_in(vcycle_frames, level=0, step="post_smooth", phase="post"),
        last_in(all_frames, level=0, step="post_smooth", phase="post") or all_frames[-1],
    ]


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("dump_dir", type=Path, help="path to vcycle_dump/ directory")
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figs"),
        help="destination directory for GIF + stills (default: figs/)",
    )
    p.add_argument(
        "--fps", type=float, default=2.0, help="GIF frame rate (default: 2 fps)"
    )
    p.add_argument(
        "--gif-name",
        default="vcycle_sol_cx_41.gif",
        help="output GIF filename (default: vcycle_sol_cx_41.gif)",
    )
    p.add_argument(
        "--still-prefix",
        default="vcycle_still_",
        help="stills filename prefix (default: vcycle_still_)",
    )
    p.add_argument(
        "--all-vcycles",
        action="store_true",
        help="include every V-cycle in the dump tree (default: first only)",
    )
    args = p.parse_args()

    if not args.dump_dir.is_dir():
        print(f"error: {args.dump_dir} is not a directory", file=sys.stderr)
        return 1
    args.output_dir.mkdir(parents=True, exist_ok=True)

    all_frames = discover_frames(args.dump_dir)
    if not all_frames:
        print(f"error: no V-cycle dumps found in {args.dump_dir}", file=sys.stderr)
        return 1
    total = len(all_frames)
    frames = all_frames if args.all_vcycles else first_vcycle_slice(all_frames)
    print(
        f"loaded {total} total snapshots from {args.dump_dir}; "
        f"rendering {len(frames)} frame(s) "
        f"({'all V-cycles' if args.all_vcycles else 'first V-cycle only'})"
    )

    # Pre-pass: compute global min/max residual magnitude so the colour
    # normalisation is consistent across frames — frame N should be
    # visually comparable to frame M.
    mags: list[tuple[np.ndarray, float, float]] = []
    for f in frames:
        mag, dx, dz = load_cell_residual(f)
        mags.append((mag, dx, dz))
    all_vals = np.concatenate([m.ravel() for m, _, _ in mags])
    vmax = float(np.nanmax(np.abs(all_vals))) or 1e-30
    vmin = float(max(1e-30, np.nanmin(np.abs(all_vals[all_vals != 0]))))

    # Render every frame and assemble the GIF.
    images: list[Image.Image] = []
    n = len(frames)
    for i, (f, (mag, dx, dz)) in enumerate(zip(frames, mags)):
        img = render_frame(f, mag, dx, dz, vmin, vmax, i, n)
        images.append(img)

    duration_ms = int(1000.0 / args.fps)
    gif_path = args.output_dir / args.gif_name
    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        duration=duration_ms,
        loop=0,
        optimize=False,
    )
    print(f"wrote {gif_path}  ({n} frames at {args.fps} fps)")

    # Six-frame still sequence for the printable PDF.
    # Stills 01-05 come from the (possibly sliced) first V-cycle; still 06
    # comes from the full dump so it captures the converged-tail residual.
    stills = choose_stills(frames, all_frames)
    for idx, still_frame in enumerate(stills, start=1):
        if still_frame is None:
            print(f"warn: still slot {idx} not available (sparse hierarchy)")
            continue
        # Use cached magnitudes when the frame is already loaded.
        if still_frame in frames:
            j = frames.index(still_frame)
            mag, dx, dz = mags[j]
            j_idx, j_n = j, n
        else:
            mag, dx, dz = load_cell_residual(still_frame)
            j_idx, j_n = all_frames.index(still_frame), total
        img = render_frame(still_frame, mag, dx, dz, vmin, vmax, j_idx, j_n)
        out = args.output_dir / f"{args.still_prefix}{idx:02d}.png"
        img.save(out, "PNG", dpi=(140, 140))
        print(f"wrote {out}  ({still_frame.label})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
