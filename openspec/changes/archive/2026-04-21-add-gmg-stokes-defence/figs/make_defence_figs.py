#!/usr/bin/env python3
"""Generate the five process diagrams for the add-gmg-stokes-defence
defence document (tasks 9.2-9.6).

Style: primer-series (matplotlib, 140 DPI, 2D, muted palette). Each
figure is a standalone vector-style flow/architecture diagram — no
physical data, no external dependencies beyond matplotlib + numpy.

Outputs (all under openspec/changes/add-gmg-stokes-defence/figs/):

  defence_vcycle_tree.png        (task 9.2)
  defence_fgmres_controlflow.png (task 9.3)
  defence_stencil_bridge.png     (task 9.4)
  defence_picard_newton.png      (task 9.5)
  defence_active_mask.png        (task 9.6)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

# Muted primer palette — same tones the primer.pdf / codes.pdf series use.
PALETTE = {
    "bg":        "#fbfaf6",
    "ink":       "#1f2a37",
    "accent":    "#4b3b8f",
    "warm":      "#b04b3b",
    "cool":      "#3b7d8a",
    "soft":      "#d9d3c5",
    "muted":     "#7a6f5a",
    "line":      "#5a5248",
    "highlight": "#d9a64b",
}

DPI = 140


def _styled_axes(ax, xlim: Tuple[float, float], ylim: Tuple[float, float]) -> None:
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.set_facecolor(PALETTE["bg"])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def _box(ax, x: float, y: float, w: float, h: float, label: str,
         fill: str = "#ffffff", edge: str = None,
         fontsize: float = 9.5, bold: bool = False) -> None:
    edge = edge or PALETTE["line"]
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.08",
        linewidth=1.1, facecolor=fill, edgecolor=edge,
    )
    ax.add_patch(patch)
    ax.text(
        x, y, label,
        ha="center", va="center",
        fontsize=fontsize, color=PALETTE["ink"],
        fontweight="bold" if bold else "normal",
    )


def _arrow(ax, xy_from: Tuple[float, float], xy_to: Tuple[float, float],
           color: str = None, style: str = "-|>", lw: float = 1.2,
           connectionstyle: str = "arc3,rad=0") -> None:
    color = color or PALETTE["line"]
    ax.add_patch(FancyArrowPatch(
        xy_from, xy_to,
        arrowstyle=style, color=color, lw=lw,
        connectionstyle=connectionstyle,
        mutation_scale=14,
    ))


def _header(fig, title: str) -> None:
    fig.suptitle(title, fontsize=12, color=PALETTE["ink"], y=0.98)


# ---------------------------------------------------------------------------
# Figure 1 — V-cycle recursion tree (task 9.2).
# ---------------------------------------------------------------------------

def fig_vcycle_tree(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 5.2), dpi=DPI)
    _styled_axes(ax, (-4.0, 11.0), (-0.2, 5.8))
    fig.patch.set_facecolor(PALETTE["bg"])

    levels = [
        (0, "Level 0  —  fine  (321²)"),
        (1, "Level 1  —  161²"),
        (2, "Level 2  —  81²"),
        (3, "Level 3  —  41²"),
        (4, "Level 4  —  coarsest  (21²)"),
    ]

    for k, (lvl, label) in enumerate(levels):
        y = 5.0 - k * 1.0
        ax.axhline(y, color=PALETTE["soft"], lw=0.5, zorder=1, xmin=0.3)
        ax.text(-3.8, y, label, ha="left", va="center",
                fontsize=8.5, color=PALETTE["muted"])

    def down_node(lvl: int, x: float, label: str, coarsest: bool = False):
        y = 5.0 - lvl
        fill = PALETTE["warm"] if coarsest else "#ffffff"
        edge = PALETTE["warm"] if coarsest else PALETTE["line"]
        _box(ax, x, y, 1.6, 0.52, label, fill=fill, edge=edge,
             fontsize=8.2, bold=coarsest)

    def up_node(lvl: int, x: float, label: str):
        y = 5.0 - lvl
        _box(ax, x, y, 1.6, 0.52, label, fill=PALETTE["soft"],
             edge=PALETTE["line"], fontsize=8.2)

    down_node(0, 1.2, r"pre-smooth  $\nu_{pre}=2$")
    down_node(1, 2.7, r"pre-smooth  $\nu_{pre}=2$")
    down_node(2, 4.2, r"pre-smooth  $\nu_{pre}=2$")
    down_node(3, 5.7, r"pre-smooth  $\nu_{pre}=2$")
    down_node(4, 7.2, "CHOLMOD  (K_coarse LU)", coarsest=True)
    up_node(3, 7.7, r"post-smooth  $\nu_{post}=2$")
    up_node(2, 8.2, r"post-smooth  $\nu_{post}=2$")
    up_node(1, 8.7, r"post-smooth  $\nu_{post}=2$")
    up_node(0, 9.2, r"post-smooth  $\nu_{post}=2$")

    _arrow(ax, (1.95, 5.0), (2.50, 4.0))
    _arrow(ax, (3.45, 4.0), (4.00, 3.0))
    _arrow(ax, (4.95, 3.0), (5.50, 2.0))
    _arrow(ax, (6.45, 2.0), (7.00, 1.0), color=PALETTE["warm"])
    _arrow(ax, (7.40, 1.0), (7.70, 2.0), color=PALETTE["accent"])
    _arrow(ax, (7.95, 2.0), (8.20, 3.0), color=PALETTE["accent"])
    _arrow(ax, (8.45, 3.0), (8.70, 4.0), color=PALETTE["accent"])
    _arrow(ax, (8.95, 4.0), (9.20, 5.0), color=PALETTE["accent"])

    ax.text(2.2, 5.55, "restrict  r→coarse", fontsize=7.5,
            color=PALETTE["muted"], ha="left")
    ax.text(7.6, 5.55, "prolongate  δ→fine", fontsize=7.5,
            color=PALETTE["accent"], ha="left")

    ax.text(7.2, 0.3,
            r"coarsest direct solve costs $O(N^{1.5})$ FLOPs  <  0.1% of V-cycle",
            fontsize=7.5, color=PALETTE["warm"], ha="center", fontstyle="italic")

    _header(fig, "V-cycle recursion — pre/post smooth + restrict/prolongate + coarse direct solve")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=DPI, facecolor=PALETTE["bg"])
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Figure 2 — FGMRES-with-GMG-preconditioner control flow (task 9.3).
# ---------------------------------------------------------------------------

def fig_fgmres_controlflow(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.2), dpi=DPI)
    _styled_axes(ax, (0, 10), (0, 6))
    fig.patch.set_facecolor(PALETTE["bg"])

    _box(ax, 1.3, 5.2, 2.1, 0.7, "init:  $x_0 = 0$\n$r_0 = b - A x_0$",
         fill=PALETTE["soft"], fontsize=8.5)
    _box(ax, 4.5, 5.2, 2.4, 0.7, "Arnoldi:  build V_k, H_k\nfrom  $A M^{-1}$  columns",
         fill="#ffffff", fontsize=8.5)
    _box(ax, 8.2, 5.2, 2.1, 0.7, "least-squares\nsolve  $\\min\\|H_k y - β e_1\\|$",
         fill="#ffffff", fontsize=8.5)

    _arrow(ax, (2.35, 5.2), (3.30, 5.2))
    _arrow(ax, (5.70, 5.2), (7.15, 5.2))
    _arrow(ax, (8.2, 4.85), (8.2, 4.25))

    _box(ax, 8.2, 3.9, 2.4, 0.7, "update  $x_k = x_0 + M^{-1} V_k y$\n$r_k = b - A x_k$",
         fill=PALETTE["soft"], fontsize=8.5)

    _box(ax, 4.9, 3.9, 3.0, 0.95,
         r"$M^{-1} z$  =  V-cycle(H, 0)"
         "\n\n applied to each Krylov column\nreturns approximate $M^{-1} v_j$",
         fill=PALETTE["accent"] + "22", edge=PALETTE["accent"],
         fontsize=8.3)
    _arrow(ax, (5.70, 5.0), (5.40, 4.40), color=PALETTE["accent"],
           connectionstyle="arc3,rad=-0.2")
    _arrow(ax, (5.40, 4.40), (5.70, 5.0), color=PALETTE["accent"],
           connectionstyle="arc3,rad=-0.2", style="<|-")
    ax.text(4.25, 4.45, "(preconditioner application)",
            fontsize=7.3, color=PALETTE["accent"], style="italic")

    _box(ax, 1.3, 3.9, 2.1, 0.7,
         r"converged? $\|r_k\| \leq \tau$",
         fill="#ffffff", edge=PALETTE["warm"], fontsize=8.5)
    _arrow(ax, (6.95, 3.9), (2.40, 3.9))

    _arrow(ax, (0.25, 3.9), (0.25, 5.2), color=PALETTE["warm"],
           connectionstyle="arc3,rad=0.3")
    _arrow(ax, (0.25, 3.9), (0.80, 3.9), color=PALETTE["warm"],
           connectionstyle="arc3,rad=0")
    _arrow(ax, (0.80, 3.9), (0.25, 3.9), color=PALETTE["warm"])
    _arrow(ax, (0.25, 5.2), (0.80, 5.2), color=PALETTE["warm"])
    ax.text(0.15, 4.55, "restart\n(if not converged\nand k == m)",
            fontsize=7.0, color=PALETTE["warm"], ha="left")

    _box(ax, 1.3, 2.5, 2.1, 0.6,
         "done  —  $x_k$ accepted",
         fill=PALETTE["warm"] + "33", edge=PALETTE["warm"],
         fontsize=8.5, bold=True)
    _arrow(ax, (1.3, 3.55), (1.3, 2.82), color=PALETTE["warm"])

    ax.text(5.0, 1.4,
            "FGMRES  ≡  flexible variant of GMRES: $M^{-1}$ may change\n"
            "across iterations. Essential for GMG preconditioning — each\n"
            "V-cycle is only approximately $A^{-1}$, so the \"preconditioner\"\n"
            "is iteration-dependent and we need FGMRES's extended basis.",
            fontsize=8.0, color=PALETTE["muted"], ha="center")

    _header(fig, "FGMRES  +  GMG-preconditioner  —  outer Krylov  ×  inner V-cycle")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=DPI, facecolor=PALETTE["bg"])
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Figure 3 — Stencil-bridge architecture (task 9.4).
# ---------------------------------------------------------------------------

def fig_stencil_bridge(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.5), dpi=DPI)
    _styled_axes(ax, (0, 10), (0, 6))
    fig.patch.set_facecolor(PALETTE["bg"])

    _box(ax, 1.7, 5.2, 2.8, 0.8,
         "MDOODZ assembler\nStokesAssemblyDecoupled.c",
         fill=PALETTE["soft"], fontsize=8.8)

    _box(ax, 5.0, 5.2, 2.8, 0.8,
         "BuildStokesOperator\nDecoupled\n(sparse matrix K)",
         fill="#ffffff", fontsize=8.5)

    _box(ax, 8.5, 5.2, 1.4, 0.8,
         "CHOLMOD\ndirect\nfactorise",
         fill=PALETTE["warm"] + "22", edge=PALETTE["warm"], fontsize=8.3)

    _arrow(ax, (3.1, 5.2), (3.60, 5.2))
    _arrow(ax, (6.40, 5.2), (7.80, 5.2))

    _box(ax, 1.7, 3.3, 2.8, 0.8,
         "stencil bundle\n(per-cell coefficients)",
         fill=PALETTE["soft"], fontsize=8.8)
    _arrow(ax, (1.7, 4.80), (1.7, 3.70))

    _box(ax, 5.0, 3.3, 2.8, 0.8,
         "ApplyStokesOperator\nMDOODZ  (matvec, GMG-side)",
         fill=PALETTE["accent"] + "22", edge=PALETTE["accent"], fontsize=8.3)
    _arrow(ax, (3.1, 3.3), (3.60, 3.3), color=PALETTE["accent"])

    _box(ax, 5.0, 1.7, 2.8, 0.8,
         "fine-level V-cycle\n(smoother + Galerkin restrict)",
         fill="#ffffff", edge=PALETTE["accent"], fontsize=8.3)
    _arrow(ax, (5.0, 2.90), (5.0, 2.10), color=PALETTE["accent"])

    _box(ax, 8.5, 1.7, 1.4, 0.8,
         "GMG\nFGMRES",
         fill=PALETTE["accent"] + "22", edge=PALETTE["accent"],
         fontsize=8.3, bold=True)
    _arrow(ax, (6.40, 1.7), (7.80, 1.7), color=PALETTE["accent"])

    _box(ax, 5.0, 0.3, 5.2, 0.55,
         r"golden cross-check test:  $\|A_{\mathrm{assemble}} x - \mathrm{matvec}(x)\|_\infty \leq 10^{-12}$",
         fill=PALETTE["highlight"] + "55", edge=PALETTE["warm"], fontsize=8.0)

    ax.text(0.3, 4.4, "the bridge:",
            fontsize=9.0, color=PALETTE["muted"], fontweight="bold")
    ax.text(0.3, 3.95,
            "assembler and\nGMG matvec both\nconsume the same\nstencil bundle",
            fontsize=7.5, color=PALETTE["muted"], va="top")

    _header(fig,
            "Decoupled Stokes + stencil-bridge  —  one coefficient source, two operators, one golden check")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out, dpi=DPI, facecolor=PALETTE["bg"])
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Figure 4 — Picard → Newton phasing (task 9.5).
# ---------------------------------------------------------------------------

def fig_picard_newton(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 5.0), dpi=DPI)
    _styled_axes(ax, (0, 12.5), (0, 6.0))
    fig.patch.set_facecolor(PALETTE["bg"])

    _box(ax, 1.0, 4.6, 1.6, 0.7,
         "Picard it.\n$i = 0$", fill=PALETTE["soft"], fontsize=8.2)

    ax.add_patch(Rectangle(
        (2.2, 3.7), 4.3, 1.8,
        facecolor=PALETTE["accent"] + "18",
        edgecolor=PALETTE["accent"], lw=1.0,
    ))
    ax.text(4.35, 5.30, "Phase 1   —   Picard mode",
            ha="center", fontsize=9.0, color=PALETTE["accent"], fontweight="bold")
    ax.text(4.35, 4.88, r"SPD operator $A_{\mathrm{Pic}}$",
            ha="center", fontsize=8.2, color=PALETTE["ink"])
    ax.text(4.35, 4.55, "GMG-FGMRES path eligible (lin_solver = 3)",
            ha="center", fontsize=8.2, color=PALETTE["ink"])
    ax.text(4.35, 4.22, "Vanka smoother + direct coarse (SPD)",
            ha="center", fontsize=8.2, color=PALETTE["muted"])
    ax.text(4.35, 3.90, "fall back → CHOLMOD if GMG diverges",
            ha="center", fontsize=7.5, color=PALETTE["warm"], fontstyle="italic")

    ax.add_patch(Rectangle(
        (7.4, 3.7), 4.8, 1.8,
        facecolor=PALETTE["warm"] + "18",
        edgecolor=PALETTE["warm"], lw=1.0,
    ))
    ax.text(9.8, 5.30, "Phase 2   —   Newton mode",
            ha="center", fontsize=9.0, color=PALETTE["warm"], fontweight="bold")
    ax.text(9.8, 4.88, r"indefinite Jacobian $J_{\mathrm{New}}$",
            ha="center", fontsize=8.2, color=PALETTE["ink"])
    ax.text(9.8, 4.55, "smoother holds the symmetric Picard block",
            ha="center", fontsize=8.2, color=PALETTE["ink"])
    ax.text(9.8, 4.22, "coarse = CHOLMOD on restricted Picard K",
            ha="center", fontsize=8.2, color=PALETTE["muted"])
    ax.text(9.8, 3.90,
            "residual reduction \"no worse\" than CHOLMOD / CHOLMOD",
            ha="center", fontsize=7.5, color=PALETTE["warm"], fontstyle="italic")

    _arrow(ax, (1.80, 4.6), (2.18, 4.6))
    _arrow(ax, (6.52, 4.6), (7.38, 4.6),
           color=PALETTE["warm"], lw=1.4)
    ax.text(6.95, 5.05, "Picard → Newton\nswitch", fontsize=7.5,
            color=PALETTE["warm"], ha="center")

    ax.plot([1.0, 1.0, 2.2, 2.2], [4.25, 2.4, 2.4, 3.70],
            color=PALETTE["muted"], lw=0.8, linestyle="--")
    ax.text(1.6, 2.1, "next Picard iteration\nuses updated viscosity",
            fontsize=7.2, color=PALETTE["muted"], ha="center")

    ax.text(6.5, 1.2,
            "rule (D7):  the V-cycle's fine-level matvec tracks the current nonlinear mode\n"
            "(Picard K or Jacobian J), but the smoother and coarse-grid operator always\n"
            "use the Picard block — this keeps the preconditioner SPD and the coarse\n"
            "solve CHOLMOD-stable.",
            fontsize=8.0, color=PALETTE["muted"], ha="center")

    _header(fig,
            "Picard  →  Newton phasing under  lin_solver = 3   —   symmetric smoother rule (D7)")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=DPI, facecolor=PALETTE["bg"])
    plt.close(fig)
    print(f"wrote {out}")


# ---------------------------------------------------------------------------
# Figure 5 — Active-mask propagation (task 9.6).
# ---------------------------------------------------------------------------

def fig_active_mask(out: Path) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.8), dpi=DPI,
                            gridspec_kw={"wspace": 0.15})
    fig.patch.set_facecolor(PALETTE["bg"])

    def draw_panel(ax, N: int, title: str, label_font: float, is_fine: bool):
        _styled_axes(ax, (-0.6, N + 0.6), (-0.6, N + 0.6))

        rng = np.random.default_rng(7)
        topo = 0.60 * N + 0.18 * N * np.sin(np.linspace(0, 2.2 * np.pi, N + 1))
        topo = topo + 0.035 * N * rng.standard_normal(N + 1)

        # Draw the active (below topo) cells.
        for ix in range(N):
            for iz in range(N):
                cx, cz = ix + 0.5, iz + 0.5
                topo_here = 0.5 * (topo[ix] + topo[ix + 1])
                active = cz < topo_here
                colour = PALETTE["cool"] + "55" if active else "#ffffff"
                ax.add_patch(Rectangle((ix, iz), 1, 1,
                                        facecolor=colour,
                                        edgecolor=PALETTE["soft"], lw=0.3))

        # Draw topo line.
        xs = np.arange(N + 1)
        ax.plot(xs, topo, color=PALETTE["warm"], lw=1.8,
                label="marker-chain surface")

        # Mark the mask.
        for ix in range(N):
            for iz in range(N):
                cx, cz = ix + 0.5, iz + 0.5
                topo_here = 0.5 * (topo[ix] + topo[ix + 1])
                if cz < topo_here:
                    ax.plot(cx, cz, "o", color=PALETTE["accent"],
                            markersize=2.2)

        if is_fine:
            # Highlight a handful of consistency pairs.
            for ix_c in (1, 3, 5):
                iz_c = 3
                ax.add_patch(Rectangle(
                    (2 * ix_c, 2 * iz_c), 2, 2,
                    facecolor="none", edgecolor=PALETTE["highlight"], lw=1.3,
                ))
        else:
            # Highlight a handful of coarse cells.
            for ix_c in (1, 3, 5):
                iz_c = 3
                ax.add_patch(Rectangle(
                    (ix_c, iz_c), 1, 1,
                    facecolor="none", edgecolor=PALETTE["highlight"], lw=1.5,
                ))

        ax.set_title(title, fontsize=9.5, color=PALETTE["ink"])
        ax.text(N * 0.5, -0.45,
                r"$\bullet$  active  =  below  z_surf   ·   " +
                r"□  inactive  =  above  z_surf",
                ha="center", fontsize=8.0, color=PALETTE["muted"])

    draw_panel(axs[0], N=12, title="fine level (Nx × Nz = 12²)",
               label_font=8.0, is_fine=True)
    draw_panel(axs[1], N=6,  title="coarse level (Nx × Nz = 6²) — aggregation of 2×2 fine cells",
               label_font=8.0, is_fine=False)

    fig.text(0.5, 0.035,
             "coarse mask = OR of the four underlying fine-cell masks  "
             r"(if any fine cell is active, the coarse cell is active).  "
             "Preserves Galerkin consistency above the marker-chain surface.",
             ha="center", fontsize=8.3, color=PALETTE["muted"])

    _header(fig,
            "Active-mask propagation above the marker-chain free surface  —  fine ↔ coarse consistency")
    fig.tight_layout(rect=(0.01, 0.08, 0.99, 0.93))
    fig.savefig(out, dpi=DPI, facecolor=PALETTE["bg"])
    plt.close(fig)
    print(f"wrote {out}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir",
                   type=Path,
                   default=Path(__file__).parent,
                   help="destination directory for PNGs (default: script dir)")
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fig_vcycle_tree       (args.output_dir / "defence_vcycle_tree.png")
    fig_fgmres_controlflow(args.output_dir / "defence_fgmres_controlflow.png")
    fig_stencil_bridge    (args.output_dir / "defence_stencil_bridge.png")
    fig_picard_newton     (args.output_dir / "defence_picard_newton.png")
    fig_active_mask       (args.output_dir / "defence_active_mask.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
