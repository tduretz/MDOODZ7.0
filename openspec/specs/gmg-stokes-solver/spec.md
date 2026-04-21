# gmg-stokes-solver Specification

## Status

**RETIRED** (2026-04-21).

The geometric multigrid Stokes solver (`lin_solver = 3`) has been **removed from the codebase**. This file is kept so links from archived OpenSpec changes still resolve. The implementation and benchmarks that defined this capability remain in git history (see below).

## Purpose (historical)

GMG preconditioned FGMRES for the 2D Stokes problem in MDOODZ, selected via `lin_solver = 3`, with a MAC hierarchy, Vanka smoothing, and coarse UMFPACK solves.

## Retirement rationale

Post–Fix F measurements (`benchmarks/ec2/results/PERFORMANCE_REPORT.md` §§11–13, where that report exists on disk) showed no speed advantage over CHOLMOD on 2D SolVi-class problems, and catastrophic stagnation on Rifting-class problems with extreme localized viscosity contrasts. Speed was the only remaining optimisation goal; CHOLMOD remains the production path.

## Historical record

Archived OpenSpec changes (full `proposal.md` / `design.md` / `tasks.md` / `specs/`):

- `openspec/changes/archive/2026-04-19-add-gmg-stokes-solver/`
- `openspec/changes/archive/2026-04-21-add-gmg-stokes-defence/`
- `openspec/changes/archive/2026-04-18-add-gmg-upleg-fix/`

Removal change: `remove-gmg-stokes-solver` (see `openspec/changes/archive/2026-04-21-remove-gmg-stokes-solver/`).

## Successor

None. Use CHOLMOD-backed Stokes (`lin_solver = 0` or `-1`, subject to existing Newton/anisotropy dispatch rules in `InputOutput.c`).

## Active requirements

None. All former requirements are **void**.
