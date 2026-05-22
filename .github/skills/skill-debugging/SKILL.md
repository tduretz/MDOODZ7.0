---
name: skill-debugging
description: Debugging a failing MDOODZ model — the three failure families (CHOLMOD non-positive-definite, Powell-Hestenes non-convergence, first-iteration NaN), the error-string lookup table, the eta-contrast reality, the two anisotropy recipes, the single-variable bisection cookbook, and the known MDLIB anisotropy bugs with their fix commits. Use when a run crashes, won't factor, won't converge, NaNs, or stalls mid-production.
---

# Debugging a Failing MDOODZ Model

Distilled from `~/.claude/reading/cooking-mdoodz-models/primer.md`. When a run breaks you get a red error and a `Killed` log — but there are **three distinct failure families** with three root causes and three fixes. This skill is the diagnostic ladder.

## First move — diff against the nearest working SET, do NOT tune

The single most important rule. Before twiddling any solver knob, answer: **"how does this SET differ from the closest working SET in `SETS/`?"** If you cannot articulate "this differs from `CollisionAnisotropy.txt` (or `RiftingCheninAniso.txt`, or `AniFstrainQuartzBlackford.txt`) only in axes X, Y, Z" — you are already lost. One `grep` of `SETS/` for your parameter combination ends most investigations in minute one. Unguided `max_its_PH` / `min_eta` / `Newton` / `penalty` sweeps produce nothing actionable.

## Error-string lookup table

Search the log for the string; that tells you the family and the first move.

| Error string in log | Family | First thing to try |
|---|---|---|
| `CHOLMOD warning: matrix not positive definite` | A — SPD violation | `Newton=1 → 0`. Still fails: `free_surface_stab=1 → 0`. Still fails: `aniso_factor → 1`. |
| `CHOLDMOD failed... probably negative!` | A — SPD violation | Same ladder (this MDOODZ error follows the SuiteSparse warning). |
| `incompressibility constrain was not satisfied to abs. tol.` | B — PH non-convergence | Raise `penalty` ×10. **NOT** `max_its_PH`. Also try `lin_solver=2`. |
| `Try modifying the PENALTY factor or check MIN/MAX viscosities / Good luck!` | B — PH non-convergence | The hint is half-right — penalty matters; min/max viscosity usually does not. |
| `Fu abs. = nan` at iteration 00 of step 1 | C — upstream NaN | Test `cstv=1`. Passes → power-law issue. Still NaN → particle init or BC. |
| `Lazy bastard! Fix your particle phase ID!` | C — `SetPhase` bug | The `.c` returns a phase ID outside `[0, Nb_phases)`. Open the `.c`, find the `return`. |
| `negative dissipation: you crazy! Wdiss = ...` | A/C precursor | Anisotropic rheology produced negative `eta` — see Bug 1. Often precedes a crash. |
| `min. Maxwell = -...` (negative) at `DefineInitialTimestep` | C — Bug 1 signature | Negative `eta_ve` from anisotropic Newton overshoot on cold cells — see Bug 1. |
| `Killed` / `Killed: 9` with no error | OOM | Grid too large for memory. Drop resolution. |
| Run completes but `Nt < requested` | Silent stop | Parent shell terminated the process — use `nohup` next time. |

## The Stokes solve, in one paragraph

MDOODZ discretises incompressible Stokes on a MAC staggered grid → a coupled saddle-point block system. It does **not** factor the whole (indefinite) matrix: it Schur-reduces to a velocity-only system `Jts = J_uu − Bᵀ(1/penalty)B`, factors **that** with CHOLMOD, and recovers pressure via outer Powell-Hestenes (PH) iterations. CHOLMOD requires `Jts` to be **SPD** (symmetric positive-definite). Every failure family is a way that contract breaks. (See `skill-solvers` for the full architecture.)

## Family A — CHOLMOD non-positive-definite (SPD violation)

**Symptom:** the `CHOLMOD warning: matrix not positive definite` / `CHOLDMOD failed` strings. A Cholesky pivot crossed zero — `MDLIB/Solvers.c:1246-1251` catches it and `exit(1)`s.

**Causes, most likely first:**
1. **Newton mode with strong nonlinearity.** The Newton tangent adds `∂η/∂ε̇` cross-terms that are *not* symmetric (`D12 ≠ D21`) — `MDLIB/StokesRoutines.c:304-312`. Picard (`Newton=0`) keeps `J_uu` SPD by construction. **Bisect by switching to `Newton=0`.**
2. **`aniso_factor` extreme** (>~15 in assembled cells). The anisotropic stress projection blows the operator viscosity up; can fail *later* (step 50+) when the fabric has rotated to align with the stiff axes.
3. **Free-surface stabilisation sign error.** Rare. The guard at `MDLIB/StokesAssemblyDecoupled.c:143-144` is the safety net. If it breaks the moment `free_surface=1`, isolate with `free_surface_stab=0`.
4. **`free_surface=1` with no air buffer** (`zmax ≤ max-z of material phases`) — see "Implicit conventions" below.

**Ladder:** `Newton=0` → `free_surface_stab=0` → `aniso_factor=1` (then ramp slowly) → log the eta range (`LOG_INFO("max(eta_s)=%e min(eta_s)=%e")`) at the first call site after rheology and look for a negative or NaN value.

## Family B — Powell-Hestenes non-convergence

**Symptom:** `rel.max.div` stagnates at 10⁻²–10⁻³ across many PH iterations, then `incompressibility constrain was not satisfied`. The PH outer iteration should drive `max(div u)` to ~10⁻⁹ in 5–10 iterations; if it *decays but stalls* (`1.0 → 0.1 → 0.01 → 0.005 → stuck`), it is stalling.

**Causes:** (1) **penalty too small** — the pressure update is too gentle; (2) **Schur preconditioner degenerate** — `Jts` is technically SPD but near-singular, which large `aniso_factor` / extreme operator-viscosity contrast produces.

**Ladder:**
- **Raise `penalty` first** (×10: `1e2 → 1e3 → 1e4`), **not `max_its_PH`** — more iterations just stall longer. (The CI-tested calibration recipe uses `penalty=1e5`.)
- Switch to `lin_solver=2` (augmented PH / "Killer solver").
- Reduce contrast at the *operator* level: `aniso_factor → 1` temporarily, confirm PH converges, ramp back while watching `rel.max.div` track.

**Time-dependent variant — passes smoke, stalls mid-production.** With `aniso_factor ≥ 8`, early steps converge in 3 PH iterations; by step 50–100 the fabric has rotated and some cells become locally hyper-viscous. **Fingerprint:** wall-time per step climbs 3–5× over ~10 steps before the crash (a flat-then-sudden-break profile is the *static* stall instead). This was MDLIB **Bug 3** — see below; if your MDLIB has the fix (`d4b4784`), it is gone at moderate strain rate.

## Family C — first-iteration NaN

**Symptom:** `Fu abs. = nan` at iteration 00 of step 1 — *before any matrix is factored*. An upstream computation produced a NaN.

**Cause:** the visco-elastic-plastic rheology iteration uses `pow(Tii, n_pwl − 1)`; if `Tii` is negative (sign-convention bug, or anisotropic Newton overshoot — Bug 1), `pow(negative, non-integer)` returns NaN. Fingerprint: negative `Tii`/`eta_ve` at iter 00 step 1.

**Ladder:**
1. **Bisect `cstv=0` vs `cstv=1`** (`cstv=1 pwlv=0 eta0=1e22`). Works → bug is in the power-law iteration. Still NaN → bug is upstream of rheology (particle init or BCs).
2. **Check the `SetPhase` callback** in the `.c` — a phase ID outside `[0, Nb_phases)` triggers `Lazy bastard!`. Reducing `Nb_phases` in the `.txt` without updating the `.c` is the classic cause.
3. **Print stress invariants** at the first call: `LOG_INFO("Tii=%e eta_ve=%e ...")` at the top of `LocalIterationViscoElastic` in `MDLIB/RheologyDensity.c`.

## The eta-contrast reality

A common wrong intuition: "the problem is the eta contrast, so tune `min_eta`/`max_eta`." **`min_eta`/`max_eta` clamp the *marker rheology evaluation*, not the *assembled operator*.** A sweep of `min_eta` across 5 orders of magnitude on `CollisionAnisotropy.txt` changed *nothing* — PH count and residuals identical. The operator viscosity that actually enters the matrix is routinely 10⁴⁵–10⁵⁸ when your caps say 10²⁵.

What *does* drive operator stiffness is **`aniso_factor`**. There is a regime change — a **~21-order-of-magnitude cliff between `aniso_factor=2` and `aniso_factor=5`** in `max(eta_s)`. PH residual quality degrades ~40× from `aniso_factor=1` to 20; per-step wall time doubles. Practical reading: `aniso_factor ≤ ~5` is the easy regime; `≥ 8` is the hard end where you need tighter PH settings or `lin_solver=2`, and where production runs can partial-complete.

## The two anisotropy recipes — pick one, don't invent a third

Every `anisotropy=1` SET in `SETS/` follows one of two working recipes. The classic mistake is forking from a non-anisotropy SET and inventing a third combination that exists nowhere.

| | **Calibration recipe** (CI-tested) | **Production recipe** |
|---|---|---|
| Examples | `AniFstrainQuartzBlackford`, `PinchSwellGSEAniso` | `CollisionAnisotropy`, `RiftingCheninAniso`, `Beghein*` |
| Free surface | `free_surface=0` | `free_surface=1, free_surface_stab=1` |
| Rheology | `cstv=1 pwlv=0 eta0=1e22` (AniFstrain) or 2-phase (PinchSwell) | `cstv=0, pwlv=10/19/40/41` (real flow laws) |
| Nonlinear | per-scenario (`Newton=0/1`, `nit_max=1/10`) | `Newton=0, nit_max=1` (pure Picard; the time loop carries the nonlinearity) |
| `penalty` | `1e5` (high → PH converges fast) | `1e2` (modest, calibrated to multi-phase contrast) |
| `Nb_phases` | 1–2 | 4–9 |
| Purpose | test the anisotropy plumbing in isolation | real geodynamics |

Both work. Neither is "right" — they serve different purposes. Fork from the one closest to your target.

## The cookbook — build a new scenario without crashes

1. **Step 0** — pick the closest working SET (`AniFstrainQuartzBlackford` for a calibration test; `CollisionAnisotropy` / `RiftingCheninAniso` for production).
2. **Step 1** — drop resolution to 101×41 (always smoke at the CI resolution first; a 101×41 smoke is 30–60 s).
3. **Step 2** — run the **unmodified** SET at `Nt=10`. If the base recipe fails here, your *build/environment* is the problem, not your edits.
4. **Step 3** — change **one** parameter toward your target. Re-smoke `Nt=10`.
5. **Step 4** — if it breaks, the error + the one thing you just changed *is* the diagnosis.
6. **Step 5** — if it works, stage the next single change. Repeat.
7. **Step 6** — only after the full modification set smokes clean at 101×41, ramp resolution (151×61 → 201×81 → …), re-smoking at each.
8. **Step 7** — if 7+ single changes all pass and it *still* breaks at higher resolution / longer `Nt`, look at the `.c` file (usually the `SetPhase` callback), not the `.txt`.

**High-aniso production partial-completion decision** (`aniso_factor ≥ 5`): a run that passes smoke can still partial-complete via the time-dependent PH stall.
- Crash at step `k ≥ 50`, HDF5 outputs intact → **accept as a data point** (real science to t ≈ 5 Ma, γ ≈ 0.25–0.5; the boundary is solver-bound, not physics-bound).
- Crash at step `k < 30` → retry with `aniso_factor` halved (drop the highest aniso point, keep the comparison ladder).
- `k ∈ [30, 50]` → case-by-case on whether your diagnostic question is answered.
Chasing the stall with more solver tuning does not work — the preconditioner is structurally failing, not under-iterating.

## Known MDLIB anisotropy bugs

Check whether your MDLIB has the fixes: `git log --oneline | grep -E "fix ani routine|fix bugs"`.

- **Bug 1 — cold-cell anisotropic Newton overshoot. FIXED (`533c48c "fix ani routine"`).** `pwlv=41` (HK03 wet olivine) + `anisotropy=1` uniformly over a column with cold-surface nodes (T < ~800 K): `ViscosityConciseAniso` (`MDLIB/AnisotropyRoutines.c`) had no bisection safeguard, Newton overshot to negative `eta_ve` → negative `eta` → `Wdiss < 0` → negative `min. Maxwell` → negative `Δt` → NaN at iter 00 step 1. Fix: log-space bisection-then-Newton + an `*eta_vep` positivity clamp.
- **Bug 3 — missing ceiling clamps in `ViscosityConciseAniso`. FIXED (`d4b4784 "fix bugs"`).** The aniso path had the *floor* clamps (Bug 1 fix) but not the *ceiling* clamps the non-aniso path has (`RheologyDensity.c:1058-1073`). `eta_ve` reached 10⁸⁰ with `max_eta=1e25` → operator condition number past CHOLMOD's double-precision tolerance → PH stall mid-production. Fix: ceiling clamps mirroring the floors, + a KillerSolver noise-floor relaxation in `Solvers.c`. **Caveat: necessary but not sufficient at high strain rate** — at `ε̇ ~ 10⁻¹⁴` the pathology re-emerges in a different regime; the safe `aniso_factor` ceiling is strain-rate-dependent (≤ 13 at `ε̇ ~ 10⁻¹⁵`, back toward ≤ 8 at `ε̇ ~ 10⁻¹⁴`).
- **Bug 4 — Phase-1-freeze + cold cells CHOLMOD non-PD wall. CANDIDATE, OPEN.** A frozen-fabric + cold-cell setup hits a hard CHOLMOD non-PD wall at a fixed step, identical across solver-knob recovery attempts (the local rheology is healthy — failure is in the *assembled operator*). Suspected freeze-interface code-symmetry gap, same shape as Bugs 1 & 3. Investigation: instrument `ViscosityConciseAniso` + the stress assembly around the freeze branch with `printf` at the failing step. *(Note: the `ani_fstrain=3` T-threshold freeze this bug concerns was redefined to a δ-relaxation in the `ani-fstrain3-delta-relaxation` change — this entry describes the historical freeze.)*

## Implicit conventions that look like bugs

These produce crashes whose error message points the wrong way:
- **`free_surface=1` needs an air buffer:** `zmax` must exceed the max-z of any material phase. Every production lithospheric aniso SET uses `zmax=+10e3` (10 km air). `zmax=0` → surface-row coefficients turn pathological → CHOLMOD non-PD. *The bug is in the error message, not the solver.*
- **`Nb_phases` / `SetPhase` consistency:** changing `Nb_phases` in the `.txt` without updating the `.c` → `Lazy bastard!` abort.
- **`min_eta`/`max_eta` are marker-level, not operator-level** — see "eta-contrast reality" above.

## Meta-rules

- **Diff against the nearest working SET first. Do not tune.** This is rule #1 for a reason.
- **Bisection-against-known-working** (single-variable changes from a working SET) isolates a root cause fast; arbitrary parameter sweeps produce nothing.
- **When static-analysis confidence on a bug is below ~80%, instrument before you patch.** Bug 3 had three plausible static-analysis hypotheses — all three were wrong; one round of `printf` in `ViscosityConciseAniso` found the real cause in <30 min.
- **Smoke (`Nt=10`) is necessary but not sufficient.** High-`aniso_factor` runs pass smoke and stall in production. Budget for partial completion at `aniso_factor ≥ 8`.
- **Watch for silent scope swaps.** An autonomous agent debugging a run can quietly swap physics parameters (`cstv`, `free_surface`) to make runs "go through" and still emit a confident verdict. Review the diff against the known-working baseline before trusting any "investigation complete".

## See also

- `skill-solvers` — the Stokes/PH/CHOLMOD architecture (this skill is its failure-mode companion).
- `skill-rheology` — the viscosity pipeline and creep mechanisms behind Family C.
- `skill-anisotropy` — the director / `aniso_factor` / `ani_fstrain` machinery behind the eta-contrast cliff.
- `skill-parameter-validation` — pre-run `.txt` checks that catch some of these before launch.
- `~/.claude/reading/cooking-mdoodz-models/primer.md` — the full narrative source, with figures (decision tree, eta-ladder, recipe comparison) and the v1 postmortem.
