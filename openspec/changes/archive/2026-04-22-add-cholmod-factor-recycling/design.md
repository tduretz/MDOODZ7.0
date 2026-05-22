## Context

`KillerSolver` (the `lin_solver = 2` path — the only production path for non-compressible Stokes since April 2026) is called once per outer Picard iteration. On hard-solve steps in `RiftingChenin 1000×800` we measure `nit = 3–4` and `solve_s ≈ 11–13 s/iter` with the solver fraction of wall time at 63% at 4 threads. Inside `KillerSolver`, three things happen:

1. `cholmod_analyze` — already cached across time steps via `DirectSolver.Analyze` (re-run only when symbolic structure changes)
2. `cholmod_factorize(Lcml, Lfact, &c)` — **runs every call**, dominates wall time
3. `kspgcr(Kcm, udum, du, Lfact, ...)` — a GCR Krylov solve that **already uses `Lfact` as a right-preconditioner**, wrapped inside the Powell-Hestenes Schur-complement loop

The existing `DirectSolver` struct in [MDLIB/mdoodz-private.h:171-187](MDLIB/mdoodz-private.h#L171-L187) persists `Lfact` across calls. So the symbolic factor is already reused, and the Krylov preconditioner path is already wired up. What is **not** reused is the numeric factor across Picard iterations — we discard it by overwriting on every re-factorization.

Between successive Picard iterations the matrix `Kcm` changes only through viscosity updates (and occasionally plasticity flips). The non-zero pattern is identical. Typical viscosity changes between Picard iterations are 1–10% once shear bands stabilise — making a stale factor a high-quality preconditioner that GCR can correct in a handful of iterations.

## Goals / Non-Goals

**Goals:**
- Reduce mean `solve_s` by ≥20% on RiftingChenin 1000×800 at 4 threads (measured over steps 22–25 from `Breakpoint00020.dat`).
- Remain bit-equivalent by default: with `cholmod_recycle = 0`, the code path produces identical factorization-count, residual sequence, and final state as today.
- Preserve convergence: when recycling would cost more than a fresh factor (GCR stagnating), transparently fall back.
- Keep the knob simple — two integer config keys, no new CLI flags or env vars.

**Non-Goals:**
- Not introducing a new outer solver (no Newton-Krylov, no FGMRES outer loop). Those are separate projects.
- Not touching the thermal or chemical CHOLMOD paths. Thermal already uses PCG via [`thermal-factorization-cache`](openspec/specs/thermal-factorization-cache/spec.md).
- Not modifying the symbolic-analysis caching. That already works.
- Not optimizing the compressible path (`lin_solver = -1`, `DirectStokesDecoupledComp`). Addressing it would duplicate work; if results are promising we will extend as a follow-up.
- Not parallelizing CHOLMOD's serial triangular solve. That is a SuiteSparse-level limitation.

## Decisions

### Decision 1 — Recycle across Picard iterations, not across time steps

**What**: Within a single time step's nonlinear loop, reuse `Lfact` from the previous Picard iteration as the preconditioner for the next KillerSolver call. At the start of each new time step, refresh the factor before the first solve.

**Why**: Within a time step the Jacobian/Picard operator changes gradually with viscosity updates — the stale factor is a good preconditioner. Across time steps the matrix changes more (new `dt`, elasticity shifts, phase changes from advection), so starting each step with a fresh factor is the safe default. The expected gain comes from iterations 2+ of the Picard loop.

**Alternatives considered**:
- Recycle across time steps too. Risk: factor quality degrades sharply after advection/phase updates. Gain is marginal (the first Picard iteration of each step is one `cholmod_factorize` either way; the gain is in the 2nd, 3rd, 4th iterations of the *same* step). Rejected for safety.
- Recycle only when viscosity Δ < threshold. Adds a global scan on every Picard iteration; the GCR convergence budget below gives the same protection automatically and is simpler.

### Decision 2 — Reuse signal = Picard iteration index, enforced by caller

**What**: Add an explicit `picard_iteration_index` (or equivalent boolean `is_first_picard_of_step`) argument that `SolveStokesDecoupled` and `SolveStokesDefectDecoupled` pass into KillerSolver. KillerSolver calls `cholmod_factorize` unconditionally when this is 0 and when `model.cholmod_recycle == 0`; otherwise it skips factorization up to `model.cholmod_recycle_max_reuse` times.

**Why**: Picard iteration count is the natural trigger. The caller already knows it — the outer loop in `Nonlinear*Iterations` drives Picard. Plumbing it through is a few-line change and avoids stateful "staleness counter" bookkeeping inside the solver.

**Alternatives considered**:
- Stash a counter in `DirectSolver`. Rejected: introduces hidden state that must be reset on every time step, every restart, and every failure.
- Auto-detect using matrix diffing. Rejected: cost of comparing matrix entries is non-trivial relative to the savings.

### Decision 3 — Correctness fallback via GCR convergence monitoring

**What**: The Powell-Hestenes loop in KillerSolver calls `kspgcr(...)` with a relative tolerance `model.rel_tol_KSP` and a max iteration budget `model.max_its_KSP`. When recycling, we watch `its_KSP` returned from `kspgcr`:

- If `its_KSP < max_its_KSP` — the stale factor preconditioned well; proceed.
- If `its_KSP == max_its_KSP` AND the residual norm `ru` did not meet tolerance — the stale factor is too stale. Abort the current Powell-Hestenes iteration, run `cholmod_factorize(Lcml, Lfact, &c)`, retry the same Picard step.

**Why**: GCR is already the production path and already reports iteration count and residual. Cost of a failed iteration is bounded: one GCR attempt at max budget (~1s at 1000×800) plus a fresh factorization (~4s), vs the normal ~4s factorization — so the worst case is a ~25% slowdown on a single iteration, amortized across `nit` iterations. Net expected gain is still substantial.

**Alternatives considered**:
- Residual-norm threshold check instead of iteration count. Equivalent but noisier to tune.
- Unconditional retry with fresh factor (no fallback, just let `kspgcr` fail loudly). Rejected because the code then diverges silently — users would see NaN solutions before seeing a warning.

### Decision 4 — `cholmod_recycle_max_reuse` default of 3

**What**: Maximum consecutive Picard iterations that reuse the same numeric factor before forcing a refactor. Default 3.

**Why**: Matches typical `nit` range (1–4 on RiftingChenin, 1–2 on simpler scenarios). Guarantees that at `nit = 4` we do at most one refactor (the first) plus three reuses, vs the baseline of four factors. If `nit = 5+`, we refactor at iteration 4 and then reuse up to 3 more times — bounds worst-case staleness.

**Alternatives considered**:
- No limit (only rely on GCR fallback). Risk: on a 10-iteration Picard step the final factor is ~10 iterations stale — GCR may still converge but slowly, eroding the gain.
- Limit of 1 (reuse once, then refactor). Too conservative; mean gain drops from ~50% to ~25%.

### Decision 5 — Factor lifecycle = managed by `DirectSolver`

**What**: `Lfact` is allocated during the first `cholmod_analyze`, kept alive across Powell-Hestenes iterations (already the case), Picard iterations (new), and time steps (already the case), and freed only at shutdown via the existing cleanup path. No extra allocation/deallocation in the recycle path — the `cholmod_factorize(L, L, &c)` CHOLMOD API does an in-place numeric update.

**Why**: CHOLMOD's `cholmod_factorize` with a pre-analyzed factor re-uses the factor's allocated storage; no extra memory is requested per call. The "more memory" that the user signed off on comes from the factor simply **not being freed** between Picard iterations — which the code already does today for the analysis. No new `malloc` paths.

**Alternative considered**: Keep two factors (current + previous) for a potential LBFGS-style blend. Rejected as over-engineering.

### Decision 6 — Config surface

Two integer keys in the `.txt` file, both parsed via `ReadInt2` in [MDLIB/InputOutput.c](MDLIB/InputOutput.c):

| Key | Default | Meaning |
|-----|---------|---------|
| `cholmod_recycle` | `0` | `0` = disabled (baseline behaviour), `1` = enabled |
| `cholmod_recycle_max_reuse` | `3` | Max consecutive Picard iterations reusing the same factor. `0` ⇒ disabled. `1` ⇒ reuse once then refactor. |

Both keys are independent of `cholmod_threads` (Decision 7 below spells out the interaction). When `cholmod_recycle = 1`, a single `LOG_INFO` at startup confirms the mode is active, mirroring the `cholmod-threading` logging pattern.

### Decision 7 — Multithreading interaction

Two independent threading levers exist: OpenMP parallelism in `cholmod_sdmult` / MDOODZ loops (driven by `OMP_NUM_THREADS`) and CHOLMOD's internal BLAS threads (driven by `cholmod_threads`). Recycling does not change either, but it changes how time is *spent*:

- **Time that recycling removes**: serial `cholmod_factorize` inside KillerSolver. This does not scale with `OMP_NUM_THREADS`.
- **Time that recycling adds**: extra GCR iterations. Each GCR iteration is dominated by `cholmod_sdmult(Kcm, ...)` (SpMV) and the triangular solve inside `cholmod_solve2`. `cholmod_sdmult` is OpenMP-parallel and scales. The triangular solve is essentially serial in CHOLMOD.

**Implication**: the gain from recycling **grows with thread count**, because we are trading serial work for partially-parallel work. On M1 with 4 P-cores the upside should be larger than on a 1-thread measurement. On workstations with `cholmod_threads > 1`, the raw factorization gets faster too, narrowing the absolute gap — we expect recycling to still win but by a smaller percentage. Benchmark the sweep at (t=1, t=4) and (cholmod_threads=1, cholmod_threads=4) to map the interaction.

**No new parallel regions** are added — CHOLMOD's existing `nthreads_max` and OpenMP region in `cholmod_sdmult` cover the hot paths.

## Risks / Trade-offs

- **[Risk]** Stale factor causes GCR to stagnate on hard steps (strong yielding, sudden viscosity collapse). → **Mitigation**: Decision 3 fallback forces a refactor and retry. Worst case is one wasted GCR attempt (~1 s at 1000×800) — small relative to a factor call (~4 s).
- **[Risk]** Line-search or Newton iterations interact with the iteration counter. In defect-correction mode (`SolveStokesDefectDecoupled`), the matrix may change substantively (switch from Picard to Jacobian operator). → **Mitigation**: Reset the iteration index to 0 — forcing a refresh — whenever the operator switches (e.g. Picard → Newton, line search rejects a step and retries). The caller is responsible for signalling this.
- **[Risk]** On an unusual scenario (`anisotropy = 1`, fresh setup, very small `dt`) the first Picard iteration already runs GCR close to its iteration limit. Any staleness tips it over. → **Mitigation**: `cholmod_recycle_max_reuse = 0` (or equivalently `cholmod_recycle = 0`) is the documented safe fallback. Parameter validation spec enforces the default.
- **[Trade-off]** Peak RSS grows by ~1–2 GB on 1000×800 because the numeric factor stays allocated between Picard iterations. User has signed off. Irrelevant on small grids. On very large grids (2000×1600) this could mean 5–8 GB extra and an OOM risk on 16 GB machines — the config flag lets users opt out.
- **[Trade-off]** Benchmark numbers now depend on both `cholmod_recycle` and `interp_mode` / `cholmod_threads`. Skill documentation must make the `recycle=1` results distinguishable from the baseline (suggest file naming convention `im3_t4_rc1.csv`).
- **[Risk]** Silent-success mode — a subtle correctness bug (e.g. wrong-matrix GCR preconditioning) could drift the solution without failing convergence checks. → **Mitigation**: Acceptance test compares velocity/pressure fields to the baseline (non-recycled) run over 5 time steps on a small scenario; any divergence above solver tolerance fails CI.

## Migration Plan

- Land the code change behind `cholmod_recycle = 0` default. No existing runs change behaviour.
- Run the benchmark sweep from `Breakpoint00020.dat` at threads ∈ {1, 2, 4} × `cholmod_recycle` ∈ {0, 1}. Publish before/after table in `benchmark-results/PERFORMANCE_REPORT_*.md`.
- If the gain is ≥20% and the correctness test passes, update `skill-stokes-benchmark` and `skill-solvers` to recommend enabling the flag for production rifting-class runs.
- No data migration, no restart format change — `.dat` checkpoints remain compatible either way (the factor is reconstructed from the restarted matrix on the first step post-restart).

## Open Questions

- Does the fallback path need to clamp the Picard iteration cap (`nit_max`) — e.g. after a fallback, should we count the failed attempt against the Picard budget or restart it? **Tentative decision**: count as a completed iteration (the caller has already advanced `k`), so we don't loop forever on a pathological step.
- Should we also expose a per-run "never refresh" debug flag to stress-test the fallback path? Probably yes as a developer-only knob guarded by `#ifdef`. Revisit during tasks.
- In defect-correction mode, the matrix fed to KillerSolver is the Jacobian (not the Picard operator). Recycling across Picard steps *within* defect correction is still valid, but if Newton-damping halves the update we may need to refresh more aggressively. Benchmark with Newton enabled once the basic Picard case works.
