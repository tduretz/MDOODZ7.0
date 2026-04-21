## Context

The previous change `add-gmg-stokes-solver` landed the GMG core with 26/26 tests green, Picard and Newton SolVi equivalence at `5.79e-11` / `1.89e-15` vs CHOLMOD respectively, and zero regression in the existing `lin_solver = 0, 1, 2, -1` paths. What it explicitly did *not* deliver: the performance and memory numbers that are the whole motivation for GMG, the integration fixtures for the remaining benchmarks (free-surface relaxation, anisotropic shear band, power-law shear zone, conjugate shear bands, high-contrast SolKz), the public-facing defence of the effort, and the docs update to the `skill-solvers` spec.

This change closes those gaps in one focused work item. The principal deliverable is the defence document — a printable narrative that ties motivation → theory → technology → decisions → measurements → limitations into a single reading artefact that lives alongside the existing primer / codes / stokes_internals / journey / anisotropy_seismology / phenomena library. The supporting engineering is an EC2-backed perf harness (so the measurements in the defence are real and reproducible), a V-cycle visualisation pipeline (so the theory has a runnable animation anchoring it), and the sweep through the remaining `DISABLED_*` fixtures.

A secondary point captured in the proposal: if the new fixtures or the perf harness surface a latent bug in the MDOODZ engine, the fix is in scope for this change. The defence document itemises every such fix with root cause and before/after receipts.

Stakeholders: Roman (primary author and reader), Duretz (upstream maintainer and probable reviewer), any future MDOODZ contributor who will cite the defence when explaining why `lin_solver = 3` exists. A small number of the benchmark scenarios touch research directions that have seismology collaborators interested (see `anisotropy_seismology.md` §6) — the defence is also a document we could hand to them.

## Goals / Non-Goals

**Goals**

- Produce a single `defence.md` / `defence.pdf` that justifies the GMG change with measurements, not rhetoric. Page budget: 25 pages printed (~550 markdown lines + embedded figures).
- Measure GMG's memory-scaling exponent `α ≤ 1.2` and the 201² wall-time speedup `≥ 3×` vs CHOLMOD on a pinned EC2 instance, with statistical discipline sufficient to defend the numbers (N = 5 repetitions, median + IQR, thermal-throttling rejection).
- Provide a reusable EC2 perf-harness infrastructure — scripted, documented, checked-in — that future MDOODZ performance work can re-invoke without redoing the AWS plumbing.
- Land a V-cycle animation that a reader (or reviewer, or seismology collaborator) can watch and immediately understand how multigrid actually works on a real MDOODZ problem.
- Enable the remaining `DISABLED_*` integration fixtures with dual-solver equivalence checks, extending correctness coverage to free-surface, anisotropic, and strongly-localising regimes.
- Carry the `skill-solvers` docs-update that was deferred from the previous change into its proper home.

**Non-Goals**

- New numerical methods. No changes to the GMG core, the FGMRES driver, or the stencil-bridge (`StokesAssemblyGMG`) beyond debug instrumentation (`gmg_dump_vcycle`).
- 3D extension. Out of scope here as it was out of scope in the previous change.
- Cross-cloud portability. The harness targets AWS via the user's existing CLI setup; Azure/GCP are not supported.
- Continuous perf regression in GitHub Actions CI. The harness runs on-demand from the MacBook agent workstation, not on every PR. A nightly or weekly cron is a future follow-up.
- Competitive comparison against LaMEM / ASPECT / StagYY. The defence compares GMG vs CHOLMOD *within MDOODZ*, not MDOODZ against other codes.
- Seismology-coupling work (shear-wave splitting forward-modelling from `anisotropy_seismology.md` §6). Out of scope; separate future change.

## Decisions

### D1: Run the perf sweep locally on the MacBook M1 inside a hard 1-hour budget

**Choice**: Every perf-sweep run for *this* change is executed on the Apple M1 MacBook (8 hardware cores, 16 GiB RAM, `hostname`-pinned, `macos-arm64`) under a hard wall-clock budget of 3600 s. The EC2 scaffolding (`provision.sh`, `run_perf_sweep.sh`, `teardown.sh`, `grids.yaml`) is kept in-tree as future-use infrastructure but is not on the critical path; the active entry point is `benchmarks/ec2/run_local.sh` driven by `grids_local.yaml`.

**Why** — originally (spec-v1) this decision pinned a `c7i.8xlarge` EC2 instance. Three things invalidated that plan during implementation:
1. AWS access was not reliably available on the timeline this change needs.
2. Every other test and fixture in this change (dual-solver SolVi, TopoRelax, SolKz, anisotropy, V-cycle animation) has been exercised repeatedly on the same MacBook without architecture-related failure — the "ARM not audited" risk from v1 D1 is empirically closed.
3. A local run sidesteps the "corrupted-by-interruption" argument against Spot pricing (you can let the laptop sit for one hour) and removes the cloud-cost blast radius entirely.

The memory wall that motivated GMG (CHOLMOD factor ≈ 40 GiB at 801², ≈ 1.3 GiB at 321²) is still demonstrated, but on the MacBook it is *half-empirical, half-analytical*: empirical up to 321² (which fits in 16 GiB for both solvers), analytical beyond (the 801² factor size is a well-known CHOLMOD property and is quoted in §2 with the formula, not the RSS number). The honest ledger in defence.md §6 carries this split.

**Alternatives considered**
- *Wait for AWS access and keep v1 D1*: rejected — blocking the defence on cloud availability is exactly the kind of dependency this change was supposed to eliminate for reviewers.
- *Rent a one-shot c7i.8xlarge via a third-party shell (e.g. a colleague's account)*: rejected — worse reproducibility story than "the author's laptop", which any reviewer can inspect directly.
- *Shrink the sweep to <15 minutes so it fits in coffee-break iteration*: rejected — four grid points for the memory fit and a 201² strong-scaling sweep are the minimum to make the defence credible; that's ~40 minutes wall-clock, not 15.
- *Keep EC2 entirely and do the MacBook sweep as a preview*: rejected — doubling the sweep matrix burns the budget twice and gives a reviewer two sets of numbers that will never exactly agree. Better to land the MacBook numbers cleanly and flag EC2 c7i replication as §8 future work (`add-gmg-perf-regression-ci`).

### D2: Measurement protocol — N=3 repetitions, median + IQR, wall-time-variance throttle proxy

**Choice**: Every (grid, solver, thread-count) tuple runs N=3 times. We report the median wall time and peak RSS plus the IQR as an uncertainty band. Thermal throttling on Apple Silicon is detected indirectly: a run whose wall time is > 1.35× the running median of its tuple is flagged `thermal_throttled=true`, and the harness retries up to 2 times before accepting an outlier.

**Why** — the v1 D1 amendment halves the available wall-clock budget (from ~4 hours on EC2 to 1 hour on the laptop), forcing the sample count down. N=3 is the bare minimum that still gives a median (robust to one outlier); the IQR narrows the confidence band honestly. Direct throttle detection on Apple Silicon requires `sudo powermetrics` streaming (the `turbostat` equivalent does not exist on macOS), which is awkward inside an unattended harness; the wall-time-variance proxy is noisier but cheap and catches the dominant symptom (a sudden 2×-3× slowdown on a mid-sweep run).

**Alternatives considered**
- *N=5 with a smaller grid sweep*: rejected — the dominant cost is the top grid (321²), so N=5 at 321² alone is 10 extra minutes and it still leaves us short of budget on the strong-scaling sweep.
- *N=1 with "warm-up"*: useless for defensible numbers.
- *`sudo powermetrics` for real throttle detection*: would work but requires the harness to prompt for the user password in the middle of a run — breaks unattended operation, rejected.
- *Run the sweep overnight with N=5*: plausible; deferred to the EC2 follow-up in §8 where there's no budget pressure.

### D3: Harness output format — CSV-per-run + top-level summary + regression plots

**Choice**: every run emits a row into a single `results.csv` with columns `(timestamp, instance_id, git_sha, grid_nx, grid_nz, lin_solver, n_threads, repetition, wall_time_s, peak_rss_kb, n_its_fgmres, n_picard, thermal_throttled)`. After all runs complete, `analyse.py` aggregates: median + IQR tables, memory-scaling fit (`log RSS = log a + α log DOF`), wall-time scaling curves, strong-scaling efficiency plots. Output: `results.csv`, `summary.md` (rendered into the defence), and `figs/perf_*.png`.

**Why**: CSV is lingua franca. Having every individual run as a row (not pre-aggregated) means future questions ("what was the IQR at 81²?") can be answered without re-running. Git SHA and instance ID in every row prevents apples-to-oranges comparisons across reruns.

**Alternatives considered**
- *HDF5 per run*: richer but overkill; nothing in perf capture needs hierarchical or binary structure.
- *Send metrics to CloudWatch*: integrates with AWS but locks the data away; defeats reproducibility.
- *Append to Postgres*: centralised but adds an operational dependency.

### D4: Defence document — narrative structure and page budget

**Choice**: the defence has a fixed eight-section structure (see below), a 25-page printable budget, primer-series tone and visual style. It embeds the V-cycle animation by reference (URL or file-path to the GIF) because PDFs can't embed GIFs — and includes a still-frame sequence for print.

```
1. Executive summary (1 page)
2. The motivation: memory wall on CHOLMOD at scale (2 pages)
3. Theory: multigrid, Vanka, FGMRES, the stencil bridge (5 pages; references primer + stokes_internals for depth)
4. The implementation journey (3 pages; references journey.md for detail)
5. Design decisions (D1–D11 + this change's D1–D10) — 1 line each in a table (2 pages)
6. Measurements: perf + memory + strong-scaling numbers with plots (6 pages)
7. Test receipts: the 26+N integration tests and what they prove (3 pages)
8. Limitations and future work (2 pages)
9. Reading list (1 page)
```

**Why**: 25 pages is the upper bound of what a reviewer will read in one sitting. The fixed structure prevents scope creep. Sections 3, 4, 5 are deliberately short because we already have the primer / journey / stokes_internals / codes PDFs — the defence points at them rather than duplicates them. The meat is sections 6 and 7 — the measurements and the tests — because that's what's *new* in this change.

**Alternatives considered**
- *One monolithic 60-page document absorbing primer content*: rejected; duplicates existing material and makes the defence harder to update.
- *A slide deck instead of a markdown document*: rejected; Roman prints and reads long-form, not slides.
- *Jupyter notebook format*: rejected; great for exploration, poor for archival-quality documents.

### D5: V-cycle animation — instrumentation + rendering pipeline

**Choice**: add a debug-only parameter `gmg_dump_vcycle = 1` to `.txt` input. When active, on the first V-cycle of the first FGMRES iteration of the first Picard iteration of the first time step, dump per-level snapshots of `u_v`, `u_p`, `r_v`, `r_p` to HDF5 immediately before and after each operator (pre-smooth, restrict, coarse solve, prolongate, post-smooth). A Python script `benchmarks/vcycle_vis/make_gif.py` reads the HDF5, produces one frame per operator-step with two panels (residual magnitude field + its radial power spectrum), and chains them into an animated GIF at 2 fps.

**Why**: animating a V-cycle on a real 41×41 SolCx problem is more pedagogically valuable than any hand-drawn diagram. The two-panel format (spatial + spectral) directly visualises the multigrid mechanism — high-frequency error dying in the smoother, low-frequency error dying in the coarse correction. The instrumentation is debug-only and adds no overhead to production runs.

**Alternatives considered**
- *Dump only the residual field, not the spectrum*: simpler, but the spectrum panel is the one that makes the multigrid "aha!" click. Keep both.
- *Live visualisation during run*: rejected; adds runtime dependencies (matplotlib in the C build loop). HDF5-then-post-process is cleaner.
- *Programmatic animation via a Jupyter widget*: rejected for the same reason — defence readers will watch a GIF, not launch a notebook.

### D6: Process diagrams — 3–5 clean figures, matplotlib + no 3D, consistent palette

**Choice**: the defence carries 3–5 process diagrams, all matplotlib-rendered at 140 DPI, 2D only, using the existing palette from the primer series (muted blues / reds / greens). Each diagram is a single concept: the V-cycle recursion tree, the FGMRES-with-GMG-preconditioner control flow, the decoupled Stokes + stencil-bridge architecture, the Picard→Newton phasing, the active-mask propagation above the marker-chain free surface. No diagram is a dashboard.

**Why**: Roman explicitly asked for "clear, not overloaded but informative and intuitive" diagrams. Matching the primer's visual language makes the defence feel like the last volume in the same library. One-concept-per-diagram keeps each figure printable at a readable size.

**Alternatives considered**
- *Mermaid.js diagrams*: good for simple flowcharts but doesn't handle the V-cycle recursion cleanly. Reject.
- *Hand-drawn / Excalidraw*: warmer but inconsistent with the existing library. Reject.
- *3D isometric illustrations*: rejected — readable on screen, illegible on print.

### D7: Repository layout — `benchmarks/ec2/` with a dual entry-point (local + EC2)

**Choice** (amended):

```
benchmarks/ec2/
  README.md              — how to run both modes, what you need on the workstation
  _common.sh             — shared helpers (logging, yaml probe, git sha capture)
  _parse_grids.py        — grids(_local).yaml -> (sweep,nx,nz,lin,thr) CSV expansion
  run_local.sh           — PRIMARY entry point for this change (amended D1).
                           Drives SolViPerf on the laptop inside a 3600-s budget.
  grids_local.yaml       — MacBook-adapted grid/solver/thread spec (4 grids, N=3, 8-thread cap)
  analyse.py             — aggregate results.csv, produce summary.md and perf_*.png
  provision.sh           — (KEPT, dormant) launch c7i.8xlarge via AWS CLI
  run_perf_sweep.sh      — (KEPT, dormant) orchestrate remote runs on an EC2 instance
  teardown.sh            — (KEPT, dormant) terminate instance, collect logs
  grids.yaml             — (KEPT, dormant) c7i.8xlarge/64-GiB grid spec (5 grids, N=5, 32-thread)
  mocks/                 — AWS CLI mock for dry-run regression tests
  results/               — gitignored; results.csv, run.log, figs/*.png land here

SETS/SolViPerf.c         — parameterised SolVi driver (argv[1] = txt path); PERF_JSON sentinel
SETS/SolViPerf.txt       — stub required by add_set() cmake rule; unused at runtime

openspec/changes/add-gmg-stokes-defence/
  defence.md / defence.pdf / figs/ / build_pdf.py
(on archive, the above move to Downloads/add-gmg-stokes-solver/ to join the reading library)
```

The "ec2" directory name stays despite no longer being cloud-only — the EC2 harness continues to live there and a future `add-gmg-perf-regression-ci` change will re-activate it. Renaming the directory mid-change would cost more than the mild naming confusion is worth; both scripts cite their mode in the header.

**Why** — v1 D7 assumed a pure EC2 entry point; the amended layout keeps the EC2 scaffold entirely intact (for the follow-up change to use verbatim) while landing a separate `run_local.sh` that shares the same CSV schema, yaml grammar, and `analyse.py` pipeline. The CSV column `instance_type` distinguishes `macbook-m1-local` rows from future `c7i.8xlarge` rows, so both sets can coexist in one aggregated history.

**Alternatives considered**
- *Rename directory to `benchmarks/perf/`*: rejected — touches every README and commit reference in `benchmarks/ec2/`; pure cosmetic churn.
- *Keep two separate directories (`benchmarks/ec2/` + `benchmarks/local/`)*: rejected — forces either code duplication or a third `benchmarks/_common/` layer. Current layout keeps one `_common.sh` and one `_parse_grids.py` shared by both modes.
- *Delete the EC2 scaffold as "unused"*: rejected — it was designed under the full v1 D1/D2 rigour, is covered by the `mocks/aws` dry-run tests, and the follow-up change needs it as-is.

### D8: Engine-fix discipline — how surfaced bugs get handled

**Choice**: when a new fixture or perf run exposes a bug in MDOODZ core, the fix lands in this change *if* the diagnosis + implementation + before/after comparison fits within ~200 LOC of core change and half a working day of debugging. Larger fixes (architectural changes, new subsystems, multi-file refactors) are deferred to targeted follow-up changes with explicit links from STATUS.md. Every in-scope fix is individually listed in the defence with:

- the originating test (name + file path),
- a one-paragraph root cause,
- before/after output receipts (same simulation, same seed, two outputs compared),
- the commit hash of the fix.

**Why**: the engineering-value of the defence depends on it being honest about what was found. Hiding in-scope engine fixes would be cheating. Caging them to a single paragraph each (with a pointer to git for detail) keeps the document readable. The 200-LOC / half-day threshold is subjective but gives a clear gate — past that, the fix is its own change.

**Alternatives considered**
- *Fix everything in scope regardless of size*: risks the change dragging on for months and the defence never shipping.
- *Fix nothing and defer everything to follow-ups*: makes the defence weaker by separating "GMG exposed bug X" from "bug X fixed".
- *Hide the fixes in a generic "improvements" bullet*: dishonest, rejected.

### D9: `gmg_dump_vcycle` parameter — default off, debug only, never on for production

**Choice**: new optional `.txt` parameter `gmg_dump_vcycle` with values `0` (default, off) and `1` (on). When `1`, the instrumentation in `MultigridStokes.c` writes per-level HDF5 snapshots of residual and solution fields at every V-cycle operator in the *first* V-cycle of the *first* FGMRES iteration of the *first* Picard iteration of the *first* time step, then never again (to avoid flooding disk). A `LOG_INFO` line identifies the output directory. Snapshots are written to `${writer_subfolder}/vcycle_dump/level_{N}_{step_name}.h5`.

**Why**: the animation only needs one V-cycle's worth of data. Dumping every V-cycle would write gigabytes of HDF5 on a long run. The one-shot semantics keep the cost bounded and the interpretation simple.

**Alternatives considered**
- *Dump every V-cycle*: too much data.
- *Dump only final residual*: loses the animation.
- *Control via compile-time `#define`*: rejected; rebuilds would be painful for an end-user debug tool.

### D10: `skill-solvers` MODIFIED delta — section scope and wording

**Choice**: the MODIFIED delta adds a new section to the existing `skill-solvers/spec.md` that (a) documents the `lin_solver = 3` value and its behaviour, (b) lists the `gmg_*` tuning parameters with defaults and valid ranges, (c) explains the dual-solver comparison pattern for future contributors, (d) points at the convergence-factor diagnostic (`ρ < 0.15` grid-independent) as a correctness canary. This is additive; no existing requirement is removed or renamed.

**Why**: the skill is where MDOODZ documents what solvers exist and how to use them. `lin_solver = 3` has shipped; the documentation has to catch up. The convergence-factor diagnostic is worth mentioning because it's the easiest sanity check a user can run if they suspect the GMG path is misbehaving.

**Alternatives considered**
- *Don't modify `skill-solvers`; rely on the defence document for user-facing docs*: rejected; the skill is the contract. Users and future contributors look there first.
- *Write a brand-new `skill-gmg-stokes` spec*: rejected; fragments the solver documentation.

## Risks / Trade-offs

- **~~AWS cost blow-up from repeated sweep runs~~** (v1 risk; obsolete under amended D1). The dormant EC2 harness retains the `$50` budget cap via `teardown.sh` for the follow-up change that re-activates it.
- **Local sweep blows the 1-hour wall-clock budget** → Mitigated by (a) the hard `wallclock_budget_s` cap in `grids_local.yaml` enforced by `run_local.sh::project_fits_budget`, (b) a pre-flight smoke test (51² GMG) that anchors the projection, (c) dropping 801² from the empirical sweep (too expensive on M1) and arguing the CHOLMOD memory wall analytically in defence §2.
- **Thermal throttling on Apple Silicon** → Mitigated by the wall-time-variance proxy in amended D2. If more than 20% of runs for a given configuration exceed the 1.35× median threshold, the harness writes a `WARNING` into `summary.md` recommending rerunning at a cooler laptop temperature (cold start, no external display, fan-mode if any).
- **Memory-scaling measurement finds `α > 1.2`** → This is acceptable. The defence reports it honestly and scopes the gap as future work. The 1.2 target came from the previous change's spec; if reality is tighter, the spec was over-ambitious rather than the implementation being broken.
- **V-cycle animation on a real problem shows high-frequency *growth* in the smoother (e.g. from the sign-convention gotcha's ghost, or from Vanka instability)** → Treated as a test finding, not a defence-document failure. If the animation reveals a bug, the fix lands per D8 and the animation is re-run.
- **Engine fixes from §D8 leak behaviour changes into production** → Mitigated by the individual before/after receipts and by the existing dual-solver fixtures catching unintended changes on benchmarks.
- **Defence page budget blown** → Mitigated by the fixed eight-section structure and the explicit page-per-section budget in D4. If a section wants to grow past its budget, content is pushed to a dedicated follow-up piece (like `journey.md` was pushed out of the primer originally).
- **`skill-solvers` delta is out of date before it's merged** (because GMG semantics evolve during this change) → Mitigated by writing the delta *last* in the task order, after the perf numbers and fixture decisions are locked.

## Migration Plan

1. Land the perf-harness scaffolding first: `benchmarks/ec2/` directory with both the v1 EC2 scripts and the amended-D1 `run_local.sh` + `grids_local.yaml` + `SETS/SolViPerf.c`. `provision.sh` + `teardown.sh` dry-run-tested against `mocks/aws`; `run_local.sh` dry-run + single-config live-tested. No MDOODZ core changes yet.
2. Land the `gmg_dump_vcycle` debug parameter plus the HDF5 dumping logic. No animation yet, just the data.
3. Add the `DISABLED_*` fixture `.txt` twins and enable the easiest ones (free-surface relaxation, SolKz). Any surfaced bug handled per D8.
4. Run the local MacBook perf sweep inside the 1-hour budget. Capture `results.csv`. Iterate on `analyse.py` until the perf plots and tables are defence-ready. *The EC2 sweep is deferred to the `add-gmg-perf-regression-ci` follow-up (§8 of defence); the dormant harness stays in-tree.*
5. Run the V-cycle dump on a 41×41 SolCx problem. Build the animation via `make_gif.py`.
6. Enable the remaining fixtures (anisotropic shear band, power-law shear zone, conjugate shear bands). Handle any surfaced bugs per D8.
7. Write the defence document, slotting in the perf plots, animation still-frames, and the test receipts. Flag MacBook-vs-EC2 provenance explicitly in §6.
8. Land the `skill-solvers` MODIFIED delta.
9. Final validation: full CI green, defence PDF builds, `openspec validate add-gmg-stokes-defence --strict` clean.

**Rollback**: purely additive; no production path changes behaviour. Disabling the change amounts to (a) removing the `gmg_dump_vcycle` parameter parse (default-off means ignoring it is safe), (b) not using `benchmarks/ec2/`, (c) deleting the defence document. Any in-scope engine fix from D8 is a separate commit and can be reverted independently if it proves problematic.

## Open Questions

- **~~Instance type choice should be re-validated once Graviton (c7g/c8g) performance has been audited on MDOODZ.~~** Resolved under amended D1: MacBook M1 (ARM) ran every fixture in this change cleanly, so ARM ports are no longer a known risk; when EC2 reruns land in the follow-up change, `c7g` is the preferred starting point (cheaper, confirmed-working architecture family).
- **~~Do we want the perf harness to support a local (non-EC2) mode for quick iteration?~~** Resolved under amended D1 + D7: local mode is now the primary entry point (`run_local.sh`) and the cloud harness is the dormant one.
- **Animation style: 2 fps playback, or step-on-click?** 2 fps is the default; if the GIF is too fast to see the spectrum panel change, drop to 1 fps. Test during D5 and pick.
- **Should the defence call out non-GMG solver improvements too (e.g. the `InputOutput.c` remap fix from the previous change's retrospective §5.2)?** Yes, one paragraph — those fixes are part of the story of what GMG's discipline surfaced. Keeps the narrative honest.
- **Retention policy on `benchmarks/ec2/results/`**: don't commit, but don't delete either. Local folder, gitignored, user-managed. A future "perf-history" change could add a simple cloud-storage sync (S3 bucket with the raw CSVs) but out of scope here.
