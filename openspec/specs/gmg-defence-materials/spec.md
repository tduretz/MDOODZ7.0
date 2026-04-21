# gmg-defence-materials Specification

## Purpose
TBD - created by archiving change add-gmg-stokes-defence. Update Purpose after archive.
## Requirements
### Requirement: Defence document follows the fixed nine-section structure

The defence SHALL be a single Markdown document `defence.md` with a companion `defence.pdf` built via the same `build_pdf.py` pipeline used for the existing reading library. The document SHALL follow the nine-section structure declared in design D4:

1. Executive summary (~1 page)
2. Motivation — the memory wall on CHOLMOD (~2 pages)
3. Theory — multigrid, Vanka, FGMRES, the stencil bridge, pointing at primer + stokes_internals for depth (~5 pages)
4. Implementation journey — pointing at journey.md for detail (~3 pages)
5. Design-decisions table — D1–D11 from `add-gmg-stokes-solver` + D1–D10 from this change, one line each (~2 pages)
6. Measurements — perf + memory + strong-scaling with plots (~6 pages)
7. Test receipts — what the integration fixtures prove (~3 pages)
8. Limitations and future work (~2 pages)
9. Reading list (~1 page)

Total printed length SHALL NOT exceed 25 pages.

#### Scenario: Structure is enforceable by review

- **WHEN** a reviewer compares the shipped `defence.md` against the nine-section contract above
- **THEN** every numbered section is present in order, and the printed PDF page count at the chosen CSS is within the 25-page budget

#### Scenario: Sections 3 and 4 reference rather than duplicate

- **WHEN** a reader reaches the theory or journey section of the defence
- **THEN** the section's prose is ≤ 5 pages and ≤ 3 pages respectively, and it includes explicit markdown links pointing at `primer.md`, `stokes_internals.md`, and `journey.md` for the deeper treatment

### Requirement: Measurements section carries defensible numbers

The Measurements section (§6) SHALL report, at minimum: the memory-scaling exponent α with its 95% confidence interval for both `lin_solver = 0` and `lin_solver = 3`, the 201² wall-time ratio `t(CHOLMOD) / t(GMG)`, the strong-scaling plot for 1–32 threads at 201×201, and a table of FGMRES iteration counts per grid size at viscosity contrast ∈ {1, 1e2, 1e4, 1e6}. Every number SHALL be accompanied by its provenance: the `results.csv` rows that produced it, the commit SHA of the MDOODZ build, and the EC2 instance type.

#### Scenario: Every claim links back to raw data

- **WHEN** a reader questions a headline number in §6 (e.g. "α = 1.12")
- **THEN** the defence either contains or links to the `results.csv` rows used to derive it and the `analyse.py` invocation that produced it, so the claim is reproducible

#### Scenario: Regression bounds are tested against

- **WHEN** the measurement produces `α > 1.2` or wall-time ratio `t(CHOLMOD) / t(GMG) < 3` at 201²
- **THEN** the defence reports the shortfall honestly in §8 Limitations, not as a success — and records the deviation from the spec bound with a hypothesis on why

### Requirement: V-cycle animation is produced and embedded

A V-cycle animation SHALL be produced from a single instrumented MDOODZ run on a reproducible small problem (default: 41×41 SolCx). The animation SHALL be a two-panel GIF showing (left) the spatial residual magnitude field, (right) the residual's radial power spectrum, with one frame per V-cycle operator step (pre-smooth, restrict, coarse solve, prolongate, post-smooth) at every level. Default playback 2 fps.

The animation SHALL be embedded in `defence.md` via a markdown image reference to a local GIF file. Because PDFs cannot embed animations, the PDF build SHALL additionally include a 6-frame still-sequence summary taken from the animation (pre-smooth, restrict-1, coarse-solve, prolongate-1, post-smooth, converged), each with its operator-step label.

#### Scenario: Animation is reproducible from the instrumented run

- **WHEN** an operator sets `gmg_dump_vcycle = 1` on a specified 41×41 SolCx config and runs `python benchmarks/vcycle_vis/make_gif.py`
- **THEN** the pipeline produces `figs/vcycle_sol_cx_41.gif` and `figs/vcycle_still_*.png` from the dumped HDF5, deterministically up to the run's floating-point ordering

#### Scenario: PDF build includes the still-sequence

- **WHEN** `python build_pdf.py defence.md` is invoked
- **THEN** the resulting `defence.pdf` includes the 6-frame still-sequence rendered as a table-of-images, and the GIF itself is referenced by path for electronic readers

### Requirement: Process diagrams are clean, single-concept, print-friendly

The defence SHALL include 3–5 process diagrams. Each diagram SHALL visualise exactly one concept. Diagrams SHALL use the primer-series matplotlib palette (muted blues / reds / greens), 140 DPI render, 2D only, and a colour-blind-safe colour choice for any colour-encoded signal. No diagram SHALL attempt to be a dashboard (multiple unrelated signals in one figure).

The minimum diagram set SHALL cover:

- the multigrid V-cycle recursion tree
- the FGMRES-with-GMG-preconditioner control flow
- the decoupled Stokes + stencil-bridge architecture (MDOODZ assembler ↔ `StokesAssemblyGMG` matvec ↔ V-cycle)
- the Picard→Newton phasing under `lin_solver = 3`
- the active-mask propagation above the marker-chain free surface

#### Scenario: Each diagram is printable at a readable size

- **WHEN** `defence.pdf` is printed on A4 with the standard primer CSS
- **THEN** every process diagram fits within a single page at a scale where all labels are readable without magnification

#### Scenario: Diagrams do not duplicate primer figures

- **WHEN** the defence's diagrams are compared against the primer-library figures under `figs/`
- **THEN** no defence diagram is a copy of an existing primer figure; each is a new figure tailored to a defence concept

### Requirement: Engine fixes surfaced during benchmarking are itemised in the defence

Every bug in the MDOODZ engine that is surfaced by the new fixtures or the perf harness AND fixed in this change SHALL receive its own subsection in the defence's §7 Test Receipts. Each subsection SHALL contain: the originating test (name + file path), a one-paragraph root cause, before/after receipts (same simulation and seed), and the commit SHA of the fix. Engine bugs that exceed the D8 threshold and are deferred to follow-up changes SHALL be listed in §8 Limitations with a pointer to the follow-up change's identifier.

#### Scenario: An in-scope fix is recorded

- **WHEN** a fixture surfaces a bug in `StokesAssemblyDecoupled.c` that is diagnosed and fixed in ≤ 200 LOC and ≤ half a day
- **THEN** §7 Test Receipts gains a subsection naming the fixture, quoting the root cause, showing before/after output on the identical run, and linking to the fix commit

#### Scenario: A deferred bug is recorded without being hidden

- **WHEN** a fixture surfaces a bug larger than the D8 threshold and is deferred to a follow-up change
- **THEN** §8 Limitations names the bug, describes its observable symptom, links to the follow-up change identifier, and the offending fixture remains `DISABLED_` in the test suite

