## Context

MDOODZ7.0 has 10 existing OpenSpec workflow skills under `.github/skills/` (openspec-new-change, openspec-apply-change, etc.) but zero domain-specific documentation skills. The codebase — 36 C modules, 100+ scenarios, and a Julia visualisation toolkit — is undocumented beyond inline code comments. The proposal defines 10 new skills covering model setup, rheology, anisotropy, solvers, visualisation, build/run, physics theory, scenario gallery, parameter validation, and a code-to-physics glossary.

Each skill is a standalone `SKILL.md` file under `.github/skills/<skill-name>/`. Copilot loads skill files based on YAML frontmatter matching (the `description` field is matched against user queries). Skills are pure documentation — no runtime code, no build changes, no dependencies.

## Goals / Non-Goals

**Goals:**
- Create 10 skill directories under `.github/skills/`, each with a `SKILL.md` containing structured domain knowledge
- Skills are self-contained: each can answer questions in its domain without requiring users to read other skills
- Skills reference specific source files, function names, struct fields, and `.txt` parameter keys so Copilot can provide grounded, code-linked answers
- The parameter-validation skill provides an actionable checklist workflow (not just reference)
- The code-glossary skill is structured as a lookup table for fast scanning

**Non-Goals:**
- No automated linting or runtime validation tools — the parameter-validation skill provides guidance, not a script
- No modifications to existing C source, headers, build system, or Julia scripts
- No API reference generation from code — skills are hand-authored with domain expertise
- No user-facing CLI commands or VS Code extensions
- No inter-skill dependencies or cross-referencing requirements (skills are independent)

## Decisions

### 1. SKILL.md file format: YAML frontmatter + structured prose

Each `SKILL.md` follows the pattern established by existing openspec skills:

```yaml
---
name: <skill-name>
description: <1-2 sentence description matching user queries>
---
```

Followed by structured markdown sections with code references, examples, and workflow guidance.

**Rationale**: This matches the existing `.github/skills/` convention (see `openspec-explore/SKILL.md`). Copilot's skill matching uses the `description` field, so descriptions must be optimised for query matching (e.g., "how to set up a simulation" should match `skill-model-setup`).

**Alternative considered**: Using `.prompt.md` files instead of skills. Rejected because skills are auto-matched by Copilot based on user queries, while prompts require explicit invocation — skills are the right mechanism for queryable documentation.

### 2. One skill per domain, flat structure

Each of the 10 skills gets its own directory: `.github/skills/<name>/SKILL.md`. No nesting, no sub-skills, no shared includes.

**Rationale**: Keeps each skill self-contained and independently maintainable. A user asking about rheology gets everything they need from one file. Copilot loads one skill at a time, so splitting a domain across multiple skills would fragment answers.

**Alternative considered**: Fewer, larger skills (e.g., combining rheology + anisotropy + physics-theory into one). Rejected because the description matching would be too broad, and the skill files would exceed practical context window sizes.

### 3. Skill description strings optimised for natural language matching

Each skill's `description` field will be written as a natural sentence covering the key query terms users would use. Examples:

- `skill-model-setup`: "Explain how to create and configure MDOODZ simulations — .c callback files, .txt parameter files, SetPhase, SetTemperature, boundary conditions, and scenario setup."
- `skill-rheology`: "Document MDOODZ rheological framework — flow laws, viscosity computation, dislocation creep, diffusion creep, Peierls, elasticity, plasticity, and material properties."
- `skill-code-glossary`: "Code-to-physics glossary for MDOODZ — variable names, struct fields, function names mapped to physical meaning, equations, and SI units."
- `skill-parameter-validation`: "Validate MDOODZ .txt configuration files — check parameter ranges, material properties, scaling values, switch compatibility, and common pitfalls."

**Rationale**: The description is Copilot's primary matching signal. Including both code terms (`SetPhase`, `.txt`) and natural language terms ("how to create", "validate") maximises the chance of correct skill activation.

### 4. Code references use relative paths and function names, not line numbers

Skills reference code via file paths (`MDLIB/AnisotropyRoutines.c`) and function/struct names (`ViscosityConciseAniso`, `mat_prop`), never via line numbers.

**Rationale**: Line numbers change as the code evolves. File paths and symbol names are more stable anchors. Copilot can search for referenced symbols if it needs more context.

### 5. Parameter-validation skill uses checklist format with ranges

The validation skill is structured as a numbered checklist with acceptable ranges and specific warnings, not as a reference document. Each check item includes: what to check → acceptable range → what to do if out of range.

**Rationale**: A checklist is actionable — users can invoke it on a `.txt` file and get step-by-step validation. A pure reference would require users to know what to look for.

### 6. Code-glossary skill uses table format for rapid lookup

The glossary skill organises mappings as categorised tables:

```
| Code variable | Physical quantity | Symbol | Units | Notes |
|---|---|---|---|---|
| sxxd | deviatoric stress xx | τ'_xx | Pa | "d" = deviatoric |
```

**Rationale**: Table format enables fast scanning. Users looking up a variable name can find it instantly in a table, whereas prose descriptions would require reading paragraphs.

**Alternative considered**: Alphabetical listing. Rejected in favour of categorical grouping (tensors, material properties, grid fields, etc.) because related variables should appear together.

### 7. Scenario-gallery skill groups by geodynamic process

Scenarios are grouped by application (rifting, subduction, collision, shear, anisotropy tests, thermal/chemical, benchmarks, magmatic, miscellaneous), not alphabetically.

**Rationale**: Users typically know what physical process they want to model, not the exact scenario name. Grouping by process lets them find relevant starting points quickly.

## Risks / Trade-offs

**[Skills become stale as code evolves]** → Skills reference function/struct names that may change. Mitigation: Keep references at the module/function level rather than implementation details. Add a note in each skill about which source files to check for the latest API.

**[Skill description matching may be imprecise]** → Users might get the wrong skill for their query. Mitigation: Make descriptions distinct with non-overlapping keyword sets. Test with representative queries during review.

**[Skills may be too long for context window]** → Some skills (scenario-gallery, code-glossary) could be large. Mitigation: Use concise table format and brief descriptions. If a skill exceeds ~300 lines, consider splitting (but defer this to feedback).

**[Physical parameter ranges may be too narrow]** → The validation skill encodes "typical" ranges that might reject valid exotic configurations (e.g., very low viscosity for asthenosphere simulations). Mitigation: Frame checks as warnings ("this is unusual") not errors, and note that advanced users may intentionally use extreme values.

**[10 skills is a large batch]** → Creating all 10 at once makes review harder. Mitigation: Implementation tasks will be ordered so foundational skills (build-and-run, model-setup) come first, enabling incremental testing.
