# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this package is

`penm` (Perturbing Elastic Network Models) builds Elastic Network Models (ENMs) of
proteins, perturbs them, and measures the difference between wild-type and mutant in
energy, structure and motion.

It was extracted in 2026-08 from the package now called `penmscan`, which kept the
scanning layer — sweeping a perturbation across all sites and tabulating
(`amrs`/`smrs`/`admrs`/`sdmrs`, the `mrs_*` functions). **Those live in `penmscan`,
not here.** If a task calls for a site-by-site scan, it belongs in that package.

The name means "perturb ENM" and is meant to cover perturbation types beyond mutation
(external forces, single-contact perturbations) as they are added.

## Architecture

### ENM creation — `R/enm.R`

`set_enm(pdb, node, model, d_max, frustrated)` builds a `prot` object with components
`param`, `nodes`, `graph`, `eij`, `kmat`, `nma`. It has **no defaults — all five
arguments are required.**

- `model`: `"anm"`, `"ming_wall"`, `"hnm"`, `"hnm0"`, `"pfanm"`, `"reach"`
- `node`: `"ca"` (alpha carbons) or `"sc"` (side chains)
- **`frustrated = TRUE` currently errors** — `R/enm.R:29` is `stopifnot(!frustrated)`.
  The `TRUE` branch is disabled, not merely untested.

The pipeline is `set_enm_*` setters calling `calculate_enm_*` workers. Both are
internal; the `calculate_*` five carry `@export` *and* `@noRd` together, so they are
exported without help pages. That contradiction is inherited from penmscan and is
known, not a discovery.

### Mutation perturbation — `R/penm.R`

`get_mutant_site(wt, site_mut, mutation, mut_model, mut_dl_sigma, mut_sd_min, seed)`
returns a mutant `prot`. Two models:

- `lfenm` — Linear Force ENM. Perturbs edge equilibrium lengths, keeps the contact map.
- `sclfenm` — Self-Consistent LFENM. Recalculates the ENM from mutant coordinates.
  **Something is off about it and it hasn't been looked into** — see below.

`mutation = 0` returns `wt` unmutated. Anything other than those two `mut_model`
values hits a `stop()`.

### Perturbation response — the `delta_*` families

Measured wild-type vs. mutant, each available **by site** (`i`) and **by mode** (`n`):

- `delta_structure_*` — `dr2` (displacement), `de2` (deformation energy), `df2` (force)
- `delta_motion_*` — `dmsf`, `dh`, `rwsip`, `bhat`
- `delta_energy` — `ddg_dv`, `ddg_tds`, `ddgact_dv`, `ddgact_tds`, `delta_energy_dvs`

`dr2i` and `dr2n` are the same displacement in two bases, so for a single mutant they
sum to the same total.

### Seeding — `R/seed.R`

All RNG seeding goes through `mut_seed(seed, site_mut, mutation, ensemble = 1L)`,
which hashes that tuple with `digest::digest(..., "xxhash32")`. This
replaced an earlier `seed + site_mut * mutation` product key that collided across
sites. `check_seeds_distinct()` errors on collision rather than warning.

### Data flow

1. PDB file → `bio3d::read.pdb()` → pdb object
2. pdb object → `set_enm()` → `prot`
3. `prot` → `get_mutant_site()` → mutant `prot`
4. (wt, mut) → `delta_*` functions → per-site or per-mode profiles

## sclfenm — something smells, unexplored

**Something is off about `sclfenm`, and what exactly is not yet known.** Julian has
not investigated it; as of 2026-08-13 it is an open question he intends to explore
and think about, not a diagnosed defect with a known fix. Treat every statement below
as an observation, not a conclusion.

What is actually observed:

- Its tests skip, with the message `"Skip sclfenm test until sclefnm is fixed"`
  (`test_penm.R`, `test_penm_sc.R`).
- The refresh scripts guard its fixtures behind `skip <- TRUE`, so
  `tests/testthat/fixtures/mut_qf.rda` predates the 2026-08 seed-key change and
  `mut_sc_qf.rda` was never created. Both are unused while the tests skip.
- Two markers left in the source: `R/penm.R:117` (`#TODO revise this: mut parameters
  are w.r.t. w0, not wt`) on the `lij` update, and `R/penm.R:226` (frustrated handling
  in `mutate_graph()`).
- sclfenm changes the *number of graph edges* (e.g. 956 → 962 for 2acy site 80),
  which follows from recalculating the contact map from mutant coordinates. This one
  looks like the model working as designed rather than part of the smell — but that
  reading has not been checked against the science either.

Whether these are one problem, several, or mostly harmless is unresolved.

**So: do not regenerate the fixtures, un-skip the tests, or "fix" the TODOs as a side
effect of other work.** Not because the model is known wrong, but because acting
would bake in an answer to a question that is still open. If a task touches sclfenm,
stop and ask.

## House conventions

- **Naming:** `snake_case`.
- **Roxygen on every function, internal ones included** — it helps development and
  maintenance. Making something internal means dropping `@export` and keeping
  `@noRd`; **never delete a roxygen block to make a function internal**. "Has no man
  page" is the intent for an internal function, not a defect.
- **Preserve `@rdname` grouping.** Ten shared pages organise most exports. Do not
  flatten them into one page per function.
- **Imports live in one place:** `R/penm-imports.R` and `R/penm-package.R`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data manipulation.
- NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen and
  run `document()`, never hand-edit them.

### Test coverage is narrow

The suite covers `set_enm`, `get_mutant_site` and `mut_seed` — 27 tests. The whole
`delta_*` / `ddg_*` / `dgact_*` family, which is most of the exports and includes
everything `msamodel` depends on, has **no test at all**. A green `devtools::test()`
therefore says little about an edit in that territory; the diff is the safety net.
See the deletion discipline in the global CLAUDE.md.

## Known state, deliberately accepted

Baseline `check()`: **0 errors, 3 warnings, 3 notes** (with `--as-cran`). All are
pre-existing documentation gaps inherited from penmscan, deferred on purpose:

- Undocumented arguments: `kmat_sqrt` (`delta_structure_by_site`), and
  `prot`/`ideal`/`pdb_site_active`/`beta` (`dgact_dv`, `dgact_tds`).
- Undocumented code objects: the `@export` + `@noRd` functions.
- Unstated test dependencies: `here`, `tidyverse` are called via `library()` in tests
  but not declared in `Suggests`.
- Several hundred NSE "no visible binding" notes, not suppressed with
  `globalVariables()`.

`get_wcn`, `get_dactive`, `delta_structure_dvmi`, `delta_structure_dvsi` are exported
with zero callers anywhere. Real implementations, kept exported pending a decision.

## Development commands

```bash
# INNER LOOP (constant) — this is the working loop
Rscript -e "devtools::load_all()"
Rscript -e "testthat::test_file('tests/testthat/test_penm.R')"   # one file
Rscript -e "devtools::test(filter='seed')"                        # name-filtered

# Document — ONLY after a roxygen / NAMESPACE / @importFrom change
Rscript -e "devtools::document()"

# Full suite — the gate before any commit touching R/, tests, or roxygen
Rscript -e "devtools::test()"

# AT A MILESTONE ONLY — never per commit
Rscript -e "devtools::check(cran = FALSE)"
```

**`check()` and the network.** The default `--as-cran` runs network checks this
machine cannot reach; each blocks until it times out. Measured 2026-08-13: **413s
wall-clock for 47s of CPU** — ~87% pure waiting, of which the world-clock call is a
60s timeout on its own. With `cran = FALSE` and `_R_CHECK_SYSTEM_CLOCK_=0` the same
check takes **39s**. Put `_R_CHECK_SYSTEM_CLOCK_=0` in `~/.Renviron`; save
`--as-cran` for actual CRAN submission. (Separately, `ideas.md` §6 in penmscan
records a *different* hang: `check()` stalls at "checking package dependencies" if
`options("repos")` is the `"@CRAN@"` placeholder.)

**Capture expensive output to a file on the first run** (`> out.txt 2>&1`), then read
that file. Re-running `check()` to see a different slice of its own output is waste.

## Test fixtures

Fixtures live in `tests/testthat/fixtures/`, with their refresh scripts beside them
(`refresh_test_enm_data.R`, `refresh_test_penm_data.R`, `refresh_test_penm_sc_data.R`).
All derive from `pdb_2acy_A.rda`; `test_seed.R` loads no fixtures at all.

Fixtures are **frozen** — regenerate intentionally, never incidentally, and never to
make a failing test pass. A fixture mismatch is a finding to investigate first.

## Relationship to penmscan and msamodel

- **`penmscan`** holds the scanning layer and may later be slimmed to `Imports: penm`.
  It is a **separate package**: read it for reference, never write to it.
- **`msamodel`** depends on `penm` and calls exactly these ten: `set_enm`,
  `get_mutant_site`, `delta_structure_dr2i`, `delta_structure_dr2n`, `ddg_dv`,
  `ddg_tds`, `ddgact_dv`, `ddgact_tds`, `get_site`, `get_pdb_site`. Because it
  re-exports `set_enm()`, **penm's help pages are what an msamodel user reads** —
  weigh that when documenting these.

Changing any of those ten signatures is a downstream break. Say so before doing it.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- **The claim may not exceed the check.** Name the scope, and put the output on screen
  in the same message as the claim about it. Anything unverified gets said explicitly.
- **Outside the plan: say it, do not do it.** Approval of a plan is not approval of
  what the plan reminded you of.
- **A question is a question.** Answer it and stop — it is not a cue to resume work
  or to append a summary.
- **An unexplained artifact is evidence.** Report it; never tidy it away.
- Push back on bad ideas rather than silently implementing them. Own contradictions
  across turns instead of smoothing them over. Terse is good.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files). Avoid
  `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger permission
  prompts here.
