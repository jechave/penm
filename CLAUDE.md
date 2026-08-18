# CLAUDE.md

Guidance for Claude Code when working in this repository.

This file holds only what the repo cannot tell you: decisions, prohibitions, and
measurements that are expensive to rediscover. Anything derivable — signatures, export
lists, test counts, `check()` results — is deliberately absent. Read the code or run
the command.

## What this package is

`penm` (Perturbing Elastic Network Models) builds Elastic Network Models (ENMs) of
proteins, perturbs them, and measures the difference between wild-type and mutant in
energy, structure and motion.

It was extracted in 2026-08 from the package now called `penmscan`, which kept the
scanning layer — sweeping a perturbation across all sites (the `mrs_*` functions).
**Those live in `penmscan`, not here.** If a task calls for a site-by-site scan, it
belongs in that package.

The name means "perturb ENM" and is meant to cover perturbation types beyond mutation
(external forces, single-contact perturbations) as they are added.

Start from `set_enm()` (`R/enm.R`) and `get_mutant_site()` (`R/penm.R`); the `delta_*`
families measure wt-vs-mutant differences, by site (`i`) or by mode (`n`).

## sclfenm — something smells, unexplored

**Something is off about `sclfenm`, and what exactly is not yet known.** Julian has
not investigated it; as of 2026-08-13 it is an open question he intends to explore
and think about, not a diagnosed defect with a known fix. Treat every statement below
as an observation, not a conclusion.

What is actually observed:

- Its tests skip, with the message `"Skip sclfenm test until sclefnm is fixed"`.
- The refresh scripts guard its fixtures behind `skip <- TRUE`, so `mut_qf.rda`
  predates the 2026-08 seed-key change and `mut_sc_qf.rda` was never created. Both are
  unused while the tests skip.
- Two `#TODO` markers in `R/penm.R`: one on the `lij` update ("mut parameters are
  w.r.t. w0, not wt"), one on frustrated handling in `mutate_graph()`.
- sclfenm changes the *number of graph edges* (e.g. 956 → 962 for 2acy site 80), which
  follows from recalculating the contact map from mutant coordinates. This one looks
  like the model working as designed rather than part of the smell — but that reading
  has not been checked against the science either.

Whether these are one problem, several, or mostly harmless is unresolved.

**So: do not regenerate the fixtures, un-skip the tests, or "fix" the TODOs as a side
effect of other work.** Not because the model is known wrong, but because acting
would bake in an answer to a question that is still open. If a task touches sclfenm,
stop and ask.

## Other standing decisions

- **`frustrated = TRUE` is disabled**, not merely untested — `set_enm()` has a
  `stopifnot(!frustrated)`. Don't enable it as a side effect of other work.
- **"No caller" is not a defect.** Many exports have no caller outside penm's own
  tests. penm exists to offer a menu of perturbation measures; one going unused in
  current work says nothing about whether it belongs in the API. Do not read the call
  graph as a usage survey or propose un-exporting on that basis.

## House conventions

- **Naming:** `snake_case`.
- **Roxygen on every function, internal ones included** — it helps development and
  maintenance. Making something internal means dropping `@export` and keeping
  `@noRd`; **never delete a roxygen block to make a function internal**. "Has no man
  page" is the intent for an internal function, not a defect.
- **Three distinct roxygen mechanisms, easily confused:**
  - `@rdname` — puts several functions on one **shared page**. Preserve this grouping;
    do not flatten to one page per function. Each group has a `@name` stub holding the
    shared title and `@param`s, with members adding `@details`.
  - `@family` — *See also* cross-links only. Not a shared page. Inert under `@noRd`.
  - `@export` + `@noRd` together — **groups nothing.** Exports to NAMESPACE while
    suppressing the man page, which is exactly what produces the "Undocumented code
    objects" WARNING. Don't reintroduce it.
- **Imports live in one place:** `R/penm-imports.R` and `R/penm-package.R`.
- **Tibbles** (not data.frames) for tabular returns; tidyverse for data manipulation.
- NAMESPACE and everything under `man/` are roxygen-generated — edit the roxygen and
  run `document()`, never hand-edit them.

## Test fixtures

Fixtures live in `tests/testthat/fixtures/`, with their refresh scripts beside them.
All derive from `pdb_2acy_A.rda`.

Fixtures are **frozen** — regenerate intentionally, never incidentally, and never to
make a failing test pass. A fixture mismatch is a finding to investigate first.

Coverage is mostly regression-style comparison against these fixtures: it pins
behaviour but does not probe edge cases. A green `devtools::test()` says little about
an edit that changes what the right answer *is* — the diff is the safety net. See the
deletion discipline in the global CLAUDE.md.

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

Note `cran = FALSE` is the weaker gate — it runs with `_R_CHECK_FORCE_SUGGESTS_:
FALSE`, so undeclared `Suggests` pass. Use `--as-cran` before an actual release.

**Capture expensive output to a file on the first run** (`> out.txt 2>&1`), then read
that file. Re-running `check()` to see a different slice of its own output is waste.

## Relationship to penmscan and msamodel

- **`penmscan`** holds the scanning layer. It is a **separate package**: read it for
  reference, never write to it. It does **not** currently depend on penm — it carries
  private copies of penm's functions — so it constrains nothing about penm's exports.
  (It may later be slimmed to `Imports: penm`; that hasn't happened.)
- **`msamodel`** does depend on penm, and re-exports `set_enm()`, so **penm's help
  pages are what an msamodel user reads.** Weigh that when documenting.

  Changing the signature of anything msamodel calls is a downstream break — say so
  before doing it. Check what it actually calls rather than trusting a list here:

  ```bash
  grep -rn "penm::\|penm:::" ../msamodel/R ../msamodel/tests ../msamodel/vignettes
  grep -n "penm" ../msamodel/NAMESPACE
  ```

  Note it reaches penm both ways — declared `importFrom(penm, ...)` and bare `penm::`
  qualification — so the NAMESPACE alone understates the coupling.

## Working style

- Don't hand-wave uncertainty — flag anything unverified. Only say
  "confirmed/verified/fixed" when you actually checked.
- **The claim may not exceed the check.** Name the scope, and put the output on screen
  in the same message as the claim about it. Anything unverified gets said explicitly.
- **A subagent's report is a lead, not a finding.** Verify before repeating it as fact.
- **Outside the plan: say it, do not do it.** Approval of a plan is not approval of
  what the plan reminded you of.
- **A question is a question.** Answer it and stop — it is not a cue to resume work
  or to append a summary.
- **An unexplained artifact is evidence.** Report it; never tidy it away.
- Push back on bad ideas rather than silently implementing them. Own contradictions
  across turns instead of smoothing them over. Terse is good.
- **Don't write derivable facts into this file.** Counts, export lists, signatures and
  `check()` results go stale and then mislead. If it can be grepped, grep it instead.
- **Tool choice:** read files with `Read` (use `offset`/`limit` for big files). Avoid
  `sed`/`head`/`cat`/`tail`/`awk`/`echo` in Bash for reading — they trigger permission
  prompts here.
