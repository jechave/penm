# CLAUDE.md

Guidance for Claude Code when working in this repository.

This file holds only what the repo cannot tell you: decisions, prohibitions, and
measurements that are expensive to rediscover. Anything derivable — signatures, export
lists, test counts, `check()` results — is deliberately absent. Read the code or run
the command.

## What this package is

`penm` ("perturbed ENM") builds Elastic Network Models of proteins and, given one,
builds the ENM of a mutant.

The mutant is itself a `prot`: its structure, energy and modes are computed exactly
as the wild type's, the two can be compared, and it can be mutated in turn — which is
what makes evolutionary trajectories possible, not just mutation scans. The `delta_*`
families are arithmetic on two `prot`s, by site (`i`) or by mode (`n`); they are a
convenience on top of the primitive, not the primitive itself.

Start from `set_enm()` (`R/enm.R`) and `get_mutant_site()` (`R/penm.R`).

## sclfenm — something smells, unexplored

**Something is off about `sclfenm`, and what exactly is not yet known.** Julian has
not investigated it; as of 2026-08-13 it is an open question he intends to explore
and think about, not a diagnosed defect with a known fix. Treat every statement below
as an observation, not a conclusion.

What is actually observed:

- Its tests skip, with the message `"Skip sclfenm test until sclefnm is fixed"`.
- The refresh scripts guard its fixtures behind `skip <- TRUE`, so `mut_qf.rda`
  predates *both* key changes (the 2026-08-13 hashing of the mutant key, and the
  2026-08-19 rename to `ensemble` that dropped a key component) and `mut_sc_qf.rda`
  was never created. Both are unused while the tests skip. Note `mut_sc_lf.rda` is a
  different file: it is not guarded, and was regenerated on 2026-08-19 with the rest.
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

- **The mutant key is `(ensemble, site_mut, mutation)`, hashed to seed the draw.**
  Keep it a pure function of that tuple — no session state — and don't add a
  second axis beside `ensemble`. `?penm_ensemble` is the canonical explanation;
  don't re-explain it elsewhere.
- **A breaking change bumps the minor version.** penm is `0.x`, so breaking
  changes go in a minor bump (`0.1.0` → `0.2.0`) and the API is not promised
  stable. **No `.9000` suffix** — that marks "a dev build after release X",
  which means nothing here because penm has no release event distinct from
  "what is in git". A dependent needs a version floor it can write.
- **`frustrated = TRUE` is disabled**, not merely untested — `set_enm()` has a
  `stopifnot(!frustrated)`. Don't enable it as a side effect of other work.

## House conventions

- **Naming:** `snake_case`.
- **Roxygen on every function, internal ones included** — it helps development and
  maintenance. Making something internal means dropping `@export` and keeping
  `@noRd`; **never delete a roxygen block to make a function internal**. "Has no man
  page" is the intent for an internal function, not a defect.
- **`@rdname` vs `@family`, easily confused:**
  - `@rdname` — puts several functions on one **shared page**. Preserve this grouping;
    do not flatten to one page per function. Each group has a `@name` stub holding the
    shared title and `@param`s, with members adding `@details`.
  - `@family` — *See also* cross-links only. Not a shared page. Inert under `@noRd`.
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
`--as-cran` for actual CRAN submission. (Separately: `check()` stalls at "checking
package dependencies" if `options("repos")` is the `"@CRAN@"` placeholder.)

Note `cran = FALSE` is the weaker gate — it runs with `_R_CHECK_FORCE_SUGGESTS_:
FALSE`, so undeclared `Suggests` pass. Use `--as-cran` before an actual release.

**Capture expensive output to a file on the first run** (`> out.txt 2>&1`), then read
that file. Re-running `check()` to see a different slice of its own output is waste.

## Downstream

**This repo is the whole world.** penm is a general-purpose library with users you
cannot see. Do not read, grep, or reason about sibling packages in `../`. A question
about penm is answered from penm.

- **Exported signatures are a contract.** Changing one is a breaking change: say so
  before doing it, and bump the minor version.
- **Never scope an API decision by counting callers.** You cannot see them, and the
  ones you could see would not be the set. Many exports have no caller outside penm's
  own tests; penm offers a menu of perturbation measures, and one going unused says
  nothing about whether it belongs in the API. Argue an export's fate from what it
  should mean, never from the call graph.
- **penm's help pages are read by people who never open penm.** Write the
  documentation for a stranger.

Accepted cost: a signature change comes with no report of what it breaks elsewhere.
If you want that, ask Julian — do not go and look.

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
