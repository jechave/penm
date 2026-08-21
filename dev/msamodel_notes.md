# msamodel: refactor needed after penm's `seed` -> `ensemble` change

Written 2026-08-19, when penm renamed `get_mutant_site()`'s `seed` argument to
`ensemble`. **Nothing in msamodel was edited.** This is the to-do list for when
you get to that package.

Note this file lives under `dev/`, which is gitignored — it is a local
reminder, and will not survive a fresh clone.

## What changed in penm

Three commits, `425029e`, `1dd60ba`, `93277b4`:

1. `get_mutant_site(..., seed = 241956)` is now `get_mutant_site(..., ensemble = 1L)`.
2. The hash key dropped its dead fourth slot: `paste(ensemble, site_mut, mutation)`.
   **Every mutant's perturbations changed.** penm's own fixtures were regenerated.
3. Malformed values now error instead of silently selecting a realization.
4. `get_mutant_site()` no longer clobbers the caller's `.Random.seed`.

The rationale: the argument was never a seed — it never reached `set.seed()`,
it was hashed with `(site_mut, mutation)` and *that* was the seed. It names
which realization of the mutational process a mutant belongs to. See
`?penm_ensemble`.

## What breaks in msamodel

### 1. `seed = NULL` default now errors — `R/spm.R:57` and `R/spm.R:204`

```r
generate_spm_core <- function(wt, ..., seed = NULL) {   # :57
  ...
      mut <- get_mutant_site(wt, j, m, ..., seed = seed)   # :75
```

`generate_spm()` at `:204` has the same default and forwards at `:207`.

Today this does **not** error. `paste()` drops NULL, so the key becomes
`"-j-m"` — verified: `paste(NULL, 1, 80, 3, sep="-")` gives `"-1-80-3"`. The
default path therefore runs under a realization nobody chose, reproducibly,
which is what makes it hard to notice.

After penm's change it raises `` `ensemble` must be a single non-missing
integer ``. That is the correct outcome, but it means **msamodel's default
invocation stops working** and needs a real default or a required argument.

### 2. Argument renamed — `R/spm.R:75`, `R/spm.R:207`

`seed = seed` must become `ensemble = <whatever msamodel calls it>`. Worth
deciding whether msamodel's own user-facing argument should also be renamed
`ensemble` for consistency, or stay `seed` and translate at the boundary. The
first is more honest; the second is less disruptive to msamodel's users.

### 3. Results will change

Independently of the rename: penm's key lost a component, so every
`(site, mutation)` now draws a different stream. **msamodel's stored fixtures
and any committed scan results are stale**, notably:

- `tests/testthat/fixtures/make-znb-fixtures.R:36,81` — `SPM_SEED <- 1024`
- `tests/testthat/test-spm-generate.R:25,40,84,95` — uses `SPM_SEED`

`test-spm-generate.R:24-26` documents the single-row reproduction trick that
this change invalidates; it will need regenerating, with the same
"is-it-only-a-reseeding" check penm used (the *set* of perturbed edges must not
move — it depends on the site and `mut_sd_min`, never on the key).

### 4. `set.seed()` at `R/spm.R:58` becomes live

```r
if (!is.null(seed)) set.seed(seed)
```

Dead twice over today: it never runs under the `seed = NULL` default, and when
a seed *was* passed it was immediately clobbered by the first
`get_mutant_site()` call. penm no longer clobbers, so it now takes effect for
non-NULL values.

**Whether msamodel's results move as a result was NOT measured** — msamodel's
suite was not run. Reading the calls the scan loop makes (`ddg_dv`, `ddg_tds`,
`ddgact_*`, `delta_structure_*`), nothing obviously draws from the RNG outside
`get_mutant_site()`, in which case the line still has no observable effect. But
that is a reading, not a check. Run the suite and find out; if nothing depends
on it, delete the line.

### 5. Stale documentation — `R/spm.R:167-169`

```r
#' @param seed Random seed, for a reproducible scan. penm seeds each mutant as
#'   `seed + site * mutation`, so scans meant to be independent need widely separated
#'   seeds, not consecutive ones.
```

Both halves are wrong. `seed + site * mutation` is the *superseded arithmetic*
key — the one that made 1710 of 2280 mutants share a stream. And the advice is
now backwards: under the hash, consecutive values are fully disjoint (measured
over 228 sites x 10 mutations: overlap 0 for 1024 vs 1025). Consecutive values
are fine; there is no need to separate them.

`R/spm.R:44` (`@param seed Optional random seed, for a reproducible scan.`)
also needs rewording — it is not optional after the change, and not a seed.

### 6. Roxygen examples and vignettes passing `seed = 1024`

- `R/spm.R:192`
- `R/model.R:40,268,330`
- `R/fitting.R:340,402`
- `R/predict.R:73,153`
- `vignettes/msamodel.Rmd.orig`, `site-analysis.Rmd.orig`,
  `mode-analysis.Rmd.orig`, `inference-methods.Rmd.orig` (plus rendered `.Rmd`)

Note the vignettes' own top-level `set.seed(1024)` calls are unrelated to
mutant generation — their comments say the fits are deterministic and the seed
is kept only for a stable committed render. Leave those alone.

### 7. Import declaration

`R/msamodel-package.R:11` imports `get_mutant_site` from penm. No change needed
to the import itself, but it confirms the coupling is direct.

## Also diverged: penmscan

Not msamodel, but worth recording. penmscan carries **private copies** of
penm's `R/seed.R` and `R/penm.R`, so it did not break — but its copies now
differ from penm's. Its `sdmrs` (`R/mutscan_sdmrs.R:69-70`) still passes
`ensemble = 1L` / `2L` as a *fourth* argument beside a scan seed. Under the new
scheme those are simply two `ensemble` values. Separate task, in that repo.
