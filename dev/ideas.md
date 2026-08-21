# Ideas / pending work — penm

Ideas, tasks and notes for things I might want to do. Nothing here is urgent;
none of it blocks current functionality.

Items 1-6 were carried over from `penmscan/dev/ideas.md` on 2026-08-14, keeping
only what applies to the code that now lives in **this** package, and re-checked
against penm's own source. Corrections found during that re-check are marked
**[corrected]**. Scanning-layer items (`smrs`, `sdmrs`, `amrs`, `mrs_all`,
`generate_mutants`, `old_name`, and the `Rd` width / stale-example notes for
`mutscan_*`) were left behind in penmscan, where those functions live.

---

## 1. Inconsistent / implicit default `seed`

`get_mutant_site()` (`R/penm.R:27`) defaults `seed = 241956`. **[corrected]** In
penmscan the complaint was that entry points *disagreed* (`1024` vs `241956` vs
none); after the split, `get_mutant_site()` is the only seeding entry point left
here, so there is no longer an inconsistency **within penm** — the open question
is just whether a scientific package should carry an implicit default at all, or
require `seed` explicitly.

Arguments for requiring it: most penmscan entry points already did, and an
implicit default silently makes two runs comparable that the user never meant to
tie together.

Not done because it is an API break, and `get_mutant_site` is one of the ten
functions `msamodel` depends on. Also, `tests/testthat/test_penm.R` and
`test_penm_sc.R` call it with no `seed` and so pin `241956`.

## 2. `set.seed()` clobbers the user's global RNG state

`R/penm.R:69` and `R/penm.R:110` both call `set.seed(mut_seed(...))`, which
permanently overwrites the caller's `.Random.seed`. A user who seeds their own
analysis and then calls `get_mutant_site()` silently loses their stream.

Fixable with `on.exit()` save/restore around those two sites, without taking on a
`withr` dependency. All seeding goes through `R/seed.R`, so the change is
localised.

## 3. Documentation gaps

These are the penm-side half of the `R CMD check` warnings. CLAUDE.md records the
baseline as **0 errors, 3 warnings, 3 notes** (`--as-cran`), all accepted on
purpose — this section is the list to work through if that ever gets tidied.

**Undocumented code objects** (WARNING) — `@export` + `@noRd` together, so
exported without a help page. **[corrected]** In penm the list is just six, not
the fifteen penmscan reported:

```
calculate_enm_nodes  calculate_enm_graph  calculate_enm_eij
calculate_enm_kmat   calculate_enm_nma    get_mutant_site
```

Decide per object: document properly, or drop `@export` and keep `@noRd` (per
CLAUDE.md, never delete the roxygen block to make something internal).

**Undocumented arguments** (WARNING) — verified present:
- `delta_energy`: `beta`, `ideal`, `pdb_site_active` (`man/delta_energy.Rd`
  documents only `wt` and `mut`, while `\usage` shows all five)
- `delta_structure_by_site`: `kmat_sqrt`
- `dgact_dv`: `prot`, `ideal`, `pdb_site_active`
- `dgact_tds`: `prot`, `ideal`, `pdb_site_active`, `beta`

**Typos** — all confirmed in the generated `.Rd`, so fix the roxygen and
re-`document()`:
- `man/delta_energy.Rd:36` — "entropic free **energ** difference"
- `man/delta_structure_by_mode.Rd:28,30,32` — "calculates **de** square"
  (three occurrences)
- `man/delta_structure_by_mode.Rd:32` — "force **vecgtor**"

**Rd line widths >100 chars** — get truncated in the PDF manual. **[corrected]**
penmscan's list named only `mutscan_*` pages, which are not here; penm's own
offenders are:

```
delta_energy.Rd (3)          delta_motion_by_mode.Rd (4)
delta_motion_by_site.Rd (1)  delta_structure_by_mode.Rd (3)
delta_structure_by_site.Rd (2)  get_umat2.Rd (1)
penm-package.Rd (5)          set_enm.Rd (2)
```

## 4. Roxygen for the ten functions `msamodel` depends on

`msamodel` calls exactly these ten, and because it re-exports `set_enm()`,
**penm's help pages are what an msamodel user reads**:

```
set_enm  get_mutant_site  delta_structure_dr2i  delta_structure_dr2n
ddg_dv   ddg_tds  ddgact_dv  ddgact_tds  get_site  get_pdb_site
```

Only `set_enm` has its own `.Rd`. The rest are exported but share *grouped* pages
via `@alias` (`delta_energy.Rd`, `delta_structure_by_site.Rd`,
`delta_structure_by_mode.Rd`, `get_prot_property.Rd`), so `?ddg_dv` resolves
while `man/ddg_dv.Rd` does not exist. That grouping is deliberate — CLAUDE.md
says preserve `@rdname` grouping, do not flatten to one page per function.

Priority order:

1. **`get_mutant_site` is exported with no documentation at all.** Confirmed: no
   `man/get_mutant_site.Rd`, and the roxygen block carries both `@export` and
   `@noRd`. A page should cover: legal `mut_model` values (`"lfenm"` /
   `"sclfenm"` — anything else hits `stop()`), that `mut_dl_sigma` is the SD of
   the per-edge equilibrium-length perturbation drawn with `rnorm` (Å), that
   `mut_sd_min` filters which edges are perturbed by sequence separation
   (`sdij >= mut_sd_min`), what `seed` does, and that `mutation = 0` returns `wt`
   unmutated.

2. **`set_enm` matters most** — it is the first call a user makes. The page
   exists but its `\item`s mostly restate the argument names (e.g. *"d_max:
   distance cutoff used to define enm contacts"*). It should say what a node is
   under each `node` value, what distinguishes the six `model` values, what units
   `d_max` is in with a workable value per node type (~10.5 Å for `ca`, ~12.5 Å
   for `sc`), and what `frustrated` changes.

   Two things that make this more pressing:
   - `set_enm()` has **no defaults on any argument** — all five are required.
   - `R/enm.R:29` is `stopifnot(!frustrated)`, so **`frustrated = TRUE` currently
     errors**. The page must say the `TRUE` branch is disabled rather than
     describe what it would do.

3. **No input validation on `set_enm()`.** A non-pdb argument fails with a
   message naming an internal generic:

   ```r
   set_enm(data.frame(x = 1), node = "ca", model = "ming_wall",
           d_max = 10.5, frustrated = FALSE)
   #> no applicable method for 'atom.select' applied to an object of class "data.frame"
   ```

   A clear error naming the expected class belongs here, not in a downstream
   wrapper. (msamodel's wrapper doing this was deleted on that side.)
   Fits the fail-loud-at-the-boundary principle.

4. **`delta_structure_by_site.Rd` / `by_mode.Rd`** already give the math (`dr2i`
   is the square of `Cf`). Add units, and state that `dr2i` and `dr2n` are the
   same displacement in two bases, so for a single mutant they sum to the same
   total — verified numerically downstream in msamodel (`all.equal` TRUE,
   relative difference 4e-14).

5. **`get_site` / `get_pdb_site` share `get_prot_property.Rd`**, which documents
   only `prot`. It does not distinguish the internal sequential index from the
   PDB numbering — a distinction that causes confusion downstream, so state it
   explicitly on that page.

Standard for this work: say what each function does and what each argument means,
with units and legal values. Verify every claim against the source; if a claim
cannot be verified, say so rather than writing something plausible.

## 5. Other `R CMD check` notes

- **`no visible binding for global variable`** — several hundred, from tidyverse
  NSE (`i`, `j`, `dr2ij`, `mij`, …). Standard fix is one `utils::globalVariables()`
  call in a single file. Currently absent from penm.
- **Unstated dependencies in tests** — `DESCRIPTION` lists only `testthat` under
  `Suggests`. **[corrected]** penmscan's note named `here`, `tictoc`, `tidyverse`;
  penm's test files actually `library()` **`here`, `tidyverse`, `bio3d`,
  `jefuns`** (no `tictoc`). `bio3d` is already in `Imports`; `jefuns` is not
  declared anywhere, and is not on CRAN.
- `CLAUDE.md` and `.claude` are already in `.Rbuildignore` — nothing to do.
  So is `dev`, so this file will not reach the built package.

## 6. `sclfenm` — open question, do not act incidentally

Something is off about `sclfenm` and **what exactly is not yet known.** As of
2026-08-13 this is an open question to explore deliberately, not a diagnosed
defect. See CLAUDE.md for the full framing.

Do **not** regenerate the fixtures, un-skip the tests, or "fix" the TODOs as a
side effect of other work — acting would bake in an answer to a question that is
still open. If a task touches sclfenm, stop and ask.

Observed markers:
- `R/penm.R:117` — `#TODO revise this: mut parameters are w.r.t. w0, not wt...`
  on the `lij` update inside the sclfenm path.
  **[corrected]** The *identical* TODO also sits at **`R/penm.R:75`**, in the
  **lfenm** path. penmscan's notes and CLAUDE.md both cite only line 117, which
  reads as sclfenm-specific. If the concern is real it touches both models —
  worth resolving what that marker actually means before assuming it is confined
  to sclfenm.
- `R/penm.R:226` — `WARNING: I'm not sure that "frustrated" case is handled well`
  on `mutate_graph()`.
- `R/penm.R:239` — a third marker, not previously listed:
  `WARNING: this is true only if edges are ordered`, on the `lij` copy in
  `mutate_graph()`.
- Tests skip with "Skip sclfenm test until sclefnm is fixed" (`test_penm.R`,
  `test_penm_sc.R`); refresh scripts guard fixtures behind `skip <- TRUE`, so
  `mut_qf.rda` predates the 2026-08 seed-key change and `mut_sc_qf.rda` was never
  created. Both unused while the tests skip. Regenerating them would freeze the
  output of a model whose correctness is unestablished, making a future correct
  implementation look like a regression.
- sclfenm changes the *number of graph edges* (e.g. 956 → 962 for 2acy site 80).
  This looks like the model working as designed — it recalculates the contact map
  from mutant coordinates — but that reading has not been checked against the
  science either.

---

## 7. ~~Unused exports awaiting a decision~~ — CLOSED 2026-08-18, not a real item

This section used to name `get_wcn`, `get_dactive`, `delta_structure_dvmi` and
`delta_structure_dvsi` as "exported with zero callers anywhere", pending a decision to
document, unexport or delete. Both halves were wrong:

- **`get_dactive` has a caller** — msamodel's `site-analysis.Rmd` vignette, via
  `penm::get_dactive()`. The "zero callers anywhere" sweep missed it because the
  vignette was not searched. (`penmscan` hits were also spurious: it does not depend
  on penm, it has private copies of everything.)
- **The other three are not special.** By the criterion "no caller outside penm's own
  tests", roughly 30 exports qualify — most of the `delta_*` family plus `get_cn`,
  `get_wcn`, `get_mlms`, `get_stress`, `get_umat2` and more. Singling out two of them
  was an artefact of an incomplete sweep, not a finding.

Julian's position (2026-08-18): the `delta_*` families are the package's reason to
exist; a measure nobody has run yet is not an unused export. So there is no decision
pending. Do not re-open this as an un-export proposal.

Note this is separate from §8 below, which is about `dvsi`'s *behaviour* under a
changed contact map — that one is still open.

---

## 8. dvsi and changed topology — revisit with sclfenm

Found 2026-08-18 while cleaning up. Relevant when sclfenm gets taken up, since
sclfenm is what changes the contact map (e.g. 956 → 962 edges for 2acy site 80).

`delta_structure_dvsi` and `delta_structure_dvsi_same_topology` compute the same
quantity two ways; verified equal to 4e-16 on 2acy site 80. The suffix is
misleading — **both** assert `stopifnot(all(gmut$edge == gwt$edge))`, so both
currently require identical topology. The real difference is speed: the
`_same_topology` one filters to non-zero edges first (~10 instead of ~960).

But their *algorithms* differ in what they could tolerate:

- `delta_structure_dvsi` uses `inner_join(... by = c("edge","i","j"))`, which
  matches by key. Verified: permute the mutant's edge rows and it returns the
  same answer. Its `stopifnot` is stricter than the algorithm needs.
- `delta_structure_dvsi_same_topology` subtracts vectors positionally, so row
  *k* of gmut must be the same edge as row *k* of gwt. Verified: on permuted
  rows it is wrong by 13903. It needs the guard absolutely.

**The open question is not about the code.** If topology changes, `inner_join`
keeps only the edges present in both, silently dropping the rest. Whether
"common contacts" is the right treatment is a modelling question — it is not
clear that an edge existing only in the mutant should contribute zero, or that
`dvsi` is even well defined across a changed contact map. That uncertainty is
probably why the `stopifnot` is there.

So: do not relax the guard on `delta_structure_dvsi` as a tidy-up. Decide the
science first, together with sclfenm (see §6).

---

## 9. `dgact_tds` accepts `ideal` but never uses it

Found 2026-08-18 while writing its roxygen (§3 gap, now closed).

`dgact_tds(prot, ideal, pdb_site_active, beta)` has `ideal` in its signature, but the
body (`R/enm_energy_activation.R`) computes only from `prot`: it diagonalises
`kmat_asite(prot, ...)` and sets `gact_ideal <- 0` with the comment "assume in the TS
the active-site is rigid". So `ideal` is inert.

Its sibling `dgact_dv` *does* use `ideal`, via `dxyz_asite(prot, ideal, ...)`.

Either the parameter is vestigial and should go, or the rigid-TS assumption is
standing in for a term that should be computed from `ideal`. The second is a science
question, not a tidy-up — do not "fix" it by deleting the argument, since that changes
the signature of an exported function (and `ddgact_tds` passes through it).

---

## 10. When sclfenm is taken up: revise the motion tests

The `delta_motion_*` measures (`dmsfi`, `dbhati`, `rwsipi`, `dhi`, `dmsfn`,
`dhn`, `rwsipn`, `nhn`) currently have **no test on inputs where they vary.**

Why: they compare two `cmat`s, and the only mutant model available to the tests
is `lfenm`, which perturbs edge equilibrium lengths and leaves `kmat` untouched.
So `cmat` is bit-identical between wt and mutant and every measure returns its
degenerate value (0, or 1 for the similarities). `test_delta.R` now asserts that
invariance in closed form, which is honest but only covers the degenerate case.

Verified 2026-08-18 that an **sclfenm mutant does exercise them properly** —
it rebuilds the contact map (960 → 962 edges for 2acy site 80), so `kmat` and
`cmat` genuinely change, and all eight functions return finite non-degenerate
values:

```
dmsfi  [-0.0153, 0.0040]     dmsfn  [-0.0211, 0.0007]
dbhati [3.6e-08, 0.0028]     dhn    [-0.0139, 0.0160]
rwsipi [0.9985, 1]           rwsipn [0.0009, 0.567]
dhi    [-0.116, 0.054]       nhn    [1, 13.31]
```

That is not being done now, deliberately: freezing sclfenm output would bake in
an answer to the open question in §6. **When sclfenm is settled, add regression
fixtures for the motion measures against an sclfenm mutant** — that is the
comparison that actually exercises them.

Related finding, already fixed: `delta_motion_nhn` computed `p*log(p)` without
guarding `p == 0`, giving `0 * -Inf = NaN` at 29 of 288 modes. It was invisible
because the fixture had frozen the NaN and `expect_equal(NaN, NaN)` passes.
Ironically the bug fired *only* in the degenerate lfenm case: identical `umat`
makes the overlap matrix exactly the identity, so the off-diagonal zeros are
exact. With sclfenm the overlaps are small-but-nonzero and no NaN arises.
