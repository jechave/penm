# Where to pick this up

State as of 2026-08-21. Read this first if resuming.

---

## What is here

| file | what |
|---|---|
| **`MODEL.md`** | the report — 14 sections, 9 figures. Also as `MODEL.html`. |
| **`REVIEW-LOG.md`** | errors found and fixed, who caught each, which checks can fail, why it took as long as it did. Also HTML. |
| `NEXT.md` | this file |
| `figs/` | the nine figures |
| `*.R` | ~60 scripts; the report's §13 table says which produces what |

Lives in `penm/dev/explorations/sclfenm/`. `dev/` is tracked in git and
`.Rbuildignore`d, so this is version-controlled but stays out of the package
build. `.rds` files are gitignored here — they are regenerable (`catalog.rds` in
3.4 s by `run_catalog.R`).

Scripts locate the package by walking up to `DESCRIPTION`, so the folder can be
moved without breaking them.

Run `Rscript validate.R` (8 checks), `Rscript derivation.R` (11 assertions),
`Rscript test_identities.R` (4) to confirm the environment still works. Scripts
that load penm take a minute; `collect_numbers.R` and the trajectory ones take
several.

---

## The state of the argument, in one paragraph

Refitting the network to the mutant structure erases the strain, so the energy
must be saved in a scalar offset (§§2–3). That works for scanning (§7). It does
*not* work for trajectories: re-referencing at every step makes the current state
the zero of energy, so every move is uphill and the walk is absorbing (§6, with a
projector proof). Julian's catalogue fixes it — reference everything to the
founder, and moves become differences of catalogued quantities, so the energy has
a fixed floor, the state space is finite, and trajectories reach a stationary
state that samples $e^{-\nu V}$ (§8). Frustration is created by every mutation
and changes the Hessian by 0.6–2.7 %, but rotates soft modes by < 0.2 % (§9).

---

## What I would do next, in order

### 1. Give the catalogue states meaning (the obvious next step)

Right now a "state" is a Gaussian draw of $\delta l$ with no identity. The
catalogue structure accommodates amino acids directly: replace *state $m$ of site
$j$* by *residue type $a$ at site $j$* and nothing else changes. Julian's 2017
appendix lists this as strategy 2 — $l_{ij}$ depending only on the identities of
the residues at $i$ and $j$, giving 210 independent parameters fitted against the
$3N-6$ equations. That turns the model from a random-perturbation toy into
something that could be compared with real substitution data.

**Concretely:** build the catalogue with 20 states per site, one per amino acid,
with $\delta l_{ij}$ drawn from a distribution conditioned on the residue pair.
The trajectory machinery in `penm_catalog.R` needs no change.

### 2. Fix penm's `sclfenm` (§11) — one line of iteration

`get_mutant_site_sclfenm()` takes one linear step, rebuilds the graph around the
result, and never re-relaxes against the new graph. The mutant sits at
$\lvert\mathbf F\rvert = 0.71$ on its own network, 33× worse than the LFENM
mutant, and the effect is driven entirely by *lost* contacts (gained ones arrive
relaxed and cost nothing).

`relax_to_min()` in `frustration.R` does the iteration in ~5 Newton steps. **I
did not touch the package.** Before changing it, note that the sclfenm fixtures
are frozen and the tests skip — see `penm/CLAUDE.md`, which says to stop and ask
before touching sclfenm.

### 3. Decide whether the catalogue needs the cross-term correction

Background-independence degrades: 9 % $\to$ 48 % (energy), 0.7 % $\to$ 17 % (structure),
between 5 and 40 substitutions. The energy error is *entirely* the cross term
(correlation 0.9999), and the cross term is computable from the current state's
strain. Correcting for it is the cheapest way to extend the catalogue's useful
range; the alternative is rebuilding the catalogue every few tens of
substitutions.

### 4. Port §§1–7 onto penm's functions

The report's first seven sections use a reimplementation (`enm_core.R`,
`sclfenm.R`) written before I read `R/penm.R`; §§8–11 use the package. They agree
to machine precision (`test_two_impls.R`), but two code paths for one piece of
physics is a defect. Note the nuance in `REVIEW-LOG.md` §4: independence is
*correct* for `sclfenm` and frustration — those are what is under investigation —
and waste for the settled parts.

---

## Open questions this exploration did not settle

- **Per-move accuracy of the catalogue.** Aggregate energy tracks to 2 %, but on
  a 30-substitution background 1 move in 12 had the *wrong sign*. Fine for
  ensemble properties; not fine for a study of individual substitutions.
- **Site-additivity at larger $\sigma$.** Measured within 2 % at $\sigma = 0.3$
  on 2acy. Untested elsewhere.
- **Whether frustration matters for anything penm computes.** Measured effect is
  small (§9) but the `frustrated = TRUE` branch has still never been exercised in
  a trajectory.
- **Path-independence under the rebuild on a real protein.** The $6\times10^{-4}$
  figure in §4 is from a synthetic 40-node helix, the one result here not
  measured on 2acy.
- **SC-LFENM's irreversibility** (`sc_lfenm.md`'s 2022 `#todo`). The catalogue
  makes the *bookkeeping* reversible, which is not the same as showing the
  underlying model is.

---

## Things to be careful of (they cost me time)

1. **`get_eij()` is not refreshed by `get_mutant_site_lfenm()`.** LFENM never
   rebuilds $\mathbf K$, so penm has no reason to. If you build a Hessian from
   penm's stored `eij` after an LFENM mutation, you mix wild-type directions with
   mutant distances. Recompute from `xyz`.
2. **Spectra must be taken at a stationary point.** The LFENM structure is the
   linear-response estimate, not a minimum. A Hessian evaluated there is not
   singular on rotations and produces three spurious modes ~$4\times10^{-4}$,
   which put $TS$ out by 8.6. `relax_to_min()` first.
3. **$\beta$ is Boltzmann's $1/kT$; selection strength is $\nu$.** Conflating
   them silently broke a figure once.
4. **The report uses `mut_sd_min = 1`**, not penm's default of 2. Every contact
   of the mutated site is perturbed, backbone neighbours included.
5. **`penm` blocks `frustrated = TRUE`** in `set_enm()`. The frustrated Hessian
   here is built directly in `frustration.R`, bypassing the package.
