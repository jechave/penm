# SC-LFENM exploration

Energy bookkeeping when the network is refit to a mutant, why naive evolutionary
trajectories fail, how to fix them, and what frustration actually changes.
2026-08-20/21.

---

## Read these three

| | what it is |
|---|---|
| **`MODEL.md`** | **The report.** 14 sections, 9 figures. Start here. |
| `REVIEW-LOG.md` | The 18 errors found while writing it, who caught each, and which checks can actually fail. |
| `NEXT.md` | Where to pick this up. Read first if resuming. |

Each also as `.html` (MathJax) and `.pdf` (printable).

## Findings in one paragraph

Refitting the network to the mutant structure erases the strain, so the energy
must be carried in a scalar offset — otherwise a mutation costing 0.6 reports 0.
That works for scanning. It does *not* work for trajectories: re-referencing at
every step makes the current state the zero of energy, so every move is uphill
(proved — the relaxation operator is an orthogonal projector) and the walk
stalls. Julian's catalogue fixes it: reference everything to the founder, and the
energy gets a fixed floor and a finite state space, so trajectories reach a
stationary state. Frustration is created by every mutation and shifts the Hessian
by 0.6–2.7 %, but rotates the soft modes by less than 0.2 %. Separately, penm's
own `sclfenm` leaves the mutant 33× further from stationarity than the LFENM one.

---

## Folders

| folder | what goes in it | run them? |
|---|---|---|
| **`R/`** | The library: files that *define functions* and are sourced by everything else. Never run directly. | no |
| **`figures/`** | One script per figure, beside the `.png` it produces. `figN_x.R` → `figN_x.png`. | yes |
| **`checks/`** | Scripts that **assert** and fail loudly. These are the evidence the report's claims rest on. | yes |
| **`analyses/`** | Scripts that **compute and print** numbers the report quotes. They do not assert; they report. | yes |
| **`attic/`** | Superseded diagnostics — steps on the way to a result now carried by something else. Kept for provenance. Nothing in the report cites them. | not needed |
| `data/` | Regenerable `.rds` artifacts. **Gitignored**; rebuilt by the scripts named in `.gitignore`. | — |

The distinction that matters: **`checks/` can fail, `analyses/` cannot.** A
script in `checks/` that always passes is a defect — see `REVIEW-LOG.md` §3,
which lists the ones that could not fail and what was done about them.

## Running anything

Every script uses the **`here`** package, anchored by the `.here` file in this
folder, so paths work from any working directory:

```r
Rscript checks/validate.R            # from the exploration root
cd checks && Rscript validate.R      # or from inside — same result
```

Scripts that load `penm` (via `devtools::load_all`) take a minute or so. The
trajectory and catalogue scripts take several. If `data/catalog.rds` is missing,
`analyses/run_catalog.R` rebuilds it in about 3 seconds.

## The checks

```
Rscript checks/validate.R          8 assertions on the offset scheme
Rscript checks/derivation.R       11 assertions walking §2's derivation
Rscript checks/test_identities.R   the 4 exact identities the catalogue rests on
Rscript checks/test_can_fail.R     breaks the code on purpose; confirms checks go red
Rscript checks/proof_projector.R   the projector proof behind §6
Rscript checks/test_stationary.R   convergence, detailed balance, reversibility
Rscript checks/test_two_impls.R    that the two implementations agree
```

All seven pass as of 2026-08-21.

## One thing to know about the code

`R/` contains **two implementations of the same physics**. `enm_core.R` and
`sclfenm.R` are a standalone reimplementation used by §§1–7; `penm_catalog.R` and
`frustration.R` use penm's own functions and are used by §§8–11. They agree to
machine precision (`checks/test_two_impls.R`).

That is deliberate in part and accidental in part. Reimplementing was *right* for
`sclfenm` and frustration — those are precisely what penm flags as unchecked
(`stopifnot(!frustrated)`, skipped tests), so an independent implementation is
what the package can be measured against, and it is what made the §11 finding
possible. It was *waste* for the settled parts. See `REVIEW-LOG.md` §4.

## Traps that cost time (all documented in `NEXT.md`)

1. `get_eij()` is not refreshed by `get_mutant_site_lfenm()` — recompute from
   `xyz` before building any Hessian.
2. Spectra must be taken at a stationary point. Off one, you get three spurious
   near-zero modes and `TS` wrong by 8.6.
3. $\beta$ is Boltzmann's $1/kT$; selection strength is $\nu$.
4. The report uses `mut_sd_min = 1`, not penm's default of 2.
5. `penm` blocks `frustrated = TRUE`; the frustrated Hessian here is built
   directly in `R/frustration.R`, bypassing the package.

## Status

The package is **not modified** by any of this. §11 identifies a real defect in
penm's `get_mutant_site_sclfenm()` (it never re-relaxes against the rebuilt
graph) with a one-line fix, deliberately not applied: the sclfenm fixtures are
frozen and its tests skip, and `penm/CLAUDE.md` says to stop and ask first.
