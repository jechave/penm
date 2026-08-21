# Review log

Errors found and fixed while producing `MODEL.md`, and what was checked. Kept
separate so the report stays clean. 2026-08-20/21.

---

## 1. Errors in my own work, and how each was caught

| # | error | caught by | fix |
|---|---|---|---|
| 1 | Wrote $V_{\text{mut}}(\mathbf r^e_{\text{mut}}) = 0$. The mutant Hamiltonian keeps the stress and relaxation terms; what is zero is $V_{\text{ref}}$, a *different* Hamiltonian. | Julian | §1 now writes all three Hamiltonians explicitly as equations (2)–(4). |
| 2 | Used a synthetic protein (a noisy helix, then a random globule) when `pdb_2acy_A` was in the package. | Julian | Everything re-run on 2acy chain A. |
| 3 | Reimplemented force and linear response instead of reading `R/penm.R` — which is why I missed that penm already linearises ($f_{ij} = -k_{ij}\delta l_{ij}$). | Julian | Catalogue built on penm's own functions. **But see §4: independence was right for `sclfenm` and frustration, which are the things under investigation; it was wrong for the settled parts.** |
| 4 | Confused *response* (marginal over mutated sites) with *influence* (marginal over response sites), and plotted the second labelled as the first. | Julian | §7 rewritten around the response matrix $R_{ij}$ and its two marginals. |
| 5 | Claimed "prestress breaks rotational invariance". Nonsense — $l_{ij}$ and $d_{ij}$ are scalars. | Julian | The three spurious modes came from evaluating the Hessian at the non-stationary LFENM structure ($\lvert\mathbf F\rvert = 0.61$). At the true minimum: exactly 6 zero modes. |
| 6 | Confused two meanings of $\beta$ (Boltzmann vs selection strength). | Julian | $\beta = 1/kT$ throughout; selection is $\nu$. |
| 7 | **Wrote hardcoded conclusions into script output, printed above numbers that contradicted them.** Three times: "SYMMETRIC by construction" above 0.49, 0.996, −0.13. | Julian | Those scripts deleted. Replacements print only what they measure. |
| 8 | RMSIP computed over *all* modes is identically 1 — two random orthonormal bases give 1. It measured nothing. | code audit | Now over the 10 softest modes: 1.000 for lfenm (frozen spectrum), 0.9937 for sclfenm. |
| 9 | §8's "exact for any state function" claim. Exact only in a *fixed* background; used mid-trajectory it defines a surrogate $V_{\text{cat}}$. | physics referee | §8 now distinguishes $V_{\text{cat}}$ from $V$, and reports the per-move error. |
| 10 | §8(c) was circular: $V_{\text{cat}}$ is separable and the proposal symmetric, so a correct sampler *must* reach the product measure. | physics referee | Now labelled a unit test of the sampler; the honest version (chain on $V_{\text{cat}}$, evaluate exact $V$) added. |
| 11 | "The error does not grow with trajectory length" — false. | physics referee | Measured uniformly: 9 % → 48 % (energy), 0.7 % → 17 % (structure), 5 → 40 substitutions. |
| 12 | The $\beta \to \nu$ rename silently broke `fig4_traj.R`: it errored on its own default, and passed `beta=` where the Metropolis rule reads `nu=`, so two "different" runs were identical. | code audit | Fixed and regenerated. Values in the report were right; the script had rotted around them. |
| 13 | Six tests that could not fail (beyond the three §12 already admits). | code audit | I2 rebuilds mutants independently; I3 reports a distribution instead of a forced 50 %; the saddle check reads the raw spectrum by rank. |
| 14 | Unseeded scripts (`run_frustration.R`, `fig3_dv_positive.R`) produced different numbers every run. | code audit | Seeded; report synced. |
| 15 | Stale numbers throughout: values quoted from one script, attributed to another. | code audit, editor | Each traced to its producing script and corrected. |
| 16 | A phantom figure panel: text cited "Panel (c)" of a two-panel figure. | editor | The number lives in fig 5(c); text corrected. |
| 17 | Algebra error in the projector proof: as written, $\mathbf f = -\mathbf B\sqrt{k}\,\mathbf u$ scales as $k^{3/2}\delta l$. | editor | Corrected, with $\mathbf B$ defined concretely. |
| 18 | Symbols used before definition ($\mathbf K$, $\mathbf C$, $\mathbf e_{ij}$, $\mathbf B$); $f$ meaning three different things in one paragraph. | editor | Symbol table extended; Step 2 rewritten. |

**The recurring failure mode is #7 and #8 and #10: producing a check that cannot
report anything but success.** Three distinct instances. Worth naming because it
is worse than an ordinary error — it manufactures confirmation.

---

## 2. Reviewers

Three, all with the code and able to run it.

**Physics referee** (twice). First round: found the projector proof I had missed
(a cleaner result than my numerical demonstration), and that my mechanism for
`keepnet` was contradicted by my own output. Second round: found errors 9, 10, 11
above, and supplied two physics gaps — where relaxedness enters the projector
proof ($\mathbf K = \mathbf B k \mathbf B^{T}$ holds only for a relaxed network:
verified, exact when relaxed, off by 0.23/12.9 when frustrated), and why exactly
three modes are spurious off-stationarity (translations are annihilated
identically; rotations only via $\mathbf F = \mathbf 0$).

**Code auditor** (twice). Verified every number against the script that produces
it. Found errors 8, 12, 13, 14, 15.

**Editor** (once). Found errors 16, 17, 18, plus structural problems: the system
description arriving in §7 after six sections of numbers drawn from it, headings
that stated topics rather than findings, and five passages of draft archaeology
that were cut.

---

## 3. Checks that exist, and which can fail

| suite | count | can fail? |
|---|---|---|
| `validate.R` V1–V5 | 5 | V1, V2 live; V3, V5 and half of V4 cannot fail — stated as a defect in §12 |
| `derivation.R` | 11 | all live: each compares an analytic step against an independent numerical evaluation |
| `test_identities.R` | 4 | I2 and I4 live after the fix; I1 is a genuine linearity test |
| `test_can_fail.R` | 2 sabotages | demonstrates V1 and V2 go red when the code is broken |
| saddle check (§9) | 1 | proved it fires: at $l = 1.5\,d$ it detects a saddle |

---

## 4. Two implementations, and whether they agree

The scripts split by when they were written:

- **§§1–7 use a reimplementation** (`enm_core.R`, `sclfenm.R`, `applications.R`)
  written before I was told to read `R/penm.R`.
- **§§8–10 use penm's own functions** (`set_enm`, `get_mutant_site`,
  `calculate_force`, `calculate_dxyz`, `ddg_dv`) via `penm_catalog.R` and
  `frustration.R`.

They were never checked against each other until asked. They agree
(`test_two_impls.R`, ANM / $d_{\max} = 10.5$ / 2acy):

| quantity | agreement |
|---|---|
| contact set | identical (960 edges, same pairs) |
| rest lengths $l_{ij}$ | exact |
| Hessian | $3.6\times10^{-15}$ (scale 12.86) |
| spectrum (288 modes) | $1.1\times10^{-14}$ |
| force from the same $\delta\mathbf l$ | $1.3\times10^{-15}$ |

So the report's numbers are consistent across the split. Two things to note
anyway:

1. **The duplication is a defect for the settled parts, and correct for the
   questioned ones.** ANM construction, force and linear response are not in
   doubt — reimplementing those was waste, and it cost the linearisation
   insight. But `sclfenm` and `frustrated` are precisely what the package flags
   as unchecked (`stopifnot(!frustrated)`, skipped tests), so building those
   independently is the point of the exercise, not a mistake. Having both is what
   made §11 possible: the independent implementation is what the package code can
   be measured *against*.
2. **The mutation model is not penm's default.** `enm_core.R` perturbs *every*
   contact of a site; penm's `generate_delta_lij` honours `mut_sd_min` (default
   2, which excludes backbone neighbours). For §§8–10 I passed
   `mut_sd_min = 1L` to match §§1–7. So the whole report uses
   $\texttt{mut\_sd\_min} = 1$, not the package default.

---

## 5. Why this took as long as it did

Mostly rework caused by errors 1–7, each of which invalidated a round of results
and required re-deriving from scratch. The reviewing was fast; the recovery was
not. Concretely:

- Building a synthetic protein instead of using `pdb_2acy_A` cost one full cycle
  of figures and numbers, all of which had to be regenerated.
- Reimplementing penm's force and response cost the linearisation insight, which
  is the foundation of the catalogue model — the central result. Reading
  `R/penm.R` first would have saved the largest single detour.
- The three hardcoded-conclusion scripts each produced a wrong answer that had to
  be found and undone.

**What would have prevented most of it:** read the existing implementation before
writing any of my own, and use the package's own test data. The nuance learned
late (§4): reuse the package for what is settled, reimplement only what is under
investigation — and say which is which from the start.
