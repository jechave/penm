# Evolutionary simulations with the generalized ENM — results and protocols

Companion to the model report. Nothing here restates the model; this records what was
run, with enough detail to repeat it, and what came out.

---

## 1. Common setup

**Structure.** 1A6M (myoglobin, oxy, 1.00 Å), chain A, CA atoms only, altloc blank or A.
N = 151 nodes. Coordinates used as deposited, no superposition or trimming.

**Network.** Springs between all pairs closer than 14 Å in the PDB structure: M = 2858.
The cutoff is a numerical truncation only, since k is negligible beyond ~13 Å.

**Force constants.** Smoothed step, evaluated at the equilibrium length of each spring:

    k(l) = (k0/2) [ 1 − tanh( (l − Rc) / W ) ],   Rc = 10.5 Å,  W = 1.0 Å

so that k = k(l_ij) and the force constants change as the parameters evolve. The
smoothing avoids the discontinuity of a hard cutoff; W = 1 Å gives k(10.5) = k0/2,
k(13) = 0.007 k0.

**Sequence → parameters.** Pseudo-sequence s of N labels, each in {0,…,19}:

    l_ij(s) = d_ij(X) + u_ij(s_i) − u_ij(s⁰_i) + v_ij(s_j) − v_ij(s⁰_j)

with u, v drawn once as iid N(0, σ²) arrays of shape (20, M), quenched for the whole
study (seed 2024 for the final runs). Both endpoints contribute, so each contact has
400 reachable lengths; s⁰ = all zeros makes the PDB structure the unfrustrated
reference. This is a function of s alone — no accumulation — which is what gives the
state-function property. (Written as differences because that is what the code does;
whether it is equivalent to drawing u_ij(0) = 0 was not verified.)

**Energy and structure.** V(r; L) = ½ Σ k_ij (d_ij(r) − l_ij)², minimized by L-BFGS-B
with the analytic gradient, warm-started from the previous structure (ftol 1e-18,
gtol 1e-12). Hessian assembled blockwise including the transverse term
g_ij = (d_ij − l_ij)/d_ij. Structures superposed by Kabsch before any RMSD.

**Evolutionary process.** Pick a site uniformly, pick a new label uniformly, rebuild L
and k, minimize, accept under the selection rule, otherwise revert the label. Threshold
rule: accept iff V_min < V_thr. Metropolis rule where stated: accept with
min(1, exp(−α ΔV_min)).

**Parameters of the final runs.** σ = 0.87 Å, k0 = 0.25 kcal mol⁻¹ Å⁻², V_thr = 120
kcal/mol. How they were fixed is in §3.

The module used for the final runs is included as `enm_evo_model.py`.

---

## 2. Verifications (run first, before any result)

**State function.** Same sequence reached by two different orders of mutation gives
identical V_min (difference exactly 0, printed to 8 decimals). A mutation followed by
its reversion returns to the original V_min with zero drift. Holds with k = k(l), not
just with frozen k.

**K evolves.** For a random sequence: |δk|/k averages 1.5% with 107 springs changing by
more than 10%; eigenvalues move 1.1% on average. Splitting the channels by recomputing
the Hessian with the ancestral k, or with g = 0: the **k channel contributes 3.5% and
the transverse channel 1.2%**. With k frozen at the reference — as in several
intermediate runs of this session — the dominant channel does not exist at all.

**Rigidity.** 6 null modes before and after mutation; the network never goes floppy at
these σ.

---

## 3. Calibration: three observables, three parameters

The scalings decouple, which is what makes the fit possible without compromise.
Verified numerically by scaling k by 0.25 and 4 (RMSD per mutation stays 0.0293 Å in
all three cases):

| observable | depends on |
|---|---|
| RMSD per substitution | σ only — invariant under uniform scaling of k |
| ΔΔG | k σ² |
| thermal MSF | k_B T / k |

**σ = 0.87 Å** from structural divergence. Diffusive models write
D² = a² d + D₀² with d in substitutions per site (Gutin & Badretdinov 1994;
Grishin 1997). Grishin's estimates of a² are 0.87–1.28 Å² across three datasets, with a
theoretical estimate near 1. A point-mutation scan (302 mutations, all sites) gives mean
RMSD 0.0872 Å, i.e. **a² = N·RMSD² = 1.15 Å²**.

**k₀ = 0.25 kcal mol⁻¹ Å⁻²** from the ΔΔG scale: the same scan gives mean ΔΔG
+1.32 kcal/mol.

**V_thr = 120** from the shape of the ΔΔG distribution (§5).

**B-factors are not usable for calibrating k.** 1A6M was collected at 100 K, so its
B-factors are mostly static disorder, and an earlier calibration against them (giving
k ≈ 2.4) was discarded.

---

## 4. Entropic contribution to ΔΔG

Protocol: from the relaxed wild type, mutate one site (76 sites sampled), minimize,
and compute ΔU = ΔV_min and TΔS = (1/2β)[Σ' ln λ_n − Σ' ln λ'_n] over the 3N−6 non-null
modes, β from T = 298 K.

Result with the final parameters: `<ΔU> = +1.37`, `<|TΔS|> = 0.113` kcal/mol, i.e. **8%**.
The mean of TΔS is +0.006 against a spread of 0.140 — symmetric noise, 47% positive —
and **no** mutation changes the sign of ΔA relative to ΔU.

With Ming & Wall parameters and frozen k (906 mutations, all sites) the fraction was
1–3%, and there TΔS was found to be uncorrelated with the site environment
(|r| < 0.04 with contact number and with MSF) while ΔU correlates strongly (+0.80 and
−0.71). The difference between 1–3% and 8% is parameterization: the fraction scales as
1/k, and TΔS grows as σ while ΔU grows as σ².

**ΔΔG ≈ ΔV_min** to about 10%. The entropic term contributes scatter between mutations,
not a shift.

Related: dropping the transverse terms from K — the approximation of rebuilding the ENM
on the relaxed structure — changes TΔS by a factor of 2 (measured with Ming & Wall
parameters at σ = 0.3 Å), hence ΔΔG by a few percent. Structural and fluctuation
observables are untouched by it: eigenvalues differ by 0.03% on average, per-residue MSF
correlates at r > 0.9999.

---

## 5. Distribution of ΔΔG

Protocol: equilibrate a threshold chain (1400 attempted mutations), then scan from the
stationary state, 2 random labels per site.

| V_thr | V_min | acceptance | mean | sd | < −0.5 | within ±0.5 | > 2 |
|---|---|---|---|---|---|---|---|
| 120 | 120 | 0.28 | +1.49 | 1.44 | 5.6% | 19% | 30% |
| 170 | 167 | 0.46 | +0.96 | 1.30 | 9.9% | 27% | 21% |
| 190 | 190 | 0.54 | +0.84 | 1.37 | 13.9% | 33% | 18% |
| 210 | 210 | 0.65 | +0.72 | 1.34 | 15.2% | 32% | 15% |

Empirical references, which disagree among themselves by an order of magnitude:

| source | nature | stabilizing |
|---|---|---|
| Tsuboyama 2023 | exhaustive, proteolysis, ~480 domains | ~4% (< −0.5); 0.2–0.6% (< −1) |
| FireProtDB | curated literature | 14% |
| ThermoMutDB / curated ProTherm | curated literature | 27–29% (mean 1.0, sd 1.6) |
| FoldX, Tokuriki 2007 | prediction, 21 proteins | ~30%; surface 0.6, core 1.4 |

The databases over-represent stabilizing mutations (engineering studies, alanine scans);
Tsuboyama measures every substitution of each domain and is the unbiased one, though its
domains are small (<72 residues) and the observable is proteolysis-derived.

**V_thr = 120** is the choice adopted: it lands between Tsuboyama and the databases on
the stabilizing tail, and matches the curated ProTherm mean and spread.

A separate neutral-drift sweep (no selection, sampling ΔΔG spectra as V_min rises)
showed the distribution translating rigidly: sd stays 1.5–1.9 throughout while the mean
falls from +1.88 at V_min = 51 to 0.00 at V_min = 300, where random sequences sit.
Stabilizing fraction rises 15% → 21% → 33% → 52%. So the position of the mean is set by
regression toward the random-sequence mean, and the stabilizing fraction measures where
the protein sits in that range — not the selection strength directly.

---

## 6. Evolutionary behaviour

**Marginal stability.** Under threshold selection V_min rises and pins to the threshold
(99.0/100, 198.4/200, 397.8/400 in earlier runs at σ = 0.3), with a stationary band of
~2 kcal/mol and 40–68% of sampled sequences within 0.5 kcal/mol of the ceiling. This
comes from the sequence entropy, not from any assumption.

**Substitution continues.** Rate falls during the approach and then holds constant
(0.09, 0.17, 0.39 per attempt for the three thresholds in the last blocks). Plot against
accepted substitutions, not attempted mutations, or the plateau looks like a stall.

**Divergence** (final parameters, V_thr = 120, measured from an equilibrated state):

| substitutions | RMSD | corr MSF | contact map kept | mean N_n |
|---|---|---|---|---|
| 10 | 0.133 Å | 0.9995 | 0.997 | 4.3 |
| 25 | 0.294 | 0.9951 | 0.988 | 13.4 |
| 50 | 0.365 | 0.9854 | 0.987 | 17.3 |
| 100 | 0.641 | 0.9662 | 0.973 | 33.5 |
| 150 | 0.748 | 0.9269 | 0.971 | 38.5 |

Structure and contact map are the most conserved; the identity of individual modes
degrades fastest. Eigenvalues barely move (0.7% at 400 substitutions in the σ = 0.3
runs) while the eigenvectors rotate substantially.

**Mode mixing.** N_n = exp(H_n) with H_n = −Σ_m |S_nm|² ln |S_nm|², S = Q₀ᵀQ′
restricted to the internal subspace, ancestor modes as reference. Rows of |S|² are
normalized, so N_n runs from 1 (mode preserved) to 3N−6 (dissolved). Profile is an
**inverted U**: soft and stiff modes relatively preserved, maximum mixing in the middle
of the spectrum, growing with evolutionary time (medians ~3 → 33 → 3 across the spectrum
at 200 substitutions).

*Missing control*: part of the inverted U may be spectral density rather than evolution,
since gaps are smallest in the middle. The control is N_n between the ancestor and a
random perturbation of the same size; the evolutionary signal is the excess over that.
`analyses/run_control.R` in this exploration already implements that kind of matched
null for mode overlaps and can be reused.

---

## 7. Open

**Substitution rate too high.** Acceptance 0.28 at V_thr = 120. Expected: under
stability selection alone every stabilizing mutation is acceptable, so the observed
4–30% stabilizing fraction is already a lower bound on what stability would let through.
Low observed rates require activity constraints.

**Diffusive accumulation.** RMSD ∝ √n in the model, which is what Gutin & Badretdinov
and Grishin predict, but it is not established empirically that real divergence grows
that way. Wood & Pearson's linearity is between z-scores, not between RMSD and
substitutions, and the z transformation is plausibly what straightens the curve. Their
differential result — 0.18–0.19 Å per 5.4 points of identity, the same in two very
different regimes — is the robust part. Measuring RMSD against substitutions per site
directly in one family (globins, core alignments, no z-scores) is what is missing.

**The ancestor.** Re-anchoring the map on an evolved state leaves that state with zero
disorder, hence atypical, and reintroduces the artifact of the unfrustrated wild type:
from there almost every mutation is destabilizing (mean +2.25, 1% below −0.5, against
+1.49 and 5.6% for the same threshold without re-anchoring). With the map anchored on
the PDB and one continuous chain the problem does not arise.

**Width of the ΔΔG distribution.** The model gives sd ≈ 1.4; matching mean and
stabilizing fraction simultaneously would need ≈ 2.5. Giving each label its own
perturbation scale (lognormal, CV 0.8) drops mean/sd from 1.14 to 0.58 and fixes the
shape, at the cost of one more parameter.

**Displacement per unit energy.** Neither the 1/√k weighting nor heterogeneity changes
the ratio RMSD²/ΔΔG by more than a factor 2: a force localized on one site spreads its
energy evenly over the spectrum, because ⟨(q_n·f)²⟩ ∝ q_nᵀ K_{k²} q_n ∝ λ_n, so every
mode absorbs the same energy. Coupling preferentially to soft modes would need a force
covariance flatter than K, which a single-site perturbation cannot produce.

---

## 8. Caveats

- One protein, one disorder realization, chains of 10²–10³ substitutions.
- Several intermediate runs of this session used **frozen k** and are void for anything
  involving dynamics; only the runs described here use k(l).
- σ = 0.87 Å is a large perturbation. It does not compromise anything, because the
  potential is exact in δl and r_e comes from minimization — the quadratic
  approximation is in r, not in δl — but it is worth keeping in view.
- The hash-based disorder used in the early runs was a homemade sine hash; the final
  runs use numpy's generator.
