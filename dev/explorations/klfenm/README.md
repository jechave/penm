# klfenm

Generalized elastic network model: the equilibrium lengths `l_ij` are the
parameters, the force constants follow from them as `k_ij = k(l_ij)`, and the
structure is obtained by minimization. Unlike LFENM, there is no linear response
and no forcing — the name is a leftover and should be changed.

## Contents

- `klfenm_report.tex` / `.pdf` — the model: potential, quadratic expansion,
  normal modes, thermal properties, conformational fluctuations, and the
  mutational changes. Frustration is in an appendix.

- `klfenm_perturbation.tex` / `.pdf` — first-order perturbation theory, kept
  separate so the main derivation stays exact. Two results worth not
  rediscovering: δK must be projected onto the internal subspace before any
  eigenvector or covariance formula is applied, or the error does not vanish as
  the mutation shrinks; and δK has two contributions of the same order, one from
  the parameter change at fixed conformation and one from the conformational
  change at fixed parameters, the second involving third derivatives of V.

- `simulations-report.md` — protocols and results of the evolutionary
  simulations: verifications, parameter calibration against structural
  divergence and ΔΔG, the entropic contribution, ΔΔG distributions, divergence
  of structure versus dynamics, mode mixing, and what is open.

- `enm_evo_model.py` — the module used for those simulations. Needs `ca.npy`,
  the CA coordinates of 1A6M.

Everything else in this directory was scaffolding and has been removed.
