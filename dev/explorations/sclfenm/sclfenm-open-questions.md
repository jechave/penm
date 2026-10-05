# Open questions about sclfenm

Notes from a working session in October 2026. Nothing here was implemented: this
records where the argument stands, so it does not have to be rebuilt from
scratch. It is written to be read cold.

---

## 0. What the models are

`penm` represents a protein as an elastic network: nodes are residues, springs
join pairs closer than `d_max` in the native structure, and the potential is

    V(r) = Σ_ij ½ k_ij (d_ij(r) − l_ij)²

with `d_ij(r)` the distance in conformation `r`, `l_ij` the spring's rest length
and `k_ij` its force constant. A mutation at a site perturbs the `l_ij` of the
springs touching it, and the structure relaxes.

Three ways of building the mutant, which differ in what is recomputed:

- **lfenm** — the structure moves by linear response, `r_e' = r_e + K⁺f`, and the
  network is untouched. `K` is therefore identical for wild type and mutant, so
  the mutation changes the structure but not the dynamics.

- **sclfenm** — the structure moves by linear response as in lfenm, and then the
  network is *rebuilt* on the new structure: the contact list, the `k_ij` and the
  unit vectors are recomputed from the relaxed coordinates. Now `K` does change,
  which is what makes mutational effects on dynamics accessible. "Self-consistent
  LFENM" is what the name means.

- **the generalized model** (developed separately, in `../klfenm/`) — the `l_ij`
  are the primitive parameters, `k_ij = k(l_ij)` follows from them, and the
  mutant's structure is obtained by *minimizing* its own potential rather than by
  linear response. Numbers quoted below as measured "in the generalized model"
  come from there.

There is also a flag `frustrated`. A network is *frustrated* when its springs
cannot all be relaxed at once, so `d_ij ≠ l_ij` at the minimum. The exact Hessian
of the potential above then has two parts per spring,

    K_ij = −k_ij [ e_ij e_ijᵀ + g_ij (I − e_ij e_ijᵀ) ],   g_ij = (d_ij − l_ij)/d_ij,

a longitudinal one and a transverse one. `frustrated = FALSE` drops the
transverse term, i.e. sets `g_ij = 0`; `frustrated = TRUE` keeps it and is
currently blocked in `set_enm()`.

---

## 1. One claim retracted: sclfenm does not need to relax to a minimum

It is tempting to object that sclfenm's mutant structure is a single linear step,
not a minimum of any potential, so the Hessian is evaluated off a stationary
point and the spectrum is meaningless.

**That objection is wrong when `frustrated = FALSE`.** With `g_ij = 0` the
Hessian is `−k_ij e_ij e_ijᵀ` built from the unit vectors of whatever structure
you give it: that is the ordinary ANM Hessian of that structure, positive
semidefinite with exactly six null modes by construction, minimum or not. No
spurious modes appear. (The earlier finding in `MODEL.md` of near-zero modes and
a badly wrong TΔS was for the *frustrated* Hessian off a stationary point, which
is a different case.)

What does remain is an inconsistency of objects rather than of numbers: `V_min`
is computed from the frustrated potential, because the rebuilt graph keeps the
mutated `l_ij`, while the modes come from the rebuilt network with the transverse
term dropped. Two Hamiltonians for one mutant. Defensible, but worth stating in
the documentation.

And "fixing" sclfenm by relaxing to the true minimum would not fix it: it would
turn it into the generalized model under another name. The linear-response
structure is sclfenm's definition, not a defect.

---

## 2. What `mutate_graph()` does with the contact map

When the network is rebuilt on the new structure, pairs can enter and leave the
contact list. The code treats three classes of edge by three different rules:

| edge | `l_ij` | `k_ij` | energy |
|---|---|---|---|
| present in both graphs | kept from the previous graph — so it carries the history of past mutations | recomputed from the new `d_ij` | kept |
| new (now inside `d_max`) | set to `d_ij`, i.e. born relaxed | from the new `d_ij` | zero |
| lost (now outside `d_max`) | gone | gone | its stored strain disappears |

Also, `calculate_enm_graph()` resets `v0ij = 0`, so any energy offset is erased
at every rebuild.

One fragility: the line

    g2[g2$edge %in% g1$edge, "lij"] <- g1[g1$edge %in% g2$edge, "lij"]

assigns by position, which is correct only because both graphs come from
`expand_grid %>% arrange(i, j)` and subsets preserve relative order. Matching on
the edge key would be safer.

---

## 3. The asymmetry, and what is really wrong

A lost contact leaves taking its accumulated strain with it, so `V_min` drops; a
new contact arrives relaxed and costs nothing. Both push the same way, so contact
turnover reads as stabilizing.

That looks like a bookkeeping bug but is not. In sclfenm the structure defines
the network, so a pair that was not a contact has no `l_ij` to inherit, and the
only local choice is `l_ij = d_ij`. Within the model's own logic, a new contact
born relaxed is correct.

The real defect is a discontinuity, and it lives in `k`, not in the accounting:
with a step cutoff, how much the potential jumps when a pair crosses `d_max` does
not depend on by how much it crossed. A pair that moved 0.001 Å produces the same
jump as one that moved 1 Å. That is what makes `V_min`, the spectrum and TΔS
discontinuous functions of the structure.

---

## 4. Three ways out, and what each costs

**Freeze the wild-type contact list.** The rebuild updates `e_ij` and `k_ij` but
not which pairs are edges. No jumps, nothing enters or leaves, `V_min` and TΔS
continuous, and wild type and mutant always comparable. The cost: it is no longer
self-consistent in the strict sense, since rebuilding the network from the new
structure is supposed to include deciding which pairs are in contact.

**Smooth the cutoff**, e.g. `k(d) = (k0/2)[1 − tanh((d − d_max)/w)]`. On its own
this only halves the jump, since a pair at the cutoff still enters with
`k ≈ k0/2`. It works only together with **widening the pair list** out to where
`k` is already negligible (with `w = 1 Å`, `k(13 Å) ≈ 0.007 k0`): then a pair that
"forms" was already in the graph, carrying the `l_ij` it had, nothing is born or
dies, and only weights move. This removes both discontinuities — but carrying
`l_ij` for non-contacts is halfway to the generalized model, the remaining
difference being that `k` still comes from `d` rather than from `l`.

**Keep lost edges with `k = 0`** instead of deleting them. Equivalent to freezing
the contact list, stated so that nothing needs deciding: the pair set is fixed
and what varies is `k`. Also implies carrying `l_ij` for inactive pairs.

Reading: the clean option that stays inside sclfenm is freezing the contact list,
accepting that it is then not strictly self-consistent. The others lead to the
generalized model.

---

## 5. Comparing on the shared graph: works for energy, not for entropy

A natural idea is to restrict comparisons to the edges present in both networks.

**For `ΔV_min` it works.** The energy is a sum over edges, so dropping the same
terms from both sides is well defined and continuous. What is then measured is
the deformation energy of the shared part of the network, which is not the same
as the difference between the two proteins' energies — if a contact is genuinely
lost, that energy exists and is being discarded.

**For TΔS it does not.** The spectrum is not additive over edges. Building both
Hessians from the shared edges diagonalizes two networks that are neither
protein, and removing edges can leave them floppy.

Integrating out the excluded degrees of freedom — `K_eff = (C_cc)⁻¹`, the Schur
complement, which is what `kmat_asite()` already does for the active site — **does
not apply here**. That marginalizes *nodes*; this is a problem of *edges*, with
the same nodes throughout, so the two matrices already have the same dimension
and there is nothing to integrate out.

---

## 6. Which functions break when the contact map changes

- `delta_structure_dvmi()` — joins the two graphs on `edge`, so unshared edges are
  silently dropped and the per-site profile no longer sums to `ddg_dv()`.
  Redundant in any case: `get_stress(mut) − get_stress(wt)` gives the same
  profile, computed on each protein's own graph.
- `delta_structure_dvsi()` — asserts that the two edge lists are identical, so it
  errors.
- `delta_structure_dvsi_same_topology()` — the same comparison by row position,
  without the assertion: wrong numbers, silently.
- `ddg_tds()` — subtracts sums over the two spectra; if the number of non-null
  modes differs it compares sums of different length. Needs a guard.
- `delta_motion_*()` — compare modes of two different structures without
  projecting onto the internal subspace. The six rigid-body null directions
  include three rotations, which are generated by the node positions and so
  differ between wild type and mutant at first order in the displacement.
  Measured in the generalized model: without the projection `P = K⁺K`, the error
  in mode overlaps and in the covariance **does not vanish** as the mutation is
  made smaller — it stays around 10%; with the projection it converges linearly.
  The entropy, being a trace, projects itself, which is why that channel never
  showed the problem.
- `calculate_vs()` — the only one that handles a changed contact map correctly,
  by walking the edges of `prot` and recomputing any missing distance from
  `ideal`'s coordinates.

---

## 7. Numbers measured in the generalized model that bear on sclfenm

On myoglobin (1A6M, 151 CA nodes), not on the package, but they transfer.

**Dropping the transverse term is harmless for dynamics and not for entropy.**
Comparing the frustrated Hessian with the rebuilt one at the same structure:
eigenvalues differ by 0.03% on average, per-residue mean-square fluctuations
correlate at r > 0.9999, but TΔS differs by a factor of two, and at the largest
mutation tested it differs in sign.

**But TΔS is a small part of ΔΔG.** It is 1–10% depending on the parameters
(1–3% with Ming & Wall constants, 8% with the values fitted in that session), its
mean is two orders of magnitude below its own spread, it changes sign in half the
mutations, and it does not correlate with the site's environment (|r| < 0.04 with
contact number and with fluctuation) while ΔU correlates strongly (+0.80 and
−0.71). So `ΔΔG ≈ ΔV_min`, and the error introduced by dropping the transverse
term is a few percent of ΔΔG.

**Mutations barely move the structure.** 0.03 Å RMSD per point mutation at
σ = 0.3 Å, with the contact map conserved at 97–99.9%. So the pairs that flip in
or out are those sitting essentially at the cutoff — numerically marginal, not
structural reorganization. This bounds how much freezing the contact list would
really restrict a trajectory, at least per mutation; over a long trajectory it was
not measured.

**The linear-response error is diffusive, not systematic.** Cycles do not close
under repeated linear response with rebuilds: after a forward-and-back cycle the
residual grows as √n, reaching 3–7% of one step's displacement after 20
substitutions, and it is **independent of the order of reversion** (last-in
first-out, first-in first-out and random agree to 6%). Most of it comes from
recomputing the forces on the moved structure rather than from rebuilding `K`:
freezing `K` alone still leaves two thirds of the error, while freezing both `K`
and the forces closes the cycle to machine precision — at the price of a strictly
linear model.

---

## 8. Energy bookkeeping: a correction

It is tempting to think that keeping `V_min` while rebuilding the network still
ratchets upward, because the rebuild zeroes the strain and the next mutation is
measured from an unstrained reference.

**That does not apply to the scheme used in the 2019 paper**, where every state's
energy is referenced to the wild type and `ΔΔG_ji = ΔΔG_j0 − ΔΔG_i0`. Both terms
are computed from state 0, so the energy is a state function by construction,
reversions refund exactly, and cycles close whatever the path. The ratchet is
real only for a scheme that applies each new `δl` to the *rebuilt* lengths
sequentially — a path construction rather than an indexed one.

Separate and still true: `mut$graph$lij <- wt$graph$lij + delta_lij` defines the
mutant's parameters relative to the *current* wild type. For scans from a fixed
reference that is fine. For trajectories the state becomes the path rather than
the sequence, unless the `δl` are attached to per-site labels so that the
parameters are a function of the sequence.

---

## 9. If sclfenm is kept as an option

The three models are nested rather than three different physics. lfenm is the
first order of the generalized model — expanding the stationarity condition of
the mutant potential gives exactly `r_e' − r_e ≃ K⁺f` — and sclfenm is the
intermediate case: structure to first order, network rebuilt around it.

So the package could expose `lfenm` and the generalized model, with sclfenm
reproducible as a configuration of the latter: linear-response structure instead
of exact minimization, `k` recomputed from the structure. One code path, and
results that depend on sclfenm remain reproducible.
