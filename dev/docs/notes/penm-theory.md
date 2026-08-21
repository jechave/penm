# Perturbed Elastic Network Models: theory

**Status: draft.** Self-contained. Develops the theory of elastic network
models, their perturbation by mutation, and the observables used to compare a
mutant with its wild type. Open questions are recorded in §8 without taking a
position.

Notation follows Echave, *Biophysical origin of long-range functional
constraints in enzyme evolution* (2017), whose §2 and appendix are the basis of
§§1–3 here.

---

## 1. The elastic network model

### 1.1 Nodes

A protein of $N$ residues is represented by $N$ point nodes. A conformation is
the column vector of Cartesian coordinates

$$
\mathbf{r} = (x_1, y_1, z_1, \ldots, x_N, y_N, z_N)^T ,
$$

of dimension $3N$, with $\mathbf{r}_i = (x_i, y_i, z_i)^T$ the position of node
$i$. The distance between nodes $i$ and $j$ is
$d_{ij} = \lVert \mathbf{r}_j - \mathbf{r}_i \rVert$, and

$$
\mathbf{e}_{ij} = \frac{\mathbf{r}_j - \mathbf{r}_i}{d_{ij}}
$$

is the unit vector along the $i \leftrightarrow j$ direction, with
$\mathbf{e}_{ij} = -\mathbf{e}_{ji}$.

The identification of a node with a residue is a modelling choice — the
$\alpha$-carbon, the $\beta$-carbon, or a side-chain centroid — and nothing in
what follows depends on which is used. More than one node per residue is
possible; then $N$ counts nodes, not residues.

### 1.2 The network

The nodes are joined by springs, and which pairs are joined is part of the
model. Write $\mathcal{E}$ for the set of joined pairs — the **contact set** —
and $i \sim j$ for $\{i,j\} \in \mathcal{E}$. Together, nodes and contacts make
the network a graph; $M = |\mathcal{E}|$ is the number of springs.

The contact set is fixed once, from a **reference conformation**, normally the
protein's experimentally known native structure. Models differ in how they
choose it, and the choice ranges over the whole spectrum:

- the **complete graph**, $\mathcal{E} = \{\{i,j\} : i \neq j\}$, every pair
  joined, $M = N(N-1)/2$;
- a **proper subset**, containing only pairs that are close in the reference
  conformation, or close in sequence, or both, with $M$ of order $N$.

Both are in use — the parameter-free ANM takes the complete graph and lets the
force constants carry the whole distance dependence, while other models restrict
$\mathcal{E}$ and let $k_{ij}$ vary less, or not at all — and the theory below is
indifferent between them. All sums
$\sum_{i<j}$ run over $\mathcal{E}$, and $\sum_{j\sim i}$ over the neighbours of
$i$ in the graph.

For any of these choices $M \gg 3N$ for a protein of appreciable size, a fact
§3.2 depends on.

### 1.3 The potential

Each contact carries a harmonic spring, and the potential energy is the sum over
contacts

$$
\boxed{\;
V(\mathbf{r}) = \sum_{i<j} \left[ V^{\min}_{ij} + \tfrac{1}{2} k_{ij}\left(d_{ij} - l_{ij}\right)^2 \right] .
\;}
\tag{1}
$$

Each spring has three parameters:

- $V^{\min}_{ij}$, the energy at the spring's own minimum — an additive offset,
  irrelevant to forces and to the Hessian, but carried because mutations may in
  principle change it;
- $k_{ij} > 0$, the force constant, normally taken to decrease with the pair's
  separation in the reference conformation, so that distant contacts are softer
  than close ones;
- $l_{ij} > 0$, the **natural length** (rest length) — the separation at which
  that spring alone is unstrained.

How $\mathcal{E}$ is chosen and how $k_{ij}$ depends on the reference structure
are what distinguish one elastic network model from another; nothing below
depends on either choice.

The essential structural point, and the one from which everything in §3 follows,
is that $d_{ij}$ and $l_{ij}$ are quantities of different kinds. The $d_{ij}$
are functions of the conformation, and are not independent of one another: $3N$
coordinates determine all $M$ of them, and for $M > 3N$ they satisfy
constraints. The $l_{ij}$ are $M$ independent parameters of the model.

## 2. Equilibrium

### 2.1 Force and Hessian

Expanding (1) to second order about an arbitrary conformation $\mathbf{r}^0$,

$$
V(\mathbf{r}) \simeq V(\mathbf{r}^0) - \mathbf{F}(\mathbf{r}^0)\left(\mathbf{r}-\mathbf{r}^0\right) + \tfrac{1}{2}\left(\mathbf{r}-\mathbf{r}^0\right)^T \mathbf{K}(\mathbf{r}^0)\left(\mathbf{r}-\mathbf{r}^0\right) ,
\tag{2}
$$

with the force row vector and the Hessian

$$
\mathbf{F}(\mathbf{r}) = -\frac{\partial V}{\partial \mathbf{r}}, \qquad
\mathbf{K}(\mathbf{r}) = \frac{\partial^2 V}{\partial \mathbf{r}^2} .
$$

The derivatives of the distances are

$$
\frac{\partial d_{ij}}{\partial \mathbf{r}_j} = \mathbf{e}_{ij}^T, \qquad
\frac{\partial^2 d_{ij}}{\partial \mathbf{r}_i \,\partial \mathbf{r}_j} = \frac{\mathbf{e}_{ij}\mathbf{e}_{ij}^T - \mathbf{I}}{d_{ij}} ,
\tag{3}
$$

whence the force on node $i$,

$$
\boxed{\;
\mathbf{F}_i(\mathbf{r}) = \sum_{j \sim i} k_{ij}\left(d_{ij} - l_{ij}\right)\mathbf{e}_{ij}^T .
\;}
\tag{4}
$$

Each term is a central force along the line joining the pair: a stretched spring
($d_{ij} > l_{ij}$) pulls $i$ toward $j$, a compressed one pushes it away. Since
$\mathbf{e}_{ij} = -\mathbf{e}_{ji}$, the forces sum to zero,
$\sum_i \mathbf{F}_i = \mathbf{0}$: the network is isolated, with no external
force.

### 2.2 The equilibrium condition

The equilibrium conformation $\mathbf{r}^e$ is the minimum of $V$, where every
node feels zero net force:

$$
\boxed{\;
\mathbf{F}_i(\mathbf{r}^e) = \sum_{j \sim i} k_{ij}\left(d^e_{ij} - l_{ij}\right)\mathbf{e}^{e\,T}_{ij} = \mathbf{0}, \qquad i = 1 \ldots N .
\;}
\tag{5}
$$

These are $3N$ scalar equations, of which $3N-6$ are independent — six
combinations are absorbed by the rigid-body freedoms, three translations and
three rotations, along which $V$ does not vary. They are non-linear in the
coordinates, since $d^e_{ij} = \lVert \mathbf{r}^e_j - \mathbf{r}^e_i\rVert$ and
$\mathbf{e}^e_{ij}$ both depend on $\mathbf{r}^e$.

Given the parameters, $\mathbf{r}^e$ may be found by iterated quadratic
minimisation: from a current conformation $\mathbf{r}^0_n$, the minimum of the
quadratic form (2) is at

$$
\mathbf{r}^0_{n+1} = \mathbf{r}^0_n + \mathbf{K}^{-1}(\mathbf{r}^0_n)\,\mathbf{F}^T(\mathbf{r}^0_n) ,
\tag{6}
$$

which converges to a conformation satisfying (5).

### 2.3 The Hessian

From (3), the Hessian decomposes into $3\times3$ blocks

$$
\boxed{\;
\mathbf{K}_{i \neq j} = -k_{ij}\left[\mathbf{e}_{ij}\mathbf{e}_{ij}^T + g_{ij}\left(\mathbf{e}_{ij}\mathbf{e}_{ij}^T - \mathbf{I}\right)\right], \qquad
\mathbf{K}_{ii} = -\sum_{j \sim i}\mathbf{K}_{ij} ,
\;}
\tag{7}
$$

with $\mathbf{K}_{ij} = \mathbf{0}$ for $i \nsim j$, and where $g_{ij}$ is a
function of the ratio $l_{ij}/d_{ij}$ that vanishes when the spring is at its
natural length. **The sign convention for $g_{ij}$ is unsettled; see §8.1.** Up
to that sign,

$$
g_{ij} = \pm\left(\frac{l_{ij}}{d_{ij}} - 1\right).
$$

The two terms have distinct geometric meanings.

- $\mathbf{e}_{ij}\mathbf{e}_{ij}^T$ is the projector onto the line joining the
  pair: a **longitudinal** stiffness, resisting stretching of the spring.
  It is rank 1, so a pair of nodes joined by this term alone resists only
  relative motion along the line joining them.
- $\mathbf{e}_{ij}\mathbf{e}_{ij}^T - \mathbf{I}$ is minus the projector onto
  the plane perpendicular to that line: a **transverse** stiffness, resisting
  the swinging of $j$ about $i$ at fixed separation. Its coefficient is
  $g_{ij}$, so it is present only for a spring that is not at its natural
  length.

The transverse term is the mechanical signature of a spring under tension. A
taut spring resists being deflected sideways, exactly as a taut string does; a
relaxed one, to this order, does not. This is the whole of the difference
between §3.1 and §3.2 below.

The Hessian is symmetric and positive semi-definite at a minimum, with six zero
eigenvalues corresponding to the rigid-body motions. Its pseudo-inverse
$\mathbf{K}^{+}$, obtained by spectral decomposition omitting those six modes,
gives the compliance.

## 3. Frustration

### 3.1 Relaxed networks

The simplest way to satisfy the equilibrium condition (5) is to make every term
vanish separately:

$$
l_{ij} = d^e_{ij} \quad \text{for all } i \sim j .
\tag{8}
$$

Every spring then sits at its natural length at equilibrium; no spring pulls or
pushes, and the network is **relaxed**, or non-frustrated. The choice
guarantees, by construction, that the reference conformation is the minimum, and
it makes $g_{ij} = 0$ throughout, so (7) reduces to the familiar rank-1 form

$$
\mathbf{K}_{i \neq j} = -k_{ij}\,\mathbf{e}_{ij}\mathbf{e}_{ij}^T .
\tag{9}
$$

Essentially all elastic network models in use are relaxed in this sense.

The minimum energy is then just the sum of the springs' own minima,

$$
V(\mathbf{r}^e) = \sum_{i<j} V^{\min}_{ij} .
$$

### 3.2 Frustrated networks

But (8) is a special solution, not the general one. Equation (5) is $3N-6$
independent equations in $M$ unknowns $l_{ij}$, and for a globular protein
$M \gg 3N$. The system is massively underdetermined: **any** set
$\{l_{ij}\}$ satisfying (5) yields an elastic network model with the same
equilibrium conformation $\mathbf{r}^e$, and is in that sense an equally valid
model of the same structure.

The general solution has the sum in (5) vanishing without its individual terms
vanishing. Some springs are then stretched ($d^e_{ij} > l_{ij}$) and others
compressed at equilibrium, their forces cancelling node by node. Such a network
is **frustrated**: it cannot simultaneously satisfy every spring's preferred
length, and settles into a compromise in which internal stresses balance.

The minimum energy acquires a second term,

$$
\boxed{\;
V(\mathbf{r}^e) = \sum_{i<j} V^{\min}_{ij} + \tfrac{1}{2}\sum_{i<j} k_{ij}\left(d^e_{ij} - l_{ij}\right)^2 ,
\;}
\tag{10}
$$

which is non-negative and vanishes only in the relaxed case. A frustrated
network therefore sits at a higher minimum energy than the relaxed network built
on the same structure: the residual stress is stored elastically.

Frustration is expected in real proteins. A protein's native structure is a
compromise among many competing interactions, not a simultaneous optimum of all
of them.

### 3.3 Frustration changes the dynamics

The consequence that matters is not energetic but dynamical.

Because $g_{ij} \neq 0$ for a frustrated network, the transverse term in (7)
survives, and the Hessian of a frustrated network **differs** from that of the
relaxed network built on the same conformation — not by an overall scale, but in
its structure, since the transverse contribution enters with a different
geometry from the longitudinal one. Two models can therefore share an
equilibrium structure exactly and still predict different fluctuations.

Since the normal modes are the eigenvectors of $\mathbf{K}$ and the mode
stiffnesses its eigenvalues, both depend on frustration, and the departure grows
with $|l_{ij}/d^e_{ij} - 1|$. Frustration is thus in principle observable in
dynamics — through fluctuation amplitudes, mode shapes and covariances — even
where it is invisible in the structure.

### 3.4 Choosing $l_{ij}$ for a frustrated model

Because (5) does not determine $\{l_{ij}\}$, a frustrated model requires further
input. Candidate strategies, from the 2017 appendix:

1. Keep $l_{ij} = d^e_{ij}$ for sequential and secondary-structure contacts, and
   allow $l_{ij} \neq d^e_{ij}$ only for tertiary contacts.
2. Let $l_{ij}$ depend only on the identities of the amino acids at $i$ and $j$,
   giving 210 independent parameters (20 amino acids, unordered pairs with
   repetition) to be fitted against the $3N-6$ equations (5) — over-determined,
   and fittable across many structures at once.
3. Impose a constraint such as a given mean energy and maximise the entropy over
   $\{l_{ij}\}$, obtaining a distribution rather than a point estimate.
4. Make all tertiary $l_{ij}$ equal and choose which contacts exist so that (5)
   holds.
5. Start relaxed, perturb all $l_{ij}$ randomly, then minimise
   $\sum_i \lVert\mathbf{F}_i\rVert$ to restore (5).
6. Generate frustration by simulated evolution: accept mutations under a
   stress-based criterion, so mutant structures stay near the wild-type while
   their $l_{ij}$ diverge.

Strategy 6 connects to §4: **mutation generates frustration automatically**. Even
a lineage that starts relaxed becomes frustrated as substitutions accumulate, so
a model restricted to relaxed networks cannot represent an evolving protein.

## 4. Mutation

### 4.1 The perturbation

A mutation at site $p$ is modelled as a perturbation of the springs that connect
$p$ to the rest of the network. In the **linearly forced** model, only the
natural lengths are perturbed:

$$
l_{pj} \to l_{pj} + \delta l_{pj}, \qquad j \sim p ,
\tag{11}
$$

with the $\delta l_{pj}$ drawn independently from a distribution with

$$
\langle \delta l_{pj}\rangle = 0, \qquad \operatorname{Var}(\delta l_{pj}) = \sigma^2 .
$$

The remaining parameters are held fixed, $\delta V^{\min}_{ij} = 0$ and
$\delta k_{ij} = 0$. In general a mutation would change all three; perturbing
only $l_{ij}$ is the model's defining simplification.

There is no amino-acid alphabet here. A mutation is a random draw of spring
perturbations, so the model describes the *statistics* of mutation at a site
rather than any particular substitution.

Substituting (11) into (1), the mutant potential is

$$
V_{\text{mut}}(\mathbf{r}) = \sum_{i<j}\left[V^{\min}_{ij} + \tfrac{1}{2}k_{ij}\left(d_{ij} - l_{ij} - \delta l_{ij}\right)^2\right] .
\tag{12}
$$

### 4.2 Stress energy and the induced force

Evaluate (12) at the *wild-type* equilibrium conformation — the mutant structure
before it has had a chance to relax. Expanding the square and using
$d^e_{ij} = l_{ij}$ for a relaxed wild type,

$$
V_{\text{mut}}(\mathbf{r}^e_{\text{wt}}) - V_{\text{wt}}(\mathbf{r}^e_{\text{wt}}) = \tfrac{1}{2}\sum_{i<j} k_{ij}\,\delta l_{ij}^2 \;\equiv\; \delta V^{\text{stress}} .
\tag{13}
$$

This **stress energy** is the energetic cost of imposing the mutation without
allowing the structure to respond. It is non-negative, quadratic in the
perturbation, and local: only springs at the mutated site contribute.

The mutation also exerts forces. From (4), the perturbed springs are no longer
at their natural lengths, so node $i$ feels

$$
f_{ij} = -k_{ij}\,\delta l_{ij}
\tag{14}
$$

directed along $\mathbf{e}_{ij}$, giving a net force vector $\mathbf{f}$ that
vanishes except at the mutated site and its neighbours. The mutation acts on the
network exactly like a set of external forces applied along the springs at site
$p$ — the sense in which the model is "linearly forced".

### 4.3 Structural response

The network relaxes under $\mathbf{f}$. To linear order, using the compliance
$\mathbf{C} = \mathbf{K}^{+}_{\text{wt}}$,

$$
\boxed{\;
\delta\mathbf{r} = \mathbf{r}^e_{\text{mut}} - \mathbf{r}^e_{\text{wt}} = \mathbf{C}\,\mathbf{f} .
\;}
\tag{15}
$$

This is the central result: the structural effect of a mutation is a linear
response, computed once the wild-type compliance is known. Because $\mathbf{C}$
is dense while $\mathbf{f}$ is local, the response extends across the whole
protein — the origin of long-range effects in this framework.

Relaxation lowers the energy from the stress value. The relaxed minimum is

$$
V_{\text{mut}}(\mathbf{r}^e_{\text{mut}}) = \delta V^{\text{stress}} - \tfrac{1}{2}\,\delta\mathbf{r}^T\mathbf{K}\,\delta\mathbf{r} ,
\tag{16}
$$

the second term being the relaxation energy recovered.

**Relaxation is partial, and this is the essential point.** The response (15)
minimises the *total* energy; it does not restore each individual spring to its
natural length, which is generally impossible. After relaxation
$d^e_{ij} \neq l_{ij}$ for the affected springs: the mutant is frustrated, by
construction and unavoidably, whatever the wild type was. This is the mechanism
behind strategy 6 of §3.4.

### 4.4 LFENM

The **linearly forced ENM** truncates at this point: it keeps the wild-type
Hessian for the mutant,

$$
V_{\text{mut}}(\mathbf{r}) = V_{\text{mut}}(\mathbf{r}^e_{\text{mut}}) + \tfrac{1}{2}\left(\mathbf{r}-\mathbf{r}^e_{\text{mut}}\right)^T \mathbf{K}_{\text{wt}} \left(\mathbf{r}-\mathbf{r}^e_{\text{mut}}\right) .
\tag{17}
$$

The justification is order counting: $\delta l$ is small, and the change in
$\mathbf{K}$ is of higher order than the change in the minimum and its location.

The consequence is sharp. The potential energy surface is **shifted**
vertically, by $V_{\text{mut}}(\mathbf{r}^e_{\text{mut}})$, and **translated**
laterally, by $\delta\mathbf{r}$ — but it is neither rotated nor deformed. The
mutant has

- the same normal modes (eigenvectors of $\mathbf{K}_{\text{wt}}$),
- the same mode stiffnesses (eigenvalues),
- and hence the same fluctuations, covariances and conformational entropy

as the wild type, exactly. In LFENM a mutation moves the protein without
changing how it moves.

### 4.5 SC-LFENM

The weakness of LFENM is one of internal consistency. An elastic network model
is *defined* by building $\mathbf{K}$ from the equilibrium conformation; but
LFENM's mutant has equilibrium conformation $\mathbf{r}^e_{\text{mut}}$ while
retaining $\mathbf{K}(\mathbf{r}^e_{\text{wt}})$. The mutant is therefore not an
ENM of its own structure:

$$
\mathbf{K}_{\text{wt}} \neq \mathbf{K}\left(\mathbf{r}^e_{\text{mut}}\right) .
$$

The **self-consistent LFENM** repairs this by recomputing the Hessian at the
relaxed mutant structure:

$$
\mathbf{K}_{\text{mut}} = \mathbf{K}\left(\mathbf{r}^e_{\text{mut}}\right),
$$

giving

$$
V_{\text{mut}}(\mathbf{r}) = V_{\text{mut}}(\mathbf{r}^e_{\text{mut}}) + \tfrac{1}{2}\left(\mathbf{r}-\mathbf{r}^e_{\text{mut}}\right)^T \mathbf{K}_{\text{mut}}\left(\mathbf{r}-\mathbf{r}^e_{\text{mut}}\right) .
\tag{18}
$$

Now the surface is rotated and deformed as well as shifted, and the mutant's
modes, stiffnesses and fluctuations all differ from the wild type's. Mutation
changes how the protein moves — which is what makes SC-LFENM, and not LFENM, the
model in which questions about mutational effects on dynamics can be posed at
all.

Recomputing $\mathbf{K}$ raises a question the wild-type theory never had to
face. The contact set $\mathcal{E}$ and the force constants $k_{ij}$ were both
fixed from a reference conformation (§1.2), and the mutant now has a different
one. Is the mutant an elastic network model of *its own* structure in these
respects too — $\mathcal{E}$ and $k_{ij}$ rebuilt at
$\mathbf{r}^e_{\text{mut}}$ — or does it inherit the wild type's network and
merely re-evaluate the Hessian on it?

The same self-consistency that motivates recomputing $\mathbf{K}$ argues for
rebuilding. What that costs depends on the network. For a complete graph
$\mathcal{E}$ cannot change, and rebuilding only adjusts the $k_{ij}$ smoothly.
For a proper subset it is more serious: a pair may cross whatever threshold
defines membership, so $\mathcal{E}_{\text{mut}} \neq \mathcal{E}_{\text{wt}}$
and the mutant is a different graph, not a re-parameterised one. Wild-type and
mutant sums then range over different index sets and are no longer comparable
term by term. The natural convention is that a contact present only in the
mutant arrives at its natural length, $\delta l_{ij} = 0$, and so carries no
stress, while one present only in the wild type simply has no mutant
counterpart. Whether to rebuild is left open (§8.8).

The transition may be read as a two-step process at fixed structure and then at
fixed parameters: a *vertical* step at $\mathbf{r}^e_{\text{wt}}$ replacing the
wild-type springs by the mutant's, costing $\delta V^{\text{stress}}$, followed
by *relaxation* to $\mathbf{r}^e_{\text{mut}}$, recovering the second term of
(16).

## 5. Statistical mechanics

At temperature $T$, with $\beta = 1/RT$, conformations follow

$$
\rho(\mathbf{r}) = \frac{e^{-\beta V(\mathbf{r})}}{Z}, \qquad Z = \int e^{-\beta V(\mathbf{r})}\,d\mathbf{r} .
$$

Within the quadratic approximation about the minimum this is a multivariate
normal,

$$
\rho(\mathbf{r}) = \frac{\exp\left[-\tfrac{1}{2}(\mathbf{r}-\mathbf{r}^e)^T\boldsymbol{\Sigma}^{-1}(\mathbf{r}-\mathbf{r}^e)\right]}{\sqrt{|2\pi\boldsymbol{\Sigma}|}}, \qquad
\boldsymbol{\Sigma} = \beta^{-1}\mathbf{K}^{+} .
\tag{19}
$$

The covariance is the compliance, scaled by temperature: soft directions
fluctuate widely, stiff ones narrowly. In the normal-mode basis
$\boldsymbol{\Sigma}$ is diagonal with
$\sigma_n^2 = (\beta\lambda_n)^{-1}$, $\lambda_n$ the $n$-th eigenvalue of
$\mathbf{K}$, and

$$
Z = e^{-\beta V(\mathbf{r}^e)}\prod_n \sqrt{\frac{2\pi}{\beta\lambda_n}} ,
\qquad
\ln Z = -\beta V(\mathbf{r}^e) + \sum_n \tfrac{1}{2}\ln\frac{2\pi}{\beta\lambda_n} .
\tag{20}
$$

Hence the thermodynamic functions

$$
U = V(\mathbf{r}^e) + \frac{3N-6}{2\beta}, \qquad
A = V(\mathbf{r}^e) - \sum_n \frac{1}{2\beta}\ln\frac{2\pi}{\beta\lambda_n}, \qquad
TS = \frac{1}{2\beta}\sum_n\left(\ln\frac{2\pi}{\beta\lambda_n} + 1\right) .
\tag{21}
$$

The free energy splits cleanly: an energetic part fixed by the minimum, and an
entropic part fixed by the spectrum. A mutation shifts the first through (16);
it shifts the second only if it changes the spectrum — that is, only under
SC-LFENM.

## 6. Comparing mutant with wild type

Write $\delta X = X_{\text{mut}} - X_{\text{wt}}$. Since $\mathbf{K}$ and
$\boldsymbol{\Sigma}$ are $3N\times3N$, the comparisons may be resolved either
**by site** — profiles over $i = 1\ldots N$, showing where in the protein the
effect falls — or **by mode** — profiles over $n$, showing which collective
motions are affected. The two are the same object in different bases, and
totals agree.

### 6.1 Structure

The displacement (15) resolved by site and by mode:

$$
\delta r^2_i = \lVert\delta\mathbf{r}_i\rVert^2, \qquad
\delta r^2_n = \left(\mathbf{u}_n^T \delta\mathbf{r}\right)^2, \qquad
\sum_i \delta r^2_i = \sum_n \delta r^2_n ,
\tag{22}
$$

with $\mathbf{u}_n$ the wild-type modes. Weighting by stiffness gives the
**deformation energy**, $\delta e^2 = \lambda_n \,\delta r^2_n$ by mode and
$\lVert\mathbf{K}^{1/2}\delta\mathbf{r}\rVert^2$ resolved by site; weighting
twice recovers the force, $\delta f^2_n = \lambda_n^2\,\delta r^2_n$. The three
measures answer different questions: how far a site moved, how much that motion
cost, and how hard it was pushed. Soft modes dominate $\delta r^2$; stiff modes
dominate $\delta f^2$.

### 6.2 Energy

From §4.2 and §5:

- $\delta V^{\text{stress}}$, equation (13): the cost before relaxation.
- $\delta V^{\min}$, from (16): the cost after relaxation, the difference of
  minimum energies. Approximately
  $\delta V^{\min} \simeq \delta V^{\text{stress}} - \tfrac{1}{2}\sum_i \delta e_i^2$.
- $\delta(TS)$, from (21): the entropic difference, non-zero only when the
  spectrum changes.

Stress energy is strictly local — only the mutated site's springs appear in
(13) — whereas the relaxation term is delocalised, since $\delta\mathbf{r}$
extends across the protein. The difference of two such terms is what produces
site-resolved energy profiles with long-range structure.

Site-resolved stress energy follows from splitting (13),

$$
\delta V^{\text{stress}} = \tfrac{1}{2}\sum_i \delta V^{\text{stress}}_i, \qquad
\delta V^{\text{stress}}_i = \tfrac{1}{2}\sum_{j\sim i} k_{ij}\,\delta l_{ij}^2 ,
\tag{23}
$$

the factors of $\tfrac{1}{2}$ preventing double counting of edges. Averaging
over mutations at a site, with $\langle\delta l^2\rangle = \sigma^2$, gives a
response matrix with entries $\tfrac{1}{2}k_{ij}\sigma^2$ for $i \sim j$ and
zero otherwise — symmetric, and computable without simulating any mutant.

### 6.3 Motion

Only SC-LFENM gives non-zero differences here. Comparisons of the fluctuation
ensembles:

- **Mean-square fluctuation**, $\delta\sigma^2_i$ by site or
  $\delta\sigma^2_n$ by mode: change in the amplitude of motion.
- **Mode overlap**, $(\mathbf{u}^{\text{wt}\,T}_m\mathbf{u}^{\text{mut}}_n)^2$:
  how far the modes have rotated. Summarised per mode by the participation
  number $n_H = e^{H}$, $H = -\sum_m p\ln p$ with $p$ the squared overlaps —
  equal to 1 for a mode that survives unchanged, larger when a mode is mixed
  across several.
- **Ensemble distances** between the covariance matrices, per site or overall:
  the Bhattacharyya distance, and the root weighted-square inner product
  (RWSIP), which weights mode overlaps by the fluctuation amplitudes so that
  soft, biologically relevant modes dominate.

### 6.4 Activation

If a subset of sites forms an **active site**, the effective stiffness for
deforming that subset alone is obtained by inverting the corresponding block of
the covariance,

$$
\mathbf{K}^{\text{act}} = \left[\boldsymbol{\Sigma}_{\text{act}}\right]^{-1} ,
\tag{24}
$$

which is not the corresponding block of $\mathbf{K}$: it accounts for the rest
of the protein relaxing in response. The energetic cost of forcing the active
site into a reference geometry $\mathbf{r}^{\ast}$ is
$\tfrac{1}{2}\,\delta\mathbf{r}_{\text{act}}^T\mathbf{K}^{\text{act}}\delta\mathbf{r}_{\text{act}}$,
with an accompanying entropic term from the spectrum of
$\mathbf{K}^{\text{act}}$. Differences between mutant and wild type give the
mutational effect on activation free energy.

## 7. Summary of the models

| | relaxed ENM | LFENM mutant | SC-LFENM mutant |
|---|---|---|---|
| natural lengths | $l_{ij} = d^e_{ij}$ | perturbed, eq. (11) | perturbed, eq. (11) |
| frustration | none | present after relaxation | present after relaxation |
| structure | reference | $\mathbf{r}^e_{\text{wt}} + \mathbf{C}\mathbf{f}$ | the same |
| network $\mathcal{E}$, $k_{ij}$ | from reference | inherited | rebuilt or inherited — open (§8.8) |
| Hessian | $\mathbf{K}(\mathbf{r}^e)$ | $\mathbf{K}_{\text{wt}}$, unchanged | $\mathbf{K}(\mathbf{r}^e_{\text{mut}})$ |
| modes, spectrum | — | unchanged | changed |
| $\delta V^{\min}$ | — | non-zero | non-zero |
| $\delta(TS)$ | — | zero | non-zero |

The two mutant models share a structural response and differ only in whether the
Hessian is rebuilt. That single difference decides whether the model can say
anything about mutational effects on dynamics.

Note the asymmetry in how frustration enters. A *relaxed* ENM by definition has
$g_{ij} = 0$. A *mutant* is frustrated by construction (§4.3). Whether the
transverse term (7) is retained when building $\mathbf{K}_{\text{mut}}$ is
therefore a real modelling decision, not an automatic consequence — see §8.2.

## 8. Open questions

Recorded without taking positions.

### 8.1 The sign of the transverse term

Two derivations of (7) give opposite signs for $g_{ij}$, namely
$1 - l_{ij}/d_{ij}$ and $l_{ij}/d_{ij} - 1$ (2017 document, where the second
appears in a comment beginning *"But beware of the sign"*). Writing
$s_{ij} = l_{ij} - d_{ij}$, one of the two gives
$\operatorname{Tr}\mathbf{K}_{ij} = -k_{ij}(1 - 2s_{ij})$, which turns positive
for $s_{ij} > 1/2$ — an apparently negative effective force constant, noted at
the time as feeling unphysical. Unresolved.

### 8.2 Is a frustrated mutant to be given a frustrated Hessian?

A mutant is frustrated after relaxation (§4.3). Two readings coexist:

- Frustration is accounted for in the *energy*, equation (10), and the mutant's
  Hessian is the ordinary relaxed-form (9) evaluated at the mutant structure.
  On this reading the transverse term is not used.
- Frustration is physical and belongs in the Hessian, since it changes the
  dynamics (§3.3), and omitting it discards exactly the effect of interest.

These give different mutant spectra. Which is intended for SC-LFENM has not been
settled.

### 8.3 SC-LFENM is not reversible

Starting from the mutant and applying the reverse perturbation does not recover
the wild-type potential. A thermodynamically consistent evolutionary model would
presumably require some form of reversibility, or an explicit account of why it
fails.

### 8.4 Negative site contributions to deformation energy

For a relaxed network any deformation increases the energy of every spring, so
site contributions are non-negative. Under frustration, relieving the stress of
one spring may increase another's, and site-resolved contributions may take
either sign even where the total is positive. Whether the quadratic form
$\tfrac{1}{2}\delta\mathbf{r}^T\mathbf{K}\delta\mathbf{r}$ correctly captures
these compensating changes is unclear; the approximation
$\delta V^{\min} \simeq \delta V^{\text{stress}} - \tfrac{1}{2}\sum_i\delta e_i^2$
of §6.2 is affected, since its left side may be negative while
$\delta e_i^2 \geq 0$ by construction.

### 8.5 Where is the minimum for a frustrated network?

For $l_{ij} = d^e_{ij}$ the minimum is known by construction. For general
$\{l_{ij}\}$ the minimum of
$V = \tfrac{1}{2}\sum k_{ij}(d_{ij}-l_{ij})^2$ must be found by solving the
non-linear system (5), and no closed form is known. This blocks analytical
averaging of $\delta V^{\min}$ over mutations, which is currently done either
numerically or through the approximation above.

### 8.6 Reference state for perturbation

When perturbations accumulate along a trajectory, it must be settled whether
each $\delta l_{ij}$ is measured relative to the immediate predecessor or to the
original wild type. The two differ once more than one substitution has occurred
at or near a site.

### 8.7 Mutations beyond natural lengths

The model perturbs $l_{ij}$ only, holding $k_{ij}$ and $V^{\min}_{ij}$ fixed. A
general perturbation would change all three. What such a generalisation would
add, beyond §4's linear-response structure, is open.

### 8.8 Does a mutant inherit the wild type's network?

The contact set and the force constants are fixed from a reference conformation
(§1.2). When SC-LFENM rebuilds the Hessian at $\mathbf{r}^e_{\text{mut}}$, it is
unsettled whether $\mathcal{E}$ and $k_{ij}$ should be rebuilt there as well, or
inherited from the wild type.

Rebuilding is the self-consistent choice: the mutant is then an elastic network
model of its own structure, which is the whole motivation for SC-LFENM.
Inheriting keeps wild type and mutant on a common graph, so that every
difference measured between them is a difference of the same quantities, at the
cost of that self-consistency.

The tension is sharp only when $\mathcal{E}$ is a proper subset, since then
rebuilding can add and remove contacts. On a complete graph the question reduces
to whether $k_{ij}$ is re-evaluated, with no change of index set. A related
question is whether, once contacts may come and go, the comparison of wild type
with mutant remains well posed for quantities defined as sums over contacts.
