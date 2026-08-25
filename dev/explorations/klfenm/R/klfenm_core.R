## klfenm core: V(r; {l}) with k = k(l), over ALL pairs.
##
## The penm package is NOT modified by anything here. Where a penm function
## would have needed changing, it is reimplemented in this file instead.
##
## Model:
##   V(r; {l}) = sum_ij [ v0_ij + 1/2 k_ij(l_ij) (d_ij(r) - l_ij)^2 ]
##
## The parameters are the rest lengths {l_ij} over ALL N(N-1)/2 pairs. Everything
## else is derived from them: k_ij = k(l_ij), and hence the active contact set
## {ij : k_ij > 0}. A contact is made or broken by MUTATION (l crossing the
## cutoff), never by relaxation. That is the whole point of the model.

## ---------------------------------------------------------------- geometry --

## All unordered pairs (i < j), in a fixed canonical order used everywhere.
all_pairs <- function(nsites) {
  ij <- t(utils::combn(nsites, 2))
  list(i = ij[, 1], j = ij[, 2])
}

## Distances of every pair, given xyz as a 3 x N matrix.
pair_dist <- function(R, pr) {
  D <- R[, pr$j, drop = FALSE] - R[, pr$i, drop = FALSE]
  sqrt(colSums(D^2))
}

## Unit vectors along each pair, i -> j. Recomputed from R every time: a stored
## eij goes stale the moment the structure moves, and for a frustrated Hessian
## (which depends on d as well as e) a stale one corrupts K silently.
pair_eij <- function(R, pr, d = NULL) {
  D <- R[, pr$j, drop = FALSE] - R[, pr$i, drop = FALSE]
  if (is.null(d)) d <- sqrt(colSums(D^2))
  sweep(D, 2, d, "/")
}

## Sequence separation |i - j|, used to keep backbone neighbours always bonded.
pair_sdij <- function(pr) abs(pr$j - pr$i)

## ------------------------------------------------------------ k(l) variants --

## The ANM step: a spring exists iff its REST LENGTH is within the cutoff.
## Backbone neighbours (|i-j| == 1) are always bonded, matching penm's kij_anm.
k_step <- function(lij, sdij, d_max = 10.5, k = 1, ...) {
  ifelse(lij <= d_max | sdij == 1, k, 0)
}

## Smooth alternative (phase 5): a sigmoid of width `w` centred on the cutoff.
## As w -> 0 this tends to k_step. Backbone neighbours still forced to k.
k_smooth <- function(lij, sdij, d_max = 10.5, k = 1, w = 0.5, ...) {
  kk <- k / (1 + exp((lij - d_max) / w))
  ifelse(sdij == 1, k, kk)
}

k_fun_of <- function(name) switch(name,
  step   = k_step,
  smooth = k_smooth,
  stop("unknown k model: ", name))

## ------------------------------------------------------------------- state --

## A klfenm state. `l` and `k` are over ALL pairs, in `pr` order.
##  l   : rest lengths, the only free parameters
##  k   : k(l), derived -- never an independent truth
##  R   : 3 x N coordinates, the equilibrium structure of THIS potential
##  v0  : per-pair constant (the V0 / Vmin_ij slot); scalar 0 by default
new_state <- function(R, l, k, pr, sdij, v0 = 0, k_model = "step",
                      k_par = list(d_max = 10.5)) {
  structure(list(R = R, l = l, k = k, pr = pr, sdij = sdij, v0 = v0,
                 k_model = k_model, k_par = k_par),
            class = "klfenm_state")
}

## Recompute k from l. This is the ONLY place k is allowed to be set from l,
## so `k = k(l)` cannot silently drift out of sync (see checks/, item 4).
refresh_k <- function(st) {
  kf <- k_fun_of(st$k_model)
  st$k <- do.call(kf, c(list(lij = st$l, sdij = st$sdij), st$k_par))
  st
}

n_active <- function(st) sum(st$k > 0)

## Build the wild-type state from a penm prot: l_ij = d_ij(r_wt) for ALL pairs.
## The active set is then exactly the standard ANM contact set, and every
## non-contact carries a latent rest length that a mutation can pull inside the
## cutoff -- which is what makes contact FORMATION possible at all.
klfenm_set_state <- function(prot, d_max = 10.5, k_model = "step", k = 1) {
  R  <- matrix(as.vector(penm::get_xyz(prot)), nrow = 3)
  N  <- ncol(R)
  pr <- all_pairs(N)
  sd <- pair_sdij(pr)
  l  <- pair_dist(R, pr)
  st <- new_state(R, l, k = NULL, pr = pr, sdij = sd, v0 = 0,
                  k_model = k_model, k_par = list(d_max = d_max, k = k))
  refresh_k(st)
}

## --------------------------------------------------------- energy & forces --

## V at an ARBITRARY structure R, for the state's parameters. Note this takes R
## as an argument: the whole discipline of this exploration is that we minimise
## V rather than estimate its minimum, so V must be evaluable off-minimum.
klfenm_energy <- function(st, R = st$R) {
  d <- pair_dist(R, st$pr)
  a <- st$k > 0
  sum(st$v0) + 0.5 * sum(st$k[a] * (d[a] - st$l[a])^2)
}

## -dV/dr at a given structure: the Newton residual, zero at a minimum.
## NOT the LFENM perturbing force f_ij = -k_ij dl_ij -- a different quantity,
## the one penm's calculate_force() computes. This model never forms it,
## because no linear-response step is ever taken.
energy_gradient <- function(st, R = st$R) {
  pr <- st$pr; N <- ncol(R)
  a  <- which(st$k > 0)
  D  <- R[, pr$j[a], drop = FALSE] - R[, pr$i[a], drop = FALSE]
  d  <- sqrt(colSums(D^2))
  E  <- sweep(D, 2, d, "/")
  co <- st$k[a] * (d - st$l[a])            # +ve when stretched
  f  <- matrix(0, 3, N)
  ii <- pr$i[a]; jj <- pr$j[a]
  for (q in seq_along(a)) {
    f[, ii[q]] <- f[, ii[q]] + co[q] * E[, q]
    f[, jj[q]] <- f[, jj[q]] - co[q] * E[, q]
  }
  as.vector(f)
}

## Hessian. frustrated = TRUE keeps the transverse term g = l/d - 1, which is
## nonzero exactly when a spring is strained. penm's set_enm() blocks this
## branch, so it is reimplemented here rather than by touching the package.
klfenm_kmat <- function(st, R = st$R, frustrated = TRUE) {
  pr <- st$pr; N <- ncol(R)
  a  <- which(st$k > 0)
  D  <- R[, pr$j[a], drop = FALSE] - R[, pr$i[a], drop = FALSE]
  d  <- sqrt(colSums(D^2))
  E  <- sweep(D, 2, d, "/")
  I3 <- diag(3)
  K  <- matrix(0, 3 * N, 3 * N)
  ii <- pr$i[a]; jj <- pr$j[a]; kk <- st$k[a]; ll <- st$l[a]
  for (q in seq_along(a)) {
    e   <- E[, q]
    ee  <- tcrossprod(e, e)
    g   <- if (frustrated) ll[q] / d[q] - 1 else 0
    Kij <- -kk[q] * (ee + g * (ee - I3))
    bi  <- (3 * (ii[q] - 2) + 4):(3 * ii[q])
    bj  <- (3 * (jj[q] - 2) + 4):(3 * jj[q])
    K[bi, bj] <- K[bi, bj] + Kij
    K[bj, bi] <- K[bj, bi] + t(Kij)
    K[bi, bi] <- K[bi, bi] - Kij
    K[bj, bj] <- K[bj, bj] - Kij
  }
  K
}

## Spectrum with the 6 rigid-body modes dropped by TOLERANCE on the eigenvalue.
## Returns the raw sorted spectrum too, so a caller can check for the negative
## eigenvalues that appear off a stationary point (see klfenm_minimise).
klfenm_nma <- function(K, tol = 1e-8) {
  e   <- eigen(K, symmetric = TRUE)
  o   <- order(e$values)
  val <- e$values[o]; vec <- e$vectors[, o, drop = FALSE]
  keep <- val > tol * max(abs(val))
  list(raw_value = val, value = val[keep], vector = vec[, keep, drop = FALSE],
       n_zero = sum(!keep))
}

## ---------------------------------------------------------------- minimise --

## Newton iteration to the TRUE minimum of the state's own potential.
##
## This is the heart of the exploration: no linear response, no C_wt f. The
## structure is found by minimising V, and only then is anything expanded.
## Uses the FRUSTRATED Hessian, because a strained network's curvature is the
## frustrated one -- using g = 0 here would be solving a different problem.
## NON-CONVERGENCE IS AN ERROR, not a value. This is the part that matters:
## the first version returned whatever structure the loop happened to reach, so
## a failed minimisation entered a trajectory as a real energy. Measured: one
## draw in 200 hit maxit with |F| = 2.4e3 and yielded dV = 1.6e5.
##
## `damp` backtracks the Newton step (halve until the energy decreases). It is
## kept because it is cheap and standard, but measured over 60 hard starts
## (sigma = 2.0) it changed nothing: 1 failure with damping, 1 without, same
## median iteration count. So it is not what fixed the bug above -- the error
## on non-convergence is.
klfenm_minimise <- function(st, tol = 1e-8, maxit = 200, frustrated = TRUE,
                            damp = TRUE) {
  R  <- st$R
  it <- 0L
  v  <- klfenm_energy(st, R)
  for (it in seq_len(maxit)) {
    f <- energy_gradient(st, R)
    if (sqrt(sum(f^2)) < tol) break
    s  <- klfenm_nma(klfenm_kmat(st, R, frustrated))
    Cm <- s$vector %*% ((1 / s$value) * t(s$vector))
    step <- as.vector(Cm %*% f)
    if (damp) {
      a <- 1
      repeat {
        Rn <- matrix(as.vector(R) + a * step, nrow = 3)
        vn <- klfenm_energy(st, Rn)
        if (vn <= v || a < 1e-6) break
        a <- a / 2
      }
      R <- Rn; v <- vn
    } else {
      R <- matrix(as.vector(R) + step, nrow = 3)
      v <- klfenm_energy(st, R)
    }
  }
  fres <- sqrt(sum(energy_gradient(st, R)^2))
  if (fres >= tol)
    stop("klfenm_minimise(): did not converge in ", maxit,
         " iterations (|F| = ", signif(fres, 4), ", tol = ", tol, ")")
  st$R <- R
  attr(st, "fres") <- fres
  attr(st, "iter") <- it
  st
}

## ---------------------------------------------------------------- mutation --

## Perturb the rest lengths of ALL of one site's pairs -- active or not -- then
## (optionally) recompute k(l), then minimise exactly.
##
## Every pair carries an l_ij and k(l) alone decides which are "on", so breaking
## and forming happen by the same rule. Restricting the perturbation to
## currently-active pairs would break that symmetry: the network could then only
## ever thin, because a switched-off pair would be frozen out of the mutational
## process for good. (Measured when it was: 14 broken, 0 formed in 25 mutations.)
##
## `radius` optionally bounds which of the site's pairs are perturbed. The
## DEFAULT IS Inf -- every pair of the site, so one mutation can reach across the
## whole protein. A finite radius keeps it local; nothing in the report uses one
## except one check, which says so.
klfenm_mutate_site <- function(st, site, sigma = 0.3, k_update = TRUE,
                       sd_min = 1L, radius = Inf, tol = 1e-8, maxit = 200) {
  pr  <- st$pr
  own <- (pr$i == site | pr$j == site) & st$sdij >= sd_min
  if (is.finite(radius)) {
    d   <- pair_dist(st$R, pr)
    own <- own & (d <= radius | st$k > 0)   # keep every active pair in play
  }
  sel <- which(own)
  if (!length(sel)) stop("site ", site, " has no pairs to mutate")

  before_active <- st$k > 0
  st$l[sel] <- st$l[sel] + stats::rnorm(length(sel), 0, sigma)
  if (k_update) st <- refresh_k(st)
  after_active <- st$k > 0

  st <- klfenm_minimise(st, tol = tol, maxit = maxit)
  attr(st, "n_broken") <- sum(before_active & !after_active)
  attr(st, "n_formed") <- sum(!before_active & after_active)
  attr(st, "mut_site") <- site
  st
}

## ------------------------------------------------------------- observables --

## Covariance = pseudo-inverse of K, at kT = 1/beta. C = (1/beta) K^+.
klfenm_cmat <- function(K, beta = 1, tol = 1e-8) {
  s <- klfenm_nma(K, tol)
  (1 / beta) * (s$vector %*% ((1 / s$value) * t(s$vector)))
}

## Mean-square fluctuation per site: the 3x3 diagonal blocks of C.
klfenm_msf_site <- function(C) {
  N <- nrow(C) / 3
  vapply(seq_len(N), function(i) {
    b <- (3 * (i - 1) + 1):(3 * i)
    sum(diag(C[b, b, drop = FALSE]))
  }, numeric(1))
}

## TS from the spectrum. beta = 1 in ANM units is a CONVENTION, not a physical
## temperature -- k_ij = 1 is dimensionless here. State it wherever TS is
## reported; do not mix this with beta_boltzmann() in the same comparison.
klfenm_entropy <- function(K, beta = 1, tol = 1e-8) {
  v <- klfenm_nma(K, tol)$value
  sum(0.5 / beta * (log(2 * pi / (beta * v)) + 1))
}

## ------------------------------------------------- energy differences (dV) --
##
## NOTATION. Both pieces below are DIFFERENCES, and the reference is always
## named. Writing them as "V_stress" and "V_relax" hides that they are deltas,
## and hides it in the one place it matters: along a trajectory the reference
## is an evolved, STRAINED state, so V_ref(r_ref) != 0 and
##
##     dV_stress(ref -> mut)  !=  V_mut(r_ref)
##
## They coincide only when the reference is relaxed (the founder). Anywhere else
## the difference is the reference's own strain energy, which grows as the walk
## proceeds. That is exactly the V-vs-dV confusion this notation exists to stop.
##
## All are computed by EVALUATING HAMILTONIANS -- no expansion, no closed form.
## A closed form would need the cross term (nonzero at a strained reference) and
## extra terms from k^mut != k^ref, and would be wrong the moment either is
## forgotten. Evaluating is exact at any reference and is simpler.
##
##     dV_min(ref->mut) = dV_stress(ref->mut) + dV_relax(ref->mut)
##
## dV_min is the difference of MINIMA.

## dV_stress(ref -> mut) = V_mut(r_ref) - V_ref(r_ref)
## Two Hamiltonians, ONE structure: the cost of changing the parameters before
## the structure is allowed to respond. Exact, whatever state `ref` is in.
klfenm_delta_v_stress <- function(ref, mut) klfenm_energy(mut, R = ref$R) - klfenm_energy(ref, R = ref$R)

## dV_relax(ref -> mut) = V_mut(r^e_mut) - V_mut(r_ref)
## ONE Hamiltonian, two structures: what relaxation gives back. <= 0 because
## r^e_mut minimises V_mut.
klfenm_delta_v_relax <- function(ref, mut)
  klfenm_energy(mut, R = mut$R) - klfenm_energy(mut, R = ref$R)

## dV_min(ref -> mut) = V_mut(r^e_mut) - V_ref(r^e_ref)
##
## THE difference of minima -- the physical energy change of the substitution,
## and the quantity a trajectory accepts or rejects on. Identically
## klfenm_delta_v_stress + klfenm_delta_v_relax, which the checks assert.
##
## Every state comes from klfenm_mutate_site(), which minimises before it
## returns, so both are at their own minimum by construction.
klfenm_delta_v_min <- function(ref, mut)
  klfenm_energy(mut, R = mut$R) - klfenm_energy(ref, R = ref$R)

## Strain per active pair, and its energy. The measure of "how frustrated".
klfenm_strain <- function(st, R = st$R) {
  a <- which(st$k > 0)
  d <- pair_dist(R, st$pr)[a]
  s <- d - st$l[a]
  list(max_abs = max(abs(s)), rms = sqrt(mean(s^2)),
       energy = 0.5 * sum(st$k[a] * s^2))
}
