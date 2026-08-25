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
state_from_prot <- function(prot, d_max = 10.5, k_model = "step", k = 1) {
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
klf_v <- function(st, R = st$R) {
  d <- pair_dist(R, st$pr)
  a <- st$k > 0
  sum(st$v0) + 0.5 * sum(st$k[a] * (d[a] - st$l[a])^2)
}

## Exact nonlinear force, -dV/dr. Only active pairs contribute.
klf_force <- function(st, R = st$R) {
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
klf_hessian <- function(st, R = st$R, frustrated = TRUE) {
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
## eigenvalues that appear off a stationary point (see klf_minimise).
klf_spectrum <- function(K, tol = 1e-8) {
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
klf_minimise <- function(st, tol = 1e-10, maxit = 200, frustrated = TRUE) {
  R <- st$R
  it <- 0L
  for (it in seq_len(maxit)) {
    f <- klf_force(st, R)
    if (sqrt(sum(f^2)) < tol) break
    s  <- klf_spectrum(klf_hessian(st, R, frustrated))
    Cm <- s$vector %*% ((1 / s$value) * t(s$vector))
    R  <- matrix(as.vector(R) + as.vector(Cm %*% f), nrow = 3)
  }
  st$R <- R
  attr(st, "fres") <- sqrt(sum(klf_force(st, R)^2))
  attr(st, "iter") <- it
  st
}

## ---------------------------------------------------------------- mutation --

## Perturb the rest lengths of one site's ACTIVE contacts, then (optionally)
## recompute k(l), then minimise exactly.
##
## Asymmetry worth stating plainly: we perturb the site's currently-active
## contacts, because an inactive pair is not "a contact of the site". But the
## perturbed l is kept for every pair, so a pair pushed below the cutoff
## activates. Contacts therefore break AND form, but only pairs that are
## already active can be selected for perturbation.
klf_mutate <- function(st, site, sigma = 0.3, k_update = TRUE,
                       sd_min = 1L, tol = 1e-10, maxit = 200) {
  pr  <- st$pr
  sel <- which((pr$i == site | pr$j == site) & st$k > 0 & st$sdij >= sd_min)
  if (!length(sel)) stop("site ", site, " has no active contacts to mutate")

  before_active <- st$k > 0
  st$l[sel] <- st$l[sel] + stats::rnorm(length(sel), 0, sigma)
  if (k_update) st <- refresh_k(st)
  after_active <- st$k > 0

  st <- klf_minimise(st, tol = tol, maxit = maxit)
  attr(st, "n_broken") <- sum(before_active & !after_active)
  attr(st, "n_formed") <- sum(!before_active & after_active)
  attr(st, "mut_site") <- site
  st
}

## ------------------------------------------------------------- observables --

## Covariance = pseudo-inverse of K, at kT = 1/beta. C = (1/beta) K^+.
klf_cmat <- function(K, beta = 1, tol = 1e-8) {
  s <- klf_spectrum(K, tol)
  (1 / beta) * (s$vector %*% ((1 / s$value) * t(s$vector)))
}

## Mean-square fluctuation per site: the 3x3 diagonal blocks of C.
klf_msf_site <- function(C) {
  N <- nrow(C) / 3
  vapply(seq_len(N), function(i) {
    b <- (3 * (i - 1) + 1):(3 * i)
    sum(diag(C[b, b, drop = FALSE]))
  }, numeric(1))
}

## TS from the spectrum. beta = 1 in ANM units is a CONVENTION, not a physical
## temperature -- k_ij = 1 is dimensionless here. State it wherever TS is
## reported; do not mix this with beta_boltzmann() in the same comparison.
klf_ts <- function(K, beta = 1, tol = 1e-8) {
  v <- klf_spectrum(K, tol)$value
  sum(0.5 / beta * (log(2 * pi / (beta * v)) + 1))
}

## Strain per active pair, and its energy. The measure of "how frustrated".
klf_strain <- function(st, R = st$R) {
  a <- which(st$k > 0)
  d <- pair_dist(R, st$pr)[a]
  s <- d - st$l[a]
  list(max_abs = max(abs(s)), rms = sqrt(mean(s^2)),
       energy = 0.5 * sum(st$k[a] * s^2))
}
