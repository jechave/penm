## ============================================================================
## SC-LFENM with energy bookkeeping
##
## Design (see MODEL.md):
##   An ENM is a FIT to a structure: given r, return a relaxed network (l = d).
##   A mutation is therefore handled in two separate stages:
##
##     (1) ENERGETICS, on the pre-rebuild (frustrated) Hamiltonian:
##           dV = V_stress - V_relax,  exact to O(dl^3)
##     (2) DYNAMICS, on a network REBUILT (relaxed) at the mutant structure.
##
##   The rebuild erases the strain, so dV is accumulated in a saved scalar
##   offset `v_off`. The offset is additive, never enters forces or K, and
##   makes the energy accounting survive the rebuild.
## ============================================================================

source("enm_core.R")

## ---- a protein state -------------------------------------------------------
## r      3N structure
## pr,k,l network (contact pairs, force constants, natural lengths)
## v_off  accumulated energy relative to the founder (0 for the founder)
## frustrated  whether K keeps the transverse term
new_state <- function(r, d_max = 10.5, k_val = 1, frustrated = FALSE) {
  st <- build_relaxed(r, d_max, k_val)
  st$v_off <- 0; st$d_max <- d_max; st$k_val <- k_val
  st$frustrated <- frustrated
  st$nmut <- 0L
  st
}

## ---- derived quantities ----------------------------------------------------
state_K <- function(st) {
  if (!is.null(st$K_frozen)) return(st$K_frozen)   # lfenm: spectrum frozen at wt
  hessian(st$r, st$pr, st$k, st$l, st$frustrated)
}
state_nma <- function(st) nma(state_K(st))

## ---- NOTATION -------------------------------------------------------------
## r      the state's EQUILIBRIUM conformation r^e (see caveat below)
## d_ij   distance in an ARBITRARY conformation r  -> use in V(r)
## d^e_ij distance at r^e                          -> use in V_min
## V(r)       = sum_ij [ v0_ij + 1/2 k_ij (d_ij - l_ij)^2 ]   a FUNCTION of r
## V_min      = V(r^e)                                        a NUMBER
##
## The total energy of a state is a MINIMUM energy, V_min, not V at some
## arbitrary conformation:
##      V_min(state) = 1/2 sum k_ij (d^e_ij - l_ij)^2  +  v_off
##
## CAVEAT. st$r is r^e exactly only when the network was rebuilt on it
## (model = "sclfenm"), where l = d^e by construction. For "lfenm"/"keepnet"
## st$r is the LINEAR-RESPONSE estimate of r^e and carries a residual force
## O(dl^2); V evaluated there exceeds the true V_min by <= 0.01% in the runs
## here. Use state_vmin(..., exact = TRUE) to relax to the true minimum first.

## minimum energy of the state: V(r^e) + offset. THIS is "the energy of a state".
state_vmin <- function(st, exact = FALSE, tol = 1e-10, maxit = 60) {
  r <- st$r
  if (exact) r <- state_req(st, tol, maxit)$r
  venm(r, st$pr, st$k, st$l) + st$v_off
}

## the state's equilibrium conformation, found by iterating the linear response
## until the net force vanishes. For "sclfenm" states this returns st$r at once.
state_req <- function(st, tol = 1e-10, maxit = 60) {
  r <- st$r
  for (i in seq_len(maxit)) {
    f <- force(r, st$pr, st$k, st$l)
    if (sqrt(sum(f^2)) < tol) break
    K <- hessian(r, st$pr, st$k, st$l, st$frustrated)
    r <- r + as.vector(cmat_of(nma(K)) %*% f)
  }
  list(r = r, iter = i, fres = sqrt(sum(force(r, st$pr, st$k, st$l)^2)))
}



## entropic free energy, TS, from the spectrum (up to an additive constant)
state_TS <- function(st, beta = 1) {
  nm <- state_nma(st)
  sum(0.5/beta * (log(2*pi/(beta*nm$evalue)) + 1))
}

## ---- one mutation ----------------------------------------------------------
## Returns the proposed mutant WITHOUT deciding acceptance.
##
## THREE models, not two. They differ only in what happens to K:
##   model = "lfenm"    K is FROZEN at the wild type. Modes never change.
##                      (K_mut = K_wt exactly; the textbook LFENM.)
##   model = "keepnet"  contact list kept, l accumulates, but K is recomputed
##                      from the moved structure. Modes change through geometry
##                      only. (This is what "don't rebuild" naively gives.)
##   model = "sclfenm"  network REFIT (relaxed) at the mutant structure.
##                      Modes change; strain is erased and carried in v_off.
mutate <- function(st, site, sigma = 0.3, model = c("sclfenm","lfenm","keepnet"),
                   sd_min = 1L) {
  model <- match.arg(model)
  ## edges of the mutated site, optionally excluding near-sequence neighbours
  mask <- (st$pr[,1] == site | st$pr[,2] == site)
  if (sd_min > 1L) mask <- mask & (abs(st$pr[,2] - st$pr[,1]) >= sd_min)
  if (!any(mask)) stop("site has no mutable contacts")

  dl <- rep(0, nrow(st$pr)); dl[mask] <- rnorm(sum(mask), 0, sigma)
  l2 <- st$l + dl

  ## --- stage 1: energetics on the pre-rebuild Hamiltonian -------------------
  K  <- state_K(st)
  nm <- nma(K)
  C  <- cmat_of(nm)
  f  <- force(st$r, st$pr, st$k, l2)     # force induced by the perturbation
  dr <- as.vector(C %*% f)
  r2 <- st$r + dr

  V_before <- venm(st$r, st$pr, st$k, st$l)
  V_after  <- venm(r2,    st$pr, st$k, l2)      # exact on the frustrated H
  dV       <- V_after - V_before

  ## decomposition (equal to dV up to O(dl^3)); reported for diagnostics
  V_stress <- venm(st$r, st$pr, st$k, l2) - V_before
  V_relax  <- 0.5 * as.numeric(crossprod(dr, K %*% dr))

  ## --- stage 2: dynamics -----------------------------------------------------
  if (model == "sclfenm") {
    new <- build_relaxed(r2, st$d_max, st$k_val)   # relaxed refit: strain erased
    new$v_off <- st$v_off + dV                     # ...so carry it in the offset
    new$K_frozen <- NULL
  } else {
    new <- list(r = r2, pr = st$pr, k = st$k, l = l2)
    new$v_off <- st$v_off                          # strain still in the network
    ## lfenm freezes the spectrum: carry the wt Hessian forward unchanged
    new$K_frozen <- if (model == "lfenm") (if (!is.null(st$K_frozen)) st$K_frozen else K) else NULL
  }
  new$d_max <- st$d_max; new$k_val <- st$k_val
  new$frustrated <- st$frustrated; new$nmut <- st$nmut + 1L
  new$model <- model

  list(state = new, dV = dV, V_stress = V_stress, V_relax = V_relax,
       dr = dr, site = site, dl = dl,
       n_edge_before = nrow(st$pr), n_edge_after = nrow(new$pr))
}

## ---- comparisons wt vs mut -------------------------------------------------
## motion: needs both spectra
delta_motion <- function(wt, mut, beta = 1, nsub = 10) {
  a <- state_nma(wt); b <- state_nma(mut)
  Ca <- cmat_of(a); Cb <- cmat_of(b)
  N <- length(wt$r)/3
  msf <- function(C) { d <- diag(C); colSums(matrix(d, nrow = 3)) }
  ## RMSIP must be taken over a SUBSET of modes. Over two complete orthonormal
  ## bases of the same subspace sum(ov^2) == ncol(ov) identically, so RMSIP = 1
  ## whatever the spectra -- it would measure nothing. Use the softest nsub modes.
  nsub <- min(nsub, ncol(a$umat), ncol(b$umat))
  soft_a <- a$umat[, order(a$evalue)[seq_len(nsub)], drop = FALSE]
  soft_b <- b$umat[, order(b$evalue)[seq_len(nsub)], drop = FALSE]
  ov <- crossprod(soft_a, soft_b)
  list(dmsf_site = msf(Cb) - msf(Ca),
       dTS = state_TS(mut, beta) - state_TS(wt, beta),
       rmsip = sqrt(sum(ov^2)/nsub), nsub = nsub,
       n_modes = c(wt = length(a$evalue), mut = length(b$evalue)))
}
