## ============================================================================
## Catalogue model, built on penm's OWN functions (no reimplementation).
##
## Model: ANM, uniform k, contacts at d < 10.5 A, all contacts of the mutated
## site perturbed (mut_sd_min = 1).
##
## penm already linearises: calculate_force() gives f_ij = -k_ij dl_ij directly
## from dl, and calculate_dxyz() gives dr = C_wt f. So dr is linear in f, hence
## exactly additive and exactly reversible along a trajectory.
##
## The catalogue: for site j and state m = 0..M, store the change relative to
## the FOUNDER. m = 0 is "don't mutate" -- penm's own convention.
## Moves are differences:  X(m0 -> m1) = X(0 -> m1) - X(0 -> m0)
##
## beta = Boltzmann 1/kT.   NU = selection strength (a separate letter).
## ============================================================================
## find the package root by walking up, so this survives being moved
.penm_root <- local({ d <- normalizePath(".")
  while (!file.exists(file.path(d, "DESCRIPTION")) && dirname(d) != d) d <- dirname(d); d })
suppressPackageStartupMessages(devtools::load_all(.penm_root, quiet = TRUE))

founder <- function(d_max = 10.5)
  set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = d_max, frustrated = FALSE)

build_catalog <- function(wt, M = 10, sigma = 0.3, ensemble = 1L, verbose = TRUE) {
  N <- get_nsites(wt); ne <- nrow(get_graph(wt))
  dl <- array(0, dim = c(ne,  N, M + 1))
  f  <- array(0, dim = c(3*N, N, M + 1))
  dr <- array(0, dim = c(3*N, N, M + 1))
  dV <- matrix(0, N, M + 1)
  for (j in seq_len(N)) {
    for (m in seq_len(M)) {
      mut <- get_mutant_site(wt, site_mut = j, mutation = m, mut_model = "lfenm",
                             mut_dl_sigma = sigma, mut_sd_min = 1L,
                             ensemble = ensemble)
      dl[, j, m+1] <- get_graph(mut)$lij - get_graph(wt)$lij
      f [, j, m+1] <- calculate_force(wt, dl[, j, m+1])
      dr[, j, m+1] <- as.vector(get_xyz(mut) - get_xyz(wt))
      dV[  j, m+1] <- ddg_dv(wt, mut)
    }
    if (verbose && j %% 20 == 0) cat("  site", j, "\n")
  }
  list(dl=dl, f=f, dr=dr, dV=dV, M=M, N=N, ne=ne, sigma=sigma,
       ensemble=ensemble, wt=wt)
}

## a sequence = one state per site
seq_force <- function(cg, s) { tot <- numeric(3*cg$N)
  for (j in seq_along(s)) if (s[j] > 0) tot <- tot + cg$f[, j, s[j]+1]; tot }
seq_dl <- function(cg, s) { tot <- numeric(cg$ne)
  for (j in seq_along(s)) if (s[j] > 0) tot <- tot + cg$dl[, j, s[j]+1]; tot }
seq_dr <- function(cg, s) as.vector(calculate_dxyz(cg$wt, seq_force(cg, s)))

seq_dV_exact <- function(cg, s) {
  wt <- cg$wt; mut <- wt
  mut$graph$lij <- wt$graph$lij + seq_dl(cg, s)
  mut$nodes$xyz <- wt$nodes$xyz + seq_dr(cg, s)
  mut$graph$dij <- dij_edge(mut$nodes$xyz, mut$graph$i, mut$graph$j)
  ddg_dv(wt, mut)
}
seq_dV_catalog <- function(cg, s) sum(cg$dV[cbind(seq_along(s), s + 1)])

## ---- evolutionary trajectory on the catalogue ------------------------------
## A sequence s[1..N] of states. One step: pick a site, propose a new state,
## accept with a Metropolis rule on the CATALOGUE energy difference.
##   nu = selection strength (NOT a temperature; beta is reserved for 1/kT)
run_catalog_trajectory <- function(cg, nstep = 500, nu = 1, seed = NULL,
                                   s0 = NULL, record_every = 1) {
  if (!is.null(seed)) set.seed(seed)
  N <- cg$N; M <- cg$M
  s <- if (is.null(s0)) integer(N) else s0
  V <- seq_dV_catalog(cg, s)
  rec <- data.frame()
  n_acc <- 0L
  for (t in seq_len(nstep)) {
    j  <- sample(N, 1)
    m1 <- sample(setdiff(0:M, s[j]), 1)
    dV <- cg$dV[j, m1+1] - cg$dV[j, s[j]+1]      # the catalogue difference
    if (dV <= 0 || runif(1) < exp(-nu * dV)) {   # Metropolis
      s[j] <- m1; V <- V + dV; n_acc <- n_acc + 1L
    }
    if (t %% record_every == 0)
      rec <- rbind(rec, data.frame(step = t, V = V, n_mut = sum(s > 0),
                                   accept = n_acc/t))
  }
  list(s = s, record = rec, n_acc = n_acc, nu = nu)
}
