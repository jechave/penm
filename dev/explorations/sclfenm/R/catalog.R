library(here)
## ============================================================================
## CATALOGUE MODEL  (Julian's fix for the trajectory problem)
##
## The trouble with re-referencing at every step: each state becomes the zero of
## its own energy, so every move is uphill and the walk is absorbing.
##
## Instead: ONE fixed reference, the founder (m = 0). For every site j and every
## state m = 1..M, catalogue the change relative to the founder:
##      X(0 -> j,m)     for X in {dV, dr, dTS, ...}
## A move between two states of the same site is then a DIFFERENCE:
##      X(j,m0 -> j,m1)  =  X(0 -> j,m1) - X(0 -> j,m0)
## which is antisymmetric by construction, so
##   * reversibility is exact: X(m1->m0) = -X(m0->m1)
##   * dV takes BOTH signs -- negative whenever m1 is a cheaper state than m0
##   * the state space is finite: (M+1)^N sequences, not an unbounded walk
##
## Notation: beta is Boltzmann's 1/kT. Selection strength is NU, a separate
## letter, because it is an evolutionary parameter and not a temperature.
## ============================================================================
source(here("R", "sclfenm.R"))

## ---- build the catalogue ---------------------------------------------------
## For each site j and state m = 1..M, apply the mutation to the FOUNDER and
## record the changes. m = 0 is the founder itself (all changes zero).
build_catalog <- function(wt, M = 8, sigma = 0.3, model = "sclfenm",
                          beta = 1, sites = NULL, verbose = TRUE) {
  N <- length(wt$r)/3
  if (is.null(sites)) sites <- seq_len(N)
  TS0 <- state_TS(wt, beta)
  dV <- dTS <- matrix(0, N, M + 1)          # column 1 is m = 0
  dr <- array(0, dim = c(3*N, N, M + 1))
  for (j in sites) {
    for (m in seq_len(M)) {
      mu <- mutate(wt, j, sigma = sigma, model = model)
      dV[j, m+1]  <- mu$dV
      dTS[j, m+1] <- state_TS(mu$state, beta) - TS0
      dr[, j, m+1] <- mu$dr
    }
    if (verbose && j %% 20 == 0) cat("  catalogued site", j, "\n")
  }
  list(dV = dV, dTS = dTS, dr = dr, M = M, N = N, sigma = sigma,
       model = model, beta = beta, wt = wt)
}

## ---- differences between two states of the same site ------------------------
## the whole point: everything is a difference of catalogued values
cat_dV  <- function(cg, j, m0, m1) cg$dV [j, m1+1] - cg$dV [j, m0+1]
cat_dTS <- function(cg, j, m0, m1) cg$dTS[j, m1+1] - cg$dTS[j, m0+1]
cat_dr  <- function(cg, j, m0, m1) cg$dr[, j, m1+1] - cg$dr[, j, m0+1]

## ---- total state of the protein --------------------------------------------
## a sequence is a vector s[1..N] of state indices; energy is additive over
## sites to the extent that mutations at different sites are independent
seq_energy <- function(cg, s) sum(cg$dV[cbind(seq_along(s), s + 1)])
seq_dTS    <- function(cg, s) sum(cg$dTS[cbind(seq_along(s), s + 1)])
seq_dr     <- function(cg, s) {
  d <- numeric(3*cg$N)
  for (j in seq_along(s)) if (s[j] > 0) d <- d + cg$dr[, j, s[j]+1]
  d
}
