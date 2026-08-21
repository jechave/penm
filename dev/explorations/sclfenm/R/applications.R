library(here)
## Two applications of the sclfenm model. source() after sclfenm.R
source(here("R", "sclfenm.R"))

## ============================================================================
## A. SINGLE-MUTATION SCANNING
## For each site, average nmut mutations; report energy, structure, motion.
## The wild type is fixed: every mutant is measured against the SAME founder.
## ============================================================================
scan_sites <- function(wt, nmut = 10, sigma = 0.3,
                       model = "sclfenm", sites = NULL, beta = 1, verbose = TRUE) {
  N <- length(wt$r)/3
  if (is.null(sites)) sites <- seq_len(N)
  ## wt quantities computed once
  nmw <- state_nma(wt); TS_wt <- state_TS(wt, beta)
  out <- data.frame()
  for (s in sites) {
    ## TWO different displacement measures, easily confused:
    ##   dr2_self  how far the MUTATED SITE itself moves  -> "buried sites move less"
    ##   dr2_total how much the WHOLE PROTEIN moves       -> influence of site s
    dV <- dr2_tot <- dr2_self <- dTS <- numeric(nmut); ne <- integer(nmut)
    for (m in seq_len(nmut)) {
      mu <- mutate(wt, s, sigma = sigma, model = model)
      dV[m]       <- mu$dV
      dr2_tot[m]  <- sum(mu$dr^2)
      dr2_self[m] <- sum(matrix(mu$dr, nrow = 3)[, s]^2)
      dTS[m]      <- state_TS(mu$state, beta) - TS_wt
      ne[m]       <- mu$n_edge_after
    }
    out <- rbind(out, data.frame(site = s,
      dV = mean(dV), dV_sd = sd(dV),
      dr2_self = mean(dr2_self), dr2_total = mean(dr2_tot),
      rmsd = sqrt(mean(dr2_tot)/N),
      dTS = mean(dTS), dTS_sd = sd(dTS),
      d_edge = mean(ne) - nrow(wt$pr),
      cn = sum(wt$pr[,1]==s | wt$pr[,2]==s)))
    if (verbose && s %% 10 == 0) cat("  site", s, "\n")
  }
  out
}

## ============================================================================
## B. EVOLUTIONARY TRAJECTORY
## Metropolis on the total energy. Each accepted mutation becomes the new
## wild type; v_off accumulates, so dV is always measured from the founder.
## NOTE ON NOTATION: beta is Boltzmann's 1/kT and is used ONLY for thermodynamic
## quantities (state_TS). Selection strength is NU -- a separate letter, because
## it is an evolutionary parameter and calling it beta invites the reader to
## think of it as a temperature.
##
## selection:
##   "neutral"    accept everything. Energy runs away (drift).
##   "step"       accept on the STEP energy, min(1, exp(-beta*dV)).
##                Under sclfenm dV > 0 always (see MODEL.md), so this only
##                slows the rise; it does NOT bound it. Kept for comparison.
##   "threshold"  accept iff the TOTAL energy from the founder stays below
##                v_cut. This is the one that gives a stationary walk: it is a
##                hard constraint on cumulative destabilisation, which is what
##                a stability criterion actually means when every step costs
##                energy.
## ============================================================================
## exact = TRUE relaxes each proposed state to its true equilibrium conformation
## before reading its energy. Only matters for the non-rebuilt models, where the
## stored structure is the linear-response estimate of r^e (error ~1e-3 in V,
## measured); sclfenm states are already exact, so it costs nothing there.
run_trajectory <- function(wt, nstep = 100, sigma = 0.3, beta = 1, nu = 1, v_cut = Inf,
                           model = "sclfenm",
                           selection = c("threshold","step","neutral"),
                           seed = NULL, verbose = FALSE, max_try = 2000,
                           exact = TRUE) {
  selection <- match.arg(selection)
  if (!is.null(seed)) set.seed(seed)
  N <- length(wt$r)/3
  st <- wt; r_founder <- wt$r
  TS0 <- state_TS(wt, beta)
  rec <- data.frame()
  n_acc <- 0L; n_try <- 0L; stalled <- FALSE
  for (s in seq_len(nstep)) {
    n_try0 <- n_try
    repeat {                                  # try until one is accepted
      n_try <- n_try + 1L
      site <- sample(N, 1)
      mu <- mutate(st, site, sigma = sigma, model = model)
      acc <- switch(selection,
        neutral   = TRUE,
        step      = runif(1) < min(1, exp(-nu * mu$dV)),
        threshold = state_vmin(mu$state, exact = exact) < v_cut)
      if (acc) break
      if (n_try - n_try0 > max_try) {
        warning(sprintf("stalled at step %d: no acceptable mutation in %d tries", s, max_try))
        stalled <- TRUE; break
      }
    }
    if (stalled) break
    n_acc <- n_acc + 1L
    st <- mu$state
    D <- matrix(st$r - r_founder, nrow = 3)
    rec <- rbind(rec, data.frame(step = s, site = site, dV_step = mu$dV,
      V_total = state_vmin(st, exact = exact), rmsd = sqrt(mean(colSums(D^2))),
      TS = state_TS(st, beta) - TS0, n_edge = nrow(st$pr),
      accept_rate = n_acc/n_try))
    if (verbose && s %% 20 == 0)
      cat(sprintf("  step %3d  V=%8.3f  rmsd=%6.3f  acc=%.2f\n",
                  s, state_vmin(st, exact = exact), rec$rmsd[s], n_acc/n_try))
  }
  list(state = st, record = rec, n_try = n_try, n_acc = n_acc, stalled = stalled)
}
