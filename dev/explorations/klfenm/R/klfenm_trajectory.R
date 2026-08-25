## Scans and trajectories.
##
## The discipline: every proposal is evaluated by MINIMISING V exactly and
## reading off the energy there. No expansion conditions the walk. Normal modes
## are computed only to OBSERVE an accepted state, never to decide anything.

## ------------------------------------------------------------------- scan --

## Mutate each site `nmut` times, always from the SAME fixed reference, and
## record what each mutation does. Nothing accumulates: every mutant is
## independent, so this is a property of the reference, not of a walk.
klfenm_scan <- function(ref, nmut = 10, sigma = 0.3, k_update = TRUE,
                        sites = NULL, beta = 1, radius = Inf,
                        seed = NULL, verbose = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  N <- ncol(ref$R)
  if (is.null(sites)) sites <- seq_len(N)
  obs_ref <- klfenm_observe(ref, beta)

  rows <- vector("list", length(sites) * nmut)
  r <- 0L
  for (s in sites) {
    for (m in seq_len(nmut)) {
      mu <- klfenm_mutate_site(ref, s, sigma = sigma, k_update = k_update,
                               radius = radius)
      obs <- klfenm_observe(mu, beta)
      dr2 <- klfenm_dr2_site(ref, mu)
      r <- r + 1L
      rows[[r]] <- data.frame(
        site      = s,
        rep       = m,
        dv_stress = klfenm_delta_v_stress(ref, mu),
        dv_relax  = klfenm_delta_v_relax(ref, mu),
        dv_min    = klfenm_delta_v_min(ref, mu),
        dr2_self  = dr2[s],
        dr2_total = sum(dr2),
        rmsd      = sqrt(mean(dr2)),
        dts       = obs$ts - obs_ref$ts,
        n_broken  = attr(mu, "n_broken"),
        n_formed  = attr(mu, "n_formed"),
        n_active  = obs$n_active,
        n_zero    = obs$n_zero,
        strain    = obs$strain$max_abs
      )
    }
    if (verbose) cat(sprintf("  site %d done\n", s))
  }
  do.call(rbind, rows)
}

## ------------------------------------------------------------- trajectory --

## A walk in parameter space. Each step proposes a mutation at a random site,
## minimises exactly, and accepts on the EXACT dV_min -- never on an estimate.
##
## selection:
##   "neutral"   accept everything
##   "metropolis" accept with prob min(1, exp(-nu dV_min))
##   "threshold" accept iff the cumulative V stays below v_cut
klfenm_trajectory <- function(founder, nstep = 100, sigma = 0.3, nu = 1,
                              v_cut = Inf, k_update = TRUE, radius = Inf,
                              selection = c("metropolis", "neutral", "threshold"),
                              beta = 1, seed = NULL, max_try = 2000,
                              observe_every = 0, verbose = FALSE) {
  selection <- match.arg(selection)
  if (!is.null(seed)) set.seed(seed)
  N  <- ncol(founder$R)
  st <- founder
  v_founder <- klfenm_energy(founder)

  rec <- vector("list", nstep)
  obs_store <- list()
  n_try <- 0L; n_acc <- 0L; stalled <- FALSE

  for (step in seq_len(nstep)) {
    accepted <- FALSE
    tries_this_step <- 0L
    while (!accepted) {
      n_try <- n_try + 1L; tries_this_step <- tries_this_step + 1L
      if (tries_this_step > max_try) {
        warning("stalled at step ", step, " after ", max_try, " tries")
        stalled <- TRUE; break
      }
      site <- sample(N, 1)
      mu   <- klfenm_mutate_site(st, site, sigma = sigma,
                                 k_update = k_update, radius = radius)
      dv   <- klfenm_delta_v_min(st, mu)          # exact, no expansion
      accepted <- switch(selection,
        neutral    = TRUE,
        metropolis = dv <= 0 || stats::runif(1) < exp(-nu * dv),
        threshold  = klfenm_energy(mu) < v_cut)
    }
    if (stalled) break
    n_acc <- n_acc + 1L
    dr2 <- klfenm_dr2_site(founder, mu)
    rec[[step]] <- data.frame(
      step      = step,
      site      = site,
      dv_step   = dv,
      v_total   = klfenm_energy(mu),
      v_from_founder = klfenm_energy(mu) - v_founder,
      rmsd      = sqrt(mean(dr2)),
      n_active  = n_active(mu),
      n_broken  = attr(mu, "n_broken"),
      n_formed  = attr(mu, "n_formed"),
      strain    = klfenm_strain(mu)$max_abs,
      strain_e  = klfenm_strain(mu)$energy,
      tries     = tries_this_step
    )
    st <- mu
    if (observe_every > 0 && step %% observe_every == 0)
      obs_store[[as.character(step)]] <- st
    if (verbose && step %% 10 == 0)
      cat(sprintf("  step %3d  V = %8.3f  strain = %.3f  active = %d\n",
                  step, klfenm_energy(st), klfenm_strain(st)$max_abs, n_active(st)))
  }

  list(state = st, record = do.call(rbind, rec[seq_len(n_acc)]),
       snapshots = obs_store, n_try = n_try, n_acc = n_acc,
       stalled = stalled, founder = founder)
}
