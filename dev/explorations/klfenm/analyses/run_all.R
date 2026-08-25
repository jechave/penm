## Every number the report quotes. Writes data/results.rds.
## Run: Rscript analyses/run_all.R
##
## Seeded throughout: rerunning must reproduce, and any drift is a bug.

suppressMessages(library(here))
suppressMessages(devtools::load_all(here::here("..", "..", ".."), quiet = TRUE))
source(here("R", "klfenm_core.R"))
source(here("R", "klfenm_profiles.R"))
source(here("R", "klfenm_trajectory.R"))

data(pdb_2acy_A, package = "penm")
wtp <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5, frustrated = FALSE)
wt  <- klfenm_set_state(wtp)
N   <- ncol(wt$R)
BETA <- 1                # convention in ANM units; k_ij = 1 is dimensionless
out <- list(nsites = N, n_contacts_wt = n_active(wt), beta = BETA)
t_start <- Sys.time()
say <- function(...) cat(sprintf(...), sep = "")

## ============================================================ A. contact flips
## How often does k(l) change the network, as a function of sigma?
say("A. contact-flip frequency vs sigma\n")
set.seed(1)
flip <- do.call(rbind, lapply(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.6), function(sg) {
  nb <- nf <- integer(200)
  for (t in 1:200) {
    site <- sample(N, 1)
    sel <- which((wt$pr$i == site | wt$pr$j == site) & wt$sdij >= 1)
    s2 <- wt; s2$l[sel] <- s2$l[sel] + rnorm(length(sel), 0, sg); s2 <- refresh_k(s2)
    nb[t] <- sum(wt$k > 0 & s2$k == 0); nf[t] <- sum(wt$k == 0 & s2$k > 0)
  }
  data.frame(sigma = sg, frac_any = mean((nb + nf) > 0),
             mean_broken = mean(nb), mean_formed = mean(nf), max_events = max(nb + nf))
}))
out$flip <- flip; print(flip)

## ================================================== B. downhill moves and strain
## MODEL.md's theorem says dV >= 0 from a RELAXED network. Does strain open a
## downhill channel here? Measured at several strain levels along one walk.
say("\nB. downhill availability vs accumulated strain\n")
set.seed(2)
walk <- klfenm_trajectory(wt, nstep = 40, sigma = 0.3, selection = "neutral",
                          seed = 2, observe_every = 10)
probe <- function(st, ndraw = 200, seed) {
  set.seed(seed)
  dv <- replicate(ndraw, klfenm_delta_v_min(st, klfenm_mutate_site(st, sample(N, 1), 0.3)))
  data.frame(strain = klfenm_strain(st)$max_abs,
             strain_e = klfenm_strain(st)$energy,
             v = klfenm_energy(st),
             mean_dv = mean(dv), min_dv = min(dv), frac_neg = mean(dv < 0))
}
downhill <- rbind(
  cbind(subs = 0, probe(wt, seed = 100)),
  do.call(rbind, lapply(names(walk$snapshots), function(k)
    cbind(subs = as.integer(k), probe(walk$snapshots[[k]], seed = 100 + as.integer(k)))))
)
out$downhill <- downhill; print(downhill)

## ================================================================== C. scans
## Site profiles from the fixed founder, k(l) live vs frozen.
say("\nC. scans (k_update TRUE and FALSE), 8 mutations x 98 sites\n")
out$scan_on  <- klfenm_scan(wt, nmut = 8, sigma = 0.3, k_update = TRUE,  seed = 11)
out$scan_off <- klfenm_scan(wt, nmut = 8, sigma = 0.3, k_update = FALSE, seed = 11)
cn <- tabulate(c(wt$pr$i[wt$k > 0], wt$pr$j[wt$k > 0]), nbins = N)
out$cn <- cn
agg <- function(d) {
  a <- aggregate(cbind(dv_min, dr2_self, dr2_total, dts, n_broken, n_formed) ~ site,
                 data = d, FUN = mean)
  a$cn <- cn[a$site]; a
}
out$scan_on_agg <- agg(out$scan_on); out$scan_off_agg <- agg(out$scan_off)
say("   k(l) live : cor(dV,cn) = %+.3f ; cor(dr2_self,cn) = %+.3f ; dTS range [%.3f, %.3f]\n",
    cor(out$scan_on_agg$dv_min, cn), cor(out$scan_on_agg$dr2_self, cn),
    min(out$scan_on_agg$dts), max(out$scan_on_agg$dts))
say("   k frozen  : cor(dV,cn) = %+.3f ; cor(dr2_self,cn) = %+.3f ; dTS range [%.3f, %.3f]\n",
    cor(out$scan_off_agg$dv_min, cn), cor(out$scan_off_agg$dr2_self, cn),
    min(out$scan_off_agg$dts), max(out$scan_off_agg$dts))
say("   dV profiles across the two: cor = %.4f ; mean rel diff = %+.1f%%\n",
    cor(out$scan_on_agg$dv_min, out$scan_off_agg$dv_min),
    100 * mean((out$scan_on_agg$dv_min - out$scan_off_agg$dv_min) / out$scan_off_agg$dv_min))

## ============================================================ D. trajectories
say("\nD. trajectories, 120 steps, three selection strengths\n")
out$traj <- lapply(c(nu0.5 = 0.5, nu1 = 1, nu4 = 4), function(nu) {
  tr <- klfenm_trajectory(wt, nstep = 120, sigma = 0.3, nu = nu,
                          selection = "metropolis", seed = 300 + round(nu * 10),
                          observe_every = 20)
  say("   nu = %-4g : V %.2f -> %.2f, strain %.3f, active %d, acceptance %.3f, stalled %s\n",
      nu, 0, klfenm_energy(tr$state), klfenm_strain(tr$state)$max_abs,
      n_active(tr$state), tr$n_acc / tr$n_try, tr$stalled)
  tr
})

## ==================================== E. frustration vs rebuild (the headline)
## At a strained state, compare the state's OWN model with what set_enm() would
## build at the same coordinates. Three causes, separated:
##   (i)   cross term only   - same graph, g != 0 vs g = 0
##   (ii)  topology only     - g = 0 both sides, k(l) graph vs k(d) graph
##   (iii) full rebuild      - what set_enm(r_mut) actually gives
say("\nE. frustration vs rebuild, along 3 seeded trajectories\n")

## The rebuild: an ENM fitted to the structure, l = d, contacts from d.
rebuild_at <- function(st) {
  rb <- st
  d  <- pair_dist(st$R, st$pr)
  rb$l <- d                       # every spring relaxed at THIS structure
  rb   <- refresh_k(rb)           # contacts re-derived, now from d (since l = d)
  rb
}

compare_one <- function(st, beta = BETA, nmodes = 20) {
  ## all three models share the SAME coordinates
  K_fr   <- klfenm_kmat(st, frustrated = TRUE)    # own params, cross term ON
  K_g0   <- klfenm_kmat(st, frustrated = FALSE)   # own params, cross term OFF
  rb     <- rebuild_at(st)
  K_rb   <- klfenm_kmat(rb, frustrated = TRUE)    # rebuilt: l = d so g = 0 anyway

  nma_fr <- klfenm_nma(K_fr); nma_g0 <- klfenm_nma(K_g0); nma_rb <- klfenm_nma(K_rb)
  rmsf   <- function(K) sqrt(klfenm_msf_site(klfenm_cmat(K, beta)))
  r_fr <- rmsf(K_fr); r_g0 <- rmsf(K_g0); r_rb <- rmsf(K_rb)

  list(
    strain      = klfenm_strain(st)$max_abs,
    strain_e    = klfenm_strain(st)$energy,
    n_active_fr = n_active(st), n_active_rb = n_active(rb),
    d_edges     = n_active(rb) - n_active(st),
    zero_fr = nma_fr$n_zero, zero_rb = nma_rb$n_zero,
    lowest_fr = min(nma_fr$raw_value), lowest_rb = min(nma_rb$raw_value),
    ## the frustrated Hessian must actually BE frustrated
    k_gap_cross = max(abs(K_fr - K_g0)),
    k_gap_topo  = max(abs(K_g0 - K_rb)),
    k_gap_full  = max(abs(K_fr - K_rb)),
    k_scale     = max(abs(K_fr)),
    ## Frobenius norms of the two steps. That they sum to the total is
    ## trivial -- A-B plus B-C is A-C for any matrices at all, and checking it
    ## verifies nothing. What is NOT trivial, and is what licenses comparing
    ## the two norms, is whether they add in QUADRATURE: that says the two
    ## perturbations are near-orthogonal rather than one being largely a
    ## re-description of the other. Measured: within 0.34%.
    fro_cross = norm(K_fr - K_g0, "F"),
    fro_topo  = norm(K_g0 - K_rb, "F"),
    fro_full  = norm(K_fr - K_rb, "F"),
    ## (i) cross term only
    cross = list(rmsf = rel_profile(r_fr, r_g0),
                 mode = klfenm_mode_comparison(nma_fr, nma_g0, nmodes),
                 rwsip = klfenm_rwsip(nma_fr, nma_g0, nmodes)),
    ## (ii) topology only (both unfrustrated)
    topo  = list(rmsf = rel_profile(r_g0, r_rb),
                 mode = klfenm_mode_comparison(nma_g0, nma_rb, nmodes),
                 rwsip = klfenm_rwsip(nma_g0, nma_rb, nmodes)),
    ## (iii) the full rebuild, which is what an ENM practitioner would do
    full  = list(rmsf = rel_profile(r_fr, r_rb),
                 mode = klfenm_mode_comparison(nma_fr, nma_rb, nmodes),
                 rwsip = klfenm_rwsip(nma_fr, nma_rb, nmodes)),
    ts_fr = klfenm_entropy(K_fr, beta), ts_rb = klfenm_entropy(K_rb, beta),
    ## single-mode overlaps are ill-conditioned where eigenvalues are close;
    ## blocks are invariant to rotation within a block. Both are reported.
    block1 = klfenm_block_overlap(nma_fr, nma_rb, nmodes, 1),
    block3 = klfenm_block_overlap(nma_fr, nma_rb, nmodes, 3),
    block5 = klfenm_block_overlap(nma_fr, nma_rb, nmodes, 5),
    eigen_gaps = klfenm_eigen_gaps(nma_fr, nmodes)
  )
}

out$rebuild <- list()
for (sd in 1:3) {
  tr <- klfenm_trajectory(wt, nstep = 60, sigma = 0.3, selection = "neutral",
                          seed = 500 + sd, observe_every = 10)
  cmp <- lapply(names(tr$snapshots), function(k) {
    z <- compare_one(tr$snapshots[[k]]); z$subs <- as.integer(k); z$seed <- sd; z
  })
  out$rebuild[[sd]] <- cmp
  for (z in cmp)
    say("   seed %d subs %2d: strain %.3f | dK cross %.4f full %.4f (scale %.1f) | edges %+d | RMSF max %+.1f%% | worst mode overlap %.3f | RMSIP %.3f\n",
        sd, z$subs, z$strain, z$k_gap_cross, z$k_gap_full, z$k_scale, z$d_edges,
        z$full$rmsf$rel[which.max(abs(z$full$rmsf$rel))],
        min(z$full$mode$overlap_best), z$full$mode$rmsip)
}

## ============================================================ F. smooth k(l)
say("\nF. smooth k(l): does the discontinuity matter?\n")
wt_sm <- klfenm_set_state(wtp, k_model = "smooth")
wt_sm$k_par$w <- 0.5
wt_sm <- refresh_k(wt_sm)
say("   smooth wt: active(k>0.5) = %d vs step %d ; |F| = %.2e\n",
    sum(wt_sm$k > 0.5), n_active(wt), sqrt(sum(energy_gradient(wt_sm)^2)))
out$smooth <- list()
for (w in c(0.25, 0.5)) {
  s0 <- klfenm_set_state(wtp, k_model = "smooth"); s0$k_par$w <- w; s0 <- refresh_k(s0)
  s0 <- klfenm_minimise(s0)
  tr <- klfenm_trajectory(s0, nstep = 40, sigma = 0.3, nu = 1,
                          selection = "metropolis", seed = 700 + w * 100)
  out$smooth[[as.character(w)]] <- tr
  say("   w = %.2f : V -> %.2f, strain %.3f, acceptance %.3f\n",
      w, klfenm_energy(tr$state), klfenm_strain(tr$state)$max_abs, tr$n_acc / tr$n_try)
}

## ===================================== G. flat summaries the report quotes
## These were once produced ad hoc, leaving .rds files with no generating
## script. Everything the report quotes must come from a script.
say("\nG. flat summaries\n")
out$flat <- do.call(rbind, lapply(out$rebuild, function(sl) do.call(rbind, lapply(sl, function(z)
  data.frame(seed = z$seed, subs = z$subs, strain_e = z$strain_e, d_edges = z$d_edges,
             cross = max(abs(z$cross$rmsf$rel)), topo = max(abs(z$topo$rmsf$rel)),
             full = max(abs(z$full$rmsf$rel)), cor_full = z$full$rmsf$cor,
             ov_full = min(z$full$mode$overlap_best), rmsip = z$full$mode$rmsip,
             reord = z$full$mode$n_reordered,
             eig = 100 * max(abs(z$full$mode$eigen_rel_diff)),
             dts = z$ts_rb - z$ts_fr,
             fro_cross = z$fro_cross, fro_topo = z$fro_topo, fro_full = z$fro_full,
             blk1 = median(z$block1), blk3 = median(z$block3), blk5 = median(z$block5))))))
q <- out$flat
q$quad <- sqrt(q$fro_cross^2 + q$fro_topo^2)
say("   Frobenius: cross %.2f topo %.2f full %.2f ; quadrature gap %.2f%% ; ratio %.2f\n",
    median(q$fro_cross), median(q$fro_topo), median(q$fro_full),
    100 * median(abs(q$fro_full - q$quad) / q$fro_full), median(q$fro_topo / q$fro_cross))
say("   blocks: b1 %.3f b3 %.3f b5 %.3f (medians) ; worst-state b5 %.3f\n",
    median(q$blk1), median(q$blk3), median(q$blk5), min(q$blk5))

## ============================ H. is the topology term a cutoff artefact?
## Every differing edge is a marginal contact; the step function counts it at
## full weight. Re-run the SAME decomposition with a smooth k(l) on both sides.
say("\nH. cutoff-sharpness sweep\n")
sweep <- list()
for (sd in seq_along(out$rebuild)) {
  tr <- klfenm_trajectory(wt, nstep = 60, sigma = 0.3, selection = "neutral",
                          seed = 500 + sd, observe_every = 20)
  for (kk in names(tr$snapshots)) {
    st0 <- tr$snapshots[[kk]]
    for (m in c("step", "0.25", "0.50", "1.00")) {
      s2 <- st0
      if (m != "step") { s2$k_model <- "smooth"; s2$k_par$w <- as.numeric(m); s2 <- refresh_k(s2) }
      Kf <- klfenm_kmat(s2, frustrated = TRUE); Kg <- klfenm_kmat(s2, frustrated = FALSE)
      rb2 <- s2; rb2$l <- pair_dist(s2$R, s2$pr); rb2 <- refresh_k(rb2)
      Kr <- klfenm_kmat(rb2, frustrated = TRUE)
      sweep[[length(sweep) + 1]] <- data.frame(seed = sd, subs = as.integer(kk), cutoff = m,
        cross = norm(Kf - Kg, "F"), topo = norm(Kg - Kr, "F"),
        ratio = norm(Kg - Kr, "F") / norm(Kf - Kg, "F"))
    }
  }
}
out$cutoff_sweep <- do.call(rbind, sweep)
print(round(tapply(out$cutoff_sweep$ratio, out$cutoff_sweep$cutoff, median), 2))

out$elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
saveRDS(out, here("data", "results.rds"))
say("\ndone in %.1f min -> data/results.rds\n", out$elapsed)
