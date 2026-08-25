## THE NULL CONTROL for section 6.
##
## Section 6 measures how much the rebuilt model differs from the frustrated one.
## On its own that number is uninterpretable: this spectrum is dense, so ANY
## perturbation of comparable size will scramble individual eigenvectors. The
## question is whether the rebuild differs MORE than a null perturbation of the
## same Hessian magnitude does.
##
## The null: jitter the spring constants of the state's OWN network by a random
## multiplicative factor, with the active set held fixed. That changes no
## topology and adds no frustration, so anything it reproduces cannot be
## evidence about either.
##
## Run: Rscript analyses/run_control.R   -> data/control.rds

suppressMessages(library(here))
suppressMessages(devtools::load_all(here::here("..", "..", ".."), quiet = TRUE))
source(here("R", "klfenm_core.R"))
source(here("R", "klfenm_profiles.R"))
source(here("R", "klfenm_trajectory.R"))

data(pdb_2acy_A, package = "penm")
wtp <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5, frustrated = FALSE)
wt  <- klfenm_set_state(wtp)
say <- function(...) cat(sprintf(...), sep = "")

rebuild_at <- function(st) { rb <- st; rb$l <- pair_dist(st$R, st$pr); refresh_k(rb) }

## Jitter k on the SAME active set, tuned so max|dK| matches the rebuild's.
jitter_matched <- function(st, target, seed, tol = 0.05, maxit = 40) {
  K0 <- klfenm_kmat(st, frustrated = TRUE)
  a  <- which(st$k > 0)
  set.seed(seed)
  z  <- rnorm(length(a))                 # fixed direction; scale is tuned
  lo <- 0; hi <- 2
  for (it in seq_len(maxit)) {
    s  <- (lo + hi) / 2
    ct <- st; ct$k[a] <- st$k[a] * pmax(1e-6, 1 + s * z)
    d  <- max(abs(klfenm_kmat(ct, frustrated = TRUE) - K0))
    if (abs(d - target) / target < tol) break
    if (d < target) lo <- s else hi <- s
  }
  attr(ct, "dK") <- d; attr(ct, "scale") <- s
  ct
}

rows <- list(); r <- 0L
for (sd in 1:3) {
  tr <- klfenm_trajectory(wt, nstep = 60, sigma = 0.3, selection = "neutral",
                          seed = 500 + sd, observe_every = 20)
  for (k in names(tr$snapshots)) {
    st  <- tr$snapshots[[k]]
    Kfr <- klfenm_kmat(st, frustrated = TRUE); nfr <- klfenm_nma(Kfr)
    rb  <- rebuild_at(st)
    Krb <- klfenm_kmat(rb, frustrated = TRUE); nrb <- klfenm_nma(Krb)
    dK  <- max(abs(Kfr - Krb))
    rmsf <- function(K) sqrt(klfenm_msf_site(klfenm_cmat(K, 1)))
    r_fr <- rmsf(Kfr)

    real <- klfenm_mode_comparison(nfr, nrb, 20)
    real_rmsf <- rel_profile(r_fr, rmsf(Krb))

    ## three independent null draws at the same |dK|
    for (rep in 1:3) {
      ct  <- jitter_matched(st, dK, seed = 7000 + 100 * sd + 10 * as.integer(k) + rep)
      Kct <- klfenm_kmat(ct, frustrated = TRUE); nct <- klfenm_nma(Kct)
      null <- klfenm_mode_comparison(nfr, nct, 20)
      null_rmsf <- rel_profile(r_fr, rmsf(Kct))
      r <- r + 1L
      rows[[r]] <- data.frame(
        seed = sd, subs = as.integer(k), rep = rep,
        dK_real = dK, dK_null = attr(ct, "dK"),
        same_active = n_active(ct) == n_active(st),
        worst_real = min(real$overlap_best), worst_null = min(null$overlap_best),
        rmsip_real = real$rmsip,            rmsip_null = null$rmsip,
        blk3_real = median(klfenm_block_overlap(nfr, nrb, 20, 3)),
        blk3_null = median(klfenm_block_overlap(nfr, nct, 20, 3)),
        rmsf_real = max(abs(real_rmsf$rel)), rmsf_null = max(abs(null_rmsf$rel)),
        eig_real = max(abs(real$eigen_rel_diff)),
        eig_null = max(abs(null$eigen_rel_diff))
      )
    }
    say("  seed %d subs %2s: |dK| %.2f | worst mode real %.3f vs null %.3f | RMSF real %.1f%% vs null %.1f%%\n",
        sd, k, dK, min(real$overlap_best),
        median(sapply(rows[(r-2):r], function(x) x$worst_null)),
        max(abs(real_rmsf$rel)),
        median(sapply(rows[(r-2):r], function(x) x$rmsf_null)))
  }
}

d <- do.call(rbind, rows)
saveRDS(d, here("data", "control.rds"))

say("\n=== summary: does the rebuild differ MORE than a matched null? ===\n")
say("worst single-mode overlap : real %.3f vs null %.3f  (higher = more similar)\n",
    median(d$worst_real), median(d$worst_null))
say("RMSIP                     : real %.3f vs null %.3f\n",
    median(d$rmsip_real), median(d$rmsip_null))
say("median block-3 overlap    : real %.3f vs null %.3f\n",
    median(d$blk3_real), median(d$blk3_null))
say("worst-site RMSF error (%%) : real %.1f vs null %.1f\n",
    median(d$rmsf_real), median(d$rmsf_null))
say("max eigenvalue rel diff   : real %.3f vs null %.3f\n",
    median(d$eig_real), median(d$eig_null))
say("\nnull draws keeping the active set unchanged: %d of %d\n",
    sum(d$same_active), nrow(d))
