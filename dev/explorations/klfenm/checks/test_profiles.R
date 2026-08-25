## Verification for the profile and trajectory layers (plan items 7-10).
## Run: Rscript checks/test_profiles.R

suppressMessages(library(here))
suppressMessages(devtools::load_all(here::here("..", "..", ".."), quiet = TRUE))
source(here("R", "klfenm_core.R"))
source(here("R", "klfenm_profiles.R"))
source(here("R", "klfenm_trajectory.R"))

ok <- function(label, pass, detail = "") {
  cat(sprintf("%-4s %-54s %s\n", if (pass) "OK" else "FAIL", label, detail))
  if (!pass) assign("FAILED", TRUE, envir = .GlobalEnv)
}
FAILED <- FALSE

data(pdb_2acy_A, package = "penm")
wtp <- set_enm(pdb_2acy_A, node = "ca", model = "anm", d_max = 10.5, frustrated = FALSE)
wt  <- klfenm_set_state(wtp)
N   <- ncol(wt$R)

cat("\n== 7. The frustrated Hessian is actually frustrated ==\n")
## The pilot's silent failure: a state built through set_enm(frustrated=FALSE)
## carries g = 0 regardless of how strained its parameters are, so the whole
## comparison measures nothing. Assert the two Hessians genuinely differ.
tr <- klfenm_trajectory(wt, nstep = 20, sigma = 0.3, selection = "neutral", seed = 99)
st <- tr$state
Kf <- klfenm_kmat(st, frustrated = TRUE)
Kg <- klfenm_kmat(st, frustrated = FALSE)
gap <- max(abs(Kf - Kg))
ok("state carries real strain", klfenm_strain(st)$max_abs > 0.1,
   sprintf("max|d-l| = %.3f", klfenm_strain(st)$max_abs))
ok("K_frustrated != K_{g=0}", gap > 1e-6,
   sprintf("max|dK| = %.4f (scale %.2f)", gap, max(abs(Kf))))
ok("at the wild type the two coincide (l == d)",
   max(abs(klfenm_kmat(wt, frustrated = TRUE) - klfenm_kmat(wt, frustrated = FALSE))) < 1e-12,
   "g = 0 identically when relaxed")

cat("\n== 8. Six zero modes wherever a spectrum is taken ==\n")
for (nm in c("wild type", "strained state")) {
  s <- if (nm == "wild type") wt else st
  sp <- klfenm_nma(klfenm_kmat(s, frustrated = TRUE))
  ok(sprintf("%s: exactly 6 zero modes, none negative", nm),
     sp$n_zero == 6 && min(sp$raw_value) > -1e-8,
     sprintf("n_zero = %d, lowest = %+.2e", sp$n_zero, min(sp$raw_value)))
}

cat("\n== 9. Both causes of the rebuild difference are live ==\n")
## NOT "the decomposition is exact" -- (A-B)+(B-C) = A-C holds for any matrices
## and testing it verifies nothing. What must be checked is that each of the two
## steps is a real, non-empty change: the rebuild really is relaxed (so the
## transverse step does something) and the active set really does move (so the
## topology step does).
rebuild_at <- function(s) { r <- s; r$l <- pair_dist(s$R, s$pr); refresh_k(r) }
rb <- rebuild_at(st)
K_rb <- klfenm_kmat(rb, frustrated = TRUE)
ok("rebuilt network is relaxed by construction",
   klfenm_strain(rb)$max_abs < 1e-12,
   sprintf("max|d-l| = %.2e", klfenm_strain(rb)$max_abs))
ok("rebuilt g == 0, so its frustrated/unfrustrated K agree",
   max(abs(K_rb - klfenm_kmat(rb, frustrated = FALSE))) < 1e-12)
say_edges <- n_active(rb) - n_active(st)
ok("rebuild changes the contact set (topology cause is live)",
   say_edges != 0, sprintf("%+d edges", say_edges))

cat("\n== 10. Mode matching is honest ==\n")
cmp <- klfenm_mode_comparison(klfenm_nma(Kf), klfenm_nma(K_rb), nmodes = 20)
## NOT "best >= index": with a genuine one-to-one assignment a mode may be given
## a worse partner than its own index so the GLOBAL matching is optimal. That
## assertion held only for the old greedy row-maxima, where it was true by
## construction and tested nothing. The meaningful global statement:
ok("assignment maximises the total overlap",
   sum(cmp$overlap_best) >= sum(cmp$overlap_index) - 1e-12,
   sprintf("sum(best) = %.3f vs sum(index) = %.3f; %d modes reordered",
           sum(cmp$overlap_best), sum(cmp$overlap_index), cmp$n_reordered))
ok("RMSIP is reported with per-mode overlaps available",
   length(cmp$overlap_best) == cmp$n_modes,
   sprintf("RMSIP = %.4f, worst per-mode overlap = %.4f",
           cmp$rmsip, min(cmp$overlap_best)))
## a mode compared with itself must give overlap 1
self <- klfenm_mode_comparison(klfenm_nma(Kf), klfenm_nma(Kf), nmodes = 10)
ok("self-comparison gives overlap 1 and RMSIP 1",
   all(abs(self$overlap_best - 1) < 1e-9) && abs(self$rmsip - 1) < 1e-9)

## row-wise which.max is NOT a matching: two modes can claim one partner, and
## each then reports its best available overlap regardless of contention, which
## inflates every number. Measured when it did: 13 of 18 states had duplicates,
## and the worst overlap at seed 501/60 read 0.434 instead of 0.288.
ok("mode matching is one-to-one", !any(duplicated(cmp$matched_to)),
   sprintf("%d duplicate assignments among %d modes",
           sum(duplicated(cmp$matched_to)), cmp$n_modes))

cat("\n== 11. Trajectory bookkeeping ==\n")
tr2 <- klfenm_trajectory(wt, nstep = 15, sigma = 0.3, nu = 1,
                         selection = "metropolis", seed = 7)
ok("records one row per accepted step", nrow(tr2$record) == tr2$n_acc,
   sprintf("%d accepted of %d tries", tr2$n_acc, tr2$n_try))
ok("energy from the founder is cumulative and consistent",
   abs(tail(tr2$record$v_from_founder, 1) -
       (klfenm_energy(tr2$state) - klfenm_energy(wt))) < 1e-10)
## Tested against the minimiser's OWN tolerance, not a hardcoded number that
## silently goes stale when the tolerance changes (it did: 1e-9 vs tol 1e-8).
MIN_TOL <- formals(klfenm_minimise)$tol
ok("every accepted state is at its own minimum",
   sqrt(sum(energy_gradient(tr2$state)^2)) < MIN_TOL,
   sprintf("|F| = %.2e (tol = %.0e)", sqrt(sum(energy_gradient(tr2$state)^2)), MIN_TOL))

cat("\n== 12. Downhill moves: absent when relaxed, present when strained ==\n")
## Not a tautology: it is MODEL.md's theorem (dV >= 0 from a relaxed network)
## on one side, and the claim that strain opens a downhill channel on the other.
set.seed(5)
dv_wt <- replicate(80, klfenm_delta_v_min(wt, klfenm_mutate_site(wt, sample(N, 1), 0.3)))
set.seed(6)
dv_st <- replicate(80, klfenm_delta_v_min(st, klfenm_mutate_site(st, sample(N, 1), 0.3)))
ok("no downhill move from the relaxed founder", all(dv_wt > 0),
   sprintf("min dV = %+.4f", min(dv_wt)))
ok("downhill moves exist from a strained state", any(dv_st < 0),
   sprintf("min dV = %+.4f, %.1f%% negative", min(dv_st), 100 * mean(dv_st < 0)))

cat("\n")
if (FAILED) { cat("SOME CHECKS FAILED\n"); quit(status = 1) } else cat("all checks passed\n")
