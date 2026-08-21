library(here)
## Is the keepnet run stationary, or just not yet at the wall?
## Extend well past the ceiling and look at the STEADY-STATE acceptance.
source(here("R", "applications.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5)
tt <- suppressWarnings(run_trajectory(wt, nstep=70, sigma=.3, v_cut=25,
        model="keepnet", selection="threshold", seed=7, max_try=1500))
r <- tt$record
## recover per-step tries from cumulative acceptance rate
cum_try <- r$step / r$accept_rate
r$tries <- c(cum_try[1], diff(cum_try))
cat(sprintf("keepnet, v_cut=25, 70 steps requested, %d completed, stalled=%s\n\n", nrow(r), tt$stalled))
cat(sprintf("%-12s %8s %10s %10s %8s\n","block","mean V","tries/acc","rmsd","n"))
for (b in seq(1, nrow(r), by=20)) {
  ix <- b:min(b+19, nrow(r))
  cat(sprintf("steps %3d-%3d %8.2f %10.2f %10.3f %8d\n",
      min(ix), max(ix), mean(r$V_total[ix]), mean(r$tries[ix]),
      mean(r$rmsd[ix]), length(ix)))
}
cat(sprintf("\nsteady-state acceptance (last 40 steps) = %.4f\n",
    1/mean(r$tries[(nrow(r)-39):nrow(r)])))
cat(sprintf("SC-LFENM acceptance at its stall            = 0.0467\n"))
saveRDS(r, here("data", "keepnet_long.rds"))
