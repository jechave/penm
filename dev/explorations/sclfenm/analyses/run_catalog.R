library(here)
source(here("R", "penm_catalog.R"))
wt <- founder()
cat("building catalogue: 98 sites x 10 states...\n")
t0 <- Sys.time()
cg <- build_catalog(wt, M = 10, sigma = 0.3, verbose = FALSE)
cat(sprintf("done in %.1f s\n\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
saveRDS(cg, here("data", "catalog.rds"))
cat("dV over the whole catalogue: range [", sprintf("%.3f", min(cg$dV)), ",",
    sprintf("%.3f", max(cg$dV)), "], all >= 0 from the founder\n\n")
for (nu in c(0, 0.5, 1, 2, 5)) {
  tr <- run_catalog_trajectory(cg, nstep = 3000, nu = nu, seed = 42)
  r <- tr$record; tail100 <- r$V[(nrow(r)-99):nrow(r)]
  cat(sprintf("nu = %4.1f : V(final) = %7.2f  <V>_last100 = %7.2f  sd = %5.2f  sites mutated = %2d/98  accept = %.3f\n",
      nu, r$V[nrow(r)], mean(tail100), sd(tail100), r$n_mut[nrow(r)], tr$n_acc/3000))
}
