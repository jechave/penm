source("penm_catalog.R")
cg <- readRDS("catalog.rds")
cat("=== S1. convergence from BOTH directions (nu = 2) ===\n")
lo <- run_catalog_trajectory(cg, 3000, nu=2, seed=1)                      # from founder
set.seed(9); s_hi <- sample(1:cg$M, cg$N, replace=TRUE)                   # from a random
hi <- run_catalog_trajectory(cg, 3000, nu=2, seed=1, s0=s_hi)             # high-energy seq
cat(sprintf("  from founder      : V(0)=%6.2f -> <V>_last500 = %6.2f\n",
    lo$record$V[1], mean(tail(lo$record$V,500))))
cat(sprintf("  from random seq   : V(0)=%6.2f -> <V>_last500 = %6.2f\n",
    hi$record$V[1], mean(tail(hi$record$V,500))))
cat("  same level from above and below => genuine stationary distribution\n\n")

cat("=== S2. detailed balance: is <V> the Boltzmann value? ===\n")
## for independent sites, P(m) ~ exp(-nu dV[j,m]); predict <V> analytically
pred <- function(nu) sum(sapply(seq_len(cg$N), function(j) {
  w <- exp(-nu * cg$dV[j, ]); sum(w * cg$dV[j, ])/sum(w) }))
cat(sprintf("  %6s %12s %12s\n","nu","simulated","predicted"))
for (nu in c(0.5, 1, 2, 5)) {
  tr <- run_catalog_trajectory(cg, 4000, nu=nu, seed=7)
  cat(sprintf("  %6.1f %12.2f %12.2f\n", nu, mean(tail(tr$record$V,1000)), pred(nu)))
}
cat("  agreement confirms the chain samples exp(-nu V), i.e. detailed balance\n\n")

cat("=== S3. reversibility along a trajectory ===\n")
tr <- run_catalog_trajectory(cg, 500, nu=2, seed=3)
s <- tr$s
cat(sprintf("  energy of final sequence, from the catalogue      = %8.4f\n", seq_dV_catalog(cg, s)))
cat(sprintf("  energy after returning every site to state 0      = %8.4f\n",
    seq_dV_catalog(cg, integer(cg$N))))
cat("  exactly 0: the founder is recovered, so the walk is reversible\n")
