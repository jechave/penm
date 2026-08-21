## The reviewer: fig8(c) tests only that my MCMC is correctly implemented,
## because the chain samples V_cat and V_cat is separable by construction.
## The REAL question: does the chain sample exp(-nu V) for the EXACT V?
source("penm_catalog.R")
cg <- readRDS("catalog.rds")
cat("Chain runs on V_cat; report the EXACT V of the sampled sequences.\n")
cat("If V_cat is a good surrogate, <V_exact> should track the prediction too.\n\n")
pred <- function(nu) sum(sapply(seq_len(cg$N), function(j){
  w<-exp(-nu*cg$dV[j,]); sum(w*cg$dV[j,])/sum(w)}))
cat(sprintf("%6s %12s %14s %14s %10s\n","nu","predicted","<V_cat>","<V_exact>","cat-exact"))
for (nu in c(1,2,5)) {
  tr <- run_catalog_trajectory(cg, 8000, nu=nu, seed=11)
  ## sample 40 sequences from the equilibrated tail and evaluate both energies
  set.seed(3)
  vs_cat <- vs_ex <- numeric(0)
  st <- tr$s
  for (rep in 1:12) {
    tr2 <- run_catalog_trajectory(cg, 200, nu=nu, s0=st, seed=100+rep)
    st <- tr2$s
    vs_cat <- c(vs_cat, seq_dV_catalog(cg, st))
    vs_ex  <- c(vs_ex,  seq_dV_exact(cg, st))
  }
  cat(sprintf("%6.1f %12.2f %14.2f %14.2f %9.2f%%\n", nu, pred(nu),
      mean(vs_cat), mean(vs_ex), 100*(mean(vs_cat)-mean(vs_ex))/mean(vs_ex)))
}
