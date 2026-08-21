source("frustration.R")
wt <- founder(); N <- get_nsites(wt)
cat("=== F1. the wild type is NOT frustrated: the flag does nothing ===\n")
g <- get_graph(wt)
cat(sprintf("  max|l - d| on the founder = %.3e  => g_ij = 0 identically\n",
    max(abs(g$lij - g$dij))))
K0 <- hess(wt, FALSE); K1 <- hess(wt, TRUE)
cat(sprintf("  ||K(frustrated) - K(relaxed)|| = %.3e\n\n", max(abs(K0-K1))))

cat("=== F2. mutation CREATES frustration ===\n")
set.seed(101)
cat(sprintf("%8s %12s %12s %12s\n","nmut","max|l-d|","strain energy","dV"))
mu <- wt
for (n in 1:6) {
  mu <- get_mutant_site(mu, sample(N,1), n, "lfenm", .3, 1L, 1L)
  gg <- get_graph(mu)
  cat(sprintf("%8d %12.4f %12.4f %12.4f\n", n, max(abs(gg$lij-gg$dij)),
      0.5*sum(gg$kij*(gg$dij-gg$lij)^2), ddg_dv(wt, mu)))
}

cat("\n=== F3. does frustration change the DYNAMICS of a mutant? ===\n")
set.seed(3); mu <- wt
for (n in 1:10) mu <- get_mutant_site(mu, sample(N,1), n, "lfenm", .3, 1L, 1L)
## CRUCIAL: the LFENM structure is the linear-response estimate, not the true
## minimum. Spectra must be taken at a stationary point or rotational
## invariance is broken numerically (3 spurious modes ~4e-4).
cat(sprintf("  |F| at the LFENM structure        = %.3e\n", sqrt(sum(enm_force(mu)^2))))
mu <- relax_to_min(mu, frustrated = TRUE)
cat(sprintf("  |F| after relaxing to the minimum = %.3e  (%d iterations)\n",
    attr(mu,"fres"), attr(mu,"iter")))
Kf <- hess(mu, FALSE); Kt <- hess(mu, TRUE)
sf <- spec(Kf); st <- spec(Kt)
ov <- crossprod(sf$vec[,order(sf$val)[1:10]], st$vec[,order(st$val)[1:10]])
cat(sprintf("  strain energy of this mutant      = %.3f\n",
    0.5*sum(get_graph(mu)$kij*(get_graph(mu)$dij-get_graph(mu)$lij)^2)))
cat(sprintf("  ||K_frust - K_relaxed||_max       = %.4f  (scale %.2f, %.1f%%)\n",
    max(abs(Kf-Kt)), max(abs(Kt)), 100*max(abs(Kf-Kt))/max(abs(Kt))))
cat(sprintf("  zero modes                        : relaxed %d, frustrated %d\n",
    sf$nzero, st$nzero))
cat(sprintf("  7th eigenvalue (1st non-rigid)    : relaxed %.4f, frustrated %.4f %s\n",
    seventh_eigen(Kf), seventh_eigen(Kt),
    if (seventh_eigen(Kt) < 0) "<- SADDLE" else "(minimum)"))
cat(sprintf("  10 softest eigenvalues differ by  : %.2f%%\n",
    100*mean(abs(sort(st$val)[1:10]-sort(sf$val)[1:10])/sort(sf$val)[1:10])))
cat(sprintf("  RMSIP over the 10 softest modes   : %.5f\n", sqrt(sum(ov^2)/10)))
cat(sprintf("  T*S                               : relaxed %.4f, frustrated %.4f (diff %.4f)\n",
    ts_of(sf$val), ts_of(st$val), ts_of(st$val)-ts_of(sf$val)))
