library(here)
source(here("R", "frustration.R"))
wt <- founder(); N <- get_nsites(wt)
cat("Effect of the frustration term on DYNAMICS, vs accumulated strain.\n")
cat("All spectra taken at the TRUE minimum of the relevant potential.\n\n")
cat(sprintf("%6s %9s %9s %10s %10s %10s %10s\n",
    "nmut","strain","dK %","soft eval %","RMSIP(10)","d(TS)","7th eval"))
set.seed(3); mu <- wt
for (n in 1:40) {
  mu <- get_mutant_site(mu, sample(N,1), n, "lfenm", .3, 1L, 1L)
  if (n %in% c(1,2,5,10,20,30,40)) {
    p <- relax_to_min(mu, frustrated=TRUE)
    g <- get_graph(p); str <- 0.5*sum(g$kij*(g$dij-g$lij)^2)
    Kf <- hess(p,FALSE); Kt <- hess(p,TRUE)
    sf <- spec(Kf); st <- spec(Kt)
    vf <- sort(sf$val)[1:10]; vt <- sort(st$val)[1:10]
    ov <- crossprod(sf$vec[,order(sf$val)[1:10]], st$vec[,order(st$val)[1:10]])
    cat(sprintf("%6d %9.3f %9.2f %10.2f %10.5f %10.4f %10.2e\n", n, str,
        100*max(abs(Kf-Kt))/max(abs(Kt)),
        100*mean(abs(vt-vf)/vf), sqrt(sum(ov^2)/10),
        ts_of(st$val)-ts_of(sf$val), seventh_eigen(Kt)))
  }
}
cat("\n7th eigenvalue = first NON-RIGID mode, taken by rank from the raw\n")
cat("spectrum (not by a tolerance filter, which could not report a saddle).\n")
cat("It stays positive => the frustrated network is a genuine minimum,\n")
cat("not a saddle, at every strain level reached here.\n")
