library(here)
## NOTES.Rmd asks: "Rigorously, should I use K_wt or K_mut to calculate dr?"
## Never tested. Test it: compare dr = C_wt f against the SELF-CONSISTENT
## displacement (iterate until the mutant's own K is used at its own minimum).
source(here("R", "frustration.R"))
wt <- founder(); N <- get_nsites(wt)
cat("K_wt or K_mut for the structural response?\n")
cat("  dr_LRA  = C_wt f            (what penm does: one linear step)\n")
cat("  dr_SC   = relax to the true minimum of the mutant potential\n\n")
cat(sprintf("%6s %11s %11s %11s %9s %9s\n",
    "site","|dr_LRA|","|dr_SC|","|difference|","rel err","cos angle"))
set.seed(5)
for (j in sample(N, 8)) {
  mu <- get_mutant_site(wt, j, 1, "lfenm", .3, 1L, 1L)
  drL <- as.vector(get_xyz(mu)) - as.vector(get_xyz(wt))
  ms  <- relax_to_min(mu, frustrated = TRUE)
  drS <- as.vector(get_xyz(ms)) - as.vector(get_xyz(wt))
  cat(sprintf("%6d %11.5f %11.5f %11.5f %8.2f%% %9.5f\n", j,
      sqrt(sum(drL^2)), sqrt(sum(drS^2)), sqrt(sum((drS-drL)^2)),
      100*sqrt(sum((drS-drL)^2))/sqrt(sum(drS^2)),
      sum(drL*drS)/sqrt(sum(drL^2)*sum(drS^2))))
}
cat("\nSame, but as the mutation gets larger (site 40):\n")
cat(sprintf("%8s %11s %11s %9s\n","sigma","|dr_LRA|","|dr_SC|","rel err"))
for (sg in c(.1,.2,.3,.5,.8)) {
  mu <- get_mutant_site(wt, 40, 1, "lfenm", sg, 1L, 1L)
  drL <- as.vector(get_xyz(mu)) - as.vector(get_xyz(wt))
  ms <- relax_to_min(mu, frustrated=TRUE)
  drS <- as.vector(get_xyz(ms)) - as.vector(get_xyz(wt))
  cat(sprintf("%8.1f %11.5f %11.5f %8.2f%%\n", sg, sqrt(sum(drL^2)), sqrt(sum(drS^2)),
      100*sqrt(sum((drS-drL)^2))/sqrt(sum(drS^2))))
}
