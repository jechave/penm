library(here)
## The reviewer's charge: the catalogue's dV for a move is exact only in a FIXED
## background. Used mid-trajectory it inherits the site-additivity error.
## Test the single-move claim directly on an evolved background.
source(here("R", "penm_catalog.R"))
cg <- readRDS(here("data", "catalog.rds")); wt <- cg$wt; N <- cg$N

## build an evolved sequence, then compare catalogue dV vs EXACT dV for a move
set.seed(21); s <- integer(N); js <- sample(N, 30); s[js] <- sample(1:cg$M, 30, TRUE)
cat("Background: 30 sites mutated.\n\n")
cat(sprintf("%6s %8s %14s %14s %12s\n","site","move","catalogue dV","exact dV","error"))
set.seed(5); bad <- 0; n <- 12
for (t in 1:n) {
  j <- sample(js,1); m0 <- s[j]; m1 <- sample(setdiff(0:cg$M,m0),1)
  cat_dV <- cg$dV[j,m1+1] - cg$dV[j,m0+1]
  s1 <- s; s1[j] <- m1
  ex_dV <- seq_dV_exact(cg, s1) - seq_dV_exact(cg, s)
  if (sign(cat_dV) != sign(ex_dV)) bad <- bad + 1
  cat(sprintf("%6d %4d->%2d %14.5f %14.5f %12.5f%s\n", j, m0, m1, cat_dV, ex_dV,
      cat_dV-ex_dV, if (sign(cat_dV)!=sign(ex_dV)) "  <- SIGN WRONG" else ""))
}
cat(sprintf("\nsign disagreements: %d/%d\n", bad, n))
