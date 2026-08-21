## The catalogue energy sums single-site dV over sites. That is a SEPARATE
## approximation from the difference identity (which is exact). Quantify it:
## compare sum_j dV(0->j,s_j) against the exact dV of the combined sequence.
source("penm_catalog.R")
cg <- readRDS("catalog.rds")
cat("Site-additivity of the catalogue energy (2acy, ANM):\n")
cat("  exact = dV of the sequence, all perturbations applied together\n")
cat("  catalogue = sum of single-site dV\n\n")
cat(sprintf("%8s %12s %12s %10s %10s\n","sites","catalogue","exact","abs err","rel err"))
set.seed(11)
for (ns in c(1,2,5,10,20,40,98)) {
  ct <- ex <- numeric(4)
  for (rep in 1:4) {
    s <- integer(cg$N)
    js <- sample(cg$N, ns); s[js] <- sample(1:cg$M, ns, replace=TRUE)
    ct[rep] <- seq_dV_catalog(cg, s); ex[rep] <- seq_dV_exact(cg, s)
  }
  cat(sprintf("%8d %12.3f %12.3f %10.3f %9.1f%%\n", ns, mean(ct), mean(ex),
      mean(ct-ex), 100*mean((ct-ex)/ex)))
}
