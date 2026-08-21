library(here)
source(here("R", "applications.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
cn <- tabulate(c(wt$pr[,1], wt$pr[,2]), nbins=N)

cat("What does scan_sites() call 'dr2'?  It is sum(mu$dr^2): the TOTAL squared\n")
cat("displacement of the WHOLE protein caused by mutating site j.\n")
cat("That is a measure of how much site j PERTURBS others -- not how much\n")
cat("site j itself MOVES. Those are different quantities.\n\n")

K <- state_K(wt); C <- cmat_of(nma(K))
set.seed(5); nmut <- 15
tot <- selfmove <- numeric(N)
for (j in 1:N) {
  a <- b <- 0
  for (m in 1:nmut) {
    mu <- mutate(wt, j, sigma=.3, model="sclfenm")
    D <- matrix(mu$dr, nrow=3)
    a <- a + sum(mu$dr^2)          # total, whole protein
    b <- b + sum(D[,j]^2)          # displacement of the mutated site itself
  }
  tot[j] <- a/nmut; selfmove[j] <- b/nmut
}
cat(sprintf("cor(total dr2 , cn) = %+.3f   <- what I plotted\n", cor(tot,cn)))
cat(sprintf("cor(site's own move, cn) = %+.3f   <- 'buried sites move less'\n", cor(selfmove,cn)))
cat("\nAlso: how compliant is each site? diag of C, per site:\n")
Cd <- diag(C); msf <- colSums(matrix(Cd, nrow=3))
cat(sprintf("cor(msf , cn) = %+.3f   <- buried sites fluctuate less (equilibrium)\n", cor(msf,cn)))
saveRDS(list(cn=cn,tot=tot,self=selfmove,msf=msf), here("data", "diag_dr2.rds"))
