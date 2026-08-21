## CLAIM (reviewer): V_relax <= V_stress ALWAYS, because M = sqrt(k) B' C B sqrt(k)
## is an orthogonal projector. If true this proves dV >= 0 on any relaxed network,
## for any structure -- no simulation needed.
source("sclfenm.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5)
N <- length(r0)/3; ne <- nrow(wt$pr)
K <- state_K(wt); C <- cmat_of(nma(K)); E <- evecs(wt$r, wt$pr)

## incidence-with-unit-vectors map B: (3N x nedge), so that f = -B (k dl)
B <- matrix(0, 3*N, ne)
for (e in seq_len(ne)) {
  i <- wt$pr[e,1]; j <- wt$pr[e,2]
  B[(3*i-2):(3*i), e] <-  E[,e]
  B[(3*j-2):(3*j), e] <- -E[,e]
}
sk <- sqrt(wt$k)
M <- (sk * t(B)) %*% C %*% (B * rep(sk, each=3*N))   # sqrt(k) B' C B sqrt(k)
M <- diag(sk) %*% t(B) %*% C %*% B %*% diag(sk)

cat("Is M a projector (M^2 = M)?\n")
cat(sprintf("  max|M^2 - M| = %.3e   (scale max|M| = %.3f)\n",
    max(abs(M %*% M - M)), max(abs(M))))
ev <- eigen(M, symmetric=TRUE)$values
cat(sprintf("  eigenvalues: %d near 1, %d near 0, rank = %d  (3N-6 = %d)\n",
    sum(abs(ev-1)<1e-8), sum(abs(ev)<1e-8), sum(abs(ev-1)<1e-8), 3*N-6))
cat(sprintf("  range of eigenvalues: [%.3e, %.6f]\n", min(ev), max(ev)))

cat("\nConsequence: with u = sqrt(k) dl,\n")
cat("  V_stress = 1/2 |u|^2   and   V_relax = 1/2 u' M u\n")
cat("  M projector => u'Mu <= |u|^2  => V_relax <= V_stress  => dV >= 0. QED\n\n")
## verify on random dl
set.seed(1); ok <- TRUE
for (t in 1:200) {
  site <- sample(N,1); mask <- wt$pr[,1]==site|wt$pr[,2]==site
  dl <- rep(0, ne); dl[mask] <- rnorm(sum(mask),0,.3)
  u <- sk*dl
  Vs <- 0.5*sum(u^2); Vr <- 0.5*as.numeric(t(u) %*% M %*% u)
  f <- force(wt$r,wt$pr,wt$k,wt$l+dl); dr <- as.vector(C%*%f)
  Vr2 <- 0.5*as.numeric(crossprod(dr, K%*%dr))
  if (Vr > Vs + 1e-10) ok <- FALSE
  if (t==1) cat(sprintf("  check: V_relax via M = %.6f ; via dr'Kdr = %.6f\n", Vr, Vr2))
}
cat(sprintf("  over 200 mutations: V_relax <= V_stress always? %s\n", ok))
