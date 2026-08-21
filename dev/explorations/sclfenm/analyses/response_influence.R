library(here)
## Response vs influence, computed analytically (exact mean over mutations).
## R[i,j] = < |dr_i|^2 >  for random mutations at site j.
##   response profile of a mutation at j : column R[,j]
##   response  at site i (marginal over j): rowSums(R)   or rowMeans(R)
##   influence of site j (marginal over i): colSums(R)
source(here("R", "applications.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
K <- state_K(wt); C <- cmat_of(nma(K)); E <- evecs(wt$r, wt$pr)
blk <- function(M,i,j) M[(3*i-2):(3*i),(3*j-2):(3*j)]
sigma <- 0.3
R <- matrix(0,N,N)
for (e in seq_len(nrow(wt$pr))) {           # each contact, perturbed independently
  a<-wt$pr[e,1]; b<-wt$pr[e,2]; ev<-E[,e]; varf<-(wt$k[e]*sigma)^2
  for (i in 1:N) {
    r <- varf*sum(((blk(C,i,b)-blk(C,i,a)) %*% ev)^2)
    R[i,a] <- R[i,a] + r                     # this edge belongs to site a...
    R[i,b] <- R[i,b] + r                     # ...and to site b
  }
}
cn <- tabulate(c(wt$pr[,1],wt$pr[,2]), nbins=N)
resp <- rowSums(R); infl <- colSums(R)
cat("2acy chain A, analytical (no sampling)\n\n")
cat(sprintf("  cor( response(i) , cn_i ) = %+.3f\n", cor(resp,cn)))
cat(sprintf("  cor( influence(j), cn_j ) = %+.3f\n", cor(infl,cn)))
cat(sprintf("  cor( response, influence ) = %+.3f\n", cor(resp,infl)))
cat(sprintf("\n  sum(R) = %.4f ; sum(resp) = %.4f ; sum(infl) = %.4f  (all equal)\n",
    sum(R), sum(resp), sum(infl)))
cat(sprintf("\n  R asymmetry: max|R-t(R)|/max(R) = %.3f\n", max(abs(R-t(R)))/max(R)))
cat("  R is NOT symmetric, and should not be: mutating j forces all cn_j of its\n")
cat("  contacts, so column sums scale with cn_j while row sums do not.\n")
saveRDS(list(R=R,cn=cn,resp=resp,infl=infl), here("data", "response_influence.rds"))
