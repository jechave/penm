## Confirm the MECHANISM: for lfenm, is the whole dV error the cross term?
source("catalog.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
j <- 40; set.seed(99); mask <- wt$pr[,1]==j | wt$pr[,2]==j
dls <- lapply(1:6, function(m){d<-rep(0,nrow(wt$pr)); d[mask]<-rnorm(sum(mask),0,.3); d})

## exact dV on a background, and its three-term decomposition
parts <- function(st, dl) {
  K <- state_K(st); C <- cmat_of(nma(K)); l2 <- st$l+dl
  f <- force(st$r,st$pr,st$k,l2); dr <- as.vector(C%*%f)
  s2 <- list(r=st$r+dr,pr=st$pr,k=st$k,l=l2,frustrated=st$frustrated)
  s2$r <- state_req(s2)$r; drt <- s2$r - st$r
  strain <- dists(st$r,st$pr) - st$l
  c(exact = venm(s2$r,s2$pr,s2$k,l2)-venm(st$r,st$pr,st$k,st$l),
    quad  = 0.5*sum(st$k*dl^2),
    cross = -sum(st$k*strain*dl),
    relax = -0.5*as.numeric(crossprod(drt, K %*% drt)))
}
b0 <- sapply(dls, function(d) parts(wt, d))
st <- wt; set.seed(7)
for (n in 1:20) st <- mutate(st, sample(setdiff(1:N,j),1), sigma=.3, model="lfenm")$state
b1 <- sapply(dls, function(d) parts(st, d))
cat("lfenm background, 20 substitutions. Per state of site 40:\n\n")
cat(sprintf("%6s %10s %10s %10s %10s %12s\n","state","dV(fnd)","dV(evolved)","difference","cross term","diff-cross"))
for (m in 1:6)
  cat(sprintf("%6d %10.4f %10.4f %10.4f %10.4f %12.2e\n", m,
      b0["exact",m], b1["exact",m], b1["exact",m]-b0["exact",m],
      b1["cross",m], (b1["exact",m]-b0["exact",m]) - b1["cross",m]))
cat(sprintf("\nquadratic term identical on both backgrounds? %s\n",
    isTRUE(all.equal(unname(b0["quad",]), unname(b1["quad",])))))
cat(sprintf("relaxation term:  founder mean %.5f, evolved mean %.5f (K frozen => equal)\n",
    mean(b0["relax",]), mean(b1["relax",])))
cat(sprintf("cor(difference, cross term) = %.4f\n",
    cor(b1["exact",]-b0["exact",], b1["cross",])))
