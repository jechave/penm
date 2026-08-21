library(here)
## Background-independence of the catalogue, for ALL THREE models.
## PREDICTION for lfenm: K is frozen, so
##   dr  = C_wt f  depends only on dl                        -> EXACT
##   dV  = 1/2 sum k dl^2 - 1/2 dr'K dr - CROSS TERM
##         cross term = -sum k (d^e - l) dl, nonzero on a STRAINED background
##         and lfenm backgrounds ARE strained                -> NOT exact
source(here("R", "catalog.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3

effect <- function(st, dl) {
  K <- state_K(st); C <- cmat_of(nma(K))
  l2 <- st$l + dl
  f <- force(st$r, st$pr, st$k, l2); dr <- as.vector(C %*% f)
  s2 <- list(r=st$r+dr, pr=st$pr, k=st$k, l=l2, frustrated=st$frustrated)
  s2$r <- state_req(s2)$r
  list(dV = venm(s2$r,s2$pr,s2$k,l2)-venm(st$r,st$pr,st$k,st$l), dr = s2$r-st$r)
}
j <- 40
set.seed(99); mask <- wt$pr[,1]==j | wt$pr[,2]==j
dls <- lapply(1:6, function(m){d<-rep(0,nrow(wt$pr)); d[mask]<-rnorm(sum(mask),0,.3); d})
base <- lapply(dls, function(d) effect(wt, d))
aV <- sapply(base,`[[`,"dV"); aR <- sapply(base, function(x) sqrt(sum(x$dr^2)))

for (md in c("lfenm","keepnet","sclfenm")) {
  cat(sprintf("\n=== %s background ===\n", md))
  cat(sprintf("%8s %10s %10s %12s %12s\n","steps","strain","cor(dV)","rel err dV","rel err |dr|"))
  st <- wt; set.seed(7)
  for (n in 1:20) {
    st <- mutate(st, sample(setdiff(1:N,j),1), sigma=.3, model=md)$state
    if (n %in% c(5,10,20)) {
      now <- lapply(dls, function(d) effect(st, d))
      bV <- sapply(now,`[[`,"dV"); bR <- sapply(now, function(x) sqrt(sum(x$dr^2)))
      strain <- venm(st$r,st$pr,st$k,st$l)
      cat(sprintf("%8d %10.3f %10.5f %11.2f%% %11.3f%%\n", n, strain, cor(aV,bV),
          100*mean(abs(bV-aV)/abs(aV)), 100*mean(abs(bR-aR)/aR)))
    }
  }
}
