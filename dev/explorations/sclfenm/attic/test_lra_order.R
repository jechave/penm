library(here)
source(here("R", "enm_core.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r_wt <- prot_2acy()$r; wt <- build_relaxed(r_wt)
K <- hessian(wt$r,wt$pr,wt$k,wt$l); C <- cmat_of(nma(K))
site <- 11; mask <- wt$pr[,1]==site | wt$pr[,2]==site
set.seed(7); base <- rnorm(sum(mask))

cat(sprintf("%8s %12s %12s %12s %10s\n","sigma","V_direct","V_decomp","abs.err","err/s^3"))
for (s in c(.025,.05,.1,.2,.4,.8)) {
  dl <- rep(0,nrow(wt$pr)); dl[mask] <- base*s
  lm_ <- wt$l+dl
  f <- force(r_wt,wt$pr,wt$k,lm_); dr <- as.vector(C%*%f); rm_ <- r_wt+dr
  Vd <- venm(rm_,wt$pr,wt$k,lm_)
  Vc <- venm(r_wt,wt$pr,wt$k,lm_) - 0.5*as.numeric(t(dr)%*%K%*%dr)
  e <- Vd-Vc
  cat(sprintf("%8.3f %12.6f %12.6f %12.3e %10.4f\n", s, Vd, Vc, e, e/s^3))
}
