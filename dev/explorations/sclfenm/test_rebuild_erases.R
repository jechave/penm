source("enm_core.R"); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r_wt <- make_prot(40); wt <- build_relaxed(r_wt)
K <- hessian(wt$r,wt$pr,wt$k,wt$l); C <- cmat_of(nma(K))
site<-11; mask <- wt$pr[,1]==site|wt$pr[,2]==site
set.seed(7); dl<-rep(0,nrow(wt$pr)); dl[mask]<-rnorm(sum(mask),0,.3)
l_mut <- wt$l+dl
f <- force(r_wt,wt$pr,wt$k,l_mut); dr<-as.vector(C%*%f); r_mut<-r_wt+dr

cat("=== the pre-rebuild (frustrated) mutant ===\n")
cat("V_min                      =", venm(r_mut,wt$pr,wt$k,l_mut), "\n")
cat("residual strain sum(k(d-l)^2)/2 =",
    0.5*sum(wt$k*(dists(r_mut,wt$pr)-l_mut)^2), "\n")
cat("max |d - l| over edges     =", max(abs(dists(r_mut,wt$pr)-l_mut)), "\n\n")

cat("=== after REBUILDING a relaxed network at r_mut ===\n")
reb <- build_relaxed(r_mut, d_max = 10.5)
cat("V_min of rebuilt network   =", venm(reb$r,reb$pr,reb$k,reb$l), "  <-- ZERO: strain erased\n")
cat("edges: wt =", nrow(wt$pr), " rebuilt =", nrow(reb$pr), "\n")
same <- nrow(wt$pr)==nrow(reb$pr) && all(wt$pr==reb$pr)
cat("contact set unchanged?      ", same, "\n\n")

cat("=== so dV would be reported as ===\n")
cat("  rebuild-only (WRONG):", venm(reb$r,reb$pr,reb$k,reb$l) - venm(r_wt,wt$pr,wt$k,wt$l), "\n")
cat("  with saved offset   :", venm(r_mut,wt$pr,wt$k,l_mut) - venm(r_wt,wt$pr,wt$k,wt$l), "\n")
