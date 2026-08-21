library(here)
source(here("R", "enm_core.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt0<-build_relaxed(r0); wt0$Voff<-0

## step WITH rebuild + saved offset
step_rebuild <- function(st, site, dlv, d_max=10.5) {
  K<-hessian(st$r,st$pr,st$k,st$l); C<-cmat_of(nma(K))
  mask<-st$pr[,1]==site|st$pr[,2]==site
  dl<-rep(0,nrow(st$pr)); dl[mask]<-dlv
  l2<-st$l+dl
  f<-force(st$r,st$pr,st$k,l2); dr<-as.vector(C%*%f); r2<-st$r+dr
  dV <- venm(r2,st$pr,st$k,l2) - venm(st$r,st$pr,st$k,st$l)   # exact on pre-rebuild H
  nw <- build_relaxed(r2, d_max)                              # REBUILD (relaxed)
  nw$Voff <- st$Voff + dV
  nw$dV <- dV
  nw
}
set.seed(99); sA<-11; sB<-27
dlA<-rnorm(sum(wt0$pr[,1]==sA|wt0$pr[,2]==sA),0,.3)
dlB<-rnorm(sum(wt0$pr[,1]==sB|wt0$pr[,2]==sB),0,.3)

s1<-step_rebuild(wt0,sA,dlA); s2<-step_rebuild(s1,sB,dlB)
t1<-step_rebuild(wt0,sB,dlB); t2<-step_rebuild(t1,sA,dlA)

cat("=== WITH REBUILD + saved offset ===\n")
cat("A then B: dV =", s1$dV, "+", s2$dV, " total Voff =", s2$Voff, "\n")
cat("B then A: dV =", t1$dV, "+", t2$dV, " total Voff =", t2$Voff, "\n")
cat("path difference in Voff =", s2$Voff-t2$Voff, "\n")
cat("final structure rmsd    =", sqrt(mean((s2$r-t2$r)^2)), "\n")
cat("edges A-then-B:", nrow(s2$pr), "  B-then-A:", nrow(t2$pr), "\n\n")

cat("Compare with no-rebuild totals (2.5789 / 2.5788):\n")
cat("  rebuild changes the energy accounting?  A->B:", s2$Voff, " vs 2.578868\n")
