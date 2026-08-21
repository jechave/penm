library(here)
source(here("R", "enm_core.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt0<-build_relaxed(r0)

## one mutation step WITHOUT rebuild (pure lfenm-style accumulation on fixed network)
step_norebuild <- function(st, site, dlv) {
  K<-hessian(st$r,st$pr,st$k,st$l); C<-cmat_of(nma(K))
  mask<-st$pr[,1]==site|st$pr[,2]==site
  dl<-rep(0,nrow(st$pr)); dl[mask]<-dlv
  l2<-st$l+dl
  f<-force(st$r,st$pr,st$k,l2); dr<-as.vector(C%*%f)
  V1<-venm(st$r,st$pr,st$k,st$l); r2<-st$r+dr
  V2<-venm(r2,st$pr,st$k,l2)
  list(r=r2,pr=st$pr,k=st$k,l=l2, dV=V2-V1)
}
set.seed(99)
sA<-11; sB<-27
dlA<-rnorm(sum(wt0$pr[,1]==sA|wt0$pr[,2]==sA),0,.3)
dlB<-rnorm(sum(wt0$pr[,1]==sB|wt0$pr[,2]==sB),0,.3)

## order A then B
s1<-step_norebuild(wt0,sA,dlA); s2<-step_norebuild(s1,sB,dlB)
## order B then A
t1<-step_norebuild(wt0,sB,dlB); t2<-step_norebuild(t1,sA,dlA)

cat("=== NO REBUILD (fixed network, l accumulates) ===\n")
cat("A then B: dV =", s1$dV, "+", s2$dV, "=", s1$dV+s2$dV, "\n")
cat("B then A: dV =", t1$dV, "+", t2$dV, "=", t1$dV+t2$dV, "\n")
cat("total dV path-independent? diff =", (s1$dV+s2$dV)-(t1$dV+t2$dV), "\n")
cat("final l identical?", isTRUE(all.equal(s2$l,t2$l)), "\n")
cat("final structure rmsd =", sqrt(mean((s2$r-t2$r)^2)), "\n")
