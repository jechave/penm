library(here)
source(here("R", "enm_core.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40)

run <- function(nstep, rebuild, seed=42, sigma=.3, d_max=10.5) {
  set.seed(seed)
  st <- build_relaxed(r0, d_max); st$Voff <- 0
  N <- length(r0)/3; Vs<-numeric(nstep); rms<-numeric(nstep); ne<-numeric(nstep)
  for (s in seq_len(nstep)) {
    K<-hessian(st$r,st$pr,st$k,st$l); C<-cmat_of(nma(K))
    site<-sample(N,1); mask<-st$pr[,1]==site|st$pr[,2]==site
    dl<-rep(0,nrow(st$pr)); dl[mask]<-rnorm(sum(mask),0,sigma)
    l2<-st$l+dl
    f<-force(st$r,st$pr,st$k,l2); dr<-as.vector(C%*%f); r2<-st$r+dr
    dV <- venm(r2,st$pr,st$k,l2) - venm(st$r,st$pr,st$k,st$l)
    Voff <- st$Voff + dV
    if (rebuild) { st<-build_relaxed(r2,d_max) } else { st<-list(r=r2,pr=st$pr,k=st$k,l=l2) }
    st$Voff<-Voff
    Vs[s]<-Voff; rms[s]<-sqrt(mean((r2-r0)^2)); ne[s]<-nrow(st$pr)
  }
  list(V=Vs,rmsd=rms,nedge=ne)
}
a <- run(30, rebuild=TRUE); b <- run(30, rebuild=FALSE)
cat(sprintf("%5s %14s %14s %10s %10s %8s\n","step","V_reb","V_norebuild","rmsd_reb","rmsd_nor","edges"))
for (s in c(1,2,5,10,20,30))
  cat(sprintf("%5d %14.4f %14.4f %10.4f %10.4f %8d\n", s, a$V[s], b$V[s], a$rmsd[s], b$rmsd[s], a$nedge[s]))
cat("\nedge count over run (rebuild):", paste(range(a$nedge),collapse="-"), "\n")
