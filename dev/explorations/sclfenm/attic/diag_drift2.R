library(here)
source(here("R", "applications.R")); set.seed(1)
make_prot<-function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt<-new_state(r0)

## Analytic: at the CURRENT minimum r*, with current l, the state is stationary.
## Perturbing l by dl and re-minimising:
##   dV = 1/2 sum k dl^2  +  sum k (d*-l) dl   -  1/2 dr'K dr
## The middle term is the CROSS TERM: it is linear in dl, zero-mean, and it is
## the ONLY term that can be negative. It vanishes when d*=l (relaxed state).
cat("Decompose dV into its three pieces on a FRUSTRATED state.\n\n")
st<-wt; set.seed(9)
for(i in 1:5) st<-mutate(st,sample(40,1),sigma=.4,model="keepnet")$state
d <- dists(st$r,st$pr); strain <- d - st$l
cat("network strain: sum k(d-l)^2/2 =", round(0.5*sum(st$k*strain^2),4),
    "  max|d-l| =", round(max(abs(strain)),4), "\n\n")

set.seed(21); n<-400
quad<-cross<-rel<-dvv<-numeric(n)
for(i in 1:n){
  site<-sample(40,1); mask<-st$pr[,1]==site|st$pr[,2]==site
  dl<-rep(0,nrow(st$pr)); dl[mask]<-rnorm(sum(mask),0,.3)
  quad[i]  <- 0.5*sum(st$k*dl^2)
  cross[i] <- -sum(st$k*strain*dl)          # d(1/2k(d-l)^2)/dl * dl
  l2<-st$l+dl
  K<-state_K(st); C<-cmat_of(nma(K))
  f<-force(st$r,st$pr,st$k,l2); dr<-as.vector(C%*%f)
  rel[i]   <- 0.5*as.numeric(crossprod(dr,K%*%dr))
  dvv[i]   <- venm(st$r+dr,st$pr,st$k,l2)-venm(st$r,st$pr,st$k,st$l)
}
cat(sprintf("quadratic  1/2 k dl^2 : mean %8.4f  (always > 0)\n",mean(quad)))
cat(sprintf("cross  -k(d-l)dl      : mean %8.4f  sd %6.4f  frac<0 %.2f\n",
            mean(cross),sd(cross),mean(cross<0)))
cat(sprintf("relaxation -1/2dr'Kdr : mean %8.4f  (always < 0)\n",-mean(rel)))
cat(sprintf("total dV              : mean %8.4f  min %7.4f  frac<0 %.3f\n",
            mean(dvv),min(dvv),mean(dvv<0)))
cat("\nsum of parts vs dV:", round(mean(quad+cross-rel),4), "vs", round(mean(dvv),4),"\n\n")
cat("So the cross term IS negative half the time, but it is too small here to\n")
cat("outweigh the quadratic term at sigma=0.3. Strain-relieving moves need\n")
cat("either larger existing strain, or smaller sigma. Check sigma dependence:\n\n")
for (sg in c(.05,.1,.2,.3)) {
  set.seed(21); dd<-numeric(200)
  for(i in 1:200){ site<-sample(40,1); mask<-st$pr[,1]==site|st$pr[,2]==site
    dl<-rep(0,nrow(st$pr)); dl[mask]<-rnorm(sum(mask),0,sg); l2<-st$l+dl
    K<-state_K(st); C<-cmat_of(nma(K)); f<-force(st$r,st$pr,st$k,l2)
    dr<-as.vector(C%*%f); dd[i]<-venm(st$r+dr,st$pr,st$k,l2)-venm(st$r,st$pr,st$k,st$l)}
  cat(sprintf("  sigma=%.2f  mean dV=%8.4f  frac(dV<0)=%.3f\n",sg,mean(dd),mean(dd<0)))
}
