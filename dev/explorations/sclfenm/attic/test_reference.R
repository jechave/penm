library(here)
source(here("R", "enm_core.R")); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r_wt <- make_prot(40); wt <- build_relaxed(r_wt)
K <- hessian(wt$r,wt$pr,wt$k,wt$l); C <- cmat_of(nma(K))
site<-11; mask<-wt$pr[,1]==site|wt$pr[,2]==site
set.seed(7); dl<-rep(0,nrow(wt$pr)); dl[mask]<-rnorm(sum(mask),0,.3)
l_mut<-wt$l+dl
f<-force(r_wt,wt$pr,wt$k,l_mut); dr<-as.vector(C%*%f); r_mut<-r_wt+dr

d_wt <- dists(r_wt,wt$pr); d_mut <- dists(r_mut,wt$pr)
cat("Both are 'the mutant energy' but w.r.t. different zero:\n\n")
cat("(a) V_mut(r_mut) with l_mut, zero at l_mut relaxed :",
    0.5*sum(wt$k*(d_mut-l_mut)^2), "\n")
cat("(b) V_stress - V_relax                             :",
    0.5*sum(wt$k*dl^2) - 0.5*as.numeric(t(dr)%*%K%*%dr), "\n\n")
cat("(a) is the strain of the MUTANT Hamiltonian at its own minimum.\n")
cat("(b) is the WORK done: cost at wt geometry, minus what relaxation gave back.\n\n")
cat("They are NOT equal. Difference:", 0.5*sum(wt$k*(d_mut-l_mut)^2) -
    (0.5*sum(wt$k*dl^2)-0.5*as.numeric(t(dr)%*%K%*%dr)), "\n\n")

cat("Which is dV = V_mut(r_mut) - V_wt(r_wt)?\n")
cat("  V_wt(r_wt) =", 0.5*sum(wt$k*(d_wt-wt$l)^2), " (zero: wt is relaxed)\n")
cat("  => dV = (a) =", 0.5*sum(wt$k*(d_mut-l_mut)^2), "\n\n")
cat("So (a) is the energy difference. (b) equals (a) only to O(dl^3)?  Check:\n")
cat("   (a) - (b) =", 0.5*sum(wt$k*(d_mut-l_mut)^2)-(0.5*sum(wt$k*dl^2)-0.5*as.numeric(t(dr)%*%K%*%dr)),"\n")
cat("   note earlier V_direct - V_decomp was ~0.0016; here it's larger =>\n")
cat("   the two 'decompositions' are NOT the same thing. Investigate.\n")
