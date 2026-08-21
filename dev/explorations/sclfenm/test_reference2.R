source("enm_core.R"); set.seed(1)
make_prot <- function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r_wt<-make_prot(40); wt<-build_relaxed(r_wt)
K<-hessian(wt$r,wt$pr,wt$k,wt$l); C<-cmat_of(nma(K))
site<-11; mask<-wt$pr[,1]==site|wt$pr[,2]==site
set.seed(7); dl<-rep(0,nrow(wt$pr)); dl[mask]<-rnorm(sum(mask),0,.3)
l_mut<-wt$l+dl
f<-force(r_wt,wt$pr,wt$k,l_mut); dr<-as.vector(C%*%f); r_mut<-r_wt+dr

Vs_exact <- venm(r_wt, wt$pr, wt$k, l_mut)      # exact: at wt geometry
Vs_form  <- 0.5*sum(wt$k*dl^2)                  # closed form
cat("V_stress exact (venm at r_wt) =", Vs_exact, "\n")
cat("V_stress closed form 1/2k dl^2=", Vs_form, "\n")
cat("equal? ", isTRUE(all.equal(Vs_exact,Vs_form)), "  diff =", Vs_exact-Vs_form, "\n")
cat("  (equal because d(r_wt)=l_wt exactly, so d-l_mut = -dl)\n\n")

Vrelax_K  <- 0.5*as.numeric(t(dr)%*%K%*%dr)
Vrelax_f  <- 0.5*as.numeric(t(dr)%*%f)          # = 1/2 f'C f, the LRA work
cat("V_relax via 1/2 dr'K dr =", Vrelax_K, "\n")
cat("V_relax via 1/2 dr'f    =", Vrelax_f, "\n")
cat("  (these agree iff K dr = f exactly; dr=Cf so K C f = f only on non-rigid space)\n")
cat("  diff =", Vrelax_K-Vrelax_f, "\n\n")

Vd <- venm(r_mut,wt$pr,wt$k,l_mut)
cat("EXACT dV = V_mut(r_mut)      =", Vd, "\n")
cat("approx  = Vs - 1/2 dr'K dr   =", Vs_exact-Vrelax_K, "  err =", Vd-(Vs_exact-Vrelax_K), "\n")
cat("approx  = Vs - 1/2 dr'f      =", Vs_exact-Vrelax_f, "  err =", Vd-(Vs_exact-Vrelax_f), "\n")
