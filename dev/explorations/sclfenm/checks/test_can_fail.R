library(here)
source(here("R", "sclfenm.R")); set.seed(1)
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
cat("Break the code on purpose; the checks must go red.\n\n")

## sabotage 1: drop the offset (the bug the whole design exists to prevent)
set.seed(11); m<-mutate(wt,11,model="sclfenm")
bad <- m$state; bad$v_off <- 0
cat("1. offset dropped -> energy reported as:", state_vmin(bad),
    " (true:", state_vmin(m$state), ")\n")
cat("   check 'energy conserved' would FAIL:", !isTRUE(all.equal(state_vmin(bad), m$dV)), "\n\n")

## sabotage 2: wrong Hessian sign (the 2017 bug)
hess_wrong <- function(r,pr,k,l,frustrated=FALSE){
  N<-length(r)/3;E<-evecs(r,pr);d<-dists(r,pr);K<-matrix(0,3*N,3*N);I3<-diag(3)
  for(e in seq_len(nrow(pr))){i<-pr[e,1];j<-pr[e,2];ee<-tcrossprod(E[,e])
    g<-if(frustrated) 1-l[e]/d[e] else 0            # WRONG SIGN
    Kij<--k[e]*(ee+g*(ee-I3));ii<-(3*i-2):(3*i);jj<-(3*j-2):(3*j)
    K[ii,jj]<-Kij;K[jj,ii]<-t(Kij)}
  for(i in seq_len(N)){ii<-(3*i-2):(3*i);b<-matrix(0,3,3)
    for(j in seq_len(N)) if(j!=i) b<-b-K[ii,(3*j-2):(3*j)];K[ii,ii]<-b};K}
rr<-r0; prr<-wt$pr; kk<-wt$k; ll<-wt$l*runif(length(wt$l),.6,1.4)  # frustrate
num <- { h<-1e-5;n<-length(rr);H<-matrix(0,n,n)
  for(a in 1:n)for(b in 1:n){pp<-rr;pp[a]<-pp[a]+h;pp[b]<-pp[b]+h
    pm<-rr;pm[a]<-pm[a]+h;pm[b]<-pm[b]-h;mp<-rr;mp[a]<-mp[a]-h;mp[b]<-mp[b]+h
    mm<-rr;mm[a]<-mm[a]-h;mm[b]<-mm[b]-h
    H[a,b]<-(venm(pp,prr,kk,ll)-venm(pm,prr,kk,ll)-venm(mp,prr,kk,ll)+venm(mm,prr,kk,ll))/(4*h*h)};H}
cat("2. Hessian check on a frustrated network:\n")
cat("   correct sign  max|err| =", max(abs(num-hessian(rr,prr,kk,ll,TRUE))), "\n")
cat("   WRONG   sign  max|err| =", max(abs(num-hess_wrong(rr,prr,kk,ll,TRUE))), " <- detected\n\n")

## sabotage 3: use K_mut instead of K_wt for the response
K<-state_K(wt);C<-cmat_of(nma(K))
cat("3. structural identity lfenm==sclfenm holds only because both use C_wt.\n")
cat("   If sclfenm relaxed with the rebuilt C, structures would differ.\n")
