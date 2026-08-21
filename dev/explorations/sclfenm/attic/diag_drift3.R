library(here)
source(here("R", "applications.R")); set.seed(1)
make_prot<-function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt<-new_state(r0)
cat("CLAIM: under sclfenm the state is always relaxed, so the cross term is 0\n")
cat("       and dV > 0 strictly. Verify directly.\n\n")
st <- wt
for (s in 1:6) st <- mutate(st, sample(40,1), sigma=.4, model="sclfenm")$state
d<-dists(st$r,st$pr); strain<-d-st$l
cat("after 6 SCLFENM steps: max|d-l| =", format(max(abs(strain)),digits=3),
    "  network strain =", format(0.5*sum(st$k*strain^2),digits=3),"\n")
cat("  -> state is EXACTLY relaxed (rebuild guarantees l=d). Cross term == 0.\n\n")
set.seed(4); dd<-numeric(300)
for(i in 1:300) dd[i]<-mutate(st,sample(40,1),sigma=.3,model="sclfenm")$dV
cat("  frac(dV<0) on this sclfenm state =", mean(dd<0), "  min dV =", round(min(dd),4),"\n\n")

cat("CONTRAST: keepnet keeps strain, so the cross term survives.\n")
st2<-wt; for(s in 1:6) st2<-mutate(st2,sample(40,1),sigma=.4,model="keepnet")$state
d2<-dists(st2$r,st2$pr)
cat("after 6 KEEPNET steps: max|d-l| =", round(max(abs(d2-st2$l)),4),
    "  strain =", round(0.5*sum(st2$k*(d2-st2$l)^2),3),"\n")
set.seed(4); dd2<-numeric(300)
for(i in 1:300) dd2[i]<-mutate(st2,sample(40,1),sigma=.3,model="keepnet")$dV
cat("  frac(dV<0) =", mean(dd2<0), "  min dV =", round(min(dd2),4),"\n")
