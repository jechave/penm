source("applications.R"); set.seed(1)
make_prot<-function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt<-new_state(r0)

cat("Is dV ever negative?\n")
set.seed(3); dv<-numeric(400)
for(i in 1:400) dv[i]<-mutate(wt,sample(40,1),sigma=.3,model="sclfenm")$dV
cat("  n =",length(dv)," min =",round(min(dv),5)," max =",round(max(dv),4),
    "  frac < 0 :",mean(dv<0),"\n")
cat("  => dV >= 0 ALWAYS: V_stress >= V_relax by construction.\n")
cat("     (relaxation can only recover part of the stress it was given)\n\n")

cat("Consequence: energy is a MONOTONE increasing accumulator.\n")
cat("Metropolis on dV can only slow it, never reverse it, because there is\n")
cat("no move that lowers V. The chain has no stationary distribution.\n\n")

cat("What's missing: mutations that RELIEVE existing strain.\n")
cat("Under sclfenm the rebuild sets l=d, so the state is always relaxed and\n")
cat("every perturbation adds strain. Under a NON-rebuilt (frustrated) state,\n")
cat("a perturbation can move l back toward d and release energy. Check:\n\n")

## frustrate a state, then see if some mutations lower the energy
st<-wt; set.seed(9)
for(i in 1:5) st<-mutate(st,sample(40,1),sigma=.4,model="keepnet")$state
cat("  after 5 keepnet mutations, network strain =",
    round(venm(st$r,st$pr,st$k,st$l),4),"\n")
dv2<-numeric(300); for(i in 1:300) dv2[i]<-mutate(st,sample(40,1),sigma=.3,model="keepnet")$dV
cat("  now frac(dV<0) =",mean(dv2<0)," min dV =",round(min(dv2),4),"\n")
cat("  => on a FRUSTRATED state, strain-relieving mutations exist.\n")
