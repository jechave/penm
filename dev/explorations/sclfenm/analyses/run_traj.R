library(here)
source(here("R", "applications.R"))
set.seed(1)
make_prot<-function(N=40){t<-seq(0,6*pi,length.out=N)
  as.vector(rbind(4*cos(t)+rnorm(N,0,.4),4*sin(t)+rnorm(N,0,.4),1.2*t/2+rnorm(N,0,.4)))}
r0<-make_prot(40); wt<-new_state(r0)

cat("=== NEUTRAL (no selection): energy must run away ===\n")
tn <- run_trajectory(wt, nstep=60, sigma=.3, model="sclfenm", selection="neutral", seed=7)
r<-tn$record
cat(sprintf("  step %3d  V=%8.3f rmsd=%6.3f edges=%d\n", r$step[c(1,10,30,60)],
    r$V_total[c(1,10,30,60)], r$rmsd[c(1,10,30,60)], r$n_edge[c(1,10,30,60)]))
cat("  V grows without bound: slope per step =", round(coef(lm(V_total~step,r))[2],4), "\n\n")

cat("=== STABILITY SELECTION, nu = 1 ===\n")
ts1 <- run_trajectory(wt, nstep=60, sigma=.3, nu=1, model="sclfenm",
                      selection="step", seed=7)
r1<-ts1$record
cat(sprintf("  step %3d  V=%8.3f rmsd=%6.3f edges=%d acc=%.3f\n", r1$step[c(1,10,30,60)],
    r1$V_total[c(1,10,30,60)], r1$rmsd[c(1,10,30,60)], r1$n_edge[c(1,10,30,60)],
    r1$accept_rate[c(1,10,30,60)]))
cat("  slope per step =", round(coef(lm(V_total~step,r1))[2],4),
    "  overall acceptance =", round(ts1$n_acc/ts1$n_try,4), "\n\n")

cat("=== STABILITY SELECTION, nu = 4 (stronger) ===\n")
ts4 <- run_trajectory(wt, nstep=60, sigma=.3, nu=4, model="sclfenm",
                      selection="step", seed=7)
r4<-ts4$record
cat(sprintf("  step %3d  V=%8.3f rmsd=%6.3f\n", r4$step[c(1,10,30,60)],
    r4$V_total[c(1,10,30,60)], r4$rmsd[c(1,10,30,60)]))
cat("  slope per step =", round(coef(lm(V_total~step,r4))[2],4),
    "  overall acceptance =", round(ts4$n_acc/ts4$n_try,4), "\n\n")

cat("=== summary: does selection bound the energy? ===\n")
cat(sprintf("  %-22s %10s %10s %10s\n","regime","V(60)","rmsd(60)","accept"))
for (nm in c("neutral","nu=1","nu=4")) {
  rr <- switch(nm, "neutral"=r, "nu=1"=r1, "nu=4"=r4)
  aa <- switch(nm, "neutral"=1, "nu=1"=ts1$n_acc/ts1$n_try, "nu=4"=ts4$n_acc/ts4$n_try)
  cat(sprintf("  %-22s %10.3f %10.3f %10.3f\n", nm, rr$V_total[60], rr$rmsd[60], aa))
}
saveRDS(list(neutral=r,b1=r1,b4=r4), here("data", "trajectories.rds"))
