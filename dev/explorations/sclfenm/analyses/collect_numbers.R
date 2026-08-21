library(here)
source(here("R", "applications.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
## V2 error at sigma=0.3
set.seed(11); mu <- mutate(wt, 40, sigma=.3, model="sclfenm")
cat(sprintf("V2 err  : %.3f%%\n", 100*abs(mu$dV-(mu$V_stress-mu$V_relax))/mu$dV))
## min dV over 400 mutations from the founder
set.seed(3); dv <- numeric(400)
for (i in 1:400) dv[i] <- mutate(wt, sample(N,1), sigma=.3, model="sclfenm")$dV
cat(sprintf("400 muts: min dV = %.3f  frac<0 = %.3f\n", min(dv), mean(dv<0)))
## after 30 sclfenm steps
st <- wt; set.seed(11)
for (i in 1:30) st <- mutate(st, sample(N,1), sigma=.3, model="sclfenm")$state
cat(sprintf("after 30 sclfenm steps: max|d-l| = %.2e\n",
  max(abs(dists(st$r,st$pr)-st$l))))
set.seed(4); dv2 <- numeric(300)
for (i in 1:300) dv2[i] <- mutate(st, sample(N,1), sigma=.3, model="sclfenm")$dV
cat(sprintf("   min dV there = %.3f\n", min(dv2)))
## slopes: read fig4's saved data if present, else say so rather than error
if (!file.exists(here("data","fig4_data.rds"))) {
  cat("\n(fig4_data.rds absent -- run fig4_traj.R first for the slope table)\n")
  quit(save="no")
}
d <- readRDS(here("data", "fig4_data.rds"))
for (lv in levels(d$run)) {
  x <- subset(d, run==lv)
  cat(sprintf("%-34s slope=%.3f  n=%d  finalV=%.2f\n", lv,
      coef(lm(V_total~step, x))[2], nrow(x), x$V_total[nrow(x)]))
}
