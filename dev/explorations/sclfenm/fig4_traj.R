## FIGURE 4: trajectories -- sclfenm walks into a wall, keepnet does not
source("applications.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5); N <- length(r0)/3
cat("N=30, edges =", nrow(wt$pr), "\n")

runs <- list()
runs$neutral <- run_trajectory(wt, nstep=50, sigma=.3, model="sclfenm",
                               selection="neutral", seed=7)
runs$step_b1 <- run_trajectory(wt, nstep=50, sigma=.3, nu=1, model="sclfenm",
                               selection="step", seed=7)
runs$step_b4 <- run_trajectory(wt, nstep=50, sigma=.3, nu=4, model="sclfenm",
                               selection="step", seed=7)
runs$thr_sc  <- suppressWarnings(run_trajectory(wt, nstep=50, sigma=.3, v_cut=25,
                    model="sclfenm", selection="threshold", seed=7, max_try=800))
runs$thr_kn  <- suppressWarnings(run_trajectory(wt, nstep=50, sigma=.3, v_cut=25,
                    model="keepnet", selection="threshold", seed=7, max_try=800))
lab <- c(neutral="SC-LFENM, no selection", step_b1="SC-LFENM, Metropolis (nu=1)",
         step_b4="SC-LFENM, Metropolis (nu=4)",
         thr_sc="SC-LFENM, threshold V<25", thr_kn="keepnet, threshold V<25")
d <- do.call(rbind, lapply(names(runs), function(n){
  r <- runs[[n]]$record; if(!nrow(r)) return(NULL)
  data.frame(r, run=lab[n], stalled=runs[[n]]$stalled)}))
d$run <- factor(d$run, levels=lab)
ends <- do.call(rbind, lapply(split(d, d$run), function(x) x[nrow(x),]))

pa <- ggplot(d, aes(step, V_total, colour=run)) +
  geom_hline(yintercept=25, linetype=2, colour="grey55") +
  annotate("text", x=2, y=26.5, label="V_cut = 25", colour="grey40", size=3, hjust=0) +
  geom_line(linewidth=.85) +
  geom_point(data=subset(ends, stalled), size=3, shape=4, stroke=1.2) +
  scale_colour_manual(values=c("grey50","#9c4f8e","#d1495b","#1b6ca8","#66a182"), name=NULL) +
  labs(title="(a)  Total energy along an evolutionary trajectory",
       subtitle="2acy chain A.  x = trajectory stalls: no acceptable mutation exists any more",
       x="substitutions accepted", y="V from the founder") +
  theme(legend.position="right")
pb <- ggplot(d, aes(step, rmsd, colour=run)) +
  geom_line(linewidth=.85) +
  scale_colour_manual(values=c("grey50","#9c4f8e","#d1495b","#1b6ca8","#66a182"), guide="none") +
  labs(title="(b)  Structural drift from the founder", x="substitutions accepted", y="RMSD")
sv(pa/pb, "fig4_trajectories.png", 9, 7)

cat("\n")
for (n in names(runs)) {
  r <- runs[[n]]$record
  cat(sprintf("%-32s steps %2d/50 %-8s finalV=%6.2f acc=%.4f\n",
      lab[n], nrow(r), if(runs[[n]]$stalled) "STALLED" else "ran ok",
      if(nrow(r)) r$V_total[nrow(r)] else NA, runs[[n]]$n_acc/max(runs[[n]]$n_try,1)))
}
saveRDS(d, "fig4_data.rds")
