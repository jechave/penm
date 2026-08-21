## FIGURE 8: the catalogue model solves the trajectory problem
source("penm_catalog.R"); source("fig_setup.R")
cg <- readRDS("catalog.rds")

## (a) dV between states takes both signs
set.seed(1); js <- sample(cg$N, 40); dd <- c()
for (j in js) for (a in 0:cg$M) for (b in 0:cg$M) if (a!=b)
  dd <- c(dd, cg$dV[j,b+1] - cg$dV[j,a+1])
pa <- ggplot(data.frame(dV=dd), aes(dV)) +
  geom_histogram(bins=60, fill="#1b6ca8", colour="white", linewidth=.15) +
  geom_vline(xintercept=0, linewidth=.5) +
  labs(title="(a)  Energy change between states of a site",
       subtitle=sprintf("%d ordered pairs, 40 sites: %.0f%% are downhill", length(dd), 100*mean(dd<0)),
       x=expression(Delta*V(m[0]%->%m[1])), y="count")

## (b) trajectories reach a stationary state
runs <- lapply(c(0,1,2,5), function(nu)
  cbind(run_catalog_trajectory(cg, 3000, nu=nu, seed=42)$record, nu=factor(nu)))
d <- do.call(rbind, runs)
pb <- ggplot(d, aes(step, V, colour=nu)) + geom_line(linewidth=.6) +
  scale_colour_viridis_d(option="viridis", end=.85, name=expression(nu)) +
  labs(title="(b)  Trajectories reach a stationary state",
       subtitle="energy plateaus; the level is set by selection strength",
       x="proposed mutations", y="V from the founder")

## (c) simulated vs analytical equilibrium
pred <- function(nu) sum(sapply(seq_len(cg$N), function(j){
  w <- exp(-nu*cg$dV[j,]); sum(w*cg$dV[j,])/sum(w)}))
nus <- c(.25,.5,1,2,3,5,8)
sim <- sapply(nus, function(nu) mean(tail(run_catalog_trajectory(cg,30000,nu=nu,seed=7)$record$V,12000)))
pc <- ggplot(data.frame(nu=nus, sim=sim, pred=sapply(nus,pred)), aes(pred, sim)) +
  geom_abline(slope=1, linetype=2, colour="grey50") +
  geom_point(size=2.6, colour="#d1495b") +
  labs(title="(c)  The chain samples exp(-nu V)",
       subtitle="30k-step runs; agreement to 0.6-4% confirms detailed balance",
       x="predicted <V>", y="simulated <V>")
sv((pa|pb)/pc + patchwork::plot_layout(heights=c(1,1)), "fig8_catalog.png", 11, 7.5)
cat(sprintf("downhill fraction %.3f ; nu range %.2f-%.2f ; max rel err %.2f%%\n",
    mean(dd<0), min(nus), max(nus), 100*max(abs(sim-sapply(nus,pred))/sapply(nus,pred))))
