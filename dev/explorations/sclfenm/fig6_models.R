## FIGURE 6: what the three models do to K, and the consequences
source("sclfenm.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
K_wt <- state_K(wt); nm_wt <- nma(K_wt)
res <- do.call(rbind, lapply(c("lfenm","keepnet","sclfenm"), function(md){
  set.seed(11); mu <- mutate(wt, 11, sigma=.3, model=md)
  Km <- state_K(mu$state); nmm <- nma(Km)
  n <- min(length(nm_wt$evalue), length(nmm$evalue))
  data.frame(model=md,
    dK = norm(Km-K_wt,"F"),
    dstruct = sqrt(mean((mu$state$r-wt$r)^2)),
    dV = mu$dV,
    d_eval = mean(abs(nmm$evalue[1:n]-nm_wt$evalue[1:n])),
    strain_in_net = venm(mu$state$r,mu$state$pr,mu$state$k,mu$state$l),
    v_off = mu$state$v_off,
    d_edges = mu$n_edge_after - mu$n_edge_before)}))
res$model <- factor(res$model, levels=c("lfenm","keepnet","sclfenm"))
long <- rbind(
  data.frame(model=res$model, q="||K_mut - K_wt||\n(does the Hessian move?)", v=res$dK),
  data.frame(model=res$model, q="RMSD of structure\n(identical in all three)", v=res$dstruct),
  data.frame(model=res$model, q="strain left in the network", v=res$strain_in_net),
  data.frame(model=res$model, q="energy in the saved offset", v=res$v_off))
p <- ggplot(long, aes(model, v, fill=model)) +
  geom_col(width=.6) +
  geom_text(aes(label=ifelse(v==0,"0",sprintf("%.3g",v))), vjust=-.4, size=3.1) +
  facet_wrap(~q, scales="free_y", ncol=4) +
  scale_fill_manual(values=c(lfenm="#d1495b", keepnet="#66a182", sclfenm="#1b6ca8"),
                    guide="none") +
  expand_limits(y=0) +
  labs(title="The three models differ only in what happens to K",
       subtitle="Same mutation (site 11) under each. Structure is identical; the Hessian and the energy bookkeeping are not.",
       x=NULL, y=NULL)
sv(p, "fig6_three_models.png", 11, 3.8)
print(res, row.names=FALSE, digits=4)
