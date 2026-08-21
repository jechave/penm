## FIGURE 5: single-mutation scanning -- what sclfenm gives that lfenm cannot
source("applications.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
cat("2acy chain A: N =", length(r0)/3, " edges =", nrow(wt$pr), "\n")
set.seed(2024); sc <- scan_sites(wt, nmut=15, sigma=.3, model="sclfenm", verbose=FALSE)
set.seed(2024); sl <- scan_sites(wt, nmut=15, sigma=.3, model="lfenm",  verbose=FALSE)
saveRDS(list(sclfenm=sc, lfenm=sl), "scan_results.rds")

pa <- ggplot(sc, aes(site, dV)) +
  geom_ribbon(aes(ymin=dV-dV_sd, ymax=dV+dV_sd), fill="#1b6ca8", alpha=.18) +
  geom_line(colour="#1b6ca8", linewidth=.8) + geom_point(colour="#1b6ca8", size=1.4) +
  labs(title="(a)  Energy cost of mutating each site",
       subtitle="2acy chain A, N=98; mean +/- sd over 15 mutations per site", x="site", y="dV")
pb <- ggplot(sc, aes(cn, dV)) +
  geom_point(colour="#1b6ca8", size=2, alpha=.8) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, colour="grey40", linewidth=.6) +
  labs(title="(b)  Buried sites cost more",
       subtitle=sprintf("r = %.3f", cor(sc$cn, sc$dV)), x="contact number", y="dV")
pc <- ggplot(sc, aes(cn, dr2_self)) +
  geom_point(colour="#7d669e", size=2, alpha=.8) +
  geom_smooth(method="lm", formula=y~x, se=FALSE, colour="grey40", linewidth=.6) +
  labs(title="(c)  ...and move less",
       subtitle=sprintf("displacement of the mutated site itself; r = %.3f", cor(sc$cn, sc$dr2_self)),
       x="contact number", y=expression("|"*delta*r[j]*"|"^2))
cmp <- rbind(data.frame(site=sl$site, dTS=sl$dTS, model="LFENM"),
             data.frame(site=sc$site, dTS=sc$dTS, model="SC-LFENM"))
pd <- ggplot(cmp, aes(site, dTS, colour=model)) +
  geom_hline(yintercept=0, colour="grey60", linewidth=.4) +
  geom_line(linewidth=.8) + geom_point(size=1.3) +
  scale_colour_manual(values=c(LFENM="#d1495b", "SC-LFENM"="#1b6ca8"), name=NULL) +
  labs(title="(d)  Change in conformational entropy",
       subtitle="LFENM is identically zero: it cannot see this at all",
       x="site", y="d(TS)") + theme(legend.position=c(.85,.15),
       legend.background=element_rect(fill="white",colour="grey80"))
sv((pa|pb)/(pc|pd), "fig5_scan.png", 11, 7.5)
cat(sprintf("cor(dV,cn)=%.3f  cor(dr2_self,cn)=%.3f  cor(dr2_total,cn)=%.3f\n",
    cor(sc$dV,sc$cn), cor(sc$dr2_self,sc$cn), cor(sc$dr2_total,sc$cn)))
cat(sprintf("dV identical lfenm/sclfenm: %s | dr2 identical: %s\n",
    isTRUE(all.equal(sc$dV,sl$dV)), isTRUE(all.equal(sc$dr2_self,sl$dr2_self))))
cat(sprintf("dTS lfenm max|.| = %.2e ; sclfenm range %.4f .. %.4f\n",
    max(abs(sl$dTS)), min(sc$dTS), max(sc$dTS)))
