library(here)
## FIGURE 3: why dV > 0 always under SC-LFENM
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5); N <- length(r0)/3

## (a) distribution of dV on a relaxed (sclfenm) state vs a strained (keepnet) state
st_sc <- wt
set.seed(31)
for (s in 1:30) st_sc <- mutate(st_sc, sample(N,1), sigma=.3, model="sclfenm")$state
st_kn <- wt
set.seed(9); for (s in 1:30) st_kn <- mutate(st_kn, sample(N,1), sigma=.3, model="keepnet")$state
## relax each state to its TRUE equilibrium r^e before measuring anything:
## the stored r is only the linear-response estimate for keepnet.
for (nm in c("st_sc","st_kn")) {
  st <- get(nm); rq <- state_req(st)
  cat(sprintf("%s: |F| %8.2e -> %8.2e after relaxing\n", nm,
      sqrt(sum(force(st$r,st$pr,st$k,st$l)^2)), rq$fres))
  st$r <- rq$r; assign(nm, st)
}
strain_sc <- max(abs(dists(st_sc$r,st_sc$pr)-st_sc$l))
strain_kn <- max(abs(dists(st_kn$r,st_kn$pr)-st_kn$l))

draw <- function(st, model, n=400, sigma=.3) {
  set.seed(4); v <- numeric(n)
  for (i in seq_len(n)) v[i] <- mutate(st, sample(N,1), sigma=sigma, model=model)$dV
  data.frame(dV=v, model=model)
}
d <- rbind(draw(st_sc,"sclfenm"), draw(st_kn,"keepnet"))
lab <- c(sclfenm=sprintf("SC-LFENM state after 30 mutations\nmax|d-l| = %.2f  (relaxed by the refit)", strain_sc),
         keepnet=sprintf("keepnet state after 30 mutations\nmax|d-l| = %.2f  (strain retained)", strain_kn))
d$facet <- lab[d$model]
pa <- ggplot(d, aes(dV, fill=model)) +
  geom_histogram(bins=45, colour="white", linewidth=.15) +
  geom_vline(xintercept=0, colour="black", linewidth=.5) +
  facet_wrap(~facet, ncol=1, scales="free_y") +
  scale_fill_manual(values=c(sclfenm="#1b6ca8", keepnet="#66a182"), guide="none") +
  labs(title="(a)  Energy change of 400 random mutations",
       subtitle="2acy chain A. No mutation lowers the energy of a relaxed network; some do on a strained one.",
       x="dV", y="count")

## (b) the three terms, and which can be negative
strain <- dists(st_kn$r,st_kn$pr) - st_kn$l   # now d^e - l, r^e exact
set.seed(21); n <- 400
q <- cr <- rl <- numeric(n)
for (i in seq_len(n)) {
  site<-sample(N,1); mask<-st_kn$pr[,1]==site|st_kn$pr[,2]==site
  dl<-rep(0,nrow(st_kn$pr)); dl[mask]<-rnorm(sum(mask),0,.3)
  q[i]  <- 0.5*sum(st_kn$k*dl^2)
  cr[i] <- -sum(st_kn$k*strain*dl)
  K<-state_K(st_kn); C<-cmat_of(nma(K)); l2<-st_kn$l+dl
  dr<-as.vector(C%*%force(st_kn$r,st_kn$pr,st_kn$k,l2))
  rl[i] <- -0.5*as.numeric(crossprod(dr,K%*%dr))
}
tt <- rbind(data.frame(term="quadratic\n+1/2 k dl^2", v=q),
            data.frame(term="cross\n-k(d-l) dl", v=cr),
            data.frame(term="relaxation\n-1/2 dr'K dr", v=rl))
tt$term <- factor(tt$term, levels=unique(tt$term))
pb <- ggplot(tt, aes(term, v, fill=term)) +
  geom_hline(yintercept=0, linewidth=.5) +
  geom_violin(colour=NA, alpha=.75, width=.85) +
  stat_summary(fun=mean, geom="point", size=2, colour="black") +
  scale_fill_manual(values=c("grey65","#e8a33d","#7d669e"), guide="none") +
  labs(title="(b)  The three contributions to dV",
       subtitle="Measured at the true equilibrium conformation. Only the cross term can be negative,\nand it is proportional to existing strain, so the SC-LFENM refit sets it to exactly zero.",
       x=NULL, y="contribution to dV")
sv(pa|pb, "fig3_dv_positive.png", 11, 5)
cat(sprintf("sclfenm: frac(dV<0)=%.3f min=%.3f | keepnet: frac(dV<0)=%.3f min=%.3f\n",
  mean(d$dV[d$model=="sclfenm"]<0), min(d$dV[d$model=="sclfenm"]),
  mean(d$dV[d$model=="keepnet"]<0), min(d$dV[d$model=="keepnet"])))
cat(sprintf("cross term: mean %.4f  sd %.4f  frac<0 %.2f\n", mean(cr), sd(cr), mean(cr<0)))
