library(here)
## FIGURE 2: the decomposition is exact to O(dl^3)
source(here("R", "sclfenm.R")); source(here("R", "fig_setup.R"))
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5)
K <- state_K(wt); C <- cmat_of(nma(K))
site <- 11; mask <- wt$pr[,1]==site | wt$pr[,2]==site
set.seed(7); base <- rnorm(sum(mask))
sig <- 10^seq(log10(.02), log10(1), length.out = 14)
res <- do.call(rbind, lapply(sig, function(s){
  dl <- rep(0,nrow(wt$pr)); dl[mask] <- base*s; l2 <- wt$l+dl
  f <- force(wt$r,wt$pr,wt$k,l2); dr <- as.vector(C%*%f)
  Vd <- venm(wt$r+dr, wt$pr, wt$k, l2)
  Vc <- venm(wt$r,wt$pr,wt$k,l2) - 0.5*as.numeric(crossprod(dr,K%*%dr))
  data.frame(sigma=s, exact=Vd, decomp=Vc, err=abs(Vd-Vc))}))
p1 <- ggplot(res, aes(sigma)) +
  geom_line(aes(y=exact, colour="exact  V(r_mut)"), linewidth=.9) +
  geom_point(aes(y=decomp, colour="decomposition"), size=1.9) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values=c("exact  V(r_mut)"="grey30","decomposition"="#1b6ca8"), name=NULL) +
  labs(title="Decomposition matches the exact energy",
       subtitle="dV = V_stress - V_relax, against direct evaluation of the perturbed potential",
       x=expression(sigma~"(size of "*delta*l*")"), y="energy") +
  theme(legend.position=c(.25,.85), legend.background=element_rect(fill="white",colour="grey80"))
p2 <- ggplot(res, aes(sigma, err)) +
  geom_line(colour="#d1495b", linewidth=.9) + geom_point(colour="#d1495b", size=1.6) +
  geom_line(aes(y = 0.127*sigma^3), linetype=2, colour="grey40") +
  annotate("text", x=.09, y=.127*.35^3, label=expression(0.127~sigma^3),
           colour="grey30", size=3.4, hjust=0) +
  scale_x_log10() + scale_y_log10() +
  labs(title="The error is exactly third order",
       subtitle="dashed line: a pure cubic. The only error is the linear-response truncation.",
       x=expression(sigma), y="|exact - decomposition|")
sv(p1|p2, "fig2_order.png", 10, 4)
cat("slope of log(err) vs log(sigma) =",
    round(coef(lm(log(err)~log(sigma), res))[2], 3), " (expect 3)\n")
