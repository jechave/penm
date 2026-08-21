## FIGURE 9: with vs without frustration
source("frustration.R"); source("fig_setup.R")
wt <- founder(); N <- get_nsites(wt)
set.seed(3); mu <- wt; rows <- list()
for (n in 1:40) {
  mu <- get_mutant_site(mu, sample(N,1), n, "lfenm", .3, 1L, 1L)
  if (n %in% c(1,2,3,5,8,12,16,20,25,30,35,40)) {
    p <- relax_to_min(mu, TRUE)
    g <- get_graph(p); str <- 0.5*sum(g$kij*(g$dij-g$lij)^2)
    Kf<-hess(p,FALSE); Kt<-hess(p,TRUE); sf<-spec(Kf); st<-spec(Kt)
    of <- order(sf$val)[1:10]; ot <- order(st$val)[1:10]
    ov <- crossprod(sf$vec[,of], st$vec[,ot])
    rows[[length(rows)+1]] <- data.frame(n=n, strain=str,
      dK=100*max(abs(Kf-Kt))/max(abs(Kt)),
      deval=100*mean(abs(sort(st$val)[1:10]-sort(sf$val)[1:10])/sort(sf$val)[1:10]),
      rmsip=sqrt(sum(ov^2)/10), dTS=ts_of(st$val)-ts_of(sf$val))
  }
}
d <- do.call(rbind, rows); saveRDS(d, "frustration_scan.rds")
pa <- ggplot(d, aes(strain, dK)) + geom_line(colour="#d1495b", linewidth=.8) +
  geom_point(colour="#d1495b", size=2) +
  labs(title="(a)  The Hessian differs, and the gap grows with strain",
       x="strain energy accumulated", y=expression("||"*K[frust]-K[relaxed]*"||  (% of scale)"))
pb <- ggplot(d, aes(strain, 1-rmsip)) + geom_line(colour="#1b6ca8", linewidth=.8) +
  geom_point(colour="#1b6ca8", size=2) +
  labs(title="(b)  But the soft modes barely rotate",
       subtitle="1 - RMSIP over the 10 softest modes; 0 would be identical",
       x="strain energy accumulated", y="1 - RMSIP")
pc <- ggplot(d, aes(strain, dTS)) + geom_hline(yintercept=0, colour="grey60") +
  geom_line(colour="#66a182", linewidth=.8) + geom_point(colour="#66a182", size=2) +
  labs(title="(c)  Entropy difference is small and unsigned",
       x="strain energy accumulated", y=expression(T*S[frust]-T*S[relaxed]))
sv(pa|pb|pc, "fig9_frustration.png", 12, 3.8)
print(d, row.names=FALSE, digits=3)
