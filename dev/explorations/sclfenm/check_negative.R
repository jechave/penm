source("sclfenm.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max=10.5); N <- length(r0)/3
cat("Do strain-relieving mutations exist on 2acy? Depends on accumulated strain.\n\n")
cat(sprintf("%6s %12s %12s %12s %10s\n","nmut","strain","min dV","mean dV","frac<0"))
st <- wt; set.seed(9)
for (n in 1:40) {
  st <- mutate(st, sample(N,1), sigma=.3, model="keepnet")$state
  if (n %in% c(5,10,20,30,40)) {
    rq <- state_req(st); st2 <- st; st2$r <- rq$r
    strain <- 0.5*sum(st2$k*(dists(st2$r,st2$pr)-st2$l)^2)
    set.seed(77); dv <- numeric(250)
    for (i in 1:250) dv[i] <- mutate(st2, sample(N,1), sigma=.3, model="keepnet")$dV
    cat(sprintf("%6d %12.3f %12.4f %12.4f %10.3f\n", n, strain, min(dv), mean(dv), mean(dv<0)))
  }
}
