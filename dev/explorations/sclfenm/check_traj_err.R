source("applications.R"); source("fig_setup.R")
r0 <- prot_2acy()$r; wt <- new_state(r0, d_max = 10.5); N <- length(r0)/3
cat("Along a keepnet trajectory, how far is the stored r from r^e,\n")
cat("and how wrong is the energy used by the threshold criterion?\n\n")
cat(sprintf("%5s %10s %14s %14s %12s\n","step","|F| at r","V(stored r)","V(r^e) exact","difference"))
st <- wt; set.seed(3)
for (n in 1:20) {
  st <- mutate(st, sample(N,1), sigma=.3, model="keepnet")$state
  if (n %% 5 == 0) {
    f <- sqrt(sum(force(st$r,st$pr,st$k,st$l)^2))
    va <- state_vmin(st); vb <- state_vmin(st, exact=TRUE)
    cat(sprintf("%5d %10.4f %14.6f %14.6f %12.2e\n", n, f, va, vb, va-vb))
  }
}
